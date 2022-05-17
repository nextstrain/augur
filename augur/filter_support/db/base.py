import abc
import os
from typing import Any, Callable, Dict, List, Set, Tuple
import sys
from tempfile import NamedTemporaryFile
import argparse
from augur.errors import AugurError
from augur.index import index_sequences, index_vcf
from augur.io import open_file, print_err, read_sequences, write_sequences
from augur.io.db import get_metadata_id_column
from augur.io import is_vcf, write_vcf
from augur.filter_support.exceptions import FilterException
from augur.filter_support.subsample import calculate_sequences_per_group, TooManyGroupsError, get_valid_group_by_cols


DUMMY_COL = 'dummy'
SEQUENCE_ONLY_FILTERS = (
    "min_length",
    "non_nucleotide",
)

FilterCallableReturn = Tuple[str, Dict[str, Any]]
FilterCallable = Callable[..., FilterCallableReturn]
FilterOption = Tuple[FilterCallable, Dict[str, Any]]


class FilterBase(abc.ABC):
    def __init__(self, args:argparse.Namespace):
        self.args = args

    def try_run(self):
        """Try running filter. If there is an error, remove the database file and exit."""
        try:
            self.run()
        except AugurError as e:
            self.db_cleanup()
            raise e

    def run(self, cleanup=True):
        self.set_metadata_columns()
        self.db_load_metadata()
        self.add_attributes()
        if self.use_sequences:
            self.create_sequence_index_table()
        self.db_create_date_table()
        self.include_exclude_filter()
        if self.do_subsample:
            self.subsample()
        self.db_create_output_table()
        self.write_outputs()
        self.write_report()
        if cleanup:
            self.db_cleanup()

    def set_metadata_columns(self):
        """Set column names for strain and date, to be used for metadata and intermediate database tables."""
        self.metadata_id_column = get_metadata_id_column(self.args.metadata, self.args.metadata_id_columns)
        self.date_column = 'date'
        self.db_set_sanitized_identifiers()

    @abc.abstractmethod
    def db_set_sanitized_identifiers(self):
        """Sets sanitized names for externally sourced identifiers."""
        pass

    @abc.abstractmethod
    def db_load_metadata(self):
        """Loads a metadata file into the database."""
        pass

    @abc.abstractmethod
    def db_load_sequence_index(self):
        """Loads a sequence index file into the database."""
        pass

    def add_attributes(self):
        """Check if there is a date column and if sequences are used."""
        self.has_date_col = self.db_has_date_col()
        self.use_sequences = bool(self.args.sequence_index or (self.args.sequences and not self.args.exclude_all))
        self.do_subsample = bool(self.args.group_by or self.args.subsample_max_sequences)

    @abc.abstractmethod
    def db_has_date_col(self) -> bool:
        """Returns a boolean indicating whether `self.date_column` is in the metadata."""
        pass

    def create_sequence_index_table(self):
        """Load sequence index into table, generating if necessary."""
        build_sequence_index = (self.args.sequence_index is None)

        if build_sequence_index:
            # Generate the sequence index on the fly, for backwards compatibility
            # with older workflows that don't generate the index ahead of time.
            # Create a temporary index using a random filename to avoid collisions
            # between multiple filter commands.
            with NamedTemporaryFile(delete=False) as sequence_index_file:
                self.args.sequence_index = sequence_index_file.name

            print_err(
                "Note: You did not provide a sequence index, so Augur will generate one.",
                "You can generate your own index ahead of time with `augur index` and pass it with `augur filter --sequence-index`."
            )

            if is_vcf(self.args.sequences):
                index_vcf(self.args.sequences, self.args.sequence_index)
            else:
                index_sequences(self.args.sequences, self.args.sequence_index)

        self.db_load_sequence_index()

        # Remove temporary index file, if it exists.
        if build_sequence_index:
            os.unlink(self.args.sequence_index)

        self.sequence_strains = self.db_get_sequence_index_strains()

    @abc.abstractmethod
    def db_get_sequence_index_strains(self):
        """Returns the set of all strains in the sequence index."""
        pass

    @abc.abstractmethod
    def db_create_date_table(self):
        """Creates an intermediate date table from the metadata table."""
        pass

    def include_exclude_filter(self):
        """Construct and apply filters to create intermediate filter reason table."""
        exclude_by, include_by = self.construct_filters()
        self.db_create_filter_reason_table(exclude_by, include_by)

    def construct_filters(self):
        """Construct lists of exclude and force-include filter expressions."""
        exclude_by:List[FilterOption] = []
        include_by:List[FilterOption] = []

        # Force include sequences specified in file(s).
        if self.args.include:
            # Collect the union of all given strains to include.
            for include_file in self.args.include:
                include_by.append((self.force_include_strains, {'include_file': include_file}))

        # Add sequences with particular metadata attributes.
        if self.args.include_where:
            for include_where in self.args.include_where:
                include_by.append((self.force_include_where, {'include_where': include_where}))

        # Exclude all strains by default.
        if self.args.exclude_all:
            exclude_by.append((self.filter_by_exclude_all, {}))

        # Filter by sequence index.
        if self.use_sequences:
            exclude_by.append((self.filter_by_sequence_index, {}))

        # Remove strains explicitly excluded by name.
        if self.args.exclude:
            for exclude_file in self.args.exclude:
                exclude_by.append((self.filter_by_exclude_strains, {'exclude_file': exclude_file}))

        # Exclude strain my metadata field like 'host=camel'.
        if self.args.exclude_where:
            for exclude_where in self.args.exclude_where:
                exclude_by.append((self.filter_by_exclude_where, {'exclude_where': exclude_where}))

        # Exclude strains by metadata, using SQL querying.
        if self.args.query:
            exclude_by.append((self.filter_by_query, {'query': self.args.query}))

        # TODO: check if no date column but filters require it
        if self.has_date_col:
            # Filter by ambiguous dates.
            if self.args.exclude_ambiguous_dates_by:
                exclude_by.append((self.filter_by_ambiguous_date, {'ambiguity': self.args.exclude_ambiguous_dates_by}))

            # Filter by date.
            if self.args.min_date:
                exclude_by.append((self.filter_by_min_date, {'min_date': self.args.min_date}))
            if self.args.max_date:
                exclude_by.append((self.filter_by_max_date, {'max_date': self.args.max_date}))

        # Filter by sequence length.
        if self.args.min_length:
            if is_vcf(self.args.sequences):
                print("WARNING: Cannot use min_length for VCF files. Ignoring...")
            else:
                exclude_by.append((self.filter_by_sequence_length, {'min_length': self.args.min_length}))

        if self.args.group_by:
            # Ambiguous year exclusions must be evaluated after month since month will capture ambiguous years as well.
            # Note that the order here is reversed later.
            if "year" in self.args.group_by:
                exclude_by.append((self.skip_group_by_with_ambiguous_year, {}))
            if "month" in self.args.group_by:
                exclude_by.append((self.skip_group_by_with_ambiguous_month, {}))

        # Exclude sequences with non-nucleotide characters.
        if self.args.non_nucleotide:
            exclude_by.append((self.filter_by_non_nucleotide, {}))

        # Reversed so that the order of precedence will be the same as the pandas engine.
        return list(reversed(exclude_by)), include_by

    @abc.abstractmethod
    def filter_by_exclude_all(self) -> FilterCallableReturn: pass

    @abc.abstractmethod
    def filter_by_exclude_strains(self, exclude_file) -> FilterCallableReturn: pass

    @abc.abstractmethod
    def parse_filter_query(self, query): pass

    @abc.abstractmethod
    def filter_by_exclude_where(self, exclude_where) -> FilterCallableReturn: pass

    @abc.abstractmethod
    def filter_by_query(self, query) -> FilterCallableReturn: pass

    @abc.abstractmethod
    def filter_by_ambiguous_date(self, ambiguity="any") -> FilterCallableReturn: pass

    @abc.abstractmethod
    def filter_by_min_date(self, min_date:float) -> FilterCallableReturn: pass

    @abc.abstractmethod
    def filter_by_max_date(self, max_date:float) -> FilterCallableReturn: pass

    @abc.abstractmethod
    def filter_by_sequence_index(self) -> FilterCallableReturn: pass

    @abc.abstractmethod
    def filter_by_sequence_length(self, min_length=0) -> FilterCallableReturn: pass

    @abc.abstractmethod
    def filter_by_non_nucleotide(self) -> FilterCallableReturn: pass

    @abc.abstractmethod
    def force_include_strains(self, include_file) -> FilterCallableReturn: pass

    @abc.abstractmethod
    def force_include_where(self, include_where) -> FilterCallableReturn: pass

    def skip_group_by_with_ambiguous_month(self) -> FilterCallableReturn:
        """Alias to filter_by_ambiguous_date(ambiguity="month") with a specific function name for filter reason."""
        return self.filter_by_ambiguous_date(ambiguity="month")

    def skip_group_by_with_ambiguous_year(self) -> FilterCallableReturn:
        """Alias to filter_by_ambiguous_date(ambiguity="year") with a specific function name for filter reason."""
        return self.filter_by_ambiguous_date(ambiguity="year")

    @abc.abstractmethod
    def db_create_filter_reason_table(self, exclude_by:List[FilterOption], include_by:List[FilterOption]):
        """Creates an intermediate table for filter reason."""
        pass

    @abc.abstractmethod
    def db_create_output_table(self, input_table:str):
        """Creates a final intermediate table to be used for output."""
        pass

    def subsample(self):
        """Apply subsampling to update filter reason table."""
        self.create_priorities_table()

        valid_group_by_cols = []
        if self.args.group_by:
            metadata_cols = self.db_get_metadata_cols()
            valid_group_by_cols = get_valid_group_by_cols(self.args.group_by, metadata_cols)
        self.db_create_extended_filtered_metadata_table(valid_group_by_cols)

        if self.args.subsample_max_sequences:
            if self.args.group_by:
                counts_per_group = self.db_get_counts_per_group(valid_group_by_cols)
            else:
                valid_group_by_cols = [DUMMY_COL]
                counts_per_group = [self.db_get_filtered_strains_count()]

            try:
                sequences_per_group, probabilistic_used = calculate_sequences_per_group(
                    self.args.subsample_max_sequences,
                    counts_per_group,
                    allow_probabilistic=self.args.probabilistic_sampling
                )
            except TooManyGroupsError as error:
                raise FilterException(str(error)) from error

            if (probabilistic_used):
                print(f"Sampling probabilistically at {sequences_per_group:0.4f} sequences per group, meaning it is possible to have more than the requested maximum of {self.args.subsample_max_sequences} sequences after filtering.")
            else:
                print(f"Sampling at {sequences_per_group} per group.")
        else:
            sequences_per_group = self.args.sequences_per_group

        self.db_create_group_sizes_table(valid_group_by_cols, sequences_per_group)
        self.db_update_filter_reason_table_with_subsampling(valid_group_by_cols)

    @abc.abstractmethod
    def db_get_counts_per_group(self, group_by_cols:List[str]) -> List[int]:
        """
        Returns
        -------
        list[int]
            List of counts per group.
        """
        pass

    @abc.abstractmethod
    def db_get_filtered_strains_count(self) -> int:
        """Returns the number of metadata strains that pass all filter rules.

        Note: this can return a different number before and after subsampling.
        """
        pass

    @abc.abstractmethod
    def db_update_filter_reason_table_with_subsampling(self, group_by_cols:List[str]):
        """Subsamples filtered metadata and updates the filter reason table."""
        pass

    def create_priorities_table(self):
        """Load or generate priorities table."""
        if self.args.priority:
            self.db_load_priorities_table()
        else:
            self.db_generate_priorities_table(self.args.subsample_seed)

    @abc.abstractmethod
    def db_load_priorities_table(self):
        """Loads a priorities file into the database."""
        pass

    @abc.abstractmethod
    def db_generate_priorities_table(self, seed:int=None):
        """Generates a priorities table with random priorities."""
        pass

    @abc.abstractmethod
    def db_get_metadata_cols(self) -> Set[str]:
        """Returns a set of metadata column names."""
        pass

    @abc.abstractmethod
    def db_create_extended_filtered_metadata_table(self, group_by_cols:List[str]):
        """Creates a new table with rows as filtered metadata, with additional columns."""
        pass

    @abc.abstractmethod
    def db_create_group_sizes_table(self, group_by_cols:List[str], sequences_per_group:float):
        """Creates an intermediate table for group sizes."""
        pass

    def write_outputs(self):
        """Write various outputs."""
        if self.args.output:
            self.read_and_output_sequences()
        if self.args.output_strains:
            self.db_output_strains()
        if self.args.output_metadata:
            self.db_output_metadata()
        if self.args.output_log:
            self.db_output_log()

    def read_and_output_sequences(self):
        """Read sequences and output all that passed filtering."""
        valid_strains = self.db_get_strains_passed()

        # Write output starting with sequences, if they've been requested. It is
        # possible for the input sequences and sequence index to be out of sync
        # (e.g., the index is a superset of the given sequences input), so we need
        # to update the set of strains to keep based on which strains are actually
        # available.
        if is_vcf(self.args.sequences):
            if self.args.output:
                # Get the samples to be deleted, not to keep, for VCF
                dropped_samps = list(self.sequence_strains - valid_strains)
                write_vcf(self.args.sequences, self.args.output, dropped_samps)
        elif self.args.sequences:
            sequences = read_sequences(self.args.sequences)

            # If the user requested sequence output, stream to disk all sequences
            # that passed all filters to avoid reading sequences into memory first.
            # Even if we aren't emitting sequences, we track the observed strain
            # names in the sequence file as part of the single pass to allow
            # comparison with the provided sequence index.
            if self.args.output:
                observed_sequence_strains = set()
                with open_file(self.args.output, "wt") as output_handle:
                    for sequence in sequences:
                        observed_sequence_strains.add(sequence.id)

                        if sequence.id in valid_strains:
                            write_sequences(sequence, output_handle, 'fasta')
            else:
                observed_sequence_strains = {sequence.id for sequence in sequences}

            if self.use_sequences and self.sequence_strains != observed_sequence_strains:
                # Warn the user if the expected strains from the sequence index are
                # not a superset of the observed strains.
                if self.sequence_strains is not None and observed_sequence_strains > self.sequence_strains:
                    print_err(
                        "WARNING: The sequence index is out of sync with the provided sequences.",
                        "Metadata and strain output may not match sequence output."
                    )

                # Update the set of available sequence strains.
                self.sequence_strains = observed_sequence_strains

    @abc.abstractmethod
    def db_output_strains(self):
        """Writes the final output strains to a file, one strain per line.
        Row order from the original metadata file is retained.
        """
        pass

    @abc.abstractmethod
    def db_output_metadata(self):
        """Writes the final output strains to a file, one strain per line,
        along with the original metadata.
        Row order from the original metadata file is retained.
        """
        pass

    @abc.abstractmethod
    def db_output_log(self):
        """Writes a file explaining the reason for excluded or force-included strains.

        This file has the following columns:
        1. Original strain column
        2. Name of the filter function responsible for inclusion/exclusion
        3. Arguments given to the filter function
        """
        pass

    def write_report(self):
        """Output a report of how many strains were dropped and reasoning."""
        total_strains_passed = self.db_get_filtered_strains_count()
        num_excluded_by_lack_of_metadata = self.get_num_excluded_by_lack_of_metadata()
        num_metadata_strains = self.db_get_num_metadata_strains()
        num_excluded_subsamp = self.db_get_num_excluded_by_subsampling()
        filter_counts = self.db_get_filter_counts()

        # Calculate the number of strains passed and filtered.
        total_strains_filtered = num_metadata_strains + num_excluded_by_lack_of_metadata - total_strains_passed

        print(f"{total_strains_filtered} strains were dropped during filtering")

        if num_excluded_by_lack_of_metadata:
            print(f"\t{num_excluded_by_lack_of_metadata} had no metadata")

        report_template_by_filter_name = {
            self.filter_by_sequence_index.__name__: "{count} had no sequence data",
            self.filter_by_exclude_all.__name__: "{count} of these were dropped by `--exclude-all`",
            self.filter_by_exclude_strains.__name__: "{count} of these were dropped because they were in {exclude_file}",
            self.filter_by_exclude_where.__name__: "{count} of these were dropped because of '{exclude_where}'",
            self.filter_by_query.__name__: "{count} of these were filtered out by the query: \"{query}\"",
            self.filter_by_ambiguous_date.__name__: "{count} of these were dropped because of their ambiguous date in {ambiguity}",
            self.filter_by_min_date.__name__: "{count} of these were dropped because they were earlier than {min_date} or missing a date",
            self.filter_by_max_date.__name__: "{count} of these were dropped because they were later than {max_date} or missing a date",
            self.filter_by_sequence_length.__name__: "{count} of these were dropped because they were shorter than minimum length of {min_length}bp",
            self.filter_by_non_nucleotide.__name__: "{count} of these were dropped because they had non-nucleotide characters",
            self.skip_group_by_with_ambiguous_year.__name__: "{count} were dropped during grouping due to ambiguous year information",
            self.skip_group_by_with_ambiguous_month.__name__: "{count} were dropped during grouping due to ambiguous month information",
            self.force_include_strains.__name__: "{count} strains were added back because they were in {include_file}",
            self.force_include_where.__name__: "{count} sequences were added back because of '{include_where}'",
        }
        for filter_name, filter_kwargs, count in filter_counts:
            if filter_kwargs:
                import json
                parameters = dict(json.loads(filter_kwargs))
            else:
                parameters = {}

            parameters["count"] = count
            print("\t" + report_template_by_filter_name[filter_name].format(**parameters))

        if self.do_subsample:
            seed_txt = ", using seed {}".format(self.args.subsample_seed) if self.args.subsample_seed else ""
            print("\t%i of these were dropped because of subsampling criteria%s" % (num_excluded_subsamp, seed_txt))

        if total_strains_passed == 0:
            raise FilterException("All samples have been dropped! Check filter rules and metadata file format.")

        print(f"{total_strains_passed} strains passed all filters")

    @abc.abstractmethod
    def db_get_metadata_strains(self) -> Set[str]:
        """Returns a set of all strain names from the original metadata."""
        pass

    @abc.abstractmethod
    def db_get_strains_passed(self) -> Set[str]:
        """Returns a set of all strain names that passed filtering and subsampling."""
        pass

    def get_num_excluded_by_lack_of_metadata(self):
        """Get number of strains excluded present in sequences but missing in metadata."""
        metadata_strains = self.db_get_metadata_strains()
        if self.use_sequences:
            return len(self.sequence_strains - metadata_strains)
        return 0

    @abc.abstractmethod
    def db_get_num_metadata_strains(self) -> int:
        """Returns the number of strains in the original metadata."""
        pass

    @abc.abstractmethod
    def db_get_num_excluded_by_subsampling(self) -> int:
        """Returns the number of strains excluded by subsampling."""
        pass

    @abc.abstractmethod
    def db_get_filter_counts(self) -> List[Tuple[str, str, int]]:
        """Returns a list of tuples for each filter function that had an effect.

        Each tuple has:
        1. Name of the filter function
        2. Arguments given to the filter function
        3. Number of strains included/excluded by the filter function
        """
        pass

    @abc.abstractmethod
    def db_cleanup(self): pass
