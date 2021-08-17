"""
Subsample an alignment
"""
from .utils import AugurException
from .filter import run as augur_filter
from .priorities import get_distance_to_focal_set, create_priorities
from .index import index_sequences
from .io import write_sequences, open_file, read_sequences, read_metadata
import yaml
from argparse import Namespace
from os import path
import pandas as pd
from tempfile import NamedTemporaryFile
import jsonschema
from pkg_resources import resource_string
import shlex

def register_arguments(parser):
    parser.add_argument('--scheme', required=True, metavar="YAML", help="subsampling scheme")
    parser.add_argument('--output-dir', required=True, metavar="PATH", help="directory to save intermediate results")
    parser.add_argument('--metadata', required=True, metavar="TSV", help="metadata")
    parser.add_argument('--alignment', required=True, metavar="FASTA", help="alignment to subsample")
    parser.add_argument('--alignment-index', required=False, metavar="INDEX", help="sequence index of alignment")
    parser.add_argument('--reference', required=True, metavar="FASTA", help="reference (which was used for alignment)")
    parser.add_argument('--include-strains-file', required=False, nargs="+", default=None, metavar="TXT", help="strains to force include")
    parser.add_argument('--exclude-strains-file', required=False, nargs="+", default=None, metavar="TXT", help="strains to force exclude")
    parser.add_argument('--output-fasta', required=True, metavar="FASTA", help="output subsampled sequences")
    parser.add_argument('--output-metadata', required=True, metavar="TSV", help="output subsampled metadata")
    parser.add_argument('--output-log', required=False, metavar="TSV", help="log file explaining why strains were excluded / included")

def run(args):
    
    config = parse_scheme(args.scheme)

    generate_sequence_index(args)

    samples = [Sample(name, data, args) for name, data in config.items()]

    graph = make_graph(samples)

    traverse_graph(
        graph,
        lambda s: s.filter()
    )

    combine_samples(args, samples)

def parse_scheme(filename):
    with open(filename) as fh:
        try:
            data = yaml.safe_load(fh)
        except yaml.YAMLError as exc:
            print(exc)
            raise AugurException(f"Error parsing subsampling scheme {filename}")
    validate_scheme(data)
    return data


def validate_scheme(scheme):
    try:
        # with open(resource_string(__package__, path.join("data", "schema-subsampling.yaml"))) as fh:
        schema = yaml.safe_load(resource_string(__package__, path.join("data", "schema-subsampling.yaml")))
    except yaml.YAMLError as err:
        raise AugurException("Subsampling schema definition is not valid YAML. Error: {}".format(err))
    # check loaded schema is itself valid -- see http://python-jsonschema.readthedocs.io/en/latest/errors/
    try:
        jsonschema.Draft6Validator.check_schema(schema)
    except jsonschema.exceptions.SchemaError as err:
        raise AugurException("Subsampling schema definition is not valid. Error: {}".format(path, err))

    try:
        jsonschema.Draft6Validator(schema).validate(scheme)
    except jsonschema.exceptions.ValidationError as err:
        print(err)
        raise AugurException("Subsampling scheme failed validation")


class Sample():
    """
    A class to hold information about a sample. A subsampling scheme will consist of multiple
    samples. Each sample may depend on the priorities based off another sample.
    """
    def __init__(self, name, config, cmd_args):
        self.name = name
        self.tmp_dir = cmd_args.output_dir
        self.alignment = cmd_args.alignment
        self.alignment_index = cmd_args.alignment_index
        self.reference = cmd_args.reference
        self.metadata = cmd_args.metadata
        self.initialise_filter_args(config, cmd_args)
        self.priorities = config.get("priorities", None) 
        print("Constructor", self.name)
    
    def initialise_filter_args(self, config, subsample_args):
        """
        Currently this method is needed as we need to call `augur filter`'s `run()` with an
        argparse instance. An improvement here would be to expose appropriate filtering
        functions and call them as needed, with the output being `return`ed rather than
        written to disk.
        """
        args = Namespace()
        args.metadata = self.metadata
        args.sequences = self.alignment
        args.sequence_index = self.alignment_index
        args.metadata_chunk_size = 100000
        args.metadata_id_columns = ["strain", "name"]
        args.min_date = None
        args.max_date = None
        args.exclude_ambiguous_dates_by = None
        args.exclude = subsample_args.exclude_strains_file
        args.exclude_all = None
        args.include = subsample_args.include_strains_file
        args.include_where = None
        args.min_length = None
        args.non_nucleotide = None
        args.probabilistic_sampling = None
        args.no_probabilistic_sampling = None
        args.priority = None
        args.subsample_seed = None
        args.output          = path.join(self.tmp_dir, f"sample.{self.name}.fasta") # filtered sequences in FASTA forma
        args.output_metadata = path.join(self.tmp_dir, f"sample.{self.name}.tsv")   # metadata for strains that passed filters
        args.output_strains  = path.join(self.tmp_dir, f"sample.{self.name}.txt")   # list of strains that passed filters (no header)
        args.output_log      = path.join(self.tmp_dir, f"sample.{self.name}.log.tsv")
        # translate the subsampling config (YAML) into arguments for `augur filter`
        # TODO - we should allow schemas to use a data structure identical to augur filter here, for instance
        # `group-by: [year month]` rather than the current ncov structure which is `group_by: "year month"`
        # right now we can use a bit of both, but this is not complete!
        if "group-by" in config:
            args.group_by = config['group-by']
        elif "group_by" in config: # ncov syntax
            args.group_by = shlex.split(config['group_by'])
        else:
            args.group_by = None
        
        if "sequences-per-group" in config:
            args.sequences_per_group = config['sequences-per-group']
        elif "seq_per_group" in config: # ncov syntax
            args.sequences_per_group = config['seq_per_group']
        else:
            args.sequences_per_group = None

        if "query" in config:
            if config['query'].startswith("--query "):
                args.query = config['query'].lstrip("--query ").strip('"').strip("'")
            else:
                args.query = config['query']
        else:
            args.query = None
        
        if "exclude-where" in config:
            args.exclude_where = config['exclude-where'] # list
        elif "exclude" in config and config['exclude'].startswith("--exclude-where"): # ncov syntax for exclude-where
            args.exclude_where = shlex.split(config['exclude'])[1:]
        else:
            args.exclude_where = None
        self.filter_args = args

        if "subsample-max-sequences" in config:
            args.subsample_max_sequences = config["subsample-max-sequences"]
        elif "max_sequences" in config: # ncov syntax
            args.subsample_max_sequences = config["max_sequences"]
        else:
            args.subsample_max_sequences = None

        if "sampling_scheme" in config:
            # ncov syntax for expressing additional filter arguments
            if config["sampling_scheme"]=="--probabilistic-sampling":
                args.probabilistic_sampling = True
            elif config["sampling_scheme"]=="--no-probabilistic-sampling":
                args.no_probabilistic_sampling = True

        # ncov syntax options still to implement:
        # include_argument = _get_specific_subsampling_setting("include", optional=True),
        # exclude_ambiguous_dates_argument = _get_specific_subsampling_setting("exclude_ambiguous_dates_by", optional=True),
        # min_date = _get_specific_subsampling_setting("min_date", optional=True),
        # max_date = _get_specific_subsampling_setting("max_date", optional=True),

    def calculate_required_priorities(self):
        """
        If computation of this sample requires priority information of another sample
        (the "focus"), then this function will compute those priorities.
        """
        if not self.priorities:
            return
        focal_sample = self.priorities.get('sample', None)
        if not focal_sample:
            raise AugurException(f"Cannot calculate priorities needed for {self.name} as the {self.get_priority_focus_name()} sample wasn't linked")
        print(f"Calculating priorities required by {self.name}")
        focal_sample.calculate_priorities()

    def calculate_priorities(self):
        """
        Calculate the priorities TSV file for the alignment in the context of this sample
        """
        print(f"Calculating proximity of {self.name}")
        get_distance_to_focal_set(
            self.alignment,
            self.reference,
            self.filter_args.output,
            path.join(self.tmp_dir, f"proximity_{self.name}.tsv"),
            ignore_seqs=["Wuhan/Hu-1/2019"]  # TODO - use the config to define this?
        )
        print(f"Calculating priorities of {self.name}")
        create_priorities(
            self.alignment_index,
            path.join(self.tmp_dir, f"proximity_{self.name}.tsv"),
            path.join(self.tmp_dir, f"priorities_{self.name}.tsv")
        )

    def get_priority_focus_name(self):
        if not self.priorities:
            return None
        return self.priorities['focus']

    def set_priority_sample(self, sample):
        if not self.priorities:
            raise AugurException(f"No priorities set for {self.name}")
        self.priorities['sample'] = sample

    def filter(self):
        print("\n---------------------------------\nCONSTRUCTING SAMPLE FOR", self.name, "\n---------------------------------")
        self.calculate_required_priorities()
        augur_filter(self.filter_args)
        # In the future, instead of `augur_filter` saving data to disk, it would return
        # data to the calling process. In lieu of that, we read the data just written.
        self.sampled_strains = set(pd.read_csv(self.filter_args.output_strains, header=None)[0])
        self.filter_log = pd.read_csv(
            self.filter_args.output_log,
            header=0,
            sep="\t",
            index_col="strain"
        )


def make_graph(samples):
    """"
    Given a config file, construct a graph of samples to perform in an iterative fashion, such that
    priorities 

    This is a DAG, however an extremely simple one which we can construct outselves rather than relying on
    extra libraries.

    Constraints:
    * Each sample can only use priorities of one other sample
    * Acyclic

    Structure:
    tuple: (sample name, list of descendent samples) where a "descendant" sample requires the linked sample to be
    created prior to it's creation. Each entry in the list has this tuple structure.
    """

    included = set() # set of samples added to graph
    graph = (None, [])

    # add all the samples which don't requre priorities to the graph
    for sample in samples:
        if not sample.get_priority_focus_name():
            graph[1].append((sample, []))
            included.add(sample.name)

    def add_descendants(level):
        parent_sample = level[0]
        descendants = level[1]
        for sample in samples:
            if sample.name in included:
                continue
            if sample.get_priority_focus_name() == parent_sample.name:
                sample.set_priority_sample(parent_sample)
                descendants.append((sample, []))
                included.add(sample.name)
        for inner_level in descendants:
            add_descendants(inner_level)

    for level in graph[1]:
        add_descendants(level)

    # from pprint import pprint
    # print("\ngraph"); pprint(graph);print("\n")

    if len(samples)!=len(included):
        AugurException("Incomplete graph construction")

    return graph

def traverse_graph(level, callback):
    this_sample, descendents = level
    # print("LEVEL", this_sample, len(descendents))
    if this_sample:
        callback(this_sample)
    for child in descendents:
        traverse_graph(child, callback)

def generate_sequence_index(args):
    if args.alignment_index:
        print("Skipping sequence index creation as an index was provided")
        return    
    print("Creating ephemeral sequence index file")
    with NamedTemporaryFile(delete=False) as sequence_index_file:
        sequence_index_path = sequence_index_file.name
        index_sequences(args.alignment, sequence_index_path)
        args.alignment_index = sequence_index_path


def combine_samples(args, samples):
    """Collect the union of strains which are included in each sample and write them to disk.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments from argparse
    samples : list[Sample]
        list of samples

    """
    print("\n\n")
    ### Form a union of each sample set, which is the subsampled strains list
    sampled_strains = set()
    for sample in samples:
        print(f"Sample \"{sample.name}\" included {len(sample.sampled_strains)} strains")
        sampled_strains.update(sample.sampled_strains)
    print(f"In total, {len(sampled_strains)} strains are included in the resulting subsampled dataset")

    ## Iterate through the input sequences, streaming a subsampled version to disk.
    sequences = read_sequences(args.alignment)
    sequences_written_to_disk = 0
    with open_file(args.output_fasta, "wt") as output_handle:
        for sequence in sequences:
            if sequence.id in sampled_strains:
                sequences_written_to_disk += 1
                write_sequences(sequence, output_handle, 'fasta')
    print(f"Writing {sequences_written_to_disk} sequences to {args.output_fasta}")

    ## Iterate through the metadata in chunks, writing out those entries which are in the subsample
    metadata_reader = read_metadata(
        args.metadata,
        id_columns=["strain", "name"], # TODO - this should be an argument
        chunk_size=10000 # TODO - argument
    )
    metadata_header = True
    metadata_mode = "w"
    metadata_written_to_disk = "???" # TODO
    for metadata in metadata_reader:
        metadata.loc[sampled_strains].to_csv(
            args.output_metadata,
            sep="\t",
            header=metadata_header,
            mode=metadata_mode,
        )
        metadata_header = False
        metadata_mode = "a"
    print(f"Writing {metadata_written_to_disk} metadata entries sequences to {args.output_metadata}")

    ## Combine the log files (from augur filter) for each sample into a larger log file
    ## Format TBD
    ## TODO
