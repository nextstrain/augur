"""
Subsample a dataset according to a provided config.

This augur functionality is _alpha_ and may change at any time.
It does not conform to the semver release standards used by Augur.
"""

import argparse
from typing import Dict, List, Optional
from .errors import AugurError
from os import path
import subprocess
import tempfile

INCOMPLETE = 'incomplete'
COMPLETE = 'complete'

def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("subsample", help=__doc__)
    parser.add_argument('--config', required=True, metavar="YAML", help="subsampling config file") # TODO: allow a string representation
    parser.add_argument('--metadata', required=True, metavar="TSV", help="sequence metadata")
    parser.add_argument('--sequences', required=True, metavar="FASTA", help="sequences in FASTA format") # TODO XXX VCF ?
    parser.add_argument('--output-metadata', required=True, metavar="TSV", help="output metadata")
    parser.add_argument('--output-sequences', required=True, metavar="FASTA", help="output sequences in FASTA format") # TODO XXX VCF ?

    optionals = parser.add_argument_group(
        title="Optional arguments",
    )
    optionals.add_argument('--dry-run', action="store_true")
    optionals.add_argument('--metadata-id-columns', metavar="NAME", nargs="+",
        help="names of possible metadata columns containing identifier information, ordered by priority. Only one ID column will be inferred.")
    optionals.add_argument('--random-seed', type=int, metavar="N",
        help="random number generator seed to allow reproducible subsampling (with same input data).")
    optionals.add_argument('--reference', metavar="FASTA", help="needed for priority calculation (but it shouldn't be!)")
    optionals.add_argument('--exclude', type=str, nargs="+", help="file(s) with list of strains to exclude")
    optionals.add_argument('--include', type=str, nargs="+", help="file(s) with list of strains to include regardless of priorities, subsampling, or absence of an entry in --sequences.")

    return parser

def parse_config(filename):
    import yaml
    with open(filename) as fh:
        try:
            config = yaml.safe_load(fh)
        except yaml.YAMLError as e:
            print(e)
            raise AugurError(f"Error parsing subsampling scheme {filename}")
    # TODO XXX - write a schema and validate against this
    if 'samples' not in config:
        raise AugurError('Config must define a "samples" key')
    return config


File = str

class Sample:
    output_strains: Optional[File]

    def __init__(self, name, dependencies):
        self.name = name
        self.dependencies = dependencies
        self.status = INCOMPLETE


# TODO: Move these classes closer to code for augur filter?
PandasQuery = str
Date = str
EqualityFilterQuery = str


class FilterSample(Sample):

    # Inputs
    # This is not optional but marking it as such for easier implementation.
    metadata: Optional[File]
    sequences: Optional[File]
    metadata_id_columns: Optional[List[str]]

    # Outputs
    output_metadata: Optional[File]
    output_sequences: Optional[File]
    output_strains: Optional[File]  # FIXME: this is redundant

    # Filters
    query: Optional[PandasQuery]
    min_date: Optional[Date]
    max_date: Optional[Date]
    exclude: Optional[List[File]]
    exclude_where: Optional[List[EqualityFilterQuery]]
    exclude_all: Optional[bool]
    include: Optional[List[File]]
    include_where: Optional[List[EqualityFilterQuery]]
    min_length: Optional[str]
    non_nucleotide: Optional[bool]

    # Sampling options
    group_by: Optional[List[str]]
    max_sequences: Optional[int]
    sequences_per_group: Optional[int]
    disable_probabilistic_sampling: Optional[bool]
    random_seed: Optional[int]

    def __init__(self, name=None, dependencies=None, **kwargs):
        super().__init__(name=name, dependencies=dependencies)

        # Initialize instance attributes
        for option in self.__annotations__.keys():
            self.__setattr__(option, None)

        # Populate instance attributes
        for option, value in kwargs.items():
            if option not in self.__annotations__:
                raise AugurError(f'Option {option!r} not allowed.')
            # TODO: Check types
            self.__setattr__(option, value)

    def args(self):
        args = ['augur', 'filter']

        # Inputs

        if self.metadata is not None:
            args.extend(['--metadata', self.metadata])

        if self.sequences is not None:
            args.extend(['--sequences', self.sequences])

        if self.metadata_id_columns is not None and len(self.metadata_id_columns) > 0:
            args.append('--metadata-id-columns')
            args.extend(self.metadata_id_columns)

        # Outputs

        if self.output_metadata is not None:
            args.extend(['--output-metadata', self.output_metadata])

        if self.output_sequences is not None:
            args.extend(['--output-sequences', self.output_sequences])

        if self.output_strains is not None:
            args.extend(['--output-strains', self.output_strains])

        # Filters

        if self.query is not None:
            args.extend(['--query', self.query])

        if self.min_date is not None:
            args.extend(['--min-date', self.min_date])

        if self.exclude is not None and len(self.exclude) > 0:
            args.append('--exclude')
            args.extend(self.exclude)

        if self.exclude_where is not None and len(self.exclude_where) > 0:
            args.append('--exclude-where')
            args.extend(self.exclude_where)

        if self.exclude_all is not None and self.exclude_all:
            args.append('--exclude-all')

        if self.include is not None and len(self.include) > 0:
            args.append('--include')
            args.extend(self.include)

        if self.include_where is not None and len(self.include_where) > 0:
            args.append('--include-where')
            args.extend(self.include_where)

        # Sampling options

        if self.group_by is not None:
            args.append('--group-by')
            args.extend(self.group_by)

        if self.max_sequences is not None:
            args.extend(['--subsample-max-sequences', str(self.max_sequences)])

        if self.sequences_per_group is not None:
            args.extend(['--sequences-per-group', str(self.sequences_per_group)])

        if self.disable_probabilistic_sampling:
            args.append('--no-probabilistic-sampling')

        if self.random_seed is not None:
            args.extend(['--subsample-seed', str(self.random_seed)])

        return args

    def exec(self, dry_run=False):
        """
        Instead of running an `augur filter` command in a subprocess like we do
        here, a nicer way would be to refactor `augur filter` to expose
        functions which can be called here so that we can get a list of returned
        strains. Doing so would provide a HUGE speedup if we could index the
        sequences+metadata on disk a single time and then use that for all
        filtering calls which subsampling performs. This is analogous to having
        an in-memory database running in a separate process we can query - not
        as fast, but much easier implementation. Using subprocess does make
        parallelisation trivial, but the above speed up would be preferable.
        """
        deps = (f"depends on {', '.join(self.dependencies)}") if self.dependencies else "no dependencies"
        print(f"Sampling for {self.name!r} ({deps})")
        for option in self.__annotations__:
            if self.__getattribute__(option) and option not in {'metadata', 'sequences', 'metadata_id_columns', 'output_metadata', 'output_sequences', 'output_strains'}:
                print(f'\t{option}: {self.__getattribute__(option)}')
        print()
        print(' '.join(str(arg) for arg in self.args()))
        print()

        if not dry_run:
            try:
                subprocess.run([str(arg) for arg in self.args()])
            except subprocess.CalledProcessError as e:
                raise AugurError(e)

        self.status = COMPLETE


class ProximitySample(Sample):
    reference: File
    focal_sequences: File
    context_sequences: File
    output_strains: File
    num_per_focal: int

    def __init__(self, name=None, dependencies=None, **kwargs):
        super().__init__(name=name, dependencies=dependencies)

        # Initialize instance attributes
        for option in self.__annotations__.keys():
            self.__setattr__(option, None)

        # Populate instance attributes
        for option, value in kwargs.items():
            if option not in self.__annotations__:
                raise AugurError(f'Option {option!r} not allowed.')
            # TODO: Check types
            self.__setattr__(option, value)

    def get_closest_sequences(self):
        """
        The function we call here (originally written for ncov) computes the minimum
        hamming distance for every sample in the (background) dataset against
        all focal strains, and then returns the _n_ closest for each strain
        """
        from augur.subsample_.get_distance_to_focal_set import get_distance_to_focal_set 
        print("Computing hamming distances & choosing closest contextual strains...")
        # the sequences are assumed to be aligned, an error will be thrown if the lengths vary
        return get_distance_to_focal_set(
            self.context_sequences,
            self.reference,
            self.focal_sequences,
            self.num_per_focal
        )

    def exec(self, dry_run=False):
        print()
        print(f"Calculate proximity for {self.name} by computing weights for {self.focal_data['sequences']} against {self.contextual_data['sequences']}")

        if not dry_run:
            closest = self.get_closest_sequences()
            with open(self.output_strains, 'w') as fh:
                print("\n".join(list(closest)), file=fh)
            print(f"\tClosest strains written to {self.output_strains} (n={len(closest)})")

        self.status = COMPLETE


# FIXME: Determine whether/where to name samples as "calls" internally. The
# concept of calls shouldn't be exposed in the config, but it's key to
# understanding the implementation. Allowing recursive samples would make the
# split more clear since samples and calls will no longer be synonymous.

def generate_calls(config: dict, args: argparse.Namespace, tmpdir):
    """
    Produce an (unordered) dictionary of calls to be made to accomplish the
    desired subsampling. Each call is either (i) a use of augur filter or (ii) a
    proximity calculation The names given to calls are config-defined, but there
    is guaranteed to be one call with the name "output".

    The separation between this function and the Filter (etc) classes is not
    quite right, but it's a WIP.
    """
    samples: Dict[str, Sample] = dict()

    total_weights = sum(sample_config['weight']
                        for name, sample_config in config['samples'].items()
                        if 'weight' in sample_config)

    for name, sample_config in config['samples'].items():
        sample: Sample = ...
        if 'priorities' in sample_config:
            # Assume all priority samples are proximity
            assert sample_config['priorities']['type'] == 'proximity'
            sample = ProximitySample(
                name=name,
                dependencies=sample_config['focal'],
                sequences=args.sequences,
                output_strains=path.join(tmpdir, f"{name}.samples.txt"),
                **sample_config,
            )
            samples[sample.name] = sample
        else:
            # Assume all non-priority samples are filter samples

            # Determine intermediate sample size
            if 'weight' in sample_config:
                max_sequences = int(config['size'] * (sample_config['weight'] / total_weights))
                del sample_config['weight']
            elif 'max_sequences' in sample_config:
                max_sequences = sample_config['max_sequences']
                del sample_config['max_sequences']
            else:
                max_sequences = None

            # Remapping keys
            if 'exclude' in sample_config:
                sample_config['exclude_where'] = sample_config['exclude']
                del sample_config['exclude']

            # Apply global exclusions at every intermediate sample
            exclude = []
            if args.exclude:
                exclude.extend(args.exclude)

            # output_strains is only necessary for inclusion in final output sample
            # output_sequences is only necessary for any downstream proximity samples
            # TODO: add these conditionally?
            sample = FilterSample(
                name=name,
                dependencies=[],
                metadata=args.metadata,
                metadata_id_columns=args.metadata_id_columns,
                sequences=args.sequences,
                max_sequences=max_sequences,
                output_strains=path.join(tmpdir, f"{name}.samples.txt"),
                output_sequences=path.join(tmpdir, f"{name}.samples.fasta"),
                random_seed=args.random_seed,
                exclude=exclude,
                **sample_config,
            )

            # TODO: augur filter optimizations, such as not including sequences if unused

            samples[sample.name] = sample

    include = [sample.output_strains for sample in samples.values()]
    if args.include:
        include.extend(args.include)

    output_sample = FilterSample(
        name='output',
        dependencies=[sample for sample in samples],
        metadata=args.metadata,
        metadata_id_columns=args.metadata_id_columns,
        sequences=args.sequences,
        exclude_all=True,
        include=include,
        output_metadata=args.output_metadata,
        output_sequences=args.output_sequences,
    )
    samples[output_sample.name] = output_sample

    # TODO XXX check acyclic
        
    # TODO XXX assert that each dependency of a call is defined as it's own call
        
    # TODO XXX prune any calls which are not themselves used in 'output' or as a dependency of another call

    return samples


def get_runnable_call(calls):
    """
    Return a call (i.e. a filter / proximity command) which can be run, either
    because it has no dependencies or because all it's dependencies have been
    computed
    """
    for name, call in calls.items():
        if call.status == COMPLETE:
            continue
        if len(call.dependencies) == 0:
            return name
        if all([calls[name].status == COMPLETE for name in call.dependencies]):
            return name
    return None

def loop(calls, dry_run):
    """
    Execute the required calls in an approprate order (i.e. taking into account
    necessary dependencies). There are plenty of ways to do this, such as making
    a proper graph, using  parallelisation, etc etc This is a nice simple
    solution however.
    """
    while name:=get_runnable_call(calls):
        calls[name].exec(dry_run)

def run(args):
    config = parse_config(args.config)
    with tempfile.TemporaryDirectory() as tmpdir:
        calls = generate_calls(config, args, tmpdir)
        loop(calls, args.dry_run)
