"""
Subsample a dataset according to a provided config.

This augur functionality is _alpha_ and may change at any time.
It does not conform to the semver release standards used by Augur.
"""

from .errors import AugurError
from os import path, mkdir
import subprocess

INCOMPLETE = 'incomplete'
COMPLETE = 'complete'

def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("subsample", help=__doc__)
    parser.add_argument('--config', required=True, metavar="YAML", help="subsampling config file") # TODO: allow a string representation
    parser.add_argument('--metadata', required=True, metavar="TSV", help="sequence metadata")
    parser.add_argument('--sequences', required=True, metavar="FASTA", help="sequences in FASTA format") # TODO XXX VCF ?
    parser.add_argument('--output-metadata', required=True, metavar="TSV", help="output metadata")
    parser.add_argument('--output-sequences', required=True, metavar="FASTA", help="output sequences in FASTA format") # TODO XXX VCF ?
    ## TODO XXX - programmatically create & remove tmp directory. But this is simpler for dev.
    parser.add_argument('--tmpdir', required=True, help="temporary directory for intermediate files")

    optionals = parser.add_argument_group(
        title="Optional arguments",
    )
    optionals.add_argument('--dry-run', action="store_true")
    optionals.add_argument('--metadata-id-columns', metavar="NAME", nargs="+",
        help="names of possible metadata columns containing identifier information, ordered by priority. Only one ID column will be inferred.")
    optionals.add_argument('--reference', metavar="FASTA", help="needed for priority calculation (but it shouldn't be!)")

    return parser

def parse_config(filename):
    import yaml
    with open(filename) as fh:
        try:
            data = yaml.safe_load(fh)
        except yaml.YAMLError as e:
            print(e)
            raise AugurError(f"Error parsing subsampling scheme {filename}")
    # TODO XXX - write a schema and validate against this
    return data

class Filter():
    def __init__(self, name, params, data_in, data_out, optional_args, depends_on, log_fname):
        self.name = name
        self.params = params
        self.data_in = data_in
        self.data_out = data_out
        self.depends_on = depends_on
        self.optional_args = optional_args
        self.log_fname = log_fname
        self.status = INCOMPLETE

    def __str__(self):
        deps = (" depends on " + ", ".join(self.depends_on)) if len(self.depends_on) else " (no dependencies)"
        return f"Augur filter for intermediate sample {self.name!r}" + deps

    def cmd(self):
        # Arg splitting would be better, but is complicated by quoting. TODO.
        cmd = f"augur filter {self.params} --metadata {self.data_in['metadata']} --sequences {self.data_in['sequences']}"
        if "metadata" in self.data_out:
            cmd += f" --output-metadata {self.data_out['metadata']}"
        if "sequences" in self.data_out:
            cmd += f" --output-sequences {self.data_out['sequences']}"
        if "strains" in self.data_out:
            cmd += f" --output-strains {self.data_out['strains']}"
        if 'metadata_id_columns' in self.optional_args:
            cmd += f" --metadata-id-columns {' '.join(self.optional_args.metadata_id_columns)}"
        return cmd

    def exec(self, dry_run=False):
        """
        Instead of running an `augur filter` command in a subprocess like we do
        here, a nicer way would be to refactor `augur filter` to expose
        functions which can be called here so that we can get a list of returned
        strains. Doing so would provide a HUGE speedup if we could load the
        sequences+metadata into memory a single time and then use that in-memory
        data for all filtering calls which subsampling performs. This is
        analogous to having an in-memory database running in a separate process
        we can query - not as fast, but much easier implementation. Using
        subprocess does make parallelisation trivial, but the above speed up
        would be preferable.
        """
        print("\n" + self.__str__())
        cmd = self.cmd()
        print("RUNNING " + cmd)
        if not dry_run:
            try:
                with open(self.log_fname, 'w') as fh:
                    subprocess.run(cmd, shell=True, check=True, text=True, stdout=fh, stderr=fh)
            except subprocess.CalledProcessError as e:
                print(e)
                raise AugurError("Check the logfile: "+self.log_fname)

        self.status = COMPLETE



class Proximity():
    def __init__(self, name, focal_data, contextual_data, num_per_focal, output, depends_on, args, logfile):
        self.name = name
        self.focal_data = focal_data
        self.contextual_data = contextual_data
        self.num_per_focal = num_per_focal
        self.output = output
        self.depends_on = depends_on
        self.args = args
        self.logfile = logfile
        self.status = INCOMPLETE

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
            self.contextual_data['sequences'],
            self.args.reference,
            self.focal_data['sequences'],
            self.num_per_focal
        )

    def __str__(self):
        return f"Calculate proximity for {self.name} by computing weights for {self.focal_data['sequences']} against {self.contextual_data['sequences']}"

    def exec(self, dry_run=False):
        print("\n" + self.__str__())
        if not dry_run:
            closest = self.get_closest_sequences()
            with open(self.output, 'w') as fh:
                print("\n".join(list(closest)), file=fh)
            print(f"\tClosest strains written to {self.output} (n={len(closest)})")
        self.status = COMPLETE


def generate_calls(config, args):
    """
    Produce an (unordered) dictionary of calls to be made to accomplish the
    desired subsampling. Each call is either (i) a use of augur filter or (ii) a
    proximity calculation The names given to calls are config-defined, but there
    is guaranteed to be one call with the name "output".

    The separation between this function and the Filter (etc) classes is not
    quite right, but it's a WIP.
    """
    calls = {}

    if 'output' not in config:
        raise AugurError('Config must define an "output" key')

    tmpdir = args.tmpdir
    if not path.isdir(args.tmpdir):
        if path.exists(args.tmpdir):
            raise AugurError(f'Provided --tmpdir {args.tmpdir!r} must be a directory')
        mkdir(args.tmpdir)

    for sample_name, sample_config in config.items():
        if sample_name == 'output':  # output is special cased
            depends_on = sample_config
            calls['output'] = Filter('output',
                '--exclude-all --include ' + ' '.join([path.join(tmpdir, f"{name}.samples.txt") for name in sample_config]),
                {'metadata': args.metadata, 'sequences': args.sequences},
                {'metadata': args.output_metadata, 'sequences': args.output_sequences},
                args,
                depends_on,
                path.join(tmpdir, "output.log.txt")
            )
            continue

        ## TODO XXX
        ## I designed this to have a 'include' parameter whereby the starting meta/seqs for this filter call could
        ## be the (joined) output of previous samples. To be implemented.
        ## Similarly, the "exclude" YAML parameter is also not implemented

        if 'include' in sample_config:
            raise AugurError("'include' subsampling functionality not yet implemented")

        depends_on = []

        if 'priorities' in sample_config:
            priority_type = sample_config['priorities'].get('type', '')
            focus_name = sample_config['priorities']['focus']
            if priority_type != 'proximity':
                raise AugurError(f"Priorities must be proximity, not {priority_type!r}")
            if sample_config.get('filter', '') != '':
                raise AugurError(f"Priorities must not be used in conjunction with filtering queries.")
            priority_fname = path.join(tmpdir, f"{focus_name}.priorities.txt")
            calls["__priorities__"+focus_name] = Proximity(focus_name,
                {'sequences': path.join(tmpdir, f"{focus_name}.fasta")},    # input (focal seqs)
                {'metadata': args.metadata, 'sequences': args.sequences},   # input (contextual seqs)
                sample_config['priorities']['num_per_focal'],
                priority_fname,                                             # output         
                [focus_name],                                               # depends on
                args,
                path.join(tmpdir, f"{sample_name}.priorities.log.txt")      # logfile
            )
            sample_config['filter'] = f"--exclude-all --include {priority_fname}"
            depends_on.append("__priorities__"+focus_name)

        calls[sample_name] = Filter(sample_name,
            sample_config['filter'],
            {'metadata': args.metadata, 'sequences': args.sequences},
            {'strains': path.join(tmpdir, sample_name+'.samples.txt')},
            args,
            depends_on,
            path.join(tmpdir, f"{sample_name}.log.txt")
        )


    # Any intermediate sample a priority calculation depends on requires FASTA input, not samples.txt
    # (This is required by distance-to-focal-set, but we could change this)
    for priority_call in [c for c in calls.values() if isinstance(c, Proximity)]:
        for dep in priority_call.depends_on:
            print("Priority calc", priority_call.name, "depends on ", dep)
            calls[dep].data_out['sequences'] = calls[dep].data_out['strains'].replace('.samples.txt', ".fasta")

    # TODO XXX check acyclic
        
    # TODO XXX assert that each dependency of a call is defined as it's own call
        
    # TODO XXX prune any calls which are not themselves used in 'output' or as a dependency of another call

    return calls


def get_runnable_call(calls):
    """
    Return a call (i.e. a filter / proximity command) which can be run, either
    because it has no dependencies or because all it's dependencies have been
    computed
    """
    for name, call in calls.items():
        if call.status==COMPLETE:
            continue
        if len(call.depends_on)==0:
            return name
        if all([calls[name].status==COMPLETE for name in call.depends_on]):
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
    calls = generate_calls(config, args)
    loop(calls, args.dry_run)
