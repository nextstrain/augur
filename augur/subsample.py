"""
Subsample a dataset according to a provided config.

This augur functionality is _alpha_ and may change at any time.
It does not conform to the semver release standards used by Augur.
"""

from .errors import AugurError

INCOMPLETE = 'incomplete'
COMPLETE = 'complete'

def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("subsample", help=__doc__)
    parser.add_argument('--config', required=True, metavar="YAML", help="subsampling config file") # TODO: allow a string representation
    parser.add_argument('--metadata', required=True, metavar="TSV", help="sequence metadata")
    parser.add_argument('--sequences', metavar="FASTA", help="sequences in FASTA format") # TODO XXX VCF ?
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
    def __init__(self, name, filter, depends_on):
        self.name = name
        self.params = filter
        self.depends_on = depends_on
        self.status = INCOMPLETE

    def __str__(self):
        if len(self.depends_on):
            return self.name + " depends on " + ", ".join(self.depends_on)
        return self.name + " (no dependencies)"

    def exec(self):
        print("RUN", self.__str__())
        self.status = COMPLETE



class Proximity():
    def __init__(self, name, focal_samples, source_seqs, depends_on):
        self.name = name
        self.focal_samples = focal_samples
        self.source_seqs = source_seqs
        self.depends_on = depends_on
        self.status = INCOMPLETE

    def __str__(self):
        return f"calculate proximity for {self.name} by computing weights for {self.source_seqs} compared with {self.focal_samples}"

    def exec(self):
        print("RUN", self.__str__())
        self.status = COMPLETE




def generate_calls(config, seq_fname, meta_fname):
    """
    Produce an (unordered) dictionary of calls to be made to accomplish the desired subsampling.
    Each call is either (i) a use of augur filter or (ii) a proximity calculation
    The names given to calls are config-defined, but there is guaranteed to be one call
    with the name "output"
    """
    calls = {}

    if 'output' not in config:
        raise AugurError('Config must define an "output" key')

    for sample_name, sample_config in config.items():
        if sample_name == 'output':
            # output is special cased
            include = ' '.join([f"{name}.samples.txt" for name in sample_config])
            calls['output'] = Filter('output',
                f"--exclude-all --include {include} --metadata {meta_fname} --sequences {seq_fname}",
                sample_config
            )
            continue

        ## TODO XXX
        ## I designed this to have a 'include' parameter whereby the starting meta/seqs for this filter call could
        ## be the (joined) output of previous samples. To be implemented.

        if 'include' in sample_config:
            raise AugurError("'include' subsampling functionality not yet implemented")

        depends_on = []

        if 'priorities' in sample_config:
            priority_type = sample_config['priorities'].get('type', '')
            focus_name = sample_config['priorities']['focus']
            if priority_type != 'proximity':
                raise AugurError(f"Priorities must be proximity, not {priority_type!r}")
            
            calls["__priorities__"+focus_name] = Proximity(focus_name,
                f"{focus_name}.samples.txt",
                f"{seq_fname}",                                         
                [focus_name]
            )
            
            depends_on.append("__priorities__"+focus_name)

        calls[sample_name] = Filter(sample_name,
            f"{sample_config['filter']} --metadata {meta_fname} --sequences {seq_fname}",
            depends_on
        )

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

def loop(calls):
    """
    Execute the required calls in an approprate order (i.e. taking into account
    necessary dependencies). There are plenty of ways to do this, such as making
    a proper graph, using  parallelisation, etc etc This is a nice simple
    solution however.
    """
    while name:=get_runnable_call(calls):
        calls[name].exec()



def run(args):
    config = parse_config(args.config)
    calls = generate_calls(config, args.metadata, args.sequences)
    loop(calls)
