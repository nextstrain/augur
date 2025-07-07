"""Subsample sequences based on a YAML configuration."""
import argparse
import yaml
from concurrent.futures import ThreadPoolExecutor
from augur import filter as augur_filter
from augur.validate import load_json_schema, validate_json


SAMPLE_OPTIONS_MAP = {
    "group_by": "--group-by",
    "group_by_weights": "--group-by-weights",
    "max_sequences": "--subsample-max-sequences",
    "sequences_per_group": "--sequences-per-group",
    "min_date": "--min-date",
    "max_date": "--max-date",
    "exclude": "--exclude",
    "exclude_where": "--exclude-where",
}

GLOBAL_OPTIONS_MAP = {
    "include": "--include",
}


# FIXME: should --include be different from --metadata-id-columns, etc.?

# FIXME: decide whether to support --priorities for v0

def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("subsample", help=__doc__)
    parser.add_argument("--metadata", required=True, help="sequence metadata")
    parser.add_argument("--sequences", required=True, help="sequence file")
    parser.add_argument("--config", required=True, help="YAML subsampling config")
    parser.add_argument("--output-metadata", required=True, help="output metadata file")
    parser.add_argument("--output-sequences", required=True, help="output sequences file")
    parser.add_argument("--nthreads", type=int, default=1, help="number of threads for parallel processing")
    # FIXME: add --subsample-seed, --metadata-id-columns, etc.
    return parser


class Sample:
    def __init__(self, name, metadata_file, sample_options, global_filter_args):
        self.name = name
        # FIXME: make this a temp file
        self.output_file = f"sample-{name}.txt"
        self.filter_args = self._build_filter_args(metadata_file, sample_options, global_filter_args)
    
    def _build_filter_args(self, metadata_file, sample_options, global_filter_args):
        filter_args = [
            "--metadata", metadata_file,
            "--output-strains", self.output_file
        ]
        
        for config_key, filter_option in SAMPLE_OPTIONS_MAP.items():
            value = sample_options.get(config_key)
            if value is None:
                continue
            if isinstance(value, list):
                filter_args.extend([filter_option, *value])
            else:
                filter_args.extend([filter_option, str(value)])

        filter_args.extend(global_filter_args)

        return filter_args
    
    def run(self):
        _run_filter(self.filter_args)


def _run_filter(args_list):
    parser = argparse.ArgumentParser()
    augur_filter.register_arguments(parser)
    args = parser.parse_args(args_list)
    augur_filter.run(args)


def parse_config(filename):
    with open(filename) as fh:
        try:
            config = yaml.safe_load(fh)
        except yaml.YAMLError as e:
            print(e)
            raise AugurError(f"Error parsing subsampling scheme {filename}")
    # FIXME: validate schema but without most of safe_load's auto conversions.
    # Some are still useful, e.g. lists used further down
    # schema = load_json_schema("schema-subsample-config.json")
    # validate_json(config, schema, args.config)
    if 'samples' not in config:
        raise AugurError('Config must define a "samples" key')
    return config


def run(args):
    config = parse_config(args.config)

    global_filter_args = []
    for config_key, filter_option in GLOBAL_OPTIONS_MAP.items():
        value = config.get(config_key)
        if value is None:
            continue
        if isinstance(value, list):
            global_filter_args.extend([filter_option, *value])

    sample_files = []
    samples = []

    for sample_name, options in config.get("samples", {}).items():
        sample = Sample(sample_name, args.metadata, options, global_filter_args)
        samples.append(sample)
        sample_files.append(sample.output_file)

    with ThreadPoolExecutor(max_workers=args.nthreads) as executor:
        executor.map(lambda sample: sample.run(), samples)

    combine_args = [
        "--metadata", args.metadata,
        "--sequences", args.sequences,
        "--exclude-all",
        "--include", *sample_files,
        "--output-metadata", args.output_metadata,
        "--output-sequences", args.output_sequences,
    ]
    _run_filter(combine_args)
