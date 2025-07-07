"""Subsample sequences based on a YAML configuration."""
import argparse
import yaml
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


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("subsample", help=__doc__)
    parser.add_argument("--metadata", required=True, help="sequence metadata")
    parser.add_argument("--sequences", required=True, help="sequence file")
    parser.add_argument("--config", required=True, help="YAML subsampling config")
    parser.add_argument("--output-metadata", required=True, help="output metadata file")
    parser.add_argument("--output-sequences", required=True, help="output sequences file")
    # FIXME: add --nthreads, --subsample-seed, --metadata-id-columns, etc.
    return parser


def _run_filter(args_list):
    parser = argparse.ArgumentParser()
    augur_filter.register_arguments(parser)
    args = parser.parse_args(args_list)
    augur_filter.run(args)


def run(args):
    with open(args.config) as fh:
        config = yaml.safe_load(fh)

    # FIXME: validate schema but without most of safe_load's auto conversions.
    # Some are still useful, e.g. lists used further down
    # schema = load_json_schema("schema-subsample-config.json")
    # validate_json(config, schema, args.config)

    sample_files = []

    for sample_name, options in config.get("samples", {}).items():
        # FIXME: make this a temp file
        sample_file = f"sample-{sample_name}.txt"
        sample_files.append(sample_file)

        sample_args = ["--metadata", args.metadata, "--output-strains", sample_file]

        for config_key, filter_option in SAMPLE_OPTIONS_MAP.items():
            value = options.get(config_key)
            if value is None:
                continue
            if isinstance(value, list):
                sample_args.extend([filter_option, *value])
            else:
                sample_args.extend([filter_option, str(value)])

        for config_key, filter_option in GLOBAL_OPTIONS_MAP.items():
            value = config.get(config_key)
            if value is None:
                continue
            if isinstance(value, list):
                sample_args.extend([filter_option, *value])

        # FIXME: parallelize with --nthreads
        _run_filter(sample_args)

    combine_args = [
        "--metadata",
        args.metadata,
        "--sequences",
        args.sequences,
        "--exclude-all",
        "--include",
        *sample_files,
        "--output-metadata",
        args.output_metadata,
        "--output-sequences",
        args.output_sequences,
    ]
    _run_filter(combine_args)
