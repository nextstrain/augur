import argparse
import shlex
from augur.filter import register_arguments


def parse_args(args: str):
    parser = argparse.ArgumentParser()
    register_arguments(parser)
    return parser.parse_args(shlex.split(args))


def write_metadata(tmpdir, metadata):
    fn = str(tmpdir / "metadata.tsv")
    with open(fn, "w") as fh:
        fh.write("\n".join(("\t".join(md) for md in metadata)))
    return fn
