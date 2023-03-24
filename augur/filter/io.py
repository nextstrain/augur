import os
import numpy as np
from collections import defaultdict

from augur.errors import AugurError


def read_priority_scores(fname):
    def constant_factory(value):
        return lambda: value

    try:
        with open(fname, encoding='utf-8') as pfile:
            return defaultdict(constant_factory(-np.inf), {
                elems[0]: float(elems[1])
                for elems in (line.strip().split('\t') if '\t' in line else line.strip().split() for line in pfile.readlines())
            })
    except Exception:
        raise AugurError(f"missing or malformed priority scores file {fname}")


def cleanup_outputs(args):
    """Remove output files. Useful when terminating midway through a loop of metadata chunks."""
    if args.output:
        _try_remove(args.output)
    if args.output_metadata:
        _try_remove(args.output_metadata)
    if args.output_strains:
        _try_remove(args.output_strains)
    if args.output_log:
        _try_remove(args.output_log)


def _try_remove(filepath):
    """Remove a file if it exists."""
    try:
        os.remove(filepath)
    except FileNotFoundError:
        pass
