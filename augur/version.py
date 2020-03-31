"""
Print the version of augur.
"""

from .__version__ import __version__


def register_arguments(parser):
    pass


def run(args):
    print("augur", _get_version())
    return 0


def _get_version():
    return __version__
