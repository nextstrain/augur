"""
Print the version of augur.
"""

from .__version__ import __version__

def register_arguments(parser):
    pass

def run(args):
    print("augur", __version__)
    return 0
