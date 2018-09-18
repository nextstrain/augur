"""
Prints the version of augur.
"""

from .__version__ import __version__

def run(args):
    print("augur", __version__)
    return 0
