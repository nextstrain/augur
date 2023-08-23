"""
Stub function and module used as a setuptools entry point.
"""

import sys
import augur
from sys import argv, exit

# Entry point for setuptools-installed script and bin/augur dev wrapper.
def main():
    # Support non-Unicode encodings by replacing Unicode characters instead of erroring.
    sys.stdout.reconfigure(errors="backslashreplace")
    sys.stderr.reconfigure(errors="backslashreplace")

    return augur.run( argv[1:] )

# Run when called as `python -m augur`, here for good measure.
if __name__ == "__main__":
    exit( main() )
