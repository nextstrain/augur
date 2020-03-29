import os
import sys

recursion_limit = os.environ.get("AUGUR_RECURSION_LIMIT")
if recursion_limit:
    sys.setrecursionlimit(int(recursion_limit))
