import sys


def print_err(*args):
    print(*args, file=sys.stderr)
