"""
Import analyses into augur pipeline from other systems
"""
from .import_beast import run_beast

def run(args):
    if "beast" in args:
        return run_beast(args)
