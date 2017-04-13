from __future__ import print_function
import os, sys
sys.path.append('..') # we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
from base.process import process

config = {
    "dir": "H7N9",
    "output": { # will move to the default config file
        "data": "processed",
        "auspice": "auspice",
    },
    "in": "prepared/H7N9_HA.json", # should be able to specify from command line
}

if __name__=="__main__":
    runner = process(config)
    runner.align()
    runner.build_tree()
