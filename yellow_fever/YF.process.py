from __future__ import print_function
import os, sys
sys.path.append('..') # we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
from base.process import process

config = {
    "dir": "yellow_fever",
    "in": "prepared/yellow-fever.json", # should be able to specify from command line
    "geo_inference": ['country'], # what traits to perform this on
    "auspice": { ## settings for auspice JSON export
        "color_options": {
            "country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete"},
        },
        "controls": {'geographic location':['country']}
    }
}

if __name__=="__main__":
    runner = process(config)
    runner.align()
    runner.build_tree()
    runner.timetree_setup_filter_run()
    runner.run_geo_inference()
    runner.save_as_nexus()
    runner.auspice_export()
