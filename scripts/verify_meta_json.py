from __future__ import print_function
import sys
sys.path.append('..') # this is an assumption and is probably wrong
import argparse
import json


def verify_defaults(defaults):
    booleans = ["mapTriplicate"];
    strings = ['geoResolution', 'colorBy', 'distanceMeasure']
    for trait in booleans:
        if (trait in defaults):
            if type(defaults[trait]) is not bool:
                print("ERROR: defaults -> {} is not boolean".format(trait))
    for trait in strings:
        if (trait in defaults):
            if type(defaults[trait]) is not str:
                print("ERROR: defaults -> {} is not a string".format(trait))

def verify_author_info(arg):
    pass

if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "Verify metadata JSONs")
    parser.add_argument('--json', required=True, type=str, help="metadata JSON")
    params = parser.parse_args()

    try:
        with open(params.json, 'r') as fh:
            data = json.load(fh)
    except:
        print("Couldn't even load the file - big problems!")


    if "defaults" in data:
        verify_defaults(data["defaults"])

    if "author_info" in data:
        verify_author_info(data["author_info"])
    else:
        print("ERROR: author_info does not exist. This build needs to be updated and will soon be incompatible with auspice.")


    # set_trace()
