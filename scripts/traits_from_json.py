from __future__ import print_function
import sys
sys.path.append('..') # this is an assumption and is probably wrong
from base.utils import parse_date
import argparse
import json
from numpy import ndarray

from augur.argparse_ import ExtendOverwriteDefault


def get_trait(attributes, trait, dateFormat):
    if trait in ["num_date", "date"]:
        try:
            date = parse_date(attributes["raw_date"], dateFormat)[1]
        except KeyError:
            return "unknown"
        if isinstance(date, (list, tuple, ndarray)):
            return (date[0] + date[1]) / 2
        return date

    else:
        try:
            return attributes[trait]
        except KeyError:
            return "unknown"


if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "Process a given JSONs")
    parser.add_argument('--json', required=True, type=str, help="prepared JSON")
    parser.add_argument('--trait', required=True, type=str, help="prepared JSON")
    parser.add_argument('--header', nargs='*', action=ExtendOverwriteDefault, type=str, help="header fields")
    parser.add_argument('--date_format', nargs='*', action=ExtendOverwriteDefault, default=["%Y-%m-%d"], type=str, help="if needed. default: [%%Y-%%m-%%d]")
    params = parser.parse_args()

    with open(params.json, 'r') as fh:
        data = json.load(fh)

    try:
        print("\t".join(params.header))
    except KeyError:
        pass

    for seq, val in data["sequences"].items():
        print("{}\t{}".format(seq, get_trait(val["attributes"], params.trait, params.date_format)))
