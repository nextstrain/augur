"""Compare JSON files line by line with deepdiff
"""
import argparse
import deepdiff
import json


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare JSON files line by line with deepdiff",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("first_json", help="first JSON to compare")
    parser.add_argument("second_json", help="second JSON to compare")
    parser.add_argument("--significant-digits", type=int, default=5, help="number of significant digits to use when comparing numeric values")

    args = parser.parse_args()

    with open(args.first_json, "r") as fh:
        first_json = json.load(fh)

    with open(args.second_json, "r") as fh:
        second_json = json.load(fh)

    print(
        deepdiff.DeepDiff(
            first_json,
            second_json,
            significant_digits=args.significant_digits
        )
    )
