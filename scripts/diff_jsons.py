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
    parser.add_argument("--exclude-paths", nargs="+", help="list of paths to exclude from consideration when performing a diff", default=["root['generated_by']['version']"])
    parser.add_argument("--exclude-regex-paths", nargs="+", help="list of path regular expressions to exclude from consideration when performing a diff")
    parser.add_argument("--pretty", action="store_true", help="use DeepDiff's pretty method when printing differences")

    args = parser.parse_args()

    with open(args.first_json, "r") as fh:
        first_json = json.load(fh)

    with open(args.second_json, "r") as fh:
        second_json = json.load(fh)

    difference = deepdiff.DeepDiff(
        first_json,
        second_json,
        significant_digits=args.significant_digits,
        exclude_paths=args.exclude_paths,
        exclude_regex_paths=args.exclude_regex_paths,
    )

    if args.pretty:
       difference = difference.pretty()

    print(difference)
