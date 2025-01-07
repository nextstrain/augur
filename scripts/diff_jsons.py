"""Compare JSON files line by line with deepdiff
"""
import argparse
import deepdiff
import json

from augur.argparse_ import ExtendOverwriteDefault


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare JSON files line by line with deepdiff",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("first_json", help="first JSON to compare")
    parser.add_argument("second_json", help="second JSON to compare")
    parser.add_argument("--significant-digits", type=int, default=5, help="number of significant digits to use when comparing numeric values")
    parser.add_argument("--exclude-paths", nargs="+", action=ExtendOverwriteDefault, help="list of paths to exclude from consideration when performing a diff", default=["root['generated_by']['version']"])
    parser.add_argument("--exclude-regex-paths", nargs="+", action=ExtendOverwriteDefault, help="list of path regular expressions to exclude from consideration when performing a diff")
    parser.add_argument("--ignore-numeric-type-changes", action="store_true", help="ignore numeric type changes in the diff (e.g., int of 1 to float of 1.0)")

    args = parser.parse_args()

    with open(args.first_json, "r") as fh:
        first_json = json.load(fh)

    with open(args.second_json, "r") as fh:
        second_json = json.load(fh)

    print(
        deepdiff.DeepDiff(
            first_json,
            second_json,
            significant_digits=args.significant_digits,
            exclude_paths=args.exclude_paths,
            exclude_regex_paths=args.exclude_regex_paths,
            ignore_numeric_type_changes=args.ignore_numeric_type_changes,
        )
    )
