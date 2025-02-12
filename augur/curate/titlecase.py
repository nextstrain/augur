"""
Applies titlecase to specified string fields
"""
import re
from typing import Optional, Set, Union

from augur.argparse_ import ExtendOverwriteDefault
from augur.errors import AugurError
from augur.io.print import print_err
from augur.types import DataErrorMethod

def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("titlecase",
    parents = [parent_subparsers.shared_parser],
    help = __doc__)

    required = parser.add_argument_group(title="REQUIRED")
    required.add_argument("--titlecase-fields", nargs="*", action=ExtendOverwriteDefault,
        help="List of fields to convert to titlecase.", required=True)

    optional = parser.add_argument_group(title="OPTIONAL")
    optional.add_argument("--articles", nargs="*", action=ExtendOverwriteDefault,
        help="List of articles that should not be converted to titlecase.")
    optional.add_argument("--abbreviations", nargs="*", action=ExtendOverwriteDefault,
        help="List of abbreviations that should not be converted to titlecase, keeps uppercase.")

    optional.add_argument("--failure-reporting",
        type=DataErrorMethod.argtype,
        choices=[ method for method in DataErrorMethod ],
        default=DataErrorMethod.ERROR_FIRST,
        help="How should failed titlecase formatting be reported.")
    return parser


def titlecase(text: Union[str, None], articles: Set[str] = set(), abbreviations: Set[str] = set()) -> Optional[str]:
    """
    Originally from nextstrain/ncov-ingest

    Returns a title cased location name from the given location name
    *tokens*. Ensures that no tokens contained in the *whitelist_tokens* are
    converted to title case.

    >>> articles = {'a', 'and', 'of', 'the', 'le'}
    >>> abbreviations = {'USA', 'DC'}

    >>> titlecase("the night OF THE LIVING DEAD", articles)
    'The Night of the Living Dead'

    >>> titlecase("BRAINE-LE-COMTE, FRANCE", articles)
    'Braine-le-Comte, France'

    >>> titlecase("auvergne-RHÔNE-alpes", articles)
    'Auvergne-Rhône-Alpes'

    >>> titlecase("washington DC, usa", articles, abbreviations)
    'Washington DC, USA'
    """
    if not isinstance(text, str):
        return None

    words = enumerate(re.split(r'\b', text))

    def changecase(index, word):
        casefold = word.casefold()
        upper = word.upper()

        if upper in abbreviations:
            return upper
        elif casefold in articles and index != 1:
            return word.lower()
        else:
            return word.title()

    return ''.join(changecase(i, w) for i, w in words)


def run(args, records):
    failures = []
    failure_reporting = args.failure_reporting

    articles = set()
    if args.articles:
        articles = set(args.articles)

    abbreviations = set()
    if args.abbreviations:
        abbreviations = set(args.abbreviations)

    for index, record in enumerate(records):
        record = record.copy()
        record_id = index

        for field in args.titlecase_fields:
            # Ignore non-existent fields but could change this to warnings if desired
            if field not in record:
                continue
            elif record[field] is None:
                continue

            titlecased_string = titlecase(record[field], articles, abbreviations)

            failure_message = f"Failed to titlecase {field!r}:{record.get(field)!r} in record {record_id!r} because the value is a {type(record.get(field)).__name__!r} and is not a string."
            if titlecased_string is None:
                if failure_reporting is DataErrorMethod.ERROR_FIRST:
                    raise AugurError(failure_message)
                if failure_reporting is DataErrorMethod.ERROR_ALL:
                    print_err(f"ERROR: {failure_message}")
                if failure_reporting is DataErrorMethod.WARN:
                    print_err(f"WARNING: {failure_message}")

                # Keep track of failures for final summary
                failures.append((record_id, field, record.get(field)))
            else:
                record[field] = titlecased_string

        yield record

    if failure_reporting is not DataErrorMethod.SILENT and failures:
        failure_message = (
            "Unable to change to titlecase for the following (record, field, field value):\n" + \
            '\n'.join(map(repr, failures))
        )
        if failure_reporting is DataErrorMethod.ERROR_ALL:
            raise AugurError(failure_message)

        elif failure_reporting is DataErrorMethod.WARN:
            print_err(f"WARNING: {failure_message}")

        else:
            raise ValueError(f"Encountered unhandled failure reporting method: {failure_reporting!r}")
