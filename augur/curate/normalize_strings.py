"""
Normalize strings to a Unicode normalization form and strip leading and trailing whitespaces.

Strings need to be normalized for predictable string comparisons, especially
in cases where strings contain diacritics (see https://unicode.org/faq/normalization.html).
"""
import unicodedata

from augur.utils import first_line


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("normalize-strings",
        parents=[parent_subparsers.shared_parser],
        help=first_line(__doc__))

    optional = parser.add_argument_group(title="OPTIONAL")
    optional.add_argument("--form", default="NFC", choices=["NFC", "NFKC", "NFD", "NFKD"],
        help="Unicode normalization form to use for normalization.")
    return parser


def normalize_strings(record, form='NFC'):
    """
    Normalizes string values in *record* to a Unicode normalization *form*
    and removes leading and trailing whitespaces from string.
    Uses `NFC` normalization form by default.

    Parameters
    ----------
    records: dict
        An input record to be normalized
    form: str, optional
        An optional Unicode normalization form

    Yields
    ------
    record: dict
        The modified record that is a shallow copy of the original record
    """
    return {
        key: (unicodedata.normalize(form, value).strip() if isinstance(value, str) else value)
        for key, value in record.items()
    }


def run(args, records):
    for record in records:
        yield normalize_strings(record, args.form)
