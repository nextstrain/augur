"""
Renames fields / columns of the input data
"""

from typing import Iterable, Literal, Union, List, Tuple
import argparse
from augur.io.print import print_err
from augur.errors import AugurError

def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("rename",
    parents = [parent_subparsers.shared_parser],
    help = __doc__)

    required = parser.add_argument_group(title="REQUIRED")
    required.add_argument("--field-map", nargs="+", required=True,
        help="Rename fields/columns via '{old_field_name}={new_field_name}'. " +
             "If the new field already exists, then the renaming of the old field will be skipped. " +
             "Multiple entries with the same '{old_field_name}' will duplicate the field/column. " +
             "Skips the field if the old field name is the same as the new field name (case-sensitive).")

    optional = parser.add_argument_group(title="OPTIONAL")
    optional.add_argument("--force", action="store_true",
        help="Force renaming of old field even if the new field already exists. " +
             "Please keep in mind this will overwrite the value of the new field.")

    return parser


def parse_field_map(field_map_arg: List[str]) -> List[Tuple[str,str]]:
    seen_new = set() # keep track of the new field names
    
    field_map = []
    for field in field_map_arg:
        fields = [n.strip() for n in field.split('=')]
        if len(fields)!=2:
            raise AugurError(f"The field-map {field!r} must contain a single '=' character.")
        old_name, new_name = fields

        # Sanity check the requests to catch typos etc
        if not old_name:
            raise AugurError(f"The field-map {field!r} doesn't specify a name for the existing field.")
        if not new_name:
            raise AugurError(f"The field-map {field!r} doesn't specify a name for the new field.")
        if new_name in seen_new:
            raise AugurError(f"Asked to rename multiple fields to {new_name!r}.")
        seen_new.add(new_name)

        field_map.append((old_name, new_name))
    return field_map


def transform_columns(existing_fields: List[str], field_map: List[Tuple[str,str]], force: bool) -> List[Tuple[str,str]]:
    """
    Calculate the mapping of old column names to new column names
    """
    # check that all columns to be renamed exist
    for idx,names in enumerate(field_map[:]):
        old_name, new_name = names
        if old_name not in existing_fields:
            print_err(f"WARNING: Asked to rename field {old_name!r} (to {new_name!r}) but it doesn't exist in the input data.")
            field_map.pop(idx)

    # iterate through field_map and remove rename requests if they would drop an existing column
    # doing this ahead-of-time allows us to preserve the order of fields using a simpler implementation
    if not force:
        for idx, fields in enumerate(field_map[:]):
            old_field, new_field = fields
            if new_field in existing_fields and new_field!=old_field:
                print_err(
                    f"WARNING: skipping rename of {old_field} because record",
                    f"already has a field named {new_field}."
                )
                field_map.pop(idx)

    names_to_change, new_names = set([f[0] for f in field_map]), set([f[1] for f in field_map])

    m = []
    for field in existing_fields:
        if field in names_to_change:
            m += [(field,new_field) for old_field, new_field in field_map if old_field==field]
        elif field in new_names:
            pass # another column is renamed to this name, so we drop it
        else:
            m.append((field, field)) # no change to field name
    return m


def run(args: argparse.Namespace, records: Iterable[dict]) -> Iterable[dict]:
    col_map: Union[Literal[False], List[Tuple[str,str]]] = False
    for record in records:
        if not col_map: # initialise using first record
            col_map = transform_columns(list(record.keys()), parse_field_map(args.field_map), args.force)
        yield({new_field:record[old_field] for old_field, new_field in col_map})
