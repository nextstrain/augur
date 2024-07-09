"""
Applies user curated geolocation rules to the geolocation fields.
"""
from collections import defaultdict
from augur.errors import AugurError
from augur.io.print import print_err
from augur.utils import first_line


class CyclicGeolocationRulesError(AugurError):
    pass


def load_geolocation_rules(geolocation_rules_file):
    """
    Loads the geolocation rules from the provided *geolocation_rules_file*.
    Returns the rules as a dict:

.. code-block:: text

    {
        regions: {
            countries: {
                divisions: {
                    locations: corrected_geolocations_tuple
                }
            }
        }
    }
    """
    geolocation_rules = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    with open(geolocation_rules_file, 'r') as rules_fh:
        for line in rules_fh:
            # ignore comments
            if line.strip()=="" or line.lstrip()[0] == '#':
                continue

            row = line.strip('\n').split('\t')
            # Skip lines that cannot be split into raw and annotated geolocations
            if len(row) != 2:
                print_err(
                    f"WARNING: Could not decode geolocation rule {line!r}.",
                    "Please make sure rules are formatted as",
                    "'region/country/division/location<tab>region/country/division/location'.")
                continue

            # remove trailing comments
            row[-1] = row[-1].partition('#')[0].rstrip()
            raw , annot = tuple( row[0].split('/') ) , tuple( row[1].split('/') )

            # Skip lines where raw or annotated geolocations cannot be split into 4 fields
            if len(raw) != 4:
                print_err(
                    f"WARNING: Could not decode the raw geolocation {row[0]!r}.",
                    "Please make sure it is formatted as 'region/country/division/location'.")
                continue

            if len(annot) != 4:
                print_err(
                    f"WARNING: Could not decode the annotated geolocation {row[1]!r}.",
                    "Please make sure it is formatted as 'region/country/division/location'.")
                continue


            geolocation_rules[raw[0]][raw[1]][raw[2]][raw[3]] = annot

    return geolocation_rules


def get_annotated_geolocation(geolocation_rules, raw_geolocation, rule_traversal = None):
    """
    Gets the annotated geolocation for the *raw_geolocation* in the provided
    *geolocation_rules*.

    Recursively traverses the *geolocation_rules* until we get the annotated
    geolocation, which must be a Tuple. Returns `None` if there are no
    applicable rules for the provided *raw_geolocation*.

    Rules are applied in the order of region, country, division, location.
    First checks the provided raw values for geolocation fields, then if there
    are not matches, tries to use general rules marked with '*'.
    """
    # Always instantiate the rule traversal as an empty list if not provided,
    # e.g. the first call of this recursive function
    if rule_traversal is None:
        rule_traversal = []

    current_rules = geolocation_rules
    # Traverse the geolocation rules based using the rule_traversal values
    for field_value in rule_traversal:
        current_rules = current_rules.get(field_value)
        # If we hit `None`, then we know there are no matching rules, so stop the rule traversal
        if current_rules is None:
            break

    # We've found the tuple of the annotated geolocation
    if isinstance(current_rules, tuple):
        return current_rules

    # We've reach the next level of geolocation rules,
    # so try to traverse the rules with the next target in raw_geolocation
    if isinstance(current_rules, dict):
        next_traversal_target = raw_geolocation[len(rule_traversal)]
        rule_traversal.append(next_traversal_target)
        return get_annotated_geolocation(geolocation_rules, raw_geolocation, rule_traversal)

    # We did not find any matching rule for the last traversal target
    if current_rules is None:
        # If we've used all general rules and we still haven't found a match,
        # then there are no applicable rules for this geolocation
        if all(value == '*' for value in rule_traversal):
            return None

        # If we failed to find matching rule with a general rule as the last
        # traversal target, then delete all trailing '*'s to reset rule_traversal
        # to end with the last index that is currently NOT a '*'
        # [A, *, B, *] => [A, *, B]
        # [A, B, *, *] => [A, B]
        # [A, *, *, *] => [A]
        if rule_traversal[-1] == '*':
            # Find the index of the first of the consecutive '*' from the
            # end of the rule_traversal
            # [A, *, B, *] => first_consecutive_general_rule_index = 3
            # [A, B, *, *] => first_consecutive_general_rule_index = 2
            # [A, *, *, *] => first_consecutive_general_rule_index = 1
            for index, field_value in reversed(list(enumerate(rule_traversal))):
                if field_value == '*':
                    first_consecutive_general_rule_index = index
                else:
                    break

            rule_traversal = rule_traversal[:first_consecutive_general_rule_index]

        # Set the final value to '*' in hopes that by moving to a general rule,
        # we can find a matching rule.
        rule_traversal[-1] = '*'

        return get_annotated_geolocation(geolocation_rules, raw_geolocation, rule_traversal)


def transform_geolocations(geolocation_rules, geolocation):
    """
    Transform the provided *geolocation* by looking it up in the provided
    *geolocation_rules*.

    This will use all rules that apply to the geolocation and rules will
    be applied in the order of region, country, division, location.

    Returns the original geolocation if no geolocation rules apply.

    Raises a `CyclicGeolocationRulesError` if more than 1000 rules have
    been applied to the raw geolocation.
    """
    transformed_values = geolocation
    rules_applied = 0
    continue_to_apply = True

    while continue_to_apply:
        annotated_values = get_annotated_geolocation(geolocation_rules, transformed_values)

        # Stop applying rules if no annotated values were found
        if annotated_values is None:
            continue_to_apply = False
        else:
            rules_applied += 1

            if rules_applied > 1000:
                raise CyclicGeolocationRulesError(
                    f"More than 1000 geolocation rules applied on the same entry {geolocation!r}."
                )

            # Create a new list of values for comparison to previous values
            new_values = list(transformed_values)
            for index, value in enumerate(annotated_values):
                # Keep original value if annotated value is '*'
                if value != '*':
                    new_values[index] = value

            # Stop applying rules if this rule did not change the values,
            # since this means we've reach rules with '*' that no longer change values
            if new_values == transformed_values:
                continue_to_apply = False

            transformed_values = new_values

    return transformed_values


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("apply-geolocation-rules",
        parents=[parent_subparsers.shared_parser],
        help=first_line(__doc__))

    parser.add_argument("--region-field", default="region",
        help="Field that contains regions in NDJSON records.")
    parser.add_argument("--country-field", default="country",
        help="Field that contains countries in NDJSON records.")
    parser.add_argument("--division-field", default="division",
        help="Field that contains divisions in NDJSON records.")
    parser.add_argument("--location-field", default="location",
        help="Field that contains location in NDJSON records.")
    parser.add_argument("--geolocation-rules", metavar="TSV", required=True,
        help="TSV file of geolocation rules with the format: " +
             "'<raw_geolocation><tab><annotated_geolocation>' where the raw and annotated geolocations " +
             "are formatted as '<region>/<country>/<division>/<location>'. " +
             "If creating a general rule, then the raw field value can be substituted with '*'." +
             "Lines starting with '#' will be ignored as comments." +
             "Trailing '#' will be ignored as comments.")

    return parser


def run(args, records):
    location_fields = [args.region_field, args.country_field, args.division_field, args.location_field]

    geolocation_rules = load_geolocation_rules(args.geolocation_rules)

    for record in records:

        annotated_values = transform_geolocations(geolocation_rules, [record.get(field, '') for field in location_fields])

        for index, field in enumerate(location_fields):
            record[field] = annotated_values[index]

        yield record
