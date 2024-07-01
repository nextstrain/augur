#!/usr/bin/env python3
"""
Parses GenBank's 'location' field of the NDJSON record from stdin to 3 separate
fields: 'country', 'division', and 'location'. Checks that a record is from
GenBank by verifying that the 'database' field has a value of "GenBank" or "RefSeq".

Outputs the modified record to stdout.
"""
import json
from sys import stdin, stderr, stdout


def parse_location(record: dict) -> dict:
    # Expected pattern for the location field is "<country_value>[:<region>][, <locality>]"
    # See GenBank docs for their "country" field:
    # https://www.ncbi.nlm.nih.gov/genbank/collab/country/
    location_field = record.get("location", "")
    if not location_field:
        print(
            "`transform-genbank-location` requires a `location` field; this record does not have one.",
            file=stderr,
        )
        # bail early because we're not gonna make any changes
        return record

    geographic_data = location_field.split(':')

    country = geographic_data[0]
    division = ''
    location = ''

    if len(geographic_data) == 2:
        division , _ , location = geographic_data[1].partition(',')

    record['country'] = country.strip()
    record['division'] = division.strip()
    record['location'] = location.strip()

    return record


if __name__ == '__main__':

    for record in stdin:
        record = json.loads(record)

        database = record.get('database', '')
        if database in {'GenBank', 'RefSeq'}:
            parse_location(record)
        else:
            if database:
                error_msg = f"""Database value of {database} not supported for `transform-genbank-location`; must be "GenBank" or "RefSeq"."""
            else:
                error_msg = "Record must contain `database` field to use `transform-genbank-location.`"

            print(error_msg, file=stderr)

        json.dump(record, stdout, allow_nan=False, indent=None, separators=',:')
        print()
