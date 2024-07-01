Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Note about using NDJSON because order will be reflected in TSV?

  $ cat >records.ndjson <<~~
  > {"accession": "record_1", "country": "country_1"}
  > {"accession": "record_2", "country": "country_2"}
  > ~~

Test that running without (rename-specific) arguments doesn't modify data

  $ $AUGUR curate rename < <(cat records.ndjson) > out1.ndjson

  $ diff records.ndjson out1.ndjson

Rename "accession" to "strain" -- the order should be preserved (i.e. strain is first column)

  $ $AUGUR curate rename --field-map "accession=strain" < <(cat records.ndjson)
  {"strain": "record_1", "country": "country_1"}
  {"strain": "record_2", "country": "country_2"}


Rename "accession" to "country" - single error message is reported as we won't overwrite without --force
and we don't change the data

  $ $AUGUR curate rename --field-map "accession=country" < <(cat records.ndjson) > out2.ndjson
  WARNING: skipping rename of accession because record already has a field named country.

  $ diff records.ndjson out2.ndjson

Rename "accession" to "country" using --force

  $ $AUGUR curate rename --field-map "accession=country" --force < <(cat records.ndjson)
  {"country": "record_1"}
  {"country": "record_2"}

Asking to rename multiple columns to the same new name is an error!

  $ AUGUR curate rename --field-map "accession=foo" "country=foo" < <(cat records.ndjson)
  ERROR: Asked to rename multiple fields to 'foo'.
  [2]

Asking to rename the same column multiple times is an error
(we may one day extend the usage to have this duplicate columns)

  $ AUGUR curate rename --field-map "country=foo" "country=bar" < <(cat records.ndjson)
  ERROR: Asked to rename field 'country' multiple times.
  [2]

Asking to rename a non-existant column is an error

  $ AUGUR curate rename --field-map "strain=foo" < <(cat records.ndjson)
  ERROR: Asked to rename field 'strain' (to 'foo') but it doesn't exist in the input data.
  [2]
