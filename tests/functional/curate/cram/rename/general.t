Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

The tests here use NDJSON I/O and don't explicitly test TSV I/O as we rely
on the general curate infrastructure to enforce that each row has the same
fields. See <https://github.com/nextstrain/augur/issues/1510> for more

  $ cat >records.ndjson <<~~
  > {"accession": "record_1", "country": "country_1"}
  > {"accession": "record_2", "country": "country_2"}
  > ~~

The --field-map argument is requried (error text not checked as it includes the entire argparse usage text)

  $ $AUGUR curate rename < <(cat records.ndjson) 2>/dev/null
  [2]

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

  $ $AUGUR curate rename --field-map "accession=foo" "country=foo" < <(cat records.ndjson)
  ERROR: Asked to rename multiple fields to 'foo'.
  [2]

Asking to rename a non-existant column raises a warning. Using a warning not an error allows the command
to be idempotent.

  $ $AUGUR curate rename --field-map "strain=foo" < <(cat records.ndjson) 1> out3.ndjson
  WARNING: Asked to rename field 'strain' (to 'foo') but it doesn't exist in the input data.

  $ diff records.ndjson out3.ndjson


Rename will re-order fields to match the first observed record (with any necessary changes applied)
This produces NDJSON output which more closely matches TSV output.

  $ cat >records.unordered.ndjson <<~~
  > {"accession": "record_1", "country": "country_1"}
  > {"country": "country_2", "accession": "record_2"}
  > ~~

  $ $AUGUR curate rename --field-map "accession=id" --force < <(cat records.unordered.ndjson)
  {"id": "record_1", "country": "country_1"}
  {"id": "record_2", "country": "country_2"}


Using --field-map without a single '=' char is an error
  $ $AUGUR curate rename --field-map "accession" < <(cat records.ndjson)
  ERROR: The field-map 'accession' must contain a single '=' character.
  [2]

Using --field-map with more than a single '=' char is an error
  $ $AUGUR curate rename --field-map "accession=id=strain" < <(cat records.ndjson)
  ERROR: The field-map 'accession=id=strain' must contain a single '=' character.
  [2]

Using --field-map with spaces surrounding field names is OK (as long as you quote the arg appropriately)
  $ $AUGUR curate rename --field-map " accession = strain " < <(cat records.ndjson)
  {"strain": "record_1", "country": "country_1"}
  {"strain": "record_2", "country": "country_2"}

Using --field-map with an empty "old field" doesn't make sense
  $ $AUGUR curate rename --field-map "=accession" < <(cat records.ndjson)
  ERROR: The field-map '=accession' doesn't specify a name for the existing field.
  [2]

Using --field-map with an empty "new field" doesn't (yet) make sense
  $ $AUGUR curate rename --field-map "accession=" < <(cat records.ndjson)
  ERROR: The field-map 'accession=' doesn't specify a name for the new field.
  [2]
