Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"


  $ cat >records.ndjson <<~~
  > {"accession": "record_1", "country": ""}
  > {"accession": "record_2", "country": "country_2"}
  > ~~

Test that --force is required if a key exists, irregardless of its associated value
(This was not the behaviour in the precursor command `transform-field-names`)

  $ $AUGUR curate rename --field-map "accession=country" < <(cat records.ndjson) 1> out1.ndjson
  WARNING: skipping rename of accession because record already has a field named country.

  $ diff records.ndjson out1.ndjson

