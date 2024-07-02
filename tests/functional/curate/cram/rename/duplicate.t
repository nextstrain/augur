Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

  $ cat >records.ndjson <<~~
  > {"accession": "record_1", "country": "country_1"}
  > {"accession": "record_2", "country": "country_2"}
  > ~~


Asking to rename the same column multiple times results in duplication of the column.
Additional columns are inserted next to the existing one, and the order of the new columns
matches the field-map

  $ $AUGUR curate rename --field-map "accession=id" "accession=genbank_accession" < <(cat records.ndjson)
  {"id": "record_1", "genbank_accession": "record_1", "country": "country_1"}
  {"id": "record_2", "genbank_accession": "record_2", "country": "country_2"}


We can use the same name to keep the original column

  $ $AUGUR curate rename --field-map "accession=id" "accession=accession" < <(cat records.ndjson)
  {"id": "record_1", "accession": "record_1", "country": "country_1"}
  {"id": "record_2", "accession": "record_2", "country": "country_2"}
