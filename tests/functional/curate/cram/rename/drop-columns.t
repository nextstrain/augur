
Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

  $ cat >records.ndjson <<~~
  > {"strain": "s_1", "country": "c_1", "accession": "a_1"}
  > {"strain": "s_2", "country": "c_2", "accession": "a_2"}
  > {"strain": "s_3", "country": "c_3", "accession": "a_3"}
  > ~~

Rename the strain column to accession -- we don't move the strain column to where accession was,
We simply rename it and drop the (old) accession column. The alternate (which I don't like) produces:
{"country": "c_1", "accession": "s_1"}
{"country": "c_2", "accession": "s_2"}
{"country": "c_3", "accession": "s_3"}

  $ $AUGUR curate rename --field-map "strain=accession" --force < <(cat records.ndjson)
  {"accession": "s_1", "country": "c_1"}
  {"accession": "s_2", "country": "c_2"}
  {"accession": "s_3", "country": "c_3"}

Similarly, renaming accession to strain re-names the column "in-place" and drops the (old) strain column

  $ $AUGUR curate rename --field-map "accession=strain" --force < <(cat records.ndjson)
  {"country": "c_1", "strain": "a_1"}
  {"country": "c_2", "strain": "a_2"}
  {"country": "c_3", "strain": "a_3"}

