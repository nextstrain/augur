Integration tests for augur parse.

  $ source "$TESTDIR"/_setup.sh
  $ pushd "$TESTDIR" > /dev/null

Try to parse Zika sequences without specifying fields.
This should fail.

  $ ${AUGUR} parse \
  >   --sequences parse/zika.fasta \
  >   --output-sequences "$TMP/sequences.fasta" \
  >   --output-metadata "$TMP/metadata.tsv"
  usage: .* (re)
  .* (re)
  .* (re)
  .* (re)
  .* (re)
  .* (re)
  augur parse: error: the following arguments are required: --fields
  [2]

Parse Zika sequences into sequences and metadata.

  $ ${AUGUR} parse \
  >   --sequences parse/zika.fasta \
  >   --output-sequences "$TMP/sequences.fasta" \
  >   --output-metadata "$TMP/metadata.tsv" \
  >   --fields strain virus accession date region country division city db segment authors url title journal paper_url \
  >   --prettify-fields region country division city \
  >   --fix-dates monthfirst

  $ diff -u "parse/sequences.fasta" "$TMP/sequences.fasta"
  $ diff -u "parse/metadata.tsv" "$TMP/metadata.tsv"
  $ rm -f "$TMP/sequences.fasta" "$TMP/metadata.tsv"

Parse Zika sequences into sequences and metadata using a different metadata field as record id (e.g. accession)

  $ ${AUGUR} parse \
  >   --sequences parse/zika.fasta \
  >   --output-sequences "$TMP/sequences.fasta" \
  >   --output-metadata "$TMP/metadata.tsv" \
  >   --output-id-field accession \
  >   --fields strain virus accession date region country division city db segment authors url title journal paper_url \
  >   --prettify-fields region country division city \
  >   --fix-dates monthfirst

  $ diff -u "parse/sequences_other.fasta" "$TMP/sequences.fasta"
  $ diff -u "parse/metadata_other.tsv" "$TMP/metadata.tsv"
  $ rm -f "$TMP/sequences.fasta" "$TMP/metadata.tsv"

Try to parse Zika sequences with a misspelled field.
This should fail.

  $ ${AUGUR} parse \
  >   --sequences parse/zika.fasta \
  >   --output-sequences "$TMP/sequences.fasta" \
  >   --output-metadata "$TMP/metadata.tsv" \
  >   --output-id-field notexist \
  >   --fields strain virus accession date region country division city db segment authors url title journal paper_url \
  >   --prettify-fields region country division city \
  >   --fix-dates monthfirst
  ERROR: Output id field 'notexist' not found in fields ['strain', 'virus', 'accession', 'date', 'region', 'country', 'division', 'city', 'db', 'segment', 'authors', 'url', 'title', 'journal', 'paper_url'].
  [2]

Parse Zika sequences into sequences and metadata, preferred default ids is 'name', then 'strain', then first field.

  $ ${AUGUR} parse \
  >   --sequences parse/zika.fasta \
  >   --output-sequences "$TMP/sequences.fasta" \
  >   --output-metadata "$TMP/metadata.tsv" \
  >   --fields strain virus name date region country division city db segment authors url title journal paper_url \
  >   --output-id-field 'name' \
  >   --prettify-fields region country division city \
  >   --fix-dates monthfirst

  $ diff -u "parse/sequences_other.fasta" "$TMP/sequences.fasta"
  $ rm -f "$TMP/sequences.fasta" "$TMP/metadata.tsv"

Parse Zika sequences into sequences and metadata when there is no 'name' field.
This should use the 2nd entry in DEFAULT_ID_COLUMNS ('name', 'strain') instead.

  $ ${AUGUR} parse \
  >   --sequences parse/zika.fasta \
  >   --output-sequences "$TMP/sequences.fasta" \
  >   --output-metadata "$TMP/metadata.tsv" \
  >   --fields col1 virus strain date region country division city db segment authors url title journal paper_url \
  >   --prettify-fields region country division city \
  >   --fix-dates monthfirst

  $ diff -u "parse/sequences_other.fasta" "$TMP/sequences.fasta"
  $ rm -f "$TMP/sequences.fasta" "$TMP/metadata.tsv"

Parse Zika sequences into sequences and metadata when no output-id-field is provided and none of the fields match DEFAULT_ID_COLUMNS (e.g. ('strain', 'name')).
This should use the first field as the id field and the metadata should not have an extra strain or name column.

  $ ${AUGUR} parse \
  >   --sequences parse/zika.fasta \
  >   --output-sequences "$TMP/sequences.fasta" \
  >   --output-metadata "$TMP/metadata.tsv" \
  >   --fields col1 virus col3 date region country division city db segment authors url title journal paper_url \
  >   --prettify-fields region country division city \
  >   --fix-dates monthfirst

  $ diff -u "parse/sequences.fasta" "$TMP/sequences.fasta"
  $ diff "parse/metadata.tsv" "$TMP/metadata.tsv" | tr '>' '+' | tr '<' '-'
  1c1
  - strain\tvirus\taccession\tdate\tregion\tcountry\tdivision\tcity\tdb\tsegment\tauthors\turl\ttitle\tjournal\tpaper_url (esc)
  ---
  + col1\tvirus\tcol3\tdate\tregion\tcountry\tdivision\tcity\tdb\tsegment\tauthors\turl\ttitle\tjournal\tpaper_url (esc)
  [1]
  $ rm -f "$TMP/sequences.fasta" "$TMP/metadata.tsv"

Parse compressed Zika sequences into sequences and metadata.

  $ ${AUGUR} parse \
  >   --sequences parse/zika.fasta.gz \
  >   --output-sequences "$TMP/sequences.fasta" \
  >   --output-metadata "$TMP/metadata.tsv" \
  >   --fields strain virus accession date region country division city db segment authors url title journal paper_url \
  >   --prettify-fields region country division city \
  >   --fix-dates monthfirst

  $ diff -u "parse/sequences.fasta" "$TMP/sequences.fasta"
  $ diff -u "parse/metadata.tsv" "$TMP/metadata.tsv"
  $ rm -f "$TMP/sequences.fasta" "$TMP/metadata.tsv"

Error on the first duplicate.

  $ cat >$TMP/data.fasta <<~~
  > >SEQ1
  > AAA
  > >SEQ1
  > AAA
  > >SEQ2
  > AAA
  > >SEQ2
  > AAA
  > ~~
  $ ${AUGUR} parse \
  >   --sequences $TMP/data.fasta \
  >   --output-sequences "$TMP/sequences.fasta" \
  >   --output-metadata "$TMP/metadata.tsv" \
  >   --fields strain
  ERROR: Duplicate found for 'SEQ1'.
  [2]

Run without --fix-dates. The date is left unchanged.

  $ cat >$TMP/data.fasta <<~~
  > >SEQ1|05/01/2020
  > AAA
  > ~~
  $ ${AUGUR} parse \
  >   --sequences $TMP/data.fasta \
  >   --output-sequences "$TMP/sequences.fasta" \
  >   --output-metadata "$TMP/metadata.tsv" \
  >   --fields strain date

  $ cat "$TMP/metadata.tsv"
  strain	date
  SEQ1	05/01/2020
  $ rm -f "$TMP/metadata.tsv"

  $ popd > /dev/null
