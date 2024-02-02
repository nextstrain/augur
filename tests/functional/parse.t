Integration tests for augur parse.

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="${AUGUR:-../../bin/augur}"

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

Parse Zika sequences into sequences and metadata, preferred default ids is 'name', then 'strain', then first field.

  $ ${AUGUR} parse \
  >   --sequences parse/zika.fasta \
  >   --output-sequences "$TMP/sequences.fasta" \
  >   --output-metadata "$TMP/metadata.tsv" \
  >   --fields strain virus name date region country division city db segment authors url title journal paper_url \
  >   --prettify-fields region country division city \
  >   --fix-dates monthfirst

  $ diff -u "parse/sequences_other.fasta" "$TMP/sequences.fasta"
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
