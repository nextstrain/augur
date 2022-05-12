Integration tests for augur parse.

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="../../bin/augur"

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

  $ popd > /dev/null
