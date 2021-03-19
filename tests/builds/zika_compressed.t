Run an example Zika build with Augur using compressed inputs and outputs where possible.

Setup test data directory, output directory, and temporarily change directories to the test data directory.
Some augur commands store the exact path of their inputs in their outputs.
Running from the test data directory allows us to use relative paths that won't differ between execution environments.

  $ TEST_DATA_DIR="$TESTDIR/zika"
  $ mkdir -p "$TMP/out"
  $ pushd "$TEST_DATA_DIR" > /dev/null
  $ export AUGUR="../../../bin/augur"

Parse a FASTA whose defline contains metadata into separate sequence and metadata files.

  $ ${AUGUR} parse \
  >   --sequences "data/zika.fasta.gz" \
  >   --output-sequences "$TMP/out/sequences.fasta.gz" \
  >   --output-metadata "$TMP/out/metadata.tsv.gz" \
  >   --fields strain virus accession date region country division city db segment authors url title journal paper_url \
  >   --prettify-fields region country division city

  $ diff -u <(gzip -c -d "results/sequences.fasta.gz") <(gzip -c -d "$TMP/out/sequences.fasta.gz")

Index sequence composition to speed up filters.

  $ ${AUGUR} index \
  >   --sequences "results/sequences.fasta.gz" \
  >   --output "$TMP/out/sequence_index.tsv.gz" \
  >   --verbose
  Analysed 12 sequences with an average length of 10598 nucleotides.

  $ diff -u <(gzip -c -d "results/sequence_index.tsv.gz") <(gzip -c -d "$TMP/out/sequence_index.tsv.gz")

Filter sequences by a minimum date and an exclusion list and only keep one sequence per country, year, and month.

  $ ${AUGUR} filter \
  >   --sequences "results/sequences.fasta.gz" \
  >   --sequence-index "results/sequence_index.tsv.gz" \
  >   --metadata "results/metadata.tsv.gz" \
  >   --exclude "config/dropped_strains.txt" \
  >   --output "$TMP/out/filtered.fasta.gz" \
  >   --group-by country year month \
  >   --sequences-per-group 1 \
  >   --subsample-seed 314159 \
  >   --no-probabilistic-sampling \
  >   --min-date 2012 > /dev/null
  $ gzip -c -d "$TMP/out/filtered.fasta.gz" | grep "^>" | wc -l
  \s*10 (re)

Switch back to the original directory where testing started.

  $ popd > /dev/null
