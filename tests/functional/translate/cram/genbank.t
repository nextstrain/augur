Setup

  $ source "$TESTDIR"/_setup.sh
  $ export ANC_DATA="$TESTDIR/../../ancestral/data/simple-genome"
  $ export DATA="$TESTDIR/../data/simple-genome"

These tests are intended to test variants of GenBank reference file formatting


Remove the mandatory source feature from the file
  $ sed '5,6d' "$DATA/reference.gb" > "reference.no-source.gb"

  $ ${AUGUR} translate \
  >  --tree $ANC_DATA/tree.nwk \
  >  --ancestral-sequences $ANC_DATA/nt_muts.ref-seq.json \
  >  --reference-sequence "reference.no-source.gb" \
  >  --output-node-data "aa_muts.json"
  ERROR: Reference .+ did not define the mandatory source feature. (re)
  [2]

Remove a nucleotide from the ORIGIN sequence so the coordinates don't match the source

  $ sed 's/TGACCATAAA/TGACCATAA/' "$DATA/reference.gb" > "reference.short-origin.gb"

  $ ${AUGUR} translate \
  >  --tree $ANC_DATA/tree.nwk \
  >  --ancestral-sequences $ANC_DATA/nt_muts.ref-seq.json \
  >  --reference-sequence "reference.short-origin.gb" \
  >  --output-node-data "aa_muts.json"
  .+ BiopythonParserWarning: .+ (re)
  .+ (re)
  ERROR: Reference .+ (re)
  [2]
