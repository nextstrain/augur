Setup

  $ source "$TESTDIR"/_setup.sh
  $ export ANC_DATA="$TESTDIR/../../ancestral/data/simple-genome"
  $ export DATA="$TESTDIR/../data/simple-genome"

Missing reference file

  $ ${AUGUR} translate \
  >  --tree $ANC_DATA/tree.nwk \
  >  --ancestral-sequences $ANC_DATA/nt_muts.ref-seq.json \
  >  --reference-sequence $DATA/reference.doesnt-exist.gff \
  >  --output-node-data "aa_muts.json" > /dev/null
  ERROR: reference sequence file '.+/reference.doesnt-exist.gff' not found (re)
  [2]
