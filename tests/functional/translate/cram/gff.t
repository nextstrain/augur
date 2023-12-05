Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../bin/augur}"
  $ export SCRIPTS="$TESTDIR/../../../../scripts"
  $ export ANC_DATA="$TESTDIR/../../ancestral/data/simple-genome"
  $ export DATA="$TESTDIR/../data/simple-genome"

These tests are intended to test variants of GFF formatting


GFF file with no valid rows

  $ head -n 3  $DATA/reference.source.gff > "reference.empty.gff"

  $ ${AUGUR} translate \
  >  --tree $ANC_DATA/tree.nwk \
  >  --ancestral-sequences $ANC_DATA/nt_muts.ref-seq.json \
  >  --reference-sequence "reference.empty.gff" \
  >  --output-node-data "aa_muts.json" > /dev/null
  ERROR: Reference 'reference.empty.gff' contains no valid data rows. .+ (re)
  [2]

GFF file with an extra record

  $ cp $DATA/reference.source.gff "reference.double.gff"

  $ echo -e "additional\tRefSeq\tsource\t1\t10\t.\t+\t.\tID=additional" >> "reference.double.gff"

  $ ${AUGUR} translate \
  >  --tree $ANC_DATA/tree.nwk \
  >  --ancestral-sequences $ANC_DATA/nt_muts.ref-seq.json \
  >  --reference-sequence "reference.double.gff" \
  >  --output-node-data "aa_muts.json"
  ERROR: Reference 'reference.double.gff' contains multiple seqids .+ (re)
  [2]
