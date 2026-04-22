Setup

  $ source "$TESTDIR"/_setup.sh


Filter using --exclude-invalid on AA sequences, which should drop strain "Seq6_invalid_chars"

  $ ${AUGUR} filter \
  >  --seq-type aa \
  >  --sequence-index "$TESTDIR/../../index/HA1_index.tsv" \
  >  --metadata "$TESTDIR/../../index/HA1.tsv" \
  >  --exclude-invalid \
  >  --output-strains filtered_strains_aa.txt
  1 strain was dropped during filtering
  	1 was dropped because it had invalid characters
  5 strains passed all filters
