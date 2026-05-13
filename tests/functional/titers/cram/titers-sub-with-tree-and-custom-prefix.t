Setup

  $ source "$TESTDIR"/_setup.sh

Test titer substitution model with alignment and tree inputs and a custom prefix for the node data attributes in the output.

  $ ${AUGUR} titers sub \
  >   --tree $TESTDIR/../data/tree.nwk \
  >   --titers $TESTDIR/../data/titers.tsv \
  >   --alignment $TESTDIR/../data/aa_seq_HA1.fasta \
  >   --gene-names HA1 \
  >   --attribute-prefix custom_prefix_ \
  >   --output $TMP/titers-sub.json > /dev/null
  Read titers from */data/titers.tsv, found: (glob)
   --- 62 strains
   --- 15 data sources
   --- 272 total measurements
  $ grep custom_prefix_cTiterSub $TMP/titers-sub.json | wc -l
  \s*120 (re)
