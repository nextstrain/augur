Setup

  $ source "$TESTDIR"/_setup.sh

Test titer substitution model with alignment and tree inputs.

  $ ${AUGUR} titers sub \
  >   --tree $TESTDIR/../data/tree.nwk \
  >   --titers $TESTDIR/../data/titers.tsv \
  >   --alignment $TESTDIR/../data/aa_seq_HA1.fasta \
  >   --gene-names HA1 \
  >   --output titers-sub.json > /dev/null
  Read titers from */data/titers.tsv, found: (glob)
   --- 62 strains
   --- 15 data sources
   --- 272 total measurements
  $ grep cTiterSub titers-sub.json | wc -l
  \s*120 (re)

Verify that the titer drops assigned per branch correspond to the expected values for this dataset.
In this example, we know that the HA1 amino acid sequence for A/Fujian/445/2003 carries a S193N substitution and that the titer model assigns a weight of 0.6 to that substitution.
The titer model assigns a higher weight of 1.22 to the opposite substitution N193S.
When we search for that sequence's per-branch titer drop, we should get the smaller value below.

  $ jq -r '.nodes["A/Fujian/445/2003"].dTiterSub' titers-sub.json
  0.6* (glob)
