Setup

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="${AUGUR:-../../../../bin/augur}"

Test titer substitution model with alignment and tree inputs.

  $ ${AUGUR} titers sub \
  >   --tree ../data/tree.nwk \
  >   --titers ../data/titers.tsv \
  >   --alignment ../data/aa_seq_HA1.fasta \
  >   --gene-names HA1 \
  >   --output $TMP/titers-sub.json > /dev/null
  Read titers from ../data/titers.tsv, found:
   --- 62 strains
   --- 15 data sources
   --- 272 total measurements
  $ grep cTiterSub $TMP/titers-sub.json | wc -l
  \s*120 (re)
