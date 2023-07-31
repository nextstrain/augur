Setup

  $ source "$TESTDIR"/_setup.sh

Infer ancestral nucleotide and amino acid sequences.

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree.nwk \
  >  --alignment $TESTDIR/../data/aligned.fasta \
  >  --annotation $TESTDIR/../data/zika_outgroup.gb \
  >  --genes ENV PRO \
  >  --translations $TESTDIR/../data/aa_sequences_%GENE.fasta \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations.json" \
  >  --output-sequences "$CRAMTMP/$TESTFILE/ancestral_sequences.fasta" > /dev/null

Check that the reference length was correctly exported as the nuc annotation

  $ grep -E "\"(ENV|PRO|nuc)\": {" "$CRAMTMP/$TESTFILE/ancestral_mutations.json"
      "ENV": {
      "PRO": {
      "nuc": {
