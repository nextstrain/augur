Setup

  $ source "$TESTDIR"/_setup.sh

Infer ancestral nucleotide and amino acid sequences.

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree.nwk \
  >  --alignment $TESTDIR/../data/aligned.fasta \
  >  --annotation $TESTDIR/../data/zika_outgroup.gb \
  >  --genes ENV PRO \
  >  --translations $TESTDIR/../data/aa_sequences_%GENE.fasta \
  >  --seed 314159 \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations.json" \
  >  --output-sequences "$CRAMTMP/$TESTFILE/ancestral_sequences.fasta" \
  >  --output-translations "$CRAMTMP/$TESTFILE/ancestral_aa_sequences_%GENE.fasta" > /dev/null

Check that the reference length was correctly exported as the nuc annotation

  $ grep -E "\"(ENV|PRO|nuc)\": {" "$CRAMTMP/$TESTFILE/ancestral_mutations.json"
      "ENV": {
      "PRO": {
      "nuc": {

Check that amino acid sequences exist for the root node of the tree.

  $ grep -A 2 "aa_sequences" "$CRAMTMP/$TESTFILE/ancestral_mutations.json"
        "aa_sequences": {
          "ENV": .* (re)
          "PRO": .* (re)

Check that internal nodes have ancestral amino acid sequences.

  $ grep "NODE" "$CRAMTMP/$TESTFILE/ancestral_aa_sequences_ENV.fasta" | wc -l
  \s*8 (re)
