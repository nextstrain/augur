Setup

  $ source "$TESTDIR"/_setup.sh

Infer ancestral nucleotide and amino acid sequences, using a genes file.

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree.nwk \
  >  --alignment $TESTDIR/../data/aligned.fasta \
  >  --annotation $TESTDIR/../data/zika_outgroup.gb \
  >  --genes $TESTDIR/../data/genes.txt \
  >  --translations $TESTDIR/../data/aa_sequences_%GENE.fasta \
  >  --seed 314159 \
  >  --output-node-data ancestral_mutations.json \
  >  --output-sequences ancestral_sequences.fasta \
  >  --output-translations ancestral_aa_sequences_%GENE.fasta > /dev/null

Check that the reference length was correctly exported as the nuc annotation

  $ grep -E "\"(ENV|PRO|nuc)\": {" ancestral_mutations.json
      "ENV": {
      "PRO": {
      "nuc": {
