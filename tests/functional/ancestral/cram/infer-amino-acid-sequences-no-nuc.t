Setup

  $ source "$TESTDIR"/_setup.sh

Infer ancestral amino acid sequences (no nucleotide reconstruction)

Firstly a single gene (ENV), using a hardcoded fasta path

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree.nwk \
  >  --annotation $TESTDIR/../data/zika_outgroup.gb \
  >  --genes ENV \
  >  --translations $TESTDIR/../data/aa_sequences_ENV.fasta \
  >  --seed 314159 \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations_1.json" \
  >  --output-translations "$CRAMTMP/$TESTFILE/ancestral_aa_sequences_ENV_1.fasta" &> /dev/null

Check that the annotations block only includes ENV, not nuc or PRO

  $ grep -E "\"(ENV|PRO|nuc)\": {" "$CRAMTMP/$TESTFILE/ancestral_mutations_1.json"
      "ENV": {

Check that amino acid sequences exist for the root node of the tree.

  $ grep -A 2 "aa_sequences" "$CRAMTMP/$TESTFILE/ancestral_mutations_1.json"
        "aa_sequences": {
          "ENV": .* (re)
        }

Check that internal nodes have ancestral amino acid sequences.

  $ grep "NODE" "$CRAMTMP/$TESTFILE/ancestral_aa_sequences_ENV_1.fasta" | wc -l
  \s*8 (re)


And the exact same, but using a %GENE pattern in the filepath

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree.nwk \
  >  --annotation $TESTDIR/../data/zika_outgroup.gb \
  >  --genes ENV \
  >  --translations $TESTDIR/../data/aa_sequences_%GENE.fasta \
  >  --seed 314159 \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations_2.json" \
  >  --output-translations "$CRAMTMP/$TESTFILE/ancestral_aa_sequences_%GENE_2.fasta" &> /dev/null

Check that the outputs are identical

  $ diff "$CRAMTMP/$TESTFILE/ancestral_mutations_1.json" "$CRAMTMP/$TESTFILE/ancestral_mutations_2.json" 

  $ diff "$CRAMTMP/$TESTFILE/ancestral_aa_sequences_ENV_1.fasta" "$CRAMTMP/$TESTFILE/ancestral_aa_sequences_ENV_2.fasta"




Finally multiple genes

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree.nwk \
  >  --annotation $TESTDIR/../data/zika_outgroup.gb \
  >  --genes ENV PRO \
  >  --translations $TESTDIR/../data/aa_sequences_%GENE.fasta \
  >  --seed 314159 \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations_3.json" \
  >  --output-translations "$CRAMTMP/$TESTFILE/ancestral_aa_sequences_%GENE_3.fasta" &> /dev/null

Check that the annotations block only includes ENV & PRO, not nuc

  $ grep -E "\"(ENV|PRO|nuc)\": {" "$CRAMTMP/$TESTFILE/ancestral_mutations_3.json"
      "ENV": {
      "PRO": {

Check that amino acid sequences exist for the root node of the tree.

  $ grep -A 2 "aa_sequences" "$CRAMTMP/$TESTFILE/ancestral_mutations_3.json"
        "aa_sequences": {
          "ENV": .* (re)
          "PRO": .* (re)

Check that internal nodes have ancestral amino acid sequences.

  $ grep "NODE" "$CRAMTMP/$TESTFILE/ancestral_aa_sequences_ENV_3.fasta" | wc -l
  \s*8 (re)

  $ grep "NODE" "$CRAMTMP/$TESTFILE/ancestral_aa_sequences_PRO_3.fasta" | wc -l
  \s*8 (re)
