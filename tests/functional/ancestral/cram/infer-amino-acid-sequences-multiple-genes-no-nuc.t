Setup

  $ source "$TESTDIR"/_setup.sh

Infer ancestral amino acid sequences (no nucleotide reconstruction)

Infer multiple genes with a provided GenBank annotation file

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree.nwk \
  >  --annotation $TESTDIR/../data/zika_outgroup.gb \
  >  --genes ENV PRO \
  >  --translations $TESTDIR/../data/aa_sequences_%GENE.fasta \
  >  --seed 314159 \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations_1.json" \
  >  --output-translations "$CRAMTMP/$TESTFILE/ancestral_aa_sequences_%GENE_1.fasta" &> /dev/null

Check that the annotations block only includes ENV & PRO, not nuc

  $ grep -E "\"(ENV|PRO|nuc)\": {" "$CRAMTMP/$TESTFILE/ancestral_mutations_1.json"
      "ENV": {
      "PRO": {

Check that amino acid sequences exist for the root node of the tree.

  $ grep -A 2 "aa_sequences" "$CRAMTMP/$TESTFILE/ancestral_mutations_1.json"
        "aa_sequences": {
          "ENV": .* (re)
          "PRO": .* (re)

Check that internal nodes have ancestral amino acid sequences.

  $ grep "NODE" "$CRAMTMP/$TESTFILE/ancestral_aa_sequences_ENV_1.fasta" | wc -l
  \s*8 (re)

  $ grep "NODE" "$CRAMTMP/$TESTFILE/ancestral_aa_sequences_PRO_1.fasta" | wc -l
  \s*8 (re)


Test that using a list of genes (ENV, PRO) and no %GENE pattern throws an error
Catches this bug <https://github.com/nextstrain/augur/pull/1958#discussion_r3034248147>

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree.nwk \
  >  --annotation $TESTDIR/../data/zika_outgroup.gb \
  >  --genes $TESTDIR/../data/genes.txt \
  >  --translations $TESTDIR/../data/aa_sequences_ENV.fasta \
  >  --seed 314159 \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations_error.json" \
  >  --output-translations "$CRAMTMP/$TESTFILE/ancestral_aa_sequences_error.fasta" > /dev/null
  ERROR: --translations must contain %GENE for multiple-gene amino acid reconstructions
  [2]
