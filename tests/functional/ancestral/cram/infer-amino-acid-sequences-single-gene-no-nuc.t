Setup

  $ source "$TESTDIR"/_setup.sh

Infer ancestral amino acid sequences (no nucleotide reconstruction)

Firstly a single gene (ENV), using a hardcoded fasta path, and an annotation file

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


For single genes the annotations file is optional. The result should be the same except the (nuc) coordinates of the CDS will differ without the annotation file

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree.nwk \
  >  --genes ENV \
  >  --translations $TESTDIR/../data/aa_sequences_ENV.fasta \
  >  --seed 314159 \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations_3.json" \
  >  --output-translations "$CRAMTMP/$TESTFILE/ancestral_aa_sequences_ENV_3.fasta" &> /dev/null


  $ diff "$CRAMTMP/$TESTFILE/ancestral_aa_sequences_ENV_1.fasta" "$CRAMTMP/$TESTFILE/ancestral_aa_sequences_ENV_3.fasta"

Nucleotide coordinates of the ENV gene are different - without the annotation file we start them at nuc_pos=1 (1-based, GFF-style)
So for this example the offset is 960 (start is 961 - 960 = 1, end is 2472 - 960 = 1512) 

  $ diff "$CRAMTMP/$TESTFILE/ancestral_mutations_1.json" "$CRAMTMP/$TESTFILE/ancestral_mutations_3.json" 
  4,6c4,5
  \<       "end": 2472, (re)
  \<       "seqid": .* (re)
  \<       "start": 961, (re)
  ---
  \>       "end": 1512, (re)
  \>       "start": 1, (re)
  [1]
