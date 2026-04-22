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


Test single gene reconstruction with an explicitly provided AA root-sequence
(See infer-amino-acid-sequences-with-root-sequence.t for using a nuc sequence as --root-sequence)

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree.nwk \
  >  --annotation $TESTDIR/../data/zika_outgroup.gb \
  >  --genes ENV \
  >  --translations $TESTDIR/../data/aa_sequences_ENV.fasta \
  >  --aa-root-sequence $TESTDIR/../data/ENV_outgroup.fasta \
  >  --seed 314159 \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations_4.json" &> /dev/null

The reference has been modified to include a leading AAA:

  $ grep -A 2 "reference" "$CRAMTMP/$TESTFILE/ancestral_mutations_4.json"
    "reference": {
      "ENV": "AAA.+ (re)
    }

Check that this results in 3 mutations between the provided root-sequence & the inferred root node

  $ grep -A 11 "NODE_0000006" "$CRAMTMP/$TESTFILE/ancestral_mutations_4.json"
      "NODE_0000006": {
        "aa_muts": {
          "ENV": [
            "A1I",
            "A2R",
            "A3C"
          ]
        },
        "aa_sequences": {
          "ENV": "IRC.+ (re)
        }
      },


Retest the above, using a %GENE placeholder in --aa-root-sequence

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree.nwk \
  >  --annotation $TESTDIR/../data/zika_outgroup.gb \
  >  --genes ENV \
  >  --translations $TESTDIR/../data/aa_sequences_ENV.fasta \
  >  --aa-root-sequence $TESTDIR/../data/%GENE_outgroup.fasta \
  >  --seed 314159 \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations_5.json" &> /dev/null


  $ diff "$CRAMTMP/$TESTFILE/ancestral_mutations_4.json" "$CRAMTMP/$TESTFILE/ancestral_mutations_5.json"


Test that accidentally providing the wrong AA root-sequence (e.g. a nuc one) results in an error

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree.nwk \
  >  --annotation $TESTDIR/../data/zika_outgroup.gb \
  >  --genes ENV \
  >  --translations $TESTDIR/../data/aa_sequences_ENV.fasta \
  >  --aa-root-sequence $TESTDIR/../data/simple-genome/reference.fasta \
  >  --seed 314159 \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations_5.json" > /dev/null
  ERROR: The provided root-sequence AA fasta for ENV has length 50 which doesn't match the length of the CDS 504 (amino acids)
  [2]
