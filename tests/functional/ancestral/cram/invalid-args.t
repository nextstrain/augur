Setup

  $ source "$TESTDIR"/_setup.sh

Input FASTA + VCF output is not possible

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree.nwk \
  >  --alignment $TESTDIR/../data/aligned.fasta \
  >  --seed 314159 \
  >  --output-vcf "output.vcf" > /dev/null
  ERROR: VCF output has been requested but the input alignment is not VCF.
  [2]

Input VCF + FASTA output is not possible (Note that the input file doesn't exist, but we exit before that's checked)

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree.nwk \
  >  --alignment $TESTDIR/../data/snps.vcf \
  >  --seed 314159 \
  >  --output-sequences "output.fasta" > /dev/null
  ERROR: Sequence (fasta) output has been requested but the input alignment is VCF.
  [2]

Output FASTA _and_ VCF is not possible

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree.nwk \
  >  --alignment $TESTDIR/../data/aligned.fasta \
  >  --seed 314159 \
  >  --output-vcf "output.vcf" \
  >  --output-sequences "output.fasta" > /dev/null
  ERROR: Both sequence (fasta) and VCF output have been requested, but these are incompatible.
  [2]


Try to infer ancestral amino acid sequences without all required arguments.
This should fail.

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree.nwk \
  >  --alignment $TESTDIR/../data/aligned.fasta \
  >  --annotation $TESTDIR/../data/zika_outgroup.gb \
  >  --genes ENV PRO \
  >  --seed 314159 \
  >  --output-node-data "ancestral_mutations.json" > /dev/null
  ERROR: For amino acid sequence reconstruction, you must provide an annotation file, a list of genes, and a template path to amino acid sequences.
  [2]

Missing tree file

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree-doesnt-exist.nwk \
  >  --alignment $TESTDIR/../data/aligned.fasta \
  >  --seed 314159 \
  >  --output-sequences "output.fasta" > /dev/null
  ERROR: The provided tree file .* doesn't exist (re)
  [2]


Attempting to use FASTA-input reference and VCF-input reference args
(The files here don't exist, but we exit before they're checked)

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree-doesnt-exist.nwk \
  >  --alignment $TESTDIR/../data/aligned.fasta \
  >  --root-sequence $TESTDIR/../data/reference.fasta \
  >  --vcf-reference $TESTDIR/../data/reference.fasta \
  >  --seed 314159 \
  >  --output-sequences "output.fasta" > /dev/null 2>"err-args.txt"
  [2]

  $ grep "augur ancestral: error: argument --vcf-reference: not allowed with argument --root-sequence" "err-args.txt"
  augur ancestral: error: argument --vcf-reference: not allowed with argument --root-sequence
