Setup

  $ source "$TESTDIR"/_setup.sh

Input FASTA + VCF output is not possible

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree.nwk \
  >  --alignment $TESTDIR/../data/aligned.fasta \
  >  --output-vcf "output.vcf" > /dev/null
  ERROR: VCF output has been requested but the input alignment is not VCF.
  [2]

Input VCF + FASTA output is not possible (Note that the input file doesn't exist, but we exit before that's checked)

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree.nwk \
  >  --alignment $TESTDIR/../data/snps.vcf \
  >  --output-sequences "output.fasta" > /dev/null
  ERROR: Sequence (fasta) output has been requested but the input alignment is VCF.
  [2]

Output FASTA _and_ VCF is not possible

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree.nwk \
  >  --alignment $TESTDIR/../data/aligned.fasta \
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
  >  --output-node-data "ancestral_mutations.json" > /dev/null
  ERROR: For amino acid sequence reconstruction, you must provide an annotation file, a list of genes, and a template path to amino acid sequences.
  [2]
