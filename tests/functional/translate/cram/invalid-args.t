Setup

  $ source "$TESTDIR"/_setup.sh
  $ export ANC_DATA="$TESTDIR/../../ancestral/data/simple-genome"
  $ export DATA="$TESTDIR/../data/simple-genome"


Input JSON + VCF output is not possible (and vice-versa)

  $ ${AUGUR} translate \
  >  --tree $ANC_DATA/tree.nwk \
  >  --ancestral-sequences $ANC_DATA/nt_muts.ref-seq.json \
  >  --reference-sequence $DATA/reference.gff \
  >  --alignment-output "translations.vcf"
  ERROR: When using a non-VCF input the --alignment-output filename must not be a VCF file
  [2]

  $ ${AUGUR} translate \
  >  --tree $ANC_DATA/tree.nwk \
  >  --ancestral-sequences input.vcf \
  >  --vcf-reference reference_in.fasta \
  >  --reference-sequence $DATA/reference.gff \
  >  --alignment-output "translations.fasta"
  ERROR: When using a VCF input the --alignment-output filename must also be a VCF file
  [2]

The arg --vcf-reference-output needs --alignment-output

  $ ${AUGUR} translate \
  >  --tree $ANC_DATA/tree.nwk \
  >  --ancestral-sequences input.vcf \
  >  --reference-sequence $DATA/reference.gff \
  >  --vcf-reference reference_in.fasta \
  >  --vcf-reference-output reference_out.fasta
  ERROR: The VCF reference output (--vcf-reference-output) needs --alignment-output
  [2]

If VCF input with requested alignment output then we require --vcf-reference-output
(This was not the case in augur 23.1.1 and earlier, where we would automatically create a filename)
  $ ${AUGUR} translate \
  >  --tree $ANC_DATA/tree.nwk \
  >  --ancestral-sequences input.vcf \
  >  --vcf-reference reference_in.fasta \
  >  --reference-sequence $DATA/reference.gff \
  >  --alignment-output "translations.vcf"
  ERROR: When using a VCF input and --alignment-output, we now require you to specify the --vcf-reference-output as well
  [2]

VCF input must have a FASTA reference provided
  $ ${AUGUR} translate \
  >  --tree $ANC_DATA/tree.nwk \
  >  --ancestral-sequences input.vcf \
  >  --reference-sequence $DATA/reference.gff \
  >  --alignment-output "translations.vcf"
  ERROR: A reference FASTA (--vcf-reference) is required with VCF-format input
  [2]
