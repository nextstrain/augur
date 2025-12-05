Integration tests for augur mask.

  $ source "$TESTDIR"/_setup.sh

Mask a VCF with a BED file and a specified output file.
Note the input VCF file is copied to the default temp dir to ensure that the
temporary variants.vcf.gz_maskTemp is created in the temp dir as well.

  $ cp "$TESTDIR/../data/variants.vcf.gz" ./
  $ ${AUGUR} mask \
  >  --sequences "variants.vcf.gz" \
  >  --mask "$TESTDIR/../data/mask_variants.bed" \
  >  --output "masked_variants.vcf" > /dev/null

  $ diff -u "$TESTDIR/../data/masked_variants.vcf" "masked_variants.vcf"
  $ rm -f "masked_variants.vcf"

Mask sequences with a BED file and a specified output file.

  $ ${AUGUR} mask \
  >  --sequences "$TESTDIR/../data/sequences.fasta" \
  >  --mask "$TESTDIR/../data/mask.bed" \
  >  --output "masked.fasta"
  3 masking sites read from .*/../data/mask.bed (re)
  Removing masked sites from FASTA file.

  $ cat "masked.fasta"
  >sequence_1
  NNGCNG
  $ rm -f "masked.fasta"
