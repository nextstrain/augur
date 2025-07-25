Integration tests for augur mask.

  $ source "$TESTDIR"/_setup.sh

Mask a VCF with a BED file and no specified output file.

  $ cp "$TESTDIR/../data/variants.vcf.gz" ./
  $ ${AUGUR} mask \
  >  --sequences "variants.vcf.gz" \
  >  --mask "$TESTDIR/../data/mask_variants.bed" > /dev/null

  $ diff -u "$TESTDIR/../data/masked_variants.vcf" <(gzip -c -d "variants.vcf.gz")
  $ rm -f "variants.vcf.gz"

Mask sequences with a BED file and no specified output file.
Since no output is provided, the input file is overridden with the masked sequences.

  $ cp "$TESTDIR/../data/sequences.fasta" "./"
  $ ${AUGUR} mask --sequences "sequences.fasta" --mask "$TESTDIR/../data/mask.bed"
  3 masking sites read from .*/../data/mask.bed (re)
  Removing masked sites from FASTA file.

  $ cat "sequences.fasta"
  >sequence_1
  NNGCNG
  $ rm -f "sequences.fasta"
