Integration tests for augur mask.

  $ source "$TESTDIR"/_setup.sh

Try masking a VCF without any specified mask.

  $ ${AUGUR} mask --sequences "$TESTDIR/../data/variants.vcf.gz"
  No masking sites provided. Must include one of --mask, --mask-from-beginning, --mask-from-end, --mask-invalid, or --mask-sites
  [1]

Mask a VCF with a BED file and no specified output file.

  $ cp "$TESTDIR/../data/variants.vcf.gz" ./
  $ ${AUGUR} mask \
  >  --sequences "variants.vcf.gz" \
  >  --mask "$TESTDIR/../data/mask_variants.bed" > /dev/null

  $ diff -u "$TESTDIR/../data/masked_variants.vcf" <(gzip -c -d "variants.vcf.gz")
  $ rm -f "variants.vcf.gz"

Mask a VCF with a BED file and a specified output file.

  $ ${AUGUR} mask \
  >  --sequences "$TESTDIR/../data/variants.vcf.gz" \
  >  --mask "$TESTDIR/../data/mask_variants.bed" \
  >  --output "masked_variants.vcf" > /dev/null

  $ diff -u "$TESTDIR/../data/masked_variants.vcf" "masked_variants.vcf"
  $ rm -f "masked_variants.vcf"

Try masking sequences without any specified mask.

  $ ${AUGUR} mask --sequences "$TESTDIR/../data/sequences.fasta"
  No masking sites provided. Must include one of --mask, --mask-from-beginning, --mask-from-end, --mask-invalid, or --mask-sites
  [1]

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

Mask one base from the beginning and the end.

  $ ${AUGUR} mask \
  >  --sequences "$TESTDIR/../data/sequences.fasta" \
  >  --mask-from-beginning 1 \
  >  --mask-from-end 1 \
  >  --output "masked.fasta"
  Removing masked sites from FASTA file.

  $ cat "masked.fasta"
  >sequence_1
  NTGCTN
  $ rm -f "masked.fasta"

Mask a specific list of sites and also mask one base from the beginning and the end.

  $ ${AUGUR} mask \
  >  --sequences "$TESTDIR/../data/sequences.fasta" \
  >  --mask-sites 3 4 \
  >  --mask-from-beginning 1 \
  >  --mask-from-end 1 \
  >  --output "masked.fasta"
  Removing masked sites from FASTA file.

  $ cat "masked.fasta"
  >sequence_1
  NTNNTN
  $ rm -f "masked.fasta"

Mask invalid nucleotides

  $ ${AUGUR} mask \
  >  --sequences "$TESTDIR/../data/invalidnucleotide.fasta" \
  >  --mask-invalid \
  >  --output "masked.fasta"
  Removing masked sites from FASTA file.

  $ cat "masked.fasta"
  >sequence_1
  ATCGNNNN
  $ rm -f "masked.fasta"
