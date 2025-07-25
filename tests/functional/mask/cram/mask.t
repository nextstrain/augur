Integration tests for augur mask.

  $ source "$TESTDIR"/_setup.sh
  $ pushd "$TESTDIR" > /dev/null

Try masking a VCF without any specified mask.

  $ ${AUGUR} mask --sequences ../data/variants.vcf.gz
  No masking sites provided. Must include one of --mask, --mask-from-beginning, --mask-from-end, --mask-invalid, or --mask-sites
  [1]

Mask a VCF with a BED file and no specified output file.

  $ cp "../data/variants.vcf.gz" "$TMP/"
  $ ${AUGUR} mask \
  >  --sequences "$TMP/variants.vcf.gz" \
  >  --mask "../data/mask_variants.bed" > /dev/null

  $ diff -u "../data/masked_variants.vcf" <(gzip -c -d "$TMP/variants.vcf.gz")
  $ rm -f "$TMP/variants.vcf.gz"

Mask a VCF with a BED file and a specified output file.

  $ ${AUGUR} mask \
  >  --sequences "../data/variants.vcf.gz" \
  >  --mask "../data/mask_variants.bed" \
  >  --output "$TMP/masked_variants.vcf" > /dev/null

  $ diff -u "../data/masked_variants.vcf" "$TMP/masked_variants.vcf"
  $ rm -f "$TMP/masked_variants.vcf"

Try masking sequences without any specified mask.

  $ ${AUGUR} mask --sequences ../data/sequences.fasta
  No masking sites provided. Must include one of --mask, --mask-from-beginning, --mask-from-end, --mask-invalid, or --mask-sites
  [1]

Mask sequences with a BED file and no specified output file.
Since no output is provided, the input file is overridden with the masked sequences.

  $ cp ../data/sequences.fasta "$TMP/"
  $ ${AUGUR} mask --sequences "$TMP/sequences.fasta" --mask ../data/mask.bed
  3 masking sites read from ../data/mask.bed
  Removing masked sites from FASTA file.

  $ cat "$TMP/sequences.fasta"
  >sequence_1
  NNGCNG
  $ rm -f "$TMP/sequences.fasta"

Mask sequences with a BED file and a specified output file.

  $ ${AUGUR} mask \
  >  --sequences ../data/sequences.fasta \
  >  --mask ../data/mask.bed \
  >  --output "$TMP/masked.fasta"
  3 masking sites read from ../data/mask.bed
  Removing masked sites from FASTA file.

  $ cat "$TMP/masked.fasta"
  >sequence_1
  NNGCNG
  $ rm -f "$TMP/masked.fasta"

Mask one base from the beginning and the end.

  $ ${AUGUR} mask \
  >  --sequences ../data/sequences.fasta \
  >  --mask-from-beginning 1 \
  >  --mask-from-end 1 \
  >  --output "$TMP/masked.fasta"
  Removing masked sites from FASTA file.

  $ cat "$TMP/masked.fasta"
  >sequence_1
  NTGCTN
  $ rm -f "$TMP/masked.fasta"

Mask a specific list of sites and also mask one base from the beginning and the end.

  $ ${AUGUR} mask \
  >  --sequences ../data/sequences.fasta \
  >  --mask-sites 3 4 \
  >  --mask-from-beginning 1 \
  >  --mask-from-end 1 \
  >  --output "$TMP/masked.fasta"
  Removing masked sites from FASTA file.

  $ cat "$TMP/masked.fasta"
  >sequence_1
  NTNNTN
  $ rm -f "$TMP/masked.fasta"

Mask invalid nucleotides

  $ ${AUGUR} mask \
  >  --sequences ../data/invalidnucleotide.fasta \
  >  --mask-invalid \
  >  --output "$TMP/masked.fasta"
  Removing masked sites from FASTA file.

  $ cat "$TMP/masked.fasta"
  >sequence_1
  ATCGNNNN
  $ rm -f "$TMP/masked.fasta"

  $ popd > /dev/null
