Integration tests for augur mask.

  $ pushd "$TESTDIR" > /dev/null

Try masking a VCF without any specified mask.

  $ augur mask --sequences mask/variants.vcf
  No masking sites provided. Must include one of --mask, --mask-from-beginning, --mask-from-end, or --mask-sites
  [1]

Mask a VCF with a BED file and no specified output file.

  $ cp "mask/variants.vcf" "$TMP/"
  $ augur mask \
  >  --sequences "$TMP/variants.vcf" \
  >  --mask "mask/mask_variants.bed" > /dev/null

  $ diff -u "mask/masked_variants.vcf" "$TMP/variants.vcf"
  $ rm -f "$TMP/variants.vcf"

Mask a VCF with a BED file and a specified output file.

  $ augur mask \
  >  --sequences "mask/variants.vcf" \
  >  --mask "mask/mask_variants.bed" \
  >  --output "$TMP/masked_variants.vcf" > /dev/null

  $ diff -u "mask/masked_variants.vcf" "$TMP/masked_variants.vcf"
  $ rm -f "$TMP/masked_variants.vcf"

Try masking sequences without any specified mask.

  $ augur mask --sequences mask/sequences.fasta
  No masking sites provided. Must include one of --mask, --mask-from-beginning, --mask-from-end, or --mask-sites
  [1]

Mask sequences with a BED file and no specified output file.
Since no output is provided, the input file is overridden with the masked sequences.

  $ cp mask/sequences.fasta "$TMP/"
  $ augur mask --sequences "$TMP/sequences.fasta" --mask mask/mask.bed
  Found 3 sites to mask in 'mask/mask.bed'
  Removing masked sites from FASTA file.

  $ cat "$TMP/sequences.fasta"
  >sequence_1
  NNGCNG
  $ rm -f "$TMP/sequences.fasta"

Mask sequences with a BED file and a specified output file.

  $ augur mask \
  >  --sequences mask/sequences.fasta \
  >  --mask mask/mask.bed \
  >  --output "$TMP/masked.fasta"
  Found 3 sites to mask in 'mask/mask.bed'
  Removing masked sites from FASTA file.

  $ cat "$TMP/masked.fasta"
  >sequence_1
  NNGCNG
  $ rm -f "$TMP/masked.fasta"

Mask one base from the beginning and the end.

  $ augur mask \
  >  --sequences mask/sequences.fasta \
  >  --mask-from-beginning 1 \
  >  --mask-from-end 1 \
  >  --output "$TMP/masked.fasta"
  Removing masked sites from FASTA file.

  $ cat "$TMP/masked.fasta"
  >sequence_1
  NTGCTN
  $ rm -f "$TMP/masked.fasta"

Mask a specific list of sites and also mask one base from the beginning and the end.

  $ augur mask \
  >  --sequences mask/sequences.fasta \
  >  --mask-sites 3 4 \
  >  --mask-from-beginning 1 \
  >  --mask-from-end 1 \
  >  --output "$TMP/masked.fasta"
  Removing masked sites from FASTA file.

  $ cat "$TMP/masked.fasta"
  >sequence_1
  NTNNTN
  $ rm -f "$TMP/masked.fasta"

  $ popd > /dev/null
