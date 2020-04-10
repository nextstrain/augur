Integration tests for augur mask.

Try masking a VCF without any specified mask.

  $ augur mask --sequences mask/variants.vcf
  usage: augur mask [-h] --sequences SEQUENCES --mask MASK [--output OUTPUT]
                    [--no-cleanup]
  augur mask: error: the following arguments are required: --mask
  [2]

Mask a VCF with a BED file and no specified output file.

  $ cp "$TESTDIR/mask/variants.vcf" "$TMP/"
  $ augur mask \
  >  --sequences "$TMP/variants.vcf" \
  >  --mask "$TESTDIR/mask/mask_variants.bed" > /dev/null

  $ diff -u "$TESTDIR/mask/masked_variants.vcf" "$TMP/variants.vcf"
  $ rm -f "$TMP/variants.vcf"

Mask a VCF with a BED file and a specified output file.

  $ augur mask \
  >  --sequences "$TESTDIR/mask/variants.vcf" \
  >  --mask "$TESTDIR/mask/mask_variants.bed" \
  >  --output "$TMP/masked_variants.vcf" > /dev/null

  $ diff -u "$TESTDIR/mask/masked_variants.vcf" "$TMP/masked_variants.vcf"
  $ rm -f "$TMP/masked_variants.vcf"

Try masking sequences without any specified mask.

  $ augur mask --sequences mask/sequences.fasta
  usage: augur mask [-h] --sequences SEQUENCES --mask MASK [--output OUTPUT]
                    [--no-cleanup]
  augur mask: error: the following arguments are required: --mask
  [2]

Mask sequences with a BED file and no specified output file.
Since no output is provided, the input file is overridden with the masked sequences.

  $ cp $TESTDIR/mask/sequences.fasta $TMP/
  $ augur mask --sequences $TMP/sequences.fasta --mask $TESTDIR/mask/mask.bed
  Found 5 sites to mask
  Removing masked sites from FASTA file.

  $ cat $TMP/sequences.fasta
  >sequence_1
  NNNCNN
  $ rm -f $TMP/sequences.fasta

Mask sequences with a BED file and a specified output file.

  $ augur mask \
  >  --sequences $TESTDIR/mask/sequences.fasta \
  >  --mask $TESTDIR/mask/mask.bed \
  >  --output $TMP/masked.fasta
  Found 5 sites to mask
  Removing masked sites from FASTA file.

  $ cat $TMP/masked.fasta
  >sequence_1
  NNNCNN
  $ rm -f $TMP/masked.fasta
