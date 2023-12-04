Setup

  $ source "$TESTDIR"/_setup.sh

  $ export DATA="$TESTDIR/../data/simple-genome"

This command mirrors the first test in general.t, however
with VCF input instead of a FASTA MSA.
The output will not have the full sequence attached to every node,
but it will have the reference sequence attached.

  $ ${AUGUR} ancestral \
  >  --tree $DATA/tree.nwk \
  >  --alignment $DATA/snps.vcf \
  >  --vcf-reference $DATA/reference.fasta \
  >  --output-node-data "nt_muts.vcf-input.ref-seq.json" \
  >  --output-vcf "nt_muts.vcf-input.ref-seq.vcf" \
  >  --inference marginal > /dev/null


  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" \
  >   "$DATA/nt_muts.ref-seq.json" \
  >   "nt_muts.vcf-input.ref-seq.json" \
  >   --exclude-regex-paths "root\['nodes'\]\['.+'\]\['sequence'\]"
  {}

