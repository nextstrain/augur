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
  >  --output-vcf nt_muts.vcf \
  >  --seed 314159 \
  >  --inference marginal > /dev/null


  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" \
  >   "$DATA/nt_muts.ref-seq.json" \
  >   "nt_muts.vcf-input.ref-seq.json" \
  >   --exclude-regex-paths "root\['nodes'\]\['.+'\]\['sequence'\]" "root\['generated_by'\]"
  {}

Here's the same mutations as in $DATA/nt_muts.ref-seq.json,
but as a VCF file with (a) strain name ordering and (b) ALT
ordering to match the way `augur ancestral` writes it.
`augur auncestral` also writes FILTER=PASS, which I include
here as it's not relevant to what I'm trying to test.
  $ cat > expected.vcf <<EOF
  > #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	node_root	sample_C	node_AB	sample_B	sample_A
  > 1	5	.	A	C	.	PASS	.	GT	1	1	1	1	1
  > 1	7	.	A	G	.	PASS	.	GT	0	0	1	1	1
  > 1	14	.	C	T	.	PASS	.	GT	1	0	1	1	1
  > 1	18	.	C	T	.	PASS	.	GT	1	1	0	0	0
  > 1	33	.	A	C,G	.	PASS	.	GT	0	0	2	2	1
  > 1	39	.	C	T	.	PASS	.	GT	0	0	0	0	1
  > 1	42	.	G	A	.	PASS	.	GT	0	0	0	1	0
  > 1	43	.	A	T	.	PASS	.	GT	0	0	1	1	1
  > EOF

Ignore the meta-lines of the output VCF
  $ grep -v '##' nt_muts.vcf > nt_muts.no-header.vcf

  $ diff expected.vcf nt_muts.no-header.vcf
