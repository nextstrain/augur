Setup

  $ source "$TESTDIR"/_setup.sh

  $ export DATA="$TESTDIR/../data/simple-genome"

We take the same `snps.vcf` file used in `general.t` but add another
allele at site 30 - sample_B has a "G". Since the root is A and this
is the only sample with G it's 'A30G'.
See <https://github.com/nextstrain/augur/issues/1380> for the bug this is testing.

  $ sed '11s/^/1\t30\t.\tA\tG,N\t.\t.\t.\tGT\t0\t1\t2\n/' \
  > "$DATA/snps.vcf" > snps2.vcf

  $ ${AUGUR} ancestral \
  >  --tree $DATA/tree.nwk \
  >  --alignment snps2.vcf \
  >  --vcf-reference $DATA/reference.fasta \
  >  --output-node-data nt_muts.json \
  >  --output-vcf nt_muts.vcf \
  >  --inference marginal > /dev/null


  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" \
  >   "$DATA/nt_muts.ref-seq.json" \
  >   nt_muts.json \
  >   --exclude-regex-paths "root\['nodes'\]\['.+'\]\['sequence'\]"
  {'iterable_item_added': {"root['nodes']['sample_B']['muts'][0]": 'A30G'}}

  $ cat > expected.vcf <<EOF
  > #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	node_root	sample_C	node_AB	sample_B	sample_A
  > 1	5	.	A	C	.	PASS	.	GT	1	1	1	1	1
  > 1	7	.	A	G	.	PASS	.	GT	0	0	1	1	1
  > 1	14	.	C	T	.	PASS	.	GT	0	0	1	1	1
  > 1	18	.	C	T	.	PASS	.	GT	1	1	0	0	0
  > 1	30	.	A	G	.	PASS	.	GT	0	0	0	1	0
  > 1	33	.	A	C,G	.	PASS	.	GT	0	0	2	2	1
  > 1	39	.	C	T	.	PASS	.	GT	0	0	0	0	1
  > 1	42	.	G	A	.	PASS	.	GT	0	0	0	1	0
  > 1	43	.	A	T	.	PASS	.	GT	0	0	1	1	1
  > EOF

Ignore the meta-lines of the output VCF
  $ grep -v '##' nt_muts.vcf > nt_muts.no-header.vcf

  $ diff expected.vcf nt_muts.no-header.vcf
