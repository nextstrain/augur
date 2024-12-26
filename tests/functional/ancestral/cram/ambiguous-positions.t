Setup

  $ source "$TESTDIR"/_setup.sh
  $ export DATA="$TESTDIR/../data/simple-genome"
  $ export SCRIPTS="$TESTDIR/../../../../scripts/"

  $ cat > reference.fasta <<EOF
  > >ref_name
  > ATGCNNNNATGC
  > EOF

  $ cat > snps.vcf <<EOF
  > #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample_A	sample_B	sample_C
  > 1	3	.	G	N	.	PASS	.	GT	1	1	1
  > 1	9	.	A	T	.	PASS	.	GT	1	1	0
  > EOF

  $ ${AUGUR} ancestral \
  >  --tree "$DATA/tree.nwk" \
  >  --alignment snps.vcf \
  >  --vcf-reference reference.fasta \
  >  --seed 314159 \
  >  --output-node-data nt_muts.json \
  >  --output-vcf nt_muts.vcf > /dev/null

  $ cat > expected.json <<EOF
  > {"nodes": {
  >   "node_root": {"muts": []},
  >   "node_AB": {"muts": ["A9T"]},
  >   "sample_A": {"muts": []},
  >   "sample_B": {"muts": []},
  >   "sample_C": {"muts": []}
  > },
  > "mask": "001011110000"
  > }
  > EOF

  $ python3 "$SCRIPTS/diff_jsons.py" \
  >   expected.json "nt_muts.json" \
  >   --exclude-regex-paths "root\['annotations'\]" "root\['generated_by'\]" "root\['reference'\]"
  {}

  $ python3 "$SCRIPTS/compare-json-vcf.py" \
  >   --tree "$DATA/tree.nwk" \
  >   --vcf nt_muts.vcf \
  >   --json nt_muts.json \
  >   --ref reference.fasta > /dev/null
