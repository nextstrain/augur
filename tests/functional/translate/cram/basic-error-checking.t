Setup

  $ source "$TESTDIR"/_setup.sh
  $ export ANC_DATA="$TESTDIR/../../ancestral/data/simple-genome"
  $ export DATA="$TESTDIR/../data/simple-genome"

Missing reference file

  $ ${AUGUR} translate \
  >  --tree $ANC_DATA/tree.nwk \
  >  --ancestral-sequences $ANC_DATA/nt_muts.ref-seq.json \
  >  --reference-sequence $DATA/reference.doesnt-exist.gff \
  >  --output-node-data "aa_muts.json" > /dev/null
  ERROR: reference sequence file '.+/reference.doesnt-exist.gff' not found (re)
  [2]

Gene length not divisible by 3

  $ cat >invalid.gff <<~~
  > ##gff-version 3
  > ##sequence-region reference_name 1 50
  > reference_name	RefSeq	region	1	50	.	+	.	ID=reference_name
  > reference_name	RefSeq	gene	10	23	.	+	.	Name=gene1;gene=gene1
  > reference_name	RefSeq	gene	36	48	.	-	.	Name=gene2;gene=gene2
  > ~~

  $ ${AUGUR} translate \
  >  --tree $ANC_DATA/tree.nwk \
  >  --ancestral-sequences $ANC_DATA/nt_muts.ref-seq.json \
  >  --reference-sequence invalid.gff \
  >  --output-node-data "aa_muts.json"
  ERROR: Reference file 'invalid.gff' has errors:
      Gene length of 'gene1' is not a multiple of 3.
      Gene length of 'gene2' is not a multiple of 3.
  [2]

Gene with compound location not divisible by 3

  $ cat >invalid.gb <<~~
  > LOCUS       Test                   100 bp    DNA     linear   VRL 01-JAN-2025
  > DEFINITION  Test genome for compound location validation
  > FEATURES             Location/Qualifiers
  >      source          1..100
  >      CDS             join(10..20,35..49)
  >                      /gene="test_gene"
  > ORIGIN
  >         1 aaaaaaaaaa tgccctgcgg gtaaaaaaaa aaaaactact tgaccataaa
  >        51 aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa
  > //
  > ~~

  $ ${AUGUR} translate \
  >  --tree $ANC_DATA/tree.nwk \
  >  --ancestral-sequences $ANC_DATA/nt_muts.ref-seq.json \
  >  --reference-sequence invalid.gb \
  >  --output-node-data "aa_muts.json"
  .* BiopythonParserWarning: Attempting to parse malformed locus line: (re)
  .* (re)
  .* (re)
  .* (re)
  .* (re)
  ERROR: Reference file 'invalid.gb' has errors:
      Gene length of 'test_gene' is not a multiple of 3.
  [2]
