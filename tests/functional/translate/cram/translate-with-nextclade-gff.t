Setup

  $ source "$TESTDIR"/_setup.sh
  $ export DATA="$TESTDIR/../data"

Translate amino acids using --use-nextclade-gff-style. CA is a gene without any
CDS so it is used as a CDS itself (named via its 'gene' attribute). PRO is a CDS
split over two rows sharing an ID which are joined into a single feature (named
via its 'Name' attribute, not the gene's). The translations are therefore the
same as for a GFF with two simple genes.

  $ cat >genemap.gff <<~~
  > ##gff-version 3
  > ##sequence-region PF13/251013_18 1 10769
  > PF13/251013_18	GenBank	gene	91	456	.	+	.	ID=gene-CA;gene=CA
  > PF13/251013_18	GenBank	gene	457	735	.	+	.	ID=gene-PRO;gene=PRO-gene
  > PF13/251013_18	GenBank	CDS	457	600	.	+	0	ID=cds-PRO;Parent=gene-PRO;Name=PRO;gene=PRO-gene
  > PF13/251013_18	GenBank	CDS	601	735	.	+	0	ID=cds-PRO;Parent=gene-PRO;Name=PRO;gene=PRO-gene
  > ~~

  $ ${AUGUR} translate \
  >   --tree "${DATA}/zika/tree.nwk" \
  >   --ancestral-sequences "${DATA}/zika/nt_muts.json" \
  >   --reference-sequence genemap.gff \
  >   --use-nextclade-gff-style \
  >   --output-node-data aa_muts.json
  Read in 3 features from reference sequence file
  Validating schema of '.+/nt_muts.json'... (re)
  Validating schema of .* (re)
  amino acid mutations written to .* (re)

  $ python3 "${SCRIPTS}/diff_jsons.py" \
  >  --exclude-regex-paths "['seqid']" "root\['annotations'\]\['PRO'\]" -- \
  >  "${DATA}/zika/aa_muts_gff.json" \
  >  aa_muts.json
  {}

  $ python3 -c 'import json; print(json.load(open("aa_muts.json"))["annotations"]["PRO"]["segments"])'
  [{'end': 600, 'start': 457}, {'end': 735, 'start': 601}]

Without --use-nextclade-gff-style only the genes are read.

  $ ${AUGUR} translate \
  >   --tree "${DATA}/zika/tree.nwk" \
  >   --ancestral-sequences "${DATA}/zika/nt_muts.json" \
  >   --reference-sequence genemap.gff \
  >   --output-node-data aa_muts_default.json > /dev/null 2>&1

  $ python3 -c 'import json; print(sorted(json.load(open("aa_muts_default.json"))["annotations"]))'
  ['CA', 'PRO-gene', 'nuc']

Requesting a subset of features via --genes.

  $ ${AUGUR} translate \
  >   --tree "${DATA}/zika/tree.nwk" \
  >   --ancestral-sequences "${DATA}/zika/nt_muts.json" \
  >   --reference-sequence genemap.gff \
  >   --use-nextclade-gff-style \
  >   --genes PRO \
  >   --output-node-data aa_muts_pro.json > /dev/null 2>&1

  $ python3 -c 'import json; print(sorted(json.load(open("aa_muts_pro.json"))["annotations"]))'
  ['PRO', 'nuc']

Feature names must be unique.

  $ cat >duplicate.gff <<~~
  > ##gff-version 3
  > ##sequence-region PF13/251013_18 1 10769
  > PF13/251013_18	GenBank	CDS	91	456	.	+	0	ID=cds-1;Name=CA
  > PF13/251013_18	GenBank	CDS	457	735	.	+	0	ID=cds-2;Name=CA
  > ~~

  $ ${AUGUR} translate \
  >   --tree "${DATA}/zika/tree.nwk" \
  >   --ancestral-sequences "${DATA}/zika/nt_muts.json" \
  >   --reference-sequence duplicate.gff \
  >   --use-nextclade-gff-style \
  >   --output-node-data aa_muts.json
  ERROR: Reference 'duplicate.gff' contains multiple genes/CDSs with the name 'CA'. Names must be unique.
  [2]

Multi-segment CDSs are not supported for VCF input.

  $ export ANC_DATA="$TESTDIR/../../ancestral/data/simple-genome"
  $ cat >simple.gff <<~~
  > ##gff-version 3
  > ##sequence-region reference_name 1 50
  > reference_name	RefSeq	CDS	10	15	.	+	0	ID=cds-1;Name=gene1
  > reference_name	RefSeq	CDS	16	24	.	+	0	ID=cds-1;Name=gene1
  > ~~

  $ ${AUGUR} translate \
  >  --tree "$ANC_DATA/tree.nwk" \
  >  --ancestral-sequences "$DATA/simple-genome/snps-inferred.vcf" \
  >  --reference-sequence simple.gff \
  >  --use-nextclade-gff-style \
  >  --output-node-data aa_muts.json \
  >  --vcf-reference "$ANC_DATA/reference.fasta"
  Read in 2 features from reference sequence file
  ERROR: 'gene1' consists of multiple segments, which is not supported for VCF input.
  [2]
