Setup

  $ source "$TESTDIR"/_setup.sh
  $ export DATA="$TESTDIR/../data"

Translate amino acids for genes using a GFF3 file where the gene names are stored in a qualifier named "gene_name".

  $ cat >genemap.gff <<~~
  > ##gff-version 3
  > ##sequence-region PF13/251013_18 1 10769
  > PF13/251013_18	GenBank	gene	91	456	.	+	.	gene_name="CA"
  > PF13/251013_18	GenBank	gene	457	735	.	+	.	gene_name="PRO"
  > ~~

  $ ${AUGUR} translate \
  >   --tree "${DATA}/zika/tree.nwk" \
  >   --ancestral-sequences "${DATA}/zika/nt_muts.json" \
  >   --reference-sequence "genemap.gff" \
  >   --output-node-data aa_muts.json
  Read in 3 features from reference sequence file
  Validating schema of '.+/nt_muts.json'... (re)
  Validating schema of .* (re)
  amino acid mutations written to .* (re)

Other than the sequence ids which will include a temporary path, the JSONs
should be identical.

  $ python3 "${SCRIPTS}/diff_jsons.py" \
  >  --exclude-regex-paths "['seqid']" -- \
  >  "${DATA}/zika/aa_muts_gff.json" \
  >  aa_muts.json
  {}
