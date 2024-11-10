Setup

  $ source "$TESTDIR"/_setup.sh
  $ export DATA="$TESTDIR/../data"

Translate amino acids for genes using a GFF3 file where the gene names are stored in a qualifier named "gene".

  $ cat >genemap.gff <<~~
  > ##gff-version 3
  > ##sequence-region PF13/251013_18 1 10769
  > PF13/251013_18	GenBank	gene	91	456	.	+	.	gene="CA"
  > PF13/251013_18	GenBank	gene	457	735	.	+	.	gene="PRO"
  > ~~

  $ ${AUGUR} translate \
  >   --tree "${DATA}/zika/tree.nwk" \
  >   --ancestral-sequences "${DATA}/zika/nt_muts.json" \
  >   --reference-sequence genemap.gff \
  >   --output-node-data aa_muts.json
  Read in 3 features from reference sequence file
  Validating schema of '.+/nt_muts.json'... (re)
  Validating schema of .* (re)
  amino acid mutations written to .* (re)

  $ python3 "${SCRIPTS}/diff_jsons.py" \
  >  --exclude-regex-paths "['seqid']" -- \
  >  "${DATA}/zika/aa_muts_gff.json" \
  >  aa_muts.json
  {}
