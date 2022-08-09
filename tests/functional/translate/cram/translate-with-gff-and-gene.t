Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Translate amino acids for genes using a GFF3 file where the gene names are stored in a qualifier named "gene".

  $ cat >$TMP/genemap.gff <<~~
  > ##gff-version 3
  > ##sequence-region PF13/251013_18 1 10769
  > PF13/251013_18	GenBank	gene	91	456	.	+	.	gene="CA"
  > PF13/251013_18	GenBank	gene	457	735	.	+	.	gene="PRO"
  > ~~

  $ ${AUGUR} translate \
  >   --tree translate/data/zika/tree.nwk \
  >   --ancestral-sequences translate/data/zika/nt_muts.json \
  >   --reference-sequence "$TMP/genemap.gff" \
  >   --output-node-data $TMP/aa_muts.json
  Validating schema of 'translate/data/zika/nt_muts.json'...
  Read in 2 features from reference sequence file
  amino acid mutations written to .* (re)
  $ python3 "../../scripts/diff_jsons.py" \
  >  --exclude-regex-paths "['seqid']" -- \
  >  translate/data/zika/aa_muts_gff.json \
  >  $TMP/aa_muts.json
  {}
