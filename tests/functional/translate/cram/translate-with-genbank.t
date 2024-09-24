Setup

  $ source "$TESTDIR"/_setup.sh
  $ export DATA="$TESTDIR/../data"

Translate amino acids for genes using a GenBank file.

  $ ${AUGUR} translate \
  >   --tree "$DATA/zika/tree.nwk" \
  >   --ancestral-sequences "$DATA/zika/nt_muts.json" \
  >   --reference-sequence "$DATA/zika/zika_outgroup.gb" \
  >   --genes CA PRO \
  >   --output-node-data aa_muts.json
  WARNING: 1 CDS features skipped as they didn't have a locus_tag or gene qualifier.
  Read in 3 features from reference sequence file
  Validating schema of '.+nt_muts.json'... (re)
  Validating schema of .* (re)
  amino acid mutations written to .* (re)

  $ python3 "$SCRIPTS/diff_jsons.py" $DATA/zika/aa_muts_genbank.json aa_muts.json \
  >  --exclude-regex-paths "root\['annotations'\]\['.+'\]\['seqid'\]" "root['meta']['updated']"
  {}
