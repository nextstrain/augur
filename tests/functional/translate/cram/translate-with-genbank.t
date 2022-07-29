Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Translate amino acids for genes using a GenBank file.

  $ ${AUGUR} translate \
  >   --tree translate/data/zika/tree.nwk \
  >   --ancestral-sequences translate/data/zika/nt_muts.json \
  >   --reference-sequence translate/data/zika/zika_outgroup.gb \
  >   --genes CA PRO \
  >   --output-node-data $TMP/aa_muts.json
  Validating schema of 'translate/data/zika/nt_muts.json'...
  Read in 3 features from reference sequence file
  amino acid mutations written to .* (re)
  $ python3 "../../scripts/diff_jsons.py" translate/data/zika/aa_muts_genbank.json $TMP/aa_muts.json
  {}
