Setup

  $ source "$TESTDIR"/_setup.sh
  $ export ANC_DATA="$TESTDIR/../../ancestral/data/simple-genome"
  $ export DATA="$TESTDIR/../data/simple-genome"

Run with a config file.

  $ cat > config.yaml <<~~
  > tree: "$ANC_DATA/tree.nwk"
  > ancestral_sequences: "$ANC_DATA/nt_muts.ref-seq.json"
  > reference_sequence: "$DATA/reference.gff"
  > output_node_data: "aa_muts.json"
  > ~~

  $ ${AUGUR} translate --config config.yaml &> /dev/null

Confirm that node data matches expected results.

  $ python3 "$SCRIPTS/diff_jsons.py" \
  >   "$DATA/aa_muts.json" \
  >   "aa_muts.json" \
  >   --exclude-regex-paths "root\['annotations'\]\['.+'\]\['seqid'\]" "root['meta']['updated']"
  {}
