Setup

  $ source "$TESTDIR"/_setup.sh

Run with a config file.

  $ cat > config.yaml <<~~
  > tree: "$TESTDIR/../data/simple-genome/tree.nwk"
  > alignment: "$TESTDIR/../data/simple-genome/sequences.fasta"
  > root_sequence: "$TESTDIR/../data/simple-genome/reference.fasta"
  > output_node_data: "nt_muts.ref-seq.json"
  > seed: 314159
  > inference: marginal
  > ~~

  $ ${AUGUR} ancestral --config config.yaml > /dev/null
  Validating schema of 'nt_muts.ref-seq.json'...

Confirm that node data matches expected results.

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" \
  >   "$TESTDIR/../data/simple-genome/nt_muts.ref-seq.json" \
  >   "nt_muts.ref-seq.json" \
  >   --exclude-paths "root['generated_by']"
  {}
