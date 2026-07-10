Setup

  $ source "$TESTDIR"/_setup.sh

Use a YAML config file for selected refine options.

  $ cat >refine.yaml <<~~
  > root:
  >   - Colombia/2016/ZC204Se
  > covariance: false
  > stochastic_resolve: false
  > max_iter: 5
  > ~~

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --output-tree tree.nwk \
  >  --output-node-data branch_lengths.json \
  >  --config refine.yaml > /dev/null

  $ python3 "$TESTDIR/../report-root.py" tree.nwk
  Tree root has a single terminal child 'Colombia/2016/ZC204Se'

Config options must not also be provided on the command line.

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" \
  >  --output-tree tree.nwk \
  >  --output-node-data branch_lengths.json \
  >  --config refine.yaml \
  >  --max-iter 2
  ERROR: --config cannot be used with these CLI options: --max-iter
  [2]

Config-backed CLI options are allowed when they do not also appear in the config file.

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" \
  >  --output-tree tree.nwk \
  >  --output-node-data branch_lengths.json \
  >  --config refine.yaml \
  >  --coalescent opt > /dev/null

Unknown config keys are rejected.

  $ cat >unknown.yaml <<~~
  > unknown: true
  > ~~

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" \
  >  --output-tree tree.nwk \
  >  --output-node-data branch_lengths.json \
  >  --config unknown.yaml
  ERROR: Invalid refine config 'unknown.yaml': unexpected config option 'unknown'
  [2]
