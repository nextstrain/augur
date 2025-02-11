Setup

  $ source "$TESTDIR"/_setup.sh

Check that rooting on an outgroup strain actually works
  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" --metadata "$TESTDIR/../data/metadata.tsv" --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --output-tree tree1.nwk --output-node-data branch_lengths1.json \
  >  --root 'Colombia/2016/ZC204Se' > /dev/null

  $ python3 "$TESTDIR/../report-root.py" tree1.nwk
  Tree root has a single terminal child 'Colombia/2016/ZC204Se'

Check that rooting on an outgroup strain with '--remove-outgroup' does remove it 
  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" --metadata "$TESTDIR/../data/metadata.tsv" --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --output-tree tree2.nwk --output-node-data branch_lengths2.json \
  >  --root 'Colombia/2016/ZC204Se' --remove-outgroup > /dev/null

  $ grep 'Colombia/2016/ZC204Se' tree2.nwk
  [1]

  $ python3 "$TESTDIR/../report-root.py" tree2.nwk
  No children of the root were terminal nodes
