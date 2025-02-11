Setup

  $ source "$TESTDIR"/_setup.sh

Attempt to root on a strain which is not in the tree is an error
  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" --metadata "$TESTDIR/../data/metadata.tsv" --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --output-tree tree.nwk --output-node-data branch_lengths.json \
  >  --root not-in-tree > /dev/null
  ERROR: Rooting with the provided strain name(s) failed.
  This is probably because the following names supplied to '--root' were not found in the tree:
    - 'not-in-tree'
  
  [2]

Check that rooting on multiple (valid) sequences + '--remove-outgroup' is an error
  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" --metadata "$TESTDIR/../data/metadata.tsv" --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --output-tree tree.nwk --output-node-data branch_lengths.json \
  >  --root 'Colombia/2016/ZC204Se' 'VEN/UF_1/2016' --remove-outgroup > /dev/null
  ERROR: --remove-outgroup is not valid with 'multiple strains' rooting - you must supply a single strain (outgroup) to root the tree
  [2]
