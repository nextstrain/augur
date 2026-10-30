Setup

  $ source "$TESTDIR"/_setup.sh

'KX369547.1', 'PRVABC59' and 'ZKC2/2016' are removed with --clock-filter-iqd 1.

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --output-tree tree.nwk \
  >  --output-node-data branch_lengths.json \
  >  --timetree \
  >  --clock-filter-iqd 1 \
  >  --seed 314159 2>&1 | grep "pruning leaf" || echo "Nothing pruned"
  pruning leaf  PRVABC59
  pruning leaf  KX369547.1
  pruning leaf  ZKC2/2016

  $ grep -q -F 'KX369547.1' tree.nwk && echo 'Present' || echo 'Pruned'
  Pruned

Use --keep-ids to force-include 'KX369547.1'.

  $ cat > include.txt <<~~
  > KX369547.1  # Keep me
  > ~~

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --output-tree tree.nwk \
  >  --output-node-data branch_lengths.json \
  >  --timetree \
  >  --clock-filter-iqd 1 \
  >  --keep-ids include.txt \
  >  --seed 314159 2>&1 | grep "pruning leaf" || echo "Nothing pruned"
  pruning leaf  PRVABC59
  pruning leaf  ZKC2/2016

  $ grep -q -F 'KX369547.1' tree.nwk && echo 'Present' || echo 'Pruned'
  Present

Multiple files can be given; ids from all of them are kept.

  $ cat > also_include.txt <<~~
  > PRVABC59
  > ZKC2/2016
  > ~~

'also_include.txt' alone does not keep 'KX369547.1'.

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --output-tree tree.nwk \
  >  --output-node-data branch_lengths.json \
  >  --timetree \
  >  --clock-filter-iqd 1 \
  >  --keep-ids also_include.txt \
  >  --seed 314159 2>&1 | grep "pruning leaf" || echo "Nothing pruned"
  pruning leaf  KX369547.1

Together, 'also_include.txt' and 'include.txt' keep all three.

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --output-tree tree.nwk \
  >  --output-node-data branch_lengths.json \
  >  --timetree \
  >  --clock-filter-iqd 1 \
  >  --keep-ids also_include.txt include.txt \
  >  --seed 314159 2>&1 | grep "pruning leaf" || echo "Nothing pruned"
  Nothing pruned

  $ for id in KX369547.1 PRVABC59 ZKC2/2016; do grep -q -F "$id" tree.nwk && echo "$id Present" || echo "$id Pruned"; done
  KX369547.1 Present
  PRVABC59 Present
  ZKC2/2016 Present
