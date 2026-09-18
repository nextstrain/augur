Integration tests for augur distance.

  $ source "$TESTDIR"/_setup.sh

Build a tree to use for distance calculations.

  $ ${AUGUR} tree \
  >   --alignment "$TESTDIR/distance/aligned.fasta" \
  >   --output "tree_raw.nwk" &> /dev/null

Refine tree to get internal node names.

  $ ${AUGUR} refine \
  >   --tree "tree_raw.nwk" \
  >   --output-tree "tree.nwk" &> /dev/null

Calculate pairwise distances and write them to an edge list.

  $ ${AUGUR} distance \
  >   --tree "tree.nwk" \
  >   --alignment "$TESTDIR/distance/aligned.fasta" \
  >   --gene-names nuc \
  >   --attribute-name snvs \
  >   --compare-to pairwise \
  >   --map "$TESTDIR/distance/distance_map_hamming.json" \
  >   --output distances.json \
  >   --output-edge-list /dev/stdout
  sequence_1	sequence_2	distance
  with_gaps	with_gaps	0
  with_gaps	some_other_seq	1
  with_gaps	no_gaps	0
  with_gaps	_R_crick_strand	0
  some_other_seq	with_gaps	1
  some_other_seq	some_other_seq	0
  some_other_seq	no_gaps	1
  some_other_seq	_R_crick_strand	1
  no_gaps	with_gaps	0
  no_gaps	some_other_seq	1
  no_gaps	no_gaps	0
  no_gaps	_R_crick_strand	0
  _R_crick_strand	with_gaps	0
  _R_crick_strand	some_other_seq	1
  _R_crick_strand	no_gaps	0
  _R_crick_strand	_R_crick_strand	0

Try to calculate distance from the root with edge list output.
This should fail because we only support edge lists for pairwise comparisons.

  $ ${AUGUR} distance \
  >   --tree "tree.nwk" \
  >   --alignment "$TESTDIR/distance/aligned.fasta" \
  >   --gene-names nuc \
  >   --attribute-name snvs \
  >   --compare-to root \
  >   --map "$TESTDIR/distance/distance_map_hamming.json" \
  >   --output distances.json \
  >   --output-edge-list distances_edge_list.tsv
  ERROR: Edge list output only works for a single pairwise comparison.
  [2]

Try to calculate distances with both root and pairwise comparisons and using edge list output.
This should fail because we only support edge list output for a single pairwise comparison.

  $ ${AUGUR} distance \
  >   --tree "tree.nwk" \
  >   --alignment "$TESTDIR/distance/aligned.fasta" \
  >   --gene-names nuc \
  >   --attribute-name snvs snvs_pairwise \
  >   --compare-to root pairwise \
  >   --map "$TESTDIR/distance/distance_map_hamming.json" \
  >   --output distances.json \
  >   --output-edge-list distances_edge_list.tsv
  ERROR: Edge list output only works for a single pairwise comparison.
  [2]

Try to calculate distances from the root with an alignment that only contains tips.
This should fail with a clear error instead of silently reporting zero distances.

  $ ${AUGUR} distance \
  >   --tree "tree.nwk" \
  >   --alignment "$TESTDIR/distance/aligned.fasta" \
  >   --gene-names nuc \
  >   --attribute-name snvs \
  >   --compare-to root \
  >   --map "$TESTDIR/distance/distance_map_hamming.json" \
  >   --output distances.json
  ERROR: Could not find a sequence for the root node 'NODE_0000000' in the given alignment(s). Comparisons to the root require the alignment to include sequences for internal nodes (e.g., ancestral sequences from `augur ancestral`).
  [2]
