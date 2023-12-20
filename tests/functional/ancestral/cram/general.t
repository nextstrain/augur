Setup

  $ source "$TESTDIR"/_setup.sh

General running of augur ancestral.
- Note that the (parsimonious) C18T on the branch of Sample_C is inferred by augur
as a root mutation of C18T + T18C reversion on Node_AB. This is reflected in the
node-data JSON we diff against.
- We supply the refererence sequence so mutations are called at the root.
- Ambiguous nucleotides (there are 3 Ns in sample_C) are inferred (default)

  $ ${AUGUR} ancestral \
  >  --tree "$TESTDIR/../data/simple-genome/tree.nwk" \
  >  --alignment "$TESTDIR/../data/simple-genome/sequences.fasta" \
  >  --root-sequence "$TESTDIR/../data/simple-genome/reference.fasta" \
  >  --output-node-data "nt_muts.ref-seq.json" \
  >  --inference marginal > /dev/null


  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" \
  >   "$TESTDIR/../data/simple-genome/nt_muts.ref-seq.json" \
  >   "nt_muts.ref-seq.json"
  {}

Same as above but without a root-sequence so no mutations inferred on the root node
(and thus the inferred reference will be different)

  $ ${AUGUR} ancestral \
  >  --tree "$TESTDIR/../data/simple-genome/tree.nwk" \
  >  --alignment "$TESTDIR/../data/simple-genome/sequences.fasta" \
  >  --output-node-data "nt_muts.no-ref-seq.json" \
  >  --inference marginal > /dev/null


  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" \
  >   "$TESTDIR/../data/simple-genome/nt_muts.ref-seq.json" \
  >   "nt_muts.no-ref-seq.json" \
  >   --exclude-paths "root['reference']['nuc']" "root['nodes']['node_root']['muts']"
  {}
