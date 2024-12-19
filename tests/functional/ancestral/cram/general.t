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
  >  --seed 314159 \
  >  --inference marginal > /dev/null


  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" \
  >   "$TESTDIR/../data/simple-genome/nt_muts.ref-seq.json" \
  >   "nt_muts.ref-seq.json" \
  >   --exclude-paths "root['generated_by']"
  {}

Same as above but without providing a `--root-sequence`. The effect of this on behaviour is:
- The JSON['reference']['nuc'] is not actually the "reference", because we
don't know what this is! Instead it's the inferred root-node sequence.
See <https://github.com/nextstrain/augur/issues/1362> for more.
- An array of 'muts' is still present on the root node, but it will never contain any
mutations (as there's nothing to compare the root node to)

  $ ${AUGUR} ancestral \
  >  --tree "$TESTDIR/../data/simple-genome/tree.nwk" \
  >  --alignment "$TESTDIR/../data/simple-genome/sequences.fasta" \
  >  --output-node-data "nt_muts.no-ref-seq.json" \
  >  --seed 314159 \
  >  --inference marginal > /dev/null


  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" \
  >   "$TESTDIR/../data/simple-genome/nt_muts.no-ref-seq.json" \
  >   "nt_muts.no-ref-seq.json" \
  >   --exclude-paths "root['generated_by']"
  {}
