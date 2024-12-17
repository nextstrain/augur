Setup

  $ source "$TESTDIR"/_setup.sh

Change the _reference_ to lowercase

  $ sed '2y/ATCGN/atcgn/' < "$TESTDIR/../data/simple-genome/reference.fasta" > ref.lower.fasta

  $ ${AUGUR} ancestral \
  >  --tree "$TESTDIR/../data/simple-genome/tree.nwk" \
  >  --alignment "$TESTDIR/../data/simple-genome/sequences.fasta" \
  >  --root-sequence ref.lower.fasta \
  >  --output-node-data "nt_muts.ref-seq.json" \
  >  --seed 314159 \
  >  --inference marginal > /dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" \
  >   "$TESTDIR/../data/simple-genome/nt_muts.ref-seq.json" \
  >   "nt_muts.ref-seq.json" \
  >   --exclude-paths "root['generated_by']"
  {}


Change the _alignment_ to (all) lowercase, so that _if_ inference is being done
without first converting to uppercase then the inferred root-sequence will also
be lowecase which will be compared against the uppercase reference

  $ cat "$TESTDIR/../data/simple-genome/sequences.fasta" | \
  >   sed '2y/ATCGN/atcgn/' | sed '4y/ATCGN/atcgn/' | sed '6y/ATCGN/atcgn/' \
  >   > sequences.lower.fasta

  $ ${AUGUR} ancestral \
  >  --tree "$TESTDIR/../data/simple-genome/tree.nwk" \
  >  --alignment sequences.lower.fasta \
  >  --root-sequence "$TESTDIR/../data/simple-genome/reference.fasta" \
  >  --output-node-data "nt_muts.ref-seq.json" \
  >  --seed 314159 \
  >  --inference marginal > /dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" \
  >   "$TESTDIR/../data/simple-genome/nt_muts.ref-seq.json" \
  >   "nt_muts.ref-seq.json" \
  >   --exclude-paths "root['generated_by']"
  {}