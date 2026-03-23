Setup

  $ source "$TESTDIR"/_setup.sh

Infer ancestral sequences for the given tree and alignment, explicitly requesting that ambiguous bases are inferred.
There should not be N bases in the inferred output sequences.

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree.nwk \
  >  --alignment $TESTDIR/../data/aligned.fasta \
  >  --infer-ambiguous \
  >  --seed 314159 \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations.json" \
  >  --output-sequences "$CRAMTMP/$TESTFILE/ancestral_sequences.fasta" &> /dev/null

  $ grep "^N" "$CRAMTMP/$TESTFILE/ancestral_sequences.fasta"
  [1]

test ambiguous bases remain in seqs and muts if we use `--keep-ambiguous`
The ambiguous nuc 'Y' (representing a 'T' or 'C') should be reported as a mutation on the sample_C branch
and  also remain in the resulting sequence.

  $ cat > aln_pos3Y.fasta <<EOF
  > >sample_A
  > AAGAA
  > >sample_B
  > AAGAA
  > >sample_C
  > AAYAA
  > EOF

  $ ${AUGUR} ancestral --seed 0 \
  >  --keep-ambiguous \
  >  --tree $TESTDIR/../data/simple-genome/tree.nwk --alignment aln_pos3Y.fasta \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations_1.json" \
  >  --output-sequences "$CRAMTMP/$TESTFILE/ancestral_sequences_1.fasta" &> /dev/null

  $ cat "$CRAMTMP/$TESTFILE/ancestral_mutations_1.json" | jq -c '.nodes.sample_C.muts'
  ["C3Y"]

  $ cat "$CRAMTMP/$TESTFILE/ancestral_mutations_1.json" | jq -c '.nodes.sample_C.sequence'
  "AAYAA"

  $ grep 'sample_C' -A 1 "$CRAMTMP/$TESTFILE/ancestral_sequences_1.fasta"
  >sample_C
  AAYAA


Same test as above but using `--infer-ambiguous` (the default) infers Y as either T or C as appropriate

  $ ${AUGUR} ancestral --seed 0 \
  >  --infer-ambiguous \
  >  --tree $TESTDIR/../data/simple-genome/tree.nwk --alignment aln_pos3Y.fasta \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations_2.json" \
  >  --output-sequences "$CRAMTMP/$TESTFILE/ancestral_sequences_2.fasta" &> /dev/null

  $ cat "$CRAMTMP/$TESTFILE/ancestral_mutations_2.json" | jq -c '.nodes.sample_C.muts'
  []

  $ cat "$CRAMTMP/$TESTFILE/ancestral_mutations_2.json" | jq -c '.nodes.sample_C.sequence'
  "AACAA"

  $ grep 'sample_C' -A 1 "$CRAMTMP/$TESTFILE/ancestral_sequences_2.fasta"
  >sample_C
  AACAA


Test the ambiguous nucleotide "X" is turned into a "N" with --keep-ambiguous
'X' has a flat profile map in treetime and thus is considered ambiguous.
Ensure it's reported as "N" regardless of whether we infer ambiguous states or not
and reported in both sequences and mutations. Note that the parent state
(i.e. the "G" in "G3N" below) is not the focus of this test!

  $ cat > aln_pos3X.fasta <<EOF
  > >sample_A
  > AAGAA
  > >sample_B
  > AAGAA
  > >sample_C
  > AAXAA
  > EOF

  $ ${AUGUR} ancestral --seed 0 \
  >  --keep-ambiguous \
  >  --tree $TESTDIR/../data/simple-genome/tree.nwk --alignment aln_pos3X.fasta \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations_3.json" \
  >  --output-sequences "$CRAMTMP/$TESTFILE/ancestral_sequences_3.fasta" &> /dev/null


  $ cat "$CRAMTMP/$TESTFILE/ancestral_mutations_3.json" | jq -c '.nodes.sample_C.muts'
  ["G3N"]

  $ cat "$CRAMTMP/$TESTFILE/ancestral_mutations_3.json" | jq -c '.nodes.sample_C.sequence'
  "AANAA"

  $ grep 'sample_C' -A 1 "$CRAMTMP/$TESTFILE/ancestral_sequences_3.fasta"
  >sample_C
  AANAA

Test the ambiguous nucleotide "X" is turned into one of A/T/G/C with --infer-ambiguous
Note the specific inferred state (TreeTime chooses 'G') is not the focus of the test here.
Note 2: TreeTime chooses the root state as G as well, so there's no mutation on sample_C

  $ ${AUGUR} ancestral --seed 0 \
  >  --infer-ambiguous \
  >  --tree $TESTDIR/../data/simple-genome/tree.nwk --alignment aln_pos3X.fasta \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations_4.json" \
  >  --output-sequences "$CRAMTMP/$TESTFILE/ancestral_sequences_4.fasta" &> /dev/null


  $ cat "$CRAMTMP/$TESTFILE/ancestral_mutations_4.json" | jq -c '.nodes.sample_C.muts'
  []

  $ cat "$CRAMTMP/$TESTFILE/ancestral_mutations_4.json" | jq -c '.nodes.sample_C.sequence'
  "AAGAA"

  $ grep 'sample_C' -A 1 "$CRAMTMP/$TESTFILE/ancestral_sequences_4.fasta"
  >sample_C
  AAGAA


Test the ambiguous nucleotide "X" on internal nodes is turned into an "N"  with --keep-ambiguous
When all descendant tips of an internal node have 'X' at a position, the internal node's sequence
should report 'N' (the standard ambiguous nucleotide) when we're not inferring ambiguous bases.

  $ cat > aln_pos3X_internal.fasta <<EOF
  > >sample_A
  > AAXAA
  > >sample_B
  > AAXAA
  > >sample_C
  > AAXAA
  > EOF

  $ ${AUGUR} ancestral --seed 0 \
  >  --keep-ambiguous \
  >  --tree $TESTDIR/../data/simple-genome/tree.nwk --alignment aln_pos3X_internal.fasta \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations_5.json" \
  >  --output-sequences "$CRAMTMP/$TESTFILE/ancestral_sequences_5.fasta" &> /dev/null


  $ cat "$CRAMTMP/$TESTFILE/ancestral_mutations_5.json" | jq -c '.nodes.node_AB.sequence'
  "AANAA"


For a fully ambiguous column (all Ns) and a reference with 'X' (also ambiguous)
we shouldn't report a mutation at the root!
(Internally, the position is masked in the reconstruction, but masked root mutations are still checked.
The reference 'X' at pos 3 should be equivalent to 'N')

  $ cat > aln_col3N.fasta <<EOF
  > >sample_A
  > AANAA
  > >sample_B
  > AANAA
  > >sample_C
  > AANAA
  > EOF

  $ cat > ref_pos3X.fasta <<EOF
  > >ref
  > AAXAA
  > EOF

  $ ${AUGUR} ancestral --seed 0 \
  >  --keep-ambiguous \
  >  --root-sequence ref_pos3X.fasta \
  >  --tree $TESTDIR/../data/simple-genome/tree.nwk --alignment aln_col3N.fasta \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations_6.json" \
  >  --output-sequences "$CRAMTMP/$TESTFILE/ancestral_sequences_6.fasta" &> /dev/null


  $ cat "$CRAMTMP/$TESTFILE/ancestral_mutations_6.json" | jq -c '[.nodes[].sequence] | unique'
  ["AANAA"]

  $ cat "$CRAMTMP/$TESTFILE/ancestral_mutations_6.json" | jq -c '[.nodes[].muts[]]'
  []

  $ cat "$CRAMTMP/$TESTFILE/ancestral_mutations_6.json" | jq -c '.reference.nuc'
  "AANAA"
