Setup

  $ source "$TESTDIR"/_setup.sh

Infer ancestral nucleotide and amino acid sequences using a root sequence that
was used for alignment. This additional argument allows the ancestral command to
assign mutations on the branch leading to the inferred root (differences between
the "root" used as an alignment reference and the inferred most recent common
ancestor).

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree.nwk \
  >  --alignment $TESTDIR/../data/aligned.fasta \
  >  --annotation $TESTDIR/../data/zika_outgroup.gb \
  >  --root-sequence $TESTDIR/../data/zika_outgroup.gb \
  >  --genes ENV PRO \
  >  --translations $TESTDIR/../data/aa_sequences_%GENE.fasta \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations.json" > /dev/null

Redo but with root sequence sequence being fasta and with mixed lower and upper case.

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree.nwk \
  >  --alignment $TESTDIR/../data/aligned.fasta \
  >  --annotation $TESTDIR/../data/zika_outgroup.gb \
  >  --root-sequence $TESTDIR/../data/zika_outgroup_lower_upper.fasta \
  >  --genes ENV PRO \
  >  --translations $TESTDIR/../data/aa_sequences_%GENE.fasta \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations_mixed_case_fasta_root.json" > /dev/null 

Check that results are identical

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" \
  >   "$CRAMTMP/$TESTFILE/ancestral_mutations.json" \
  >   "$CRAMTMP/$TESTFILE/ancestral_mutations_mixed_case_fasta_root.json"
  {}


Check that the reference length was correctly exported as the nuc annotation

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" \
  >   --exclude-regex-paths "['seqid']" -- \
  >   "$TESTDIR/../data/ancestral_mutations_with_root_sequence.json" \
  >   "$CRAMTMP/$TESTFILE/ancestral_mutations.json"
  {}
