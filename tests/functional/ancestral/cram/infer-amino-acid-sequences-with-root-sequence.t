Setup

  $ source "$TESTDIR"/_setup.sh

Infer ancestral nucleotide and amino acid sequences using a root sequence that
was used for alignment. This additional argument allows the ancestral command to
assign mutations on the branch leading to the inferred root (differences between
the "root" used as an alignment reference and the inferred most recent common
ancestor).

TODO: enable this test after Biopython 1.85 warning is properly addressed
<https://github.com/nextstrain/augur/issues/1727>

$ ${AUGUR} ancestral \
>  --tree $TESTDIR/../data/tree.nwk \
>  --alignment $TESTDIR/../data/aligned.fasta \
>  --annotation $TESTDIR/../data/zika_outgroup.gb \
>  --root-sequence $TESTDIR/../data/zika_outgroup.gb \
>  --genes ENV PRO \
>  --translations $TESTDIR/../data/aa_sequences_%GENE.fasta \
>  --seed 314159 \
>  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations.json" > /dev/null

Check that the reference length was correctly exported as the nuc annotation

$ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" \
>   --exclude-regex-paths "['seqid']" -- \
>   "$TESTDIR/../data/ancestral_mutations_with_root_sequence.json" \
>   "$CRAMTMP/$TESTFILE/ancestral_mutations.json"
{}
