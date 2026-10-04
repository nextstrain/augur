Setup

  $ source "$TESTDIR"/_setup.sh

Infer ancestral nucleotide and amino acid sequences using a Nextclade-style GFF.
The annotation contains CDSs nested under mRNAs, and multi-row CDSs (GP, ssGP)
which are joined into single features.

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/ebola/tree.nwk \
  >  --alignment $TESTDIR/../data/ebola/masked.fasta \
  >  --annotation $TESTDIR/../data/ebola/genome_annotation.gff3 \
  >  --genes GP L NP sGP ssGP VP24 VP30 VP35 VP40 \
  >  --use-nextclade-gff-style \
  >  --translations $TESTDIR/../data/ebola/translations/%GENE.fasta \
  >  --infer-ambiguous \
  >  --inference joint \
  >  --seed 314159 \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations.json" \
  >  &> /dev/null

Check that output is as expected

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" \
  >   --exclude-regex-paths "\['seqid'\]" -- \
  >   "$TESTDIR/../data/ebola/nt_muts.json" \
  >   "$CRAMTMP/$TESTFILE/ancestral_mutations.json"
  {}

The Nextclade GFF parsing flag requires an annotation file.

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/ebola/tree.nwk \
  >  --alignment $TESTDIR/../data/ebola/masked.fasta \
  >  --use-nextclade-gff-style \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations.json"
  ERROR: --use-nextclade-gff-style requires an --annotation file.
  [2]
