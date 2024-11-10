Setup

  $ source "$TESTDIR"/_setup.sh

Infer ancestral nucleotide and amino acid sequences using Nextclade GFF annotations.

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/ebola/tree.nwk \
  >  --alignment $TESTDIR/../data/ebola/masked.fasta \
  >  --annotation $TESTDIR/../data/ebola/genome_annotation.gff3 \
  >  --use-nextclade-gff-parsing \
  >  --translations $TESTDIR/../data/ebola/translations/%GENE.fasta \
  >  --infer-ambiguous \
  >  --inference joint \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations.json" \
  >  > /dev/null

Check that output is as expected

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" \
  >   --exclude-regex-paths "\['seqid'\]" -- \
  >   "$TESTDIR/../data/ebola/nt_muts.json" \
  >   "$CRAMTMP/$TESTFILE/ancestral_mutations.json"
  {}
