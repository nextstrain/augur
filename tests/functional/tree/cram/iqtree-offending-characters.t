Setup

  $ source "$TESTDIR"/_setup.sh

Build a tree with offending characters, and ensure characters are not overwritten and warning is printed.

  $ ${AUGUR} tree \
  >  --alignment "$TESTDIR/../data/aligned_offending_characters.fasta" \
  >  --method iqtree \
  >  --output tree_raw.nwk \
  >  --nthreads 1 > /dev/null
  "WARNING: Potentially offending character ''' detected in taxon name KX369547.1/Cote_d'Ivore. We recommend replacing offending characters with '_' in the alignment file to avoid issues downstream."
  [1]

  $ cat "$TESTDIR/tree_raw.nwk" | grep -q "Cote_d'Ivore" && echo "Substring found"
    Substring found

  $ cat "$TESTDIR/tree_raw.nwk" | grep -q "United-Arab-Emirates" && echo "Substring found"
    Substring found