Setup

  $ source "$TESTDIR"/_setup.sh
  $ export ANC_DATA="$TESTDIR/../../ancestral/data/simple-genome"
  $ export DATA="$TESTDIR/../data/simple-genome"

Similar tests to those in `general.t` but here testing the --genes argument.
Note that the output is a little misleading, as it's counting the 'source' GFF ID
as a feature ('nuc' in this case)

  $ ${AUGUR} translate \
  >  --tree "$ANC_DATA/tree.nwk" \
  >  --ancestral-sequences "$ANC_DATA/nt_muts.ref-seq.json" \
  >  --reference-sequence "$DATA/reference.gff" \
  >  --genes gene2 gene3 \
  >  --output-node-data "aa_muts.genes-args.json"
  Couldn't find gene gene3 in GFF or GenBank file
  Read in 2 features from reference sequence file
  Validating schema of .+ (re)
  Validating schema of .+ (re)
  amino acid mutations written to .+ (re)

  $ python3 "$SCRIPTS/diff_jsons.py" \
  >   "$DATA/aa_muts.json" \
  >   "aa_muts.genes-args.json" \
  >   --exclude-regex-paths "seqid" "gene1" "root['meta']['updated']"
  {}
Using a text file rather than command line arguments

  $ echo -e "#comment\ngene2\ngene3"> "genes.txt"

  $ ${AUGUR} translate \
  >  --tree "$ANC_DATA/tree.nwk" \
  >  --ancestral-sequences "$ANC_DATA/nt_muts.ref-seq.json" \
  >  --reference-sequence "$DATA/reference.gff" \
  >  --genes "genes.txt" \
  >  --output-node-data "aa_muts.genes-txt.json"
  Read in 2 specified genes to translate.
  Couldn't find gene gene3 in GFF or GenBank file
  Read in 2 features from reference sequence file
  Validating schema of .+ (re)
  Validating schema of .+ (re)
  amino acid mutations written to .+ (re)

  $ python3 "$SCRIPTS/diff_jsons.py" \
  >   "aa_muts.genes-args.json" \
  >   "aa_muts.genes-txt.json" \
  > --exclude-paths "root['meta']['updated']"
  {}
