Setup

  $ source "$TESTDIR"/_setup.sh
  $ export ANC_DATA="$TESTDIR/../../ancestral/data/simple-genome"
  $ export DATA="$TESTDIR/../data/simple-genome"

General running of augur translate. See the cram test general.t for `augur ancestral`
which uses many of the same files.
NOTE: The GFF reference here contains a 'source' ID because without this downstream commands
which validate the output will fail as it's missing a 'nuc' annotation.

  $ ${AUGUR} translate \
  >  --tree "$ANC_DATA/tree.nwk" \
  >  --ancestral-sequences "$ANC_DATA/nt_muts.ref-seq.json" \
  >  --reference-sequence "$DATA/reference.gff" \
  >  --output-node-data "aa_muts.json" > /dev/null

  $ python3 "$SCRIPTS/diff_jsons.py" \
  >   "$DATA/aa_muts.json" \
  >   "aa_muts.json" \
  >   --exclude-regex-paths "root\['annotations'\]\['.+'\]\['seqid'\]" "root['meta']['updated']" 
  {}

Same as above but using a GenBank file. This changes the 'type' of the annotations,
but this is irrelevant for Auspice's use and simply reflects the reference source.

  $ ${AUGUR} translate \
  >  --tree "$ANC_DATA/tree.nwk" \
  >  --ancestral-sequences "$ANC_DATA/nt_muts.ref-seq.json" \
  >  --reference-sequence "$DATA/reference.gb" \
  >  --output-node-data "aa_muts.genbank.json" > /dev/null

  $ python3 "$SCRIPTS/diff_jsons.py" \
  >   "$DATA/aa_muts.json" \
  >   "aa_muts.genbank.json" \
  >   --exclude-regex-paths "root\['annotations'\]\['.+'\]\['seqid'\]" "root\['annotations'\]\['.+'\]\['type'\]" "root['meta']['updated']" 
  {}
