Setup

  $ source "$TESTDIR"/_setup.sh
  $ export ANC_DATA="$TESTDIR/../../ancestral/data/simple-genome"
  $ export DATA="$TESTDIR/../data/simple-genome"

These tests are intended to test variants of GFF formatting


GFF file with no valid rows

  $ head -n 3  $DATA/reference.gff > "reference.empty.gff"

  $ ${AUGUR} translate \
  >  --tree $ANC_DATA/tree.nwk \
  >  --ancestral-sequences $ANC_DATA/nt_muts.ref-seq.json \
  >  --reference-sequence "reference.empty.gff" \
  >  --output-node-data "aa_muts.json" > /dev/null
  ERROR: Reference 'reference.empty.gff' contains no valid data rows. .+ (re)
  [2]

GFF file with an extra record

  $ cp $DATA/reference.gff "reference.double.gff"

  $ echo -e "additional\tRefSeq\tsource\t1\t10\t.\t+\t.\tID=additional" >> "reference.double.gff"

  $ ${AUGUR} translate \
  >  --tree $ANC_DATA/tree.nwk \
  >  --ancestral-sequences $ANC_DATA/nt_muts.ref-seq.json \
  >  --reference-sequence "reference.double.gff" \
  >  --output-node-data "aa_muts.json"
  ERROR: Reference 'reference.double.gff' contains multiple seqids .+ (re)
  [2]


GFF file with data row GFF type 'region' replaced by 'source' _and_ the
##sequence-region pragma removed. This essentially mimics the information
augur 23.1.1 and earlier would use, before augur started parsing region and/or
the ##sequence-region pragma.
  $ grep -v '##sequence-region' "$DATA/reference.gff" |
  >  sed 's/\tregion\t/\tsource\t/' > "reference-only.gff"

  $ ${AUGUR} translate \
  >  --tree $ANC_DATA/tree.nwk \
  >  --ancestral-sequences $ANC_DATA/nt_muts.ref-seq.json \
  >  --reference-sequence "reference-only.gff" \
  >  --output-node-data "aa_muts-only.json" > /dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" \
  >   "$DATA/aa_muts.json" \
  >   "aa_muts-only.json" \
  >   --exclude-regex-paths "root\['annotations'\]\['.+'\]\['seqid'\]"
  {'values_changed': {"root['annotations']['nuc']['type']": {'new_value': 'source', 'old_value': '##sequence-region pragma'}}}

GFF file with data row added with GFF type 'source' with coordinates which don't match
  $ sed '5s/^/reference_name\tRefSeq\tsource\t1\t70\t.\t+\t.\tID=reference_name\n/' \
  >   "$DATA/reference.gff" > "reference-contradicts.gff"

  $ ${AUGUR} translate \
  >  --tree $ANC_DATA/tree.nwk \
  >  --ancestral-sequences $ANC_DATA/nt_muts.ref-seq.json \
  >  --reference-sequence "reference-contradicts.gff" \
  >  --output-node-data "aa_muts.json"
  ERROR: Reference .+ contained contradictory coordinates .+ (re)
  [2]

GFF file with 'region' removed, so the only genome information is the ##sequence-region pragma
  $ egrep -v '\tregion\t' "$DATA/reference.gff" > "reference.pragma-only.gff"

  $ ${AUGUR} translate \
  >  --tree $ANC_DATA/tree.nwk \
  >  --ancestral-sequences $ANC_DATA/nt_muts.ref-seq.json \
  >  --reference-sequence "reference.pragma-only.gff" \
  >  --output-node-data "aa_muts.pragma-only.json" > /dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" \
  >   "$DATA/aa_muts.json" \
  >   "aa_muts.pragma-only.json" \
  >   --exclude-regex-paths "root\['annotations'\]\['.+'\]\['seqid'\]" "root['meta']['updated']"
  {}

GFF file with no genome coordinate information
  $ egrep -v 'region' "$DATA/reference.gff" > "reference.no-nuc-info.gff"

  $ ${AUGUR} translate \
  >  --tree $ANC_DATA/tree.nwk \
  >  --ancestral-sequences $ANC_DATA/nt_muts.ref-seq.json \
  >  --reference-sequence "reference.no-nuc-info.gff" \
  >  --output-node-data "aa_muts.json"
  ERROR: Reference .+ didn't define any information we can use to create the 'nuc' annotation. .+ (re)
  [2]