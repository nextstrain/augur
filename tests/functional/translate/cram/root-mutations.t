Setup

  $ source "$TESTDIR"/_setup.sh
  $ export ANC_DATA="$TESTDIR/../../ancestral/data/simple-genome"
  $ export DATA="$TESTDIR/../data/simple-genome"

This is the same as the "general.t" test, but we are modifying the input data
such that the reference sequence contains "G" at pos 20 (1-based), and include a
compensating mutation G20A on the root node. (This manipulation would be much
nicer using `jq` but it's not (yet) available in all nextstrain environments.)
This results in the reference translation of gene1 to be MPCE* not MPCG*. (Note
that the compensating nuc mutation doesn't actually need to be present in the
JSON, `augur translate` just looks at the sequence attached to each node.)

  $ sed '24s/^/        "G20A",\n/' "$ANC_DATA/nt_muts.ref-seq.json" | 
  > sed 's/"nuc": "AAAAAAAAAATGCCCTGCGGG/"nuc": "AAAAAAAAAATGCCCTGCGAG/' > nt_muts.json

  $ ${AUGUR} translate \
  >  --tree "$ANC_DATA/tree.nwk" \
  >  --ancestral-sequences nt_muts.json \
  >  --reference-sequence "$DATA/reference.gff" \
  >  --output-node-data "aa_muts.json" > /dev/null

The output should be a gene1 reference of MPCE* (not MPCG*). The root-sequence
is unchanged (MPCG*). There is also a mutation E4G at the root node to compensate.

  $ python3 "$SCRIPTS/diff_jsons.py" \
  >   "$DATA/aa_muts.json" \
  >   "aa_muts.json" \
  >   --exclude-regex-paths "root\['annotations'\]\['.+'\]\['seqid'\]" "root['meta']['updated']"
  {'values_changed': {"root['reference']['gene1']": {'new_value': 'MPCE*', 'old_value': 'MPCG*'}}, 'iterable_item_added': {"root['nodes']['node_root']['aa_muts']['gene1'][1]": 'E4G'}}
