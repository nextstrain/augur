Setup

  $ source "$TESTDIR"/_setup.sh

Create a copy of tests/functional/refine/data/metadata.tsv, adding partial ambiguity on the century-level (20XX) for the first strain PAN/CDC_259359_V1_V3/2015.

  $ cat >metadata.tsv <<~~
  > strain	date
  > PAN/CDC_259359_V1_V3/2015	20XX-XX-XX
  > COL/FLR_00024/2015	2015-12-XX
  > PRVABC59	2015-12-XX
  > COL/FLR_00008/2015	2015-12-XX
  > Colombia/2016/ZC204Se	2016-01-06
  > ZKC2/2016	2016-02-16
  > VEN/UF_1/2016	2016-03-25
  > DOM/2016/BB_0059	2016-04-04
  > BRA/2016/FC_6706	2016-04-08
  > DOM/2016/BB_0183	2016-04-18
  > EcEs062_16	2016-04-XX
  > HND/2016/HU_ME59	2016-05-13
  > ~~

Limit ambiguous dates to be within (2000, 2020).

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --metadata metadata.tsv \
  >  --output-tree tree.nwk \
  >  --output-node-data branch_lengths.json \
  >  --timetree \
  >  --year-bounds 2000 2020 \
  >  --coalescent opt \
  >  --date-confidence \
  >  --date-inference marginal \
  >  --clock-filter-iqd 4 \
  >  --seed 314159 \
  >  --divergence-units mutations &> /dev/null

Check that the inferred date is 2020-12-31.
TODO: Switch to use jq once it's available in a well-defined test environment.
<https://github.com/nextstrain/augur/issues/1557>

  $ python3 -c 'import json, sys; print(json.load(sys.stdin)["nodes"]["PAN/CDC_259359_V1_V3/2015"]["date"])' < branch_lengths.json
  2020-12-31

Reverse the order to check that order does not matter.

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --metadata metadata.tsv \
  >  --output-tree tree.nwk \
  >  --output-node-data branch_lengths.json \
  >  --timetree \
  >  --year-bounds 2020 2000 \
  >  --coalescent opt \
  >  --date-confidence \
  >  --date-inference marginal \
  >  --clock-filter-iqd 4 \
  >  --seed 314159 \
  >  --divergence-units mutations &> /dev/null

Run the same check as above.

  $ python3 -c 'import json, sys; print(json.load(sys.stdin)["nodes"]["PAN/CDC_259359_V1_V3/2015"]["date"])' < branch_lengths.json
  2020-12-31
