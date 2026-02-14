Setup

  $ source "$TESTDIR"/_setup.sh

Update metadata with some missing/invalid dates.

  $ cat > metadata.tsv <<~~
  > strain	date
  > PAN/CDC_259359_V1_V3/2015	
  > COL/FLR_00024/2015	unused
  > PRVABC59	NA
  > COL/FLR_00008/2015	invalid
  > Colombia/2016/ZC204Se	XXXX-XX-XX
  > ZKC2/2016	2016-02-16
  > VEN/UF_1/2016	2016-03-25
  > DOM/2016/BB_0059	2016-04-04
  > BRA/2016/FC_6706	2016-04-08
  > DOM/2016/BB_0183	2016-04-18
  > EcEs062_16	2016-04-XX
  > HND/2016/HU_ME59	2016-05-13
  > ~~

Error on missing dates.

  $ ${AUGUR} frequencies \
  >  --method kde \
  >  --tree "$TESTDIR/../data/tree.nwk" \
  >  --metadata metadata.tsv \
  >  --pivot-interval 3 \
  >  --output tip-frequencies.json > /dev/null
  ERROR: The following sequence ids are missing valid dates:
  
    'COL/FLR_00008/2015' has date 'invalid'
    'PAN/CDC_259359_V1_V3/2015' has date ''
    'PRVABC59' has date 'NA'
  
  [2]
