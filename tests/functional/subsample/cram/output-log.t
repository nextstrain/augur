Setup

  $ source "$TESTDIR"/_setup.sh

Test two overlapping samples with an output log.

  $ cat >samples.yaml <<~~
  > samples:
  >   A:
  >     max_sequences: 4
  >   B:
  >     max_sequences: 5
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/sequences.fasta \
  >   --config samples.yaml \
  >   --output-metadata subsampled.tsv \
  >   --output-sequences subsampled.fasta \
  >   --output-log log.tsv \
  >   --seed 0
  Validating schema of 'samples.yaml'...
  [A] 8 strains were dropped during filtering
  [A] 	8 were dropped because of subsampling criteria
  [A] 4 strains passed all filters
  [B] 7 strains were dropped during filtering
  [B] 	7 were dropped because of subsampling criteria
  [B] 5 strains passed all filters
  [collect samples] combining outputs from 2 samples
  [collect samples] 8 strains were dropped during filtering
  [collect samples] 	1 had no metadata
  [collect samples] 	12 were dropped by `--exclude-all`
  \[collect samples\] \\t4 were added back because they were in .*sample_A.* (re)
  \[collect samples\] \\t5 were added back because they were in .*sample_B.* (re)
  [collect samples] 5 strains passed all filters

Show the output log.

  $ tail -n +2 log.tsv | sort
  BRA/2016/FC_6706	filter_by_exclude_all	[]
  COL/FLR_00008/2015	filter_by_exclude_all	[]
  COL/FLR_00024/2015	filter_by_exclude_all	[]
  Colombia/2016/ZC204Se	filter_by_exclude_all	[]
  DOM/2016/BB_0059	filter_by_exclude_all	[]
  DOM/2016/BB_0059\\tforce_include_strains\\t"\[\[""include_file"", .*sample_A.* (re)
  DOM/2016/BB_0059\\tforce_include_strains\\t"\[\[""include_file"", .*sample_B.* (re)
  DOM/2016/BB_0183	filter_by_exclude_all	[]
  DOM/2016/BB_0183\\tforce_include_strains\\t"\[\[""include_file"", .*sample_A.* (re)
  DOM/2016/BB_0183\\tforce_include_strains\\t"\[\[""include_file"", .*sample_B.* (re)
  EcEs062_16	filter_by_exclude_all	[]
  HND/2016/HU_ME59	filter_by_exclude_all	[]
  HND/2016/HU_ME59\\tforce_include_strains\\t"\[\[""include_file"", .*sample_B.* (re)
  PRVABC59	filter_by_exclude_all	[]
  SG_018	filter_by_exclude_all	[]
  SG_018\\tforce_include_strains\\t"\[\[""include_file"", .*sample_A.* (re)
  SG_018\\tforce_include_strains\\t"\[\[""include_file"", .*sample_B.* (re)
  VEN/UF_1/2016	filter_by_exclude_all	[]
  VEN/UF_1/2016\\tforce_include_strains\\t"\[\[""include_file"", .*sample_A.* (re)
  VEN/UF_1/2016\\tforce_include_strains\\t"\[\[""include_file"", .*sample_B.* (re)
  ZKC2/2016	filter_by_exclude_all	[]
