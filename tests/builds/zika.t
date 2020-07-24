Run an example Zika build with augur.

Setup test data directory, output directory, and temporarily change directories to the test data directory.
Some augur commands store the exact path of their inputs in their outputs.
Running from the test data directory allows us to use relative paths that won't differ between execution environments.

  $ TEST_DATA_DIR="$TESTDIR/zika"
  $ mkdir -p "$TMP/out"
  $ pushd "$TEST_DATA_DIR" > /dev/null
  $ export AUGUR="../../../bin/augur"

Parse a FASTA whose defline contains metadata into separate sequence and metadata files.

  $ ${AUGUR} parse \
  >   --sequences "data/zika.fasta" \
  >   --output-sequences "$TMP/out/sequences.fasta" \
  >   --output-metadata "$TMP/out/metadata.tsv" \
  >   --fields strain virus accession date region country division city db segment authors url title journal paper_url \
  >   --prettify-fields region country division city

  $ diff -u "results/sequences.fasta" "$TMP/out/sequences.fasta"

Filter sequences by a minimum date and an exclusion list and only keep one sequence per country, year, and month.

  $ ${AUGUR} filter \
  >   --sequences "results/sequences.fasta" \
  >   --metadata "results/metadata.tsv" \
  >   --exclude "config/dropped_strains.txt" \
  >   --output "$TMP/out/filtered.fasta" \
  >   --group-by country year month \
  >   --sequences-per-group 1 \
  >   --subsample-seed 314159 \
  >   --min-date 2012 > /dev/null

  $ diff -u "results/filtered.fasta" "$TMP/out/filtered.fasta"

Align filtered sequences to a specific reference sequence and fill any gaps.

  $ ${AUGUR} align \
  >  --sequences "results/filtered.fasta" \
  >  --reference-sequence "config/zika_outgroup.gb" \
  >  --output "$TMP/out/aligned.fasta" \
  >  --fill-gaps > /dev/null

  $ diff -u "results/aligned.fasta" "$TMP/out/aligned.fasta"

Build a tree from the multiple sequence alignment.

  $ ${AUGUR} tree \
  >  --alignment "results/aligned.fasta" \
  >  --output "$TMP/out/tree_raw.nwk" \
  >  --method iqtree \
  >  --tree-builder-args "-seed 314159" > /dev/null

  $ python3 "$TESTDIR/../../scripts/diff_trees.py" "results/tree_raw.nwk" "$TMP/out/tree_raw.nwk" --significant-digits 5
  {}

Confirm that tree log file exists.

  $ test -e "results/aligned-delim.iqtree.log"
  $ echo $?
  0

Build a time tree from the existing tree topology, the multiple sequence alignment, and the strain metadata.

  $ ${AUGUR} refine \
  >  --tree "results/tree_raw.nwk" \
  >  --alignment "results/aligned.fasta" \
  >  --metadata "results/metadata.tsv" \
  >  --output-tree "$TMP/out/tree.nwk" \
  >  --output-node-data "$TMP/out/branch_lengths.json" \
  >  --timetree \
  >  --coalescent opt \
  >  --date-confidence \
  >  --date-inference marginal \
  >  --clock-filter-iqd 4 \
  >  --seed 314159 > /dev/null

Confirm that TreeTime trees match expected topology and branch lengths.

  $ python3 "$TESTDIR/../../scripts/diff_trees.py" "results/tree.nwk" "$TMP/out/tree.nwk" --significant-digits 2
  {}

Branch lengths and other annotations like dates are too stochastic across runs to consistently compare with reasonable precision.

#$ python3 "$TESTDIR/../../scripts/diff_jsons.py" "results/branch_lengths.json" "$TMP/out/branch_lengths.json" --significant-digits 0
#{}

Calculate tip frequencies from the tree.

  $ ${AUGUR} frequencies \
  >  --method kde \
  >  --tree "results/tree.nwk" \
  >  --metadata "results/metadata.tsv" \
  >  --pivot-interval 3 \
  >  --output "$TMP/out/zika_tip-frequencies.json" > /dev/null

  $ diff -u --ignore-matching-lines version "auspice/zika_tip-frequencies.json" "$TMP/out/zika_tip-frequencies.json"

Infer ancestral sequences from the tree.

  $ ${AUGUR} ancestral \
  >  --tree "results/tree.nwk" \
  >  --alignment "results/aligned.fasta" \
  >  --infer-ambiguous \
  >  --output-node-data "$TMP/out/nt_muts.json" \
  >  --inference joint > /dev/null

  $ diff -u --ignore-matching-lines version "results/nt_muts.json" "$TMP/out/nt_muts.json"

Infer ancestral traits from the tree.

  $ ${AUGUR} traits \
  >  --tree "results/tree.nwk" \
  >  --weights "config/trait_weights.csv" \
  >  --metadata "results/metadata.tsv" \
  >  --output-node-data "$TMP/out/traits.json" \
  >  --columns country region \
  >  --sampling-bias-correction 3 \
  >  --confidence > /dev/null

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" "results/traits.json" "$TMP/out/traits.json" --significant-digits 5
  {}

Implicit mugration model outputs are not written to the same directory as the traits output, so we cannot test for matching mugration models here.
See augur issue 541 (https://github.com/nextstrain/augur/issues/541) for more details.

#$ diff -u "results/treecountry.mugration_model.txt" "$TMP/out/treecountry.mugration_model.txt"
#$ diff -u "results/treeregion.mugration_model.txt" "$TMP/out/treeregion.mugration_model.txt"

Translate inferred ancestral and observed nucleotide sequences to amino acid mutations.

  $ ${AUGUR} translate \
  >  --tree "results/tree.nwk" \
  >  --ancestral-sequences "$TMP/out/nt_muts.json" \
  >  --reference-sequence "config/zika_outgroup.gb" \
  >  --output-node-data "$TMP/out/aa_muts.json" > /dev/null

  $ diff -u --ignore-matching-lines version "results/aa_muts.json" "$TMP/out/aa_muts.json"

Export JSON files as v1 auspice outputs.

  $ ${AUGUR} export v1 \
  >  --tree "results/tree.nwk" \
  >  --metadata "results/metadata.tsv" \
  >  --node-data "$TMP/out/branch_lengths.json" \
  >              "$TMP/out/traits.json" \
  >              "$TMP/out/nt_muts.json" \
  >              "$TMP/out/aa_muts.json" \
  >  --colors "config/colors.tsv" \
  >  --auspice-config "config/auspice_config_v1.json" \
  >  --output-tree "$TMP/out/v1_zika_tree.json" \
  >  --output-meta "$TMP/out/v1_zika_meta.json" \
  >  --output-sequence "$TMP/out/v1_zika_seq.json" > /dev/null

  $ ${AUGUR} validate export-v1 "$TMP/out/v1_zika_meta.json" "$TMP/out/v1_zika_tree.json" > /dev/null

Compare auspice metadata files, but ignore the "updated" field since this changes with the date the export command is run.

  $ diff -u --ignore-matching-lines updated "auspice/v1_zika_meta.json" "$TMP/out/v1_zika_meta.json"

Export JSON files as v2 auspice outputs.

  $ ${AUGUR} export v2 \
  >  --tree "results/tree.nwk" \
  >  --metadata "results/metadata.tsv" \
  >  --node-data "$TMP/out/branch_lengths.json" \
  >              "$TMP/out/traits.json" \
  >              "$TMP/out/nt_muts.json" \
  >              "$TMP/out/aa_muts.json" \
  >  --colors "config/colors.tsv" \
  >  --auspice-config "config/auspice_config_v2.json" \
  >  --output "$TMP/out/v2_zika.json" \
  >  --title 'Real-time tracking of Zika virus evolution -- v2 JSON' \
  >  --panels tree map entropy frequencies > /dev/null

  $ ${AUGUR} validate export-v2 "$TMP/out/v2_zika.json" > /dev/null

Switch back to the original directory where testing started.

  $ popd > /dev/null