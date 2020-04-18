
Unlike most of the CLI tests expected/ files, this tests expected output does *not* precisely match that of the original build.  expected/aa_muts.json has been 
updated with filepaths that work for this tests directory structure (e.g. replacing "data/Mtb_H37Rv_NCBI_Annot.gff" with "in/Mtb_H37Rv_NCBI_Annot.gff").

  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

Explicitly cd into TEST_DATA_DIR, so we can pass a relative path to --reference-sequence.  
The value for this param gets written into the output file.  So if we just passed the full path, e.g. $TEST_DATA_DIR/in/foobar.xyz,
this test would fail on other machines/environments.

  $ cd $TEST_DATA_DIR
  $ augur translate --tree $TEST_DATA_DIR/in/tree.nwk --genes $TEST_DATA_DIR/in/genes.txt --vcf-reference $TEST_DATA_DIR/in/ref.fasta --ancestral-sequences $TEST_DATA_DIR/in/nt_muts.vcf --output-node-data $TMP/out/aa_muts.json --reference-sequence in/Mtb_H37Rv_NCBI_Annot.gff --alignment-output $TMP/out/translations.vcf --vcf-reference-output $TMP/out/translations_reference.fasta >/dev/null
  Gene length of rrs_Rvnr01 is not a multiple of 3. will pad with N
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/aa_muts.json $TMP/out/aa_muts.json
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/translations.vcf $TMP/out/translations.vcf
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/translations_reference.fasta $TMP/out/translations_reference.fasta
  $ echo $?
  0
