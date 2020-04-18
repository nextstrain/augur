
The 'augur validate' tests do not map to their corresponding snakemake builds in quite the same way as other CLI tests do.
In the snakemake builds, ‘augur validate’ calls implicitly depend on a preceding ‘augur export’ call in the same Snakefile rule. 
As such, the input files (under ./in) correspond to *output* from the corresponding snakemake rule.

Further, the only output for 'augur validate' is the logging emitted to stdout.  For the purposes of this test, we write that logging
to $TMP/out and compare it to the results in expected/.

Additionally, the output from 'augur validate' contains the filepath you pass into the command.  If we passed in the standard param values
using $TEST_DATA_DIR, the test-time output would thus contain the fully-qualified path to the file - which would cause the diff to fail on different machines.
So, instead, here we explicitly cd into the test directory and pass in relative paths.

  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"
  $ cd $TEST_DATA_DIR
  $ augur validate export-v1 ./in/v1_zika_meta.json ./in/v1_zika_tree.json > $TMP/out/validate_output.txt
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/validate_output.txt $TMP/out/validate_output.txt
  $ echo $?
  0
