
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur filter --sequences $TEST_DATA_DIR/in/drm.vcf.gz --metadata $TEST_DATA_DIR/in/meta.tsv --output $TMP/out/filtered.vcf.gz --exclude $TEST_DATA_DIR/in/dropped_strains.txt >/dev/null
  $ echo $?
  0

A normal 'diff' will not work here.  This is because the gzip format embeds the files compression datetime in the header.  
As such, the newly-created $TMP/out/filtered.vcf.gz will always diff against the checked-in $TEST_DATA_DIR/expected/filtered.vcf.gz.
To handle this, we unzip each file, calculate the md5sum on its contents, and compare those results.

Note that this one-liner works in bash, but not in Ubuntu-default dash.  This is why we have 'shell = /bin/bash' in .cramrc.  

  $ diff -q <(zcat $TEST_DATA_DIR/expected/filtered.vcf.gz | md5sum | cut -f1 -d' ') <(zcat $TMP/out/filtered.vcf.gz | md5sum | cut -f1 -d' ') && echo same || echo not_same
  same
