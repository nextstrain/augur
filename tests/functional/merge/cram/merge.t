SETUP

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../bin/augur}"


BASIC USAGE

Full outer join like behaviour, column coalescing/overwriting, and recording of
each row's source file(s) in extra columns.

  $ cat >x.tsv <<~~
  > strain	a	b	c
  > one	X1a	X1b	X1c
  > two	X2a	X2b	X2c
  > ~~

  $ cat >y.tsv <<~~
  > strain	b	c	f	e	d
  > two		Y2c	Y2f	Y2e	Y2d
  > three			Y3f	Y3e	Y3d
  > ~~

  $ ${AUGUR} merge \
  >   --metadata X=x.tsv Y=y.tsv \
  >   --output-metadata - --quiet | csv2tsv --csv-delim $'\t' | tsv-pretty
  strain  a    b    c    f    e    d    __source_metadata_X  __source_metadata_Y
  one     X1a  X1b  X1c                                   1                    0
  two     X2a  X2b  Y2c  Y2f  Y2e  Y2d                    1                    1
  three                  Y3f  Y3e  Y3d                    0                    1

More than two inputs.

  $ cat >z.tsv <<~~
  > strain	g	c
  > one	Z1g	
  > two	Z2g	Z2c
  > three	Z3g	
  > ~~

  $ ${AUGUR} merge \
  >   --metadata X=x.tsv Y=y.tsv Z=z.tsv \
  >   --output-metadata - --quiet | csv2tsv --csv-delim $'\t' | tsv-pretty
  strain  a    b    c    f    e    d    g    __source_metadata_X  __source_metadata_Y  __source_metadata_Z
  one     X1a  X1b  X1c                 Z1g                    1                    0                    1
  two     X2a  X2b  Z2c  Y2f  Y2e  Y2d  Z2g                    1                    1                    1
  three                  Y3f  Y3e  Y3d  Z3g                    0                    1                    1

Supports Augur's standard id column detection.  Note that the first file's id
column name (e.g. "name" here) is used as the output id column name, per
Augur's convention of preserving the input id column name.

  $ sed '1s/^strain/name/g' < x.tsv > x-name-column.tsv
  $ ${AUGUR} merge \
  >   --metadata X=x-name-column.tsv Y=y.tsv \
  >   --output-metadata - --quiet | csv2tsv --csv-delim $'\t' | tsv-pretty
  name   a    b    c    f    e    d    __source_metadata_X  __source_metadata_Y
  one    X1a  X1b  X1c                                   1                    0
  two    X2a  X2b  Y2c  Y2f  Y2e  Y2d                    1                    1
  three                 Y3f  Y3e  Y3d                    0                    1

  $ sed '1s/^strain/name/g' < y.tsv > y-name-column.tsv
  $ ${AUGUR} merge \
  >   --metadata X=x.tsv Y=y-name-column.tsv \
  >   --output-metadata - --quiet | csv2tsv --csv-delim $'\t' | tsv-pretty
  strain  a    b    c    f    e    d    __source_metadata_X  __source_metadata_Y
  one     X1a  X1b  X1c                                   1                    0
  two     X2a  X2b  Y2c  Y2f  Y2e  Y2d                    1                    1
  three                  Y3f  Y3e  Y3d                    0                    1

Supports --metadata-id-columns.

  $ sed '1s/^strain/id/g' < x.tsv > x-id-column.tsv
  $ ${AUGUR} merge \
  >   --metadata X=x-id-column.tsv Y=y.tsv \
  >   --metadata-id-columns id strain \
  >   --output-metadata - --quiet | csv2tsv --csv-delim $'\t' | tsv-pretty
  id     a    b    c    f    e    d    __source_metadata_X  __source_metadata_Y
  one    X1a  X1b  X1c                                   1                    0
  two    X2a  X2b  Y2c  Y2f  Y2e  Y2d                    1                    1
  three                 Y3f  Y3e  Y3d                    0                    1

Supports Augur's standard delimiter detection.

  $ sed 's/\t/,/g' < x.tsv > x.csv
  $ ${AUGUR} merge \
  >   --metadata X=x.csv Y=y.tsv \
  >   --output-metadata - --quiet | csv2tsv --csv-delim $'\t' | tsv-pretty
  strain  a    b    c    f    e    d    __source_metadata_X  __source_metadata_Y
  one     X1a  X1b  X1c                                   1                    0
  two     X2a  X2b  Y2c  Y2f  Y2e  Y2d                    1                    1
  three                  Y3f  Y3e  Y3d                    0                    1

Supports --metadata-delimiters.

  $ sed 's/\t/|/g' < x.tsv > x.txt
  $ ${AUGUR} merge \
  >   --metadata X=x.txt Y=y.tsv \
  >   --metadata-delimiters '|' $'\t' \
  >   --output-metadata - --quiet | csv2tsv --csv-delim $'\t' | tsv-pretty
  strain  a    b    c    f    e    d    __source_metadata_X  __source_metadata_Y
  one     X1a  X1b  X1c                                   1                    0
  two     X2a  X2b  Y2c  Y2f  Y2e  Y2d                    1                    1
  three                  Y3f  Y3e  Y3d                    0                    1

Supports Augur's standard accepted compression formats.

  $ xz   < x.tsv > x.tsv.xz
  $ zstd < y.tsv > y.tsv.zst
  $ ${AUGUR} merge \
  >   --metadata X=x.tsv.xz Y=y.tsv.zst \
  >   --output-metadata - --quiet | csv2tsv --csv-delim $'\t' | tsv-pretty
  strain  a    b    c    f    e    d    __source_metadata_X  __source_metadata_Y
  one     X1a  X1b  X1c                                   1                    0
  two     X2a  X2b  Y2c  Y2f  Y2e  Y2d                    1                    1
  three                  Y3f  Y3e  Y3d                    0                    1


OFF THE BEATEN PATH

Metadata names are only the part before the first '='.

  $ cp x.tsv x=first.tsv
  $ ${AUGUR} merge \
  >   --metadata X=x=first.tsv Y=y.tsv \
  >   --output-metadata /dev/null
  Reading 'X' metadata from 'x=first.tsv'…
  Reading 'Y' metadata from 'y.tsv'…
  Merging metadata and writing to '/dev/null'…

Metadata field values with metachars (field or record delimiters) are handled properly.

  $ cat >metachars.csv <<~~
  > strain,comma,tab,newline
  > one,"x,x","x	x","x
  > x"
  > two,"","",""
  > ~~
  $ ${AUGUR} merge \
  >   --metadata X=x.tsv metachars=metachars.csv \
  >   --output-metadata - --quiet
  strain	a	b	c	comma	tab	newline	__source_metadata_X	__source_metadata_metachars
  one	X1a	X1b	X1c	x,x	"x	x"	"x
  x"	1	1
  two	X2a	X2b	X2c				1	1


ERROR HANDLING

At least two metadata inputs are required.

  $ ${AUGUR} merge \
  >   --metadata X=x.tsv \
  >   --output-metadata -
  ERROR: At least two metadata inputs are required for merging.
  [2]

Metadata names are required.

  $ ${AUGUR} merge \
  >   --metadata x.tsv =y.tsv \
  >   --output-metadata -
  ERROR: All metadata inputs must be assigned a name, e.g. with NAME=FILE.
  
  The following inputs were missing a name:
  
    'x.tsv'
    '=y.tsv'
  
  [2]

Metadata names must be unique.

  $ ${AUGUR} merge \
  >   --metadata data=x.tsv data=y.tsv \
  >   --output-metadata -
  ERROR: Metadata input names must be unique.
  
  The following names were used more than once:
  
    'data'
  
  [2]

Duplicates.

  $ cat >dups.tsv <<~~
  > strain	a	b	c
  > one	1a	1b	1c
  > one	2a	2b	2c
  > ~~
  $ ${AUGUR} merge \
  >   --metadata dups=dups.tsv Y=y.tsv \
  >   --output-metadata /dev/null
  Reading 'dups' metadata from 'dups.tsv'…
  Error: stepping, UNIQUE constraint failed: metadata_dups.strain (19)
  WARNING: Skipped deletion of */augur-merge-*.sqlite due to error, but you may want to clean it up yourself (e.g. if it's large). (glob)
  ERROR: sqlite3 invocation failed
  [2]

SQLITE3 env var can be used to override `sqlite3` location (and failure is
handled).

  $ cat >sqlite3 <<~~
  > #!/bin/bash
  > echo "This isn't the program you're looking for." >&2
  > exit 1
  > ~~
  $ chmod +x sqlite3
  $ SQLITE3="$PWD/sqlite3" \
  > augur merge \
  >   --metadata X=x.tsv Y=y.tsv \
  >   --output-metadata /dev/null --quiet
  This isn't the program you're looking for.
  ERROR: sqlite3 invocation failed
  [2]

No `sqlite3`.

  $ SQLITE3= \
  > augur merge \
  >   --metadata X=x.tsv Y=y.tsv \
  >   --output-metadata /dev/null --quiet
  ERROR: Unable to find the program `sqlite3`.  Is it installed?
  
  In order to use `augur merge`, the SQLite 3 CLI (version ≥3.39)
  must be installed separately.  It is typically provided by a
  Nextstrain runtime.
  
  [2]
