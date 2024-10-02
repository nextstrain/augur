Setup

  $ source "$TESTDIR"/_setup.sh

Testing records are validated appropriately by augur curate.
Create NDJSON file for testing validation catches records with different fields.

  $ cat >records.ndjson <<~~
  > {"string": "string_1"}
  > {"string": "string_2"}
  > {"string": "string_3"}
  > {"string": "string_4", "number": 123}
  > ~~

This will always pass thru the records that pass validation but should raise an
error when it encounters the record with mismatched fields.

  $ cat records.ndjson | ${AUGUR} curate passthru
  ERROR: Records do not have the same fields! Please check your input data has the same fields.
  {"string": "string_1"}
  {"string": "string_2"}
  {"string": "string_3"}
  [2]

Passing the records through multiple augur curate commands should raise the
same error when it encounters the record with mismatched fields.

  $ cat records.ndjson \
  >   | ${AUGUR} curate passthru \
  >   | ${AUGUR} curate passthru \
  >   | ${AUGUR} curate passthru
  ERROR: Records do not have the same fields! Please check your input data has the same fields.
  {"string": "string_1"}
  {"string": "string_2"}
  {"string": "string_3"}
  [2]
