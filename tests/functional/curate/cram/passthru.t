Setup

  $ source "$TESTDIR"/_setup.sh

Create NDJSON file for testing all valid JSON data types.

  $ cat >records.ndjson <<~~
  > {"string": "string", "number": 123, "object": {"string": "string"}, "array": ["string0", "string1", "string2"], "boolean1": true, "boolean2": false, "null": null}
  > ~~

Output should be exactly the same as the input.

  $ cat records.ndjson | ${AUGUR} curate passthru
  {"string": "string", "number": 123, "object": {"string": "string"}, "array": ["string0", "string1", "string2"], "boolean1": true, "boolean2": false, "null": null}
