Setup

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="${AUGUR:-../../../../bin/augur}"

Create NDJSON file for testing all valid JSON data types.

  $ cat >$TMP/records.ndjson <<~~
  > {"string": "string", "number": 123, "object": {"string": "string"}, "array": ["string0", "string1", "string2"], "boolean1": true, "boolean2": false, "null": null}
  > ~~

Output should be exactly the same as the input.

  $ cat $TMP/records.ndjson | ${AUGUR} curate passthru
  {"string": "string", "number": 123, "object": {"string": "string"}, "array": ["string0", "string1", "string2"], "boolean1": true, "boolean2": false, "null": null}
