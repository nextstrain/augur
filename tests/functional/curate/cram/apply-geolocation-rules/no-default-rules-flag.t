Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Test the default rules are ignored when using the `--no-default-rules` flag.

Create the raw location NDJSON.

  $ cat >raw_geolocation.ndjson <<~~
  > {"region": "North America", "country": "USA", "division": "Wa", "location": "Seattle"}
  > {"region": "North America", "country": "USA", "division": "Ca", "location": "Los Angeles"}
  > ~~

Create a custom rule for a specific location.

  $ cat >custom-rules.tsv <<~~
  > North America/USA/Wa/Seattle	North America/United States of America/Washington/Seattle
  > ~~

Run the command with default rules and the custom rules to keep the
result of annotations for comparison.

  $ cat raw_geolocation.ndjson \
  >   |  ${AUGUR} curate apply-geolocation-rules \
  >       --geolocation-rules custom-rules.tsv > default-and-custom-annotations.ndjson

Run the command with the custom rules and the `--no-default-rules` flag.

  $ cat raw_geolocation.ndjson \
  >   |  ${AUGUR} curate apply-geolocation-rules \
  >       --geolocation-rules custom-rules.tsv \
  >       --no-default-rules > custom-annotations-only.ndjson

Show the difference between annotations.
Note the location that is not included in the custom rules still has it's original raw value.

  $ diff -u default-and-custom-annotations.ndjson custom-annotations-only.ndjson
  --- default-and-custom-annotations.ndjson.* (re)
  \+\+\+ custom-annotations-only.ndjson.* (re)
  @@ -1,2 +1,2 @@
   {"region": "North America", "country": "United States of America", "division": "Washington", "location": "Seattle"}
  -{"region": "North America", "country": "USA", "division": "California", "location": "Los Angeles County CA"}
  +{"region": "North America", "country": "USA", "division": "Ca", "location": "Los Angeles"}
  [1]
