Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Test the custom rules have precedence over the default rules.

Create the raw location NDJSON.

  $ cat >raw_geolocation.ndjson <<~~
  > {"region": "North America", "country": "USA", "division": "Wa", "location": "Seattle"}
  > {"region": "North America", "country": "USA", "division": "Ca", "location": "Los Angeles"}
  > ~~

Create a custom rule for a specific location that is different
than the default rule in augur/data/geolocation_rules.tsv

  $ cat >custom-rules.tsv <<~~
  > North America/USA/Wa/Seattle	North America/United States of America/Washington/Seattle
  > ~~

Run the command without the custom rules to keep the
result of annotations from default rules for comparison.

  $ cat raw_geolocation.ndjson \
  >   |  ${AUGUR} curate apply-geolocation-rules > default-annotations.ndjson

Run the command with the custom rules.

  $ cat raw_geolocation.ndjson \
  >   |  ${AUGUR} curate apply-geolocation-rules \
  >       --geolocation-rules custom-rules.tsv > custom-annotations.ndjson

Show the difference between the default annotations and the custom annotations.
Note the the custom rules only affected the single location.

  $ diff -u default-annotations.ndjson custom-annotations.ndjson
  --- default-annotations.ndjson.* (re)
  \+\+\+ custom-annotations.ndjson.* (re)
  @@ -1,2 +1,2 @@
  -{"region": "North America", "country": "USA", "division": "Washington", "location": "King County WA"}
  +{"region": "North America", "country": "United States of America", "division": "Washington", "location": "Seattle"}
   {"region": "North America", "country": "USA", "division": "California", "location": "Los Angeles County CA"}
  [1]

Create custom rule with wildcards below the country level that is different
than the default rules in augur/data/geolocation_rules.tsv

  $ cat >wildcard-custom-rules.tsv <<~~
  > North America/USA/*/*	North America/United States of America/*/*
  > ~~

Run the command with the wildcard custom rules.

  $ cat raw_geolocation.ndjson \
  >   |  ${AUGUR} curate apply-geolocation-rules \
  >       --geolocation-rules wildcard-custom-rules.tsv > wildcard-custom-annotations.ndjson

Show the difference between the default annotations and the wildcard custom annotations.
The wildcard custom rules will be applied to all records that matched and
they will no longer have the default annotations below the country level.

  $ diff -u default-annotations.ndjson wildcard-custom-annotations.ndjson
  --- default-annotations.ndjson.* (re)
  \+\+\+ wildcard-custom-annotations.ndjson.* (re)
  @@ -1,2 +1,2 @@
  -{"region": "North America", "country": "USA", "division": "Washington", "location": "King County WA"}
  -{"region": "North America", "country": "USA", "division": "California", "location": "Los Angeles County CA"}
  +{"region": "North America", "country": "United States of America", "division": "Wa", "location": "Seattle"}
  +{"region": "North America", "country": "United States of America", "division": "Ca", "location": "Los Angeles"}
  [1]
