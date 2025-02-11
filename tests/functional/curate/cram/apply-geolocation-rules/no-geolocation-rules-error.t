Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Running the command without default rules or custom rules will result in an error.

  $ echo '{"region": "North America", "country": "USA", "division": "Wa", "location": "Seattle"}' \
  >   |  ${AUGUR} curate apply-geolocation-rules \
  >      --no-default-rules
  ERROR: No geolocation rules were provided! Either remove the `--no-default-rules` flag to use Augur's default geolocation rules or provide custom rules via `--geolocation-rules`.
  [2]
