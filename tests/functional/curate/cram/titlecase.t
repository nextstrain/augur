Setup

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="${AUGUR:-../../../../bin/augur}"


Test output with articles and a mixture of lower and uppercase letters.

  $ echo '{"title": "the night OF THE LIVING DEAD"}' \
  >   | ${AUGUR} curate titlecase --titlecase-fields "title" --articles "a" "and" "of" "the" "le"
  {"title": "The Night of the Living Dead"}

Test output with hyphenated location.

  $ echo '{"location": "BRAINE-LE-COMTE, FRANCE"}' \
  >   | ${AUGUR} curate titlecase --titlecase-fields "location" --articles "a" "and" "of" "the" "le"
  {"location": "Braine-le-Comte, France"}

Test output with unicode characters

  $ echo '{"location": "Auvergne-Rhône-Alpes" }' \
  >   | ${AUGUR} curate titlecase --titlecase-fields "location"
  {"location": "Auvergne-Rh\u00f4ne-Alpes"}

Test output with abbreviations

  $ echo '{"city": "Washington DC, USA"}' \
  >   | ${AUGUR} curate titlecase --titlecase-fields "city" --abbreviations "USA" "DC"
  {"city": "Washington DC, USA"}
