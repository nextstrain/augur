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

Test output with numbers

  $ echo '{"title": "2021 SARS-CoV"}' \
  >   | ${AUGUR} curate titlecase --titlecase-fields "title" --abbreviations "SARS"
  {"title": "2021 SARS-Cov"}

Test output with only numbers

  $ echo '{"int": "2021", "float": "2021.10", "address": "2021.20.30" }' \
  >   | ${AUGUR} curate titlecase --titlecase-fields "int" "float" "address"
  {"int": "2021", "float": "2021.10", "address": "2021.20.30"}

Test case that passes on empty or null values

  $ echo '{"empty": "", "null_entry":null  }' \
  >   | ${AUGUR} curate titlecase --titlecase-fields "empty" "null_entry"
  {"empty": "", "null_entry": null}


Test case that fails on a non-string int

  $ echo '{"bare_int": 2021}' \
  >   | ${AUGUR} curate titlecase --titlecase-fields "bare_int"
  ERROR: Failed to titlecase 'bare_int':2021 in record 0 because the value is a 'int' and is not a string.
  [2]

Test case that fails on complex types (e.g. arrays)

  $ echo '{"an_array": ["hello", "world"]}' \
  >   | ${AUGUR} curate titlecase --titlecase-fields "an_array"
  ERROR: Failed to titlecase 'an_array':['hello', 'world'] in record 0 because the value is a 'list' and is not a string.
  [2]

Test cases when fields do not exist, decide if this should error out and may affect ingest pipelines

  $ echo '{"region":"europe", "country":"france" }' \
  >   | ${AUGUR} curate titlecase --titlecase-fields "region" "country" "division" "location" "not exist"
  {"region": "Europe", "country": "France"}

Test output with non-string value input with `ERROR_ALL` failure reporting.
This reports a collection of all titlecase failures which is especially beneficial for automated pipelines.

  $ echo '{"bare_int": 2021, "bare_float": 1.2}' \
  >   | ${AUGUR} curate titlecase --titlecase-fields "bare_int" "bare_float" \
  >   --failure-reporting "error_all" 1> /dev/null
  ERROR: Failed to titlecase 'bare_int':2021 in record 0 because the value is a 'int' and is not a string.
  ERROR: Failed to titlecase 'bare_float':1.2 in record 0 because the value is a 'float' and is not a string.
  ERROR: Unable to change to titlecase for the following (record, field, field value):
  (0, 'bare_int', 2021)
  (0, 'bare_float', 1.2)
  [2]

Test warning on failures such as when encountering a non-string value.

  $ echo '{"bare_int": 2021}' \
  >   | ${AUGUR} curate titlecase --titlecase-fields "bare_int" \
  >   --failure-reporting "warn"
  WARNING: Failed to titlecase 'bare_int':2021 in record 0 because the value is a 'int' and is not a string.
  WARNING: Unable to change to titlecase for the following (record, field, field value):
  (0, 'bare_int', 2021)
  {"bare_int": 2021}

Test silencing on failures such as when encountering a non-string value

  $ echo '{"bare_int": 2021}' \
  >   | ${AUGUR} curate titlecase --titlecase-fields "bare_int" \
  >   --failure-reporting "silent"
  {"bare_int": 2021}

