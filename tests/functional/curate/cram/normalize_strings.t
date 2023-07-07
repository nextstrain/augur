Setup

  $ source "$TESTDIR"/_setup.sh

Test two versions of C-cedilla that look the same visually but
have different code points, therefore are considered "Not equal".

  $ export DIACRITIC_1="Ç"
  $ export DIACRITIC_2="Ç"
  $ [[ "${DIACRITIC_1}" == "${DIACRITIC_2}" ]] && echo "Equal" || echo "Not equal"
  Not equal

Create NDJSON file for testing normalize-strings with different forms

  $ cat >records.ndjson <<~~
  > {"record": 1, "diacritic_1": "${DIACRITIC_1}", "diacritic_2": "${DIACRITIC_2}"}
  > ~~

Test output with default Unicode normalization form "NFC".

  $ cat records.ndjson \
  >   | ${AUGUR} curate normalize-strings
  {"record": 1, "diacritic_1": "\u00c7", "diacritic_2": "\u00c7"}

Test output with Unicode normalization form "NFKC".

  $ cat records.ndjson \
  >   | ${AUGUR} curate normalize-strings --form NFKC
  {"record": 1, "diacritic_1": "\u00c7", "diacritic_2": "\u00c7"}

Test output with Unicode normalization form "NFD".

  $ cat records.ndjson \
  >   | ${AUGUR} curate normalize-strings --form NFD
  {"record": 1, "diacritic_1": "C\u0327", "diacritic_2": "C\u0327"}

Test output with Unicode normalization form "NFKD".

  $ cat records.ndjson \
  >   | ${AUGUR} curate normalize-strings --form NFKD
  {"record": 1, "diacritic_1": "C\u0327", "diacritic_2": "C\u0327"}
