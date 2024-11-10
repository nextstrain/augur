Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Running the command with no arguments produces the expected output

  $ echo '{"strain":"OC43"}' \
  >   | ${AUGUR} curate transform-strain-name
  {"strain": "OC43"}

Providing a strain regex to the command produces the expected output when the strain matches

  $ echo '{"strain":"OC43"}' \
  >   | ${AUGUR} curate transform-strain-name --strain-regex '^\w{2}\d{2}$'
  {"strain": "OC43"}

Providing a strain regex to the command produces an empty field and a warning when the strain doesn't match

  $ echo '{"strain":"OC43"}' \
  >   | ${AUGUR} curate transform-strain-name --strain-regex '^\d{2}\w{2}$'
  WARNING: Record number 0 has an empty string as the strain name.
  {"strain": ""}

Providing a backup field produces the expected output

  $ echo '{"potential-strain":"OC43"}' \
  >   | ${AUGUR} curate transform-strain-name --backup-fields potential-strain
  {"potential-strain": "OC43", "strain": "OC43"}


Multiple backup fields produce the expected output

  $ echo '{"potential-strain2":"OC43"}' \
  >   | ${AUGUR} curate transform-strain-name --backup-fields potential-strain potential-strain2
  {"potential-strain2": "OC43", "strain": "OC43"}
