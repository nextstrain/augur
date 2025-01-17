Setup

  $ source "$TESTDIR"/_setup.sh

Specify a warning as text.

  $ ${AUGUR} export v2 \
  >   --tree "$TESTDIR/../data/tree.nwk" \
  >   --warning 'A warning with "quotes"' \
  >   --skip-validation \
  >   --output dataset.json &>/dev/null

TODO: Switch to use jq once it's available in a well-defined test environment.
<https://github.com/nextstrain/augur/issues/1557>

  $ python3 -c 'import json, sys; print(json.load(sys.stdin)["meta"]["warning"])' < dataset.json
  A warning with "quotes"

Add a warning from a markdown file.

  $ cat >warning.md <<~~
  > A warning with "quotes".
  > ~~

  $ ${AUGUR} export v2 \
  >   --tree "$TESTDIR/../data/tree.nwk" \
  >   --warning warning.md \
  >   --skip-validation \
  >   --output dataset.json &>/dev/null

  $ python3 -c 'import json, sys; print(json.load(sys.stdin)["meta"]["warning"])' < dataset.json
  A warning with "quotes".
  

Note: there is an extra newline compared to the original contents. This extra
newline is harmless and can be explained by `augur export` embedding the
original trailing newline explicitly, followed by an additional trailing newline
when writing to extracted-warning.md.
