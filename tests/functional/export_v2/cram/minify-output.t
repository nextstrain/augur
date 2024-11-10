Setup

  $ source "$TESTDIR"/_setup.sh

Define a function to generate a Newick tree of a given size. Use lengthy node
names to create larger file sizes with less recursion.

  $ generate_newick() {
  >   local n=$1
  >   local prefix=$(printf 'N%.0s' {1..1000}})
  >   if [ $n -eq 1 ]; then
  >        echo "$prefix$n"
  >    else
  >        echo "($(generate_newick $((n-1))),$prefix$n)"
  >    fi
  > }

A small tree is not automatically minified.
The unminified output is 13KB which is considered acceptable.

  $ echo "$(generate_newick 10);" > small_tree.nwk

  $ ${AUGUR} export v2 \
  >   --tree small_tree.nwk \
  >   --skip-validation \
  >   --output output.json &>/dev/null

  $ head -c 20 output.json
  {
    "version": "v2", (no-eol)

  $ ls -l output.json | awk '{print $5}'
  13813

It can be forcefully minified with an argument.

  $ ${AUGUR} export v2 \
  >   --tree small_tree.nwk \
  >   --skip-validation \
  >   --minify-json \
  >   --output output.json &>/dev/null

  $ head -c 20 output.json
  {"version": "v2", "m (no-eol)

It can also be forcefully minified by setting AUGUR_MINIFY_JSON to any value,
even if it may seem "falsey".

  $ AUGUR_MINIFY_JSON=0 ${AUGUR} export v2 \
  >   --tree small_tree.nwk \
  >   --skip-validation \
  >   --output output.json &>/dev/null

  $ head -c 20 output.json
  {"version": "v2", "m (no-eol)


A large tree, when forcefully not minified, has an output size of 6MB which is
considered large.

  $ echo "$(generate_newick 500);" > big_tree.nwk

  $ ${AUGUR} export v2 \
  >   --tree big_tree.nwk \
  >   --skip-validation \
  >   --no-minify-json \
  >   --output output.json &>/dev/null

  $ head -c 20 output.json
  {
    "version": "v2", (no-eol)

  $ ls -l output.json | awk '{print $5}'
  6568454

This means it is automatically minified.

  $ ${AUGUR} export v2 \
  >   --tree big_tree.nwk \
  >   --skip-validation \
  >   --output output.json &>/dev/null

  $ head -c 20 output.json
  {"version": "v2", "m (no-eol)

  $ ls -l output.json | awk '{print $5}'
  561436
