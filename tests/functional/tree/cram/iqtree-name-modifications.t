Setup

  $ source "$TESTDIR"/_setup.sh

BACKGROUND: IQ tree modifies strain names to remove certain characters. We
attempt to prevent this within augur tree (by replacing them on the way into
IQ-TREE and switching them back on the way out). Beyond this certain characters
can't be written by Bio.Phylo as newick strain names (without escaping).

# Start with a trivial test

  $ echo -e ">AAA\nATGC\n>BBB\nATGC\n>CCC\nATGC\n" > trivial.mfa

  $ ${AUGUR} tree \
  >  --method iqtree \
  >  --alignment trivial.mfa \
  >  --output trivial.new 1>/dev/null

Test single quotes, which were the errors behind <https://github.com/nextstrain/augur/issues/1084>,
<https://github.com/nextstrain/avian-flu/issues/122> 

  $ echo -e ">Coted'Ivoire/IPCI-DVE-GR1901/2022\nATGC" > single-quotes.mfa
  $ echo -e ">A/Geoffroy'sCat/WA/24-037410-004-original/2024\nATGC" >> single-quotes.mfa
  $ echo -e ">need-three-taxa\nATGC" >> single-quotes.mfa

  $ ${AUGUR} tree \
  >  --method iqtree \
  >  --alignment single-quotes.mfa \
  >  --output single-quotes.new 1>/dev/null
  ERROR: Certain strain names have characters which cannot be written in newick format (by Bio.Python, at least).
  You should ensure these strain names are changed as early as possible in your analysis. The following
  names are unescaped and surrounded by double quotes.
  
  Invalid strain names:
    - .+ (re)
    - .+ (re)
  
  Invalid characters: .* (re)
  
  [2]

Test some other observed strain names to ensure they're not modified
  $ echo -e ">A/PETFOOD/USA:OR/24-037325-011/2024\nATGC" > reported.mfa
  $ echo -e ">[name>$!%]:1234\nATGC" >> reported.mfa
  $ echo -e ">simple\nATGC" >> reported.mfa

  $ ${AUGUR} tree \
  >  --method iqtree \
  >  --alignment reported.mfa \
  >  --output reported.new 1>/dev/null


Backslashes are problematic, but only some of the time. For instance (and assuming we replaced the
character for the IQ-TREE building) the first three of these strains work and the others are
escaped by `Bio.Phylo.write`. Erring on the side of caution we don't allow any backslashes,
similar to single quotes.

  $ echo '>sing\le' > backslashes.mfa
  $ echo 'ATGC' >> backslashes.mfa
  $ echo '>dou\\ble' >> backslashes.mfa
  $ echo 'ATGC' >> backslashes.mfa
  $ echo '>three\\\ee' >> backslashes.mfa
  $ echo 'ATGC' >> backslashes.mfa
  $ echo '>v~`,\' >> backslashes.mfa
  $ echo 'ATGC' >> backslashes.mfa
  $ echo '>wb)\o' >> backslashes.mfa
  $ echo 'ATGC' >> backslashes.mfa

  $ ${AUGUR} tree \
  >  --method iqtree \
  >  --alignment backslashes.mfa \
  >  --output backslashes.new 1>/dev/null
  ERROR: Certain strain names have characters which cannot be written in newick format (by Bio.Python, at least).
  You should ensure these strain names are changed as early as possible in your analysis. The following
  names are unescaped and surrounded by double quotes.
  
  Invalid strain names:
    - .+ (re)
    - .+ (re)
    - .+ (re)
    - .+ (re)
    - .+ (re)
  
  Invalid characters: .* (re)
  
  [2]

This test generates random ASCII names. It's disabled for CI as it's both
stochastic and slow but you can easily toggle it back on (by uncommenting the
function call) if you want to better test strain name handling in `augur tree`

  $ random_ascii_names() {
  >   python3 "$TESTDIR"/generate-fasta.py > random_${1}.mfa
  >  
  >   ${AUGUR} tree \
  >    --method iqtree \
  >    --alignment random_${1}.mfa \
  >    --output random_${1}.new > /dev/null
  > }

  $ for iteration in $(seq 1 100); do
  >    # random_ascii_names $iteration
  >    :
  > done