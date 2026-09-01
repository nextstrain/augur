Setup

  $ source "$TESTDIR"/_setup.sh

Run with a config file.

  $ cat > config.yaml <<~~
  > sequences:
  >   - "$TESTDIR/../../../data/align/test_unaligned_sequences.fasta"
  > output: "aligned.fasta"
  > nthreads: 1
  > ~~

  $ ${AUGUR} align --config config.yaml
  
  using mafft to align via:
  	mafft --reorder --anysymbol --nomemsave --adjustdirection --thread 1 aligned.fasta.to_align.fasta 1> aligned.fasta 2> aligned.fasta.log 
  
  	Katoh et al, Nucleic Acid Research, vol 30, issue 14
  	https://doi.org/10.1093%2Fnar%2Fgkf436
  
  Sequence "crick_strand" was reverse-complemented by the alignment program.

  $ cat "aligned.fasta"
  >with_gaps
  ---ATATA---
  >no_gaps
  GGGATATAGGG
  >some_other_seq
  --GATCTAGGG
  >crick_strand some description
  -GGATATAGG-
