Setup

  $ source "$TESTDIR"/_setup.sh

The purpose of this test file is to check format and consistency among the
3 output file types.

  $ ${AUGUR} filter \
  >  --sequences "$TESTDIR/../data/sequences.fasta" \
  >  --sequence-index "$TESTDIR/../data/sequence_index.tsv" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --min-date 2012 \
  >  --subsample-max-sequences 5 \
  >  --subsample-seed 314159 \
  >  --no-probabilistic-sampling \
  >  --output-metadata filtered_metadata.tsv \
  >  --output-strains filtered_strains.txt \
  >  --output filtered.fasta 2>/dev/null

Check that the header row is identical between input and output metadata.

  $ diff \
  >   <(head -n 1 "$TESTDIR/../data/metadata.tsv") \
  >   <(head -n 1 filtered_metadata.tsv)

Check that the row for a strain is identical between input and output metadata.

  $ strain=Colombia/2016/ZC204Se
  $ diff \
  >   <(grep -F "$strain" "$TESTDIR/../data/metadata.tsv") \
  >   <(grep -F "$strain" filtered_metadata.tsv)

Check the order of strains in the filtered strains file.

  $ cat filtered_strains.txt
  Colombia/2016/ZC204Se
  ZKC2/2016
  DOM/2016/BB_0059
  BRA/2016/FC_6706
  EcEs062_16

Check that the order of strains in the metadata is the same as above.

  $ diff \
  >   <(cat filtered_strains.txt) \
  >   <(tail -n+2 filtered_metadata.tsv | cut -f 1)

Check the order of strains in the FASTA sequence output.

  $ grep ">" filtered.fasta
  >Colombia/2016/ZC204Se
  >ZKC2/2016
  >DOM/2016/BB_0059
  >BRA/2016/FC_6706
  >EcEs062_16

Check the order of strains in the FASTA sequence output.

  $ grep ">" filtered.fasta
  >Colombia/2016/ZC204Se
  >ZKC2/2016
  >DOM/2016/BB_0059
  >BRA/2016/FC_6706
  >EcEs062_16

Check the first 10 bases of a particular sequence in the FASTA sequence output.

  $ grep -A 1 Colombia/2016/ZC204Se filtered.fasta
  >Colombia/2016/ZC204Se
  gacagttcga.* (re)
