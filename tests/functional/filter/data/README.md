# Test files for `augur filter`

## Overview

These files are provided for functional testing of `augur filter`.

## Files

- Data on 12 Zika sequences:
    - `metadata.tsv`: File to be used with `augur filter --metadata`.
    - `sequences.fasta`: File to be used with `augur filter --sequences`. If updated, please regenerate `sequence_index.tsv`.
    - `priorities.tsv`: File to be used with `augur filter --priority`.
    - `include.txt`: File to be used with `augur filter --include`.
    - `sequence_index.tsv`: File to be used with `augur filter --sequence-index`. Generated using:

        ```
        augur index --sequences sequences.fasta --output sequence_index.tsv
        ```

- Data on 165 Tuberculosis sequences:
    - `tb_metadata.tsv`: File to be used with `augur filter --metadata`.
    - `tb.vcf.gz`: File to be used with `augur filter --sequences`.
