"""
Convert a prepared JSON file from augur into a FASTA file.
"""
import argparse
import Bio
import json
import logging
import pandas as pd
import sys

sys.path.append('..')

from base.sequences_process import sequence_set


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a prepared JSON file from augur into a FASTA file.")
    parser.add_argument("json", help="prepared JSON from augur")
    parser.add_argument("--metadata", help="tab-delimited file to dump prepared JSON metadata to")

    args = parser.parse_args()

    # Setup the logger.
    logger = logging.getLogger(__name__)

    # Load the JSON data.
    with open(args.json, "r") as fh:
        data = json.load(fh)

    # Prepare a sequence set.
    sequences = sequence_set(
        logger,
        data["sequences"],
        data["reference"],
        data["info"]["date_format"]
    )

    # Add the reference to output sequences if it isn't already included.
    output_sequences = sequences.seqs.values()
    if not sequences.reference_in_dataset:
        output_sequences.append(sequences.reference_seq)

    # Write sequences to standard out.
    Bio.SeqIO.write(output_sequences, sys.stdout, "fasta")

    # Prepare metadata if it has been requested.
    if args.metadata:
        metadata = [sequences.seqs[seq].attributes for seq in sequences.seqs]
        metadata_df = pd.DataFrame(metadata)
        metadata_df.to_csv(args.metadata, sep="\t", index=False)
