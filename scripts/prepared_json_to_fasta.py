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

    # Write sequences to standard out.
    output_sequences = sequences.seqs.values()
    Bio.SeqIO.write(output_sequences, sys.stdout, "fasta")

    # Prepare metadata if it has been requested.
    if args.metadata:
        metadata = [sequences.seqs[seq].attributes for seq in sequences.seqs]
        metadata_df = pd.DataFrame(metadata)
        metadata_df = metadata_df.rename(columns={"num_date": "prepared_num_date"})
        metadata_df.to_csv(args.metadata, sep="\t", index=False)
