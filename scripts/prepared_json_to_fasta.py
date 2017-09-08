"""
Convert a prepared JSON file from augur into a FASTA file.
"""
import argparse
import Bio
import json
import logging
import sys

sys.path.append('..')

from base.sequences_process import sequence_set


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("json", help="prepared JSON from augur")
    parser.add_argument("fasta", help="FASTA output for sequences in JSON file")

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

    # Write sequences to disk.
    Bio.SeqIO.write(output_sequences, args.fasta, "fasta")
