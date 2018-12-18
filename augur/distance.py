"""Calculate the distance between amino acid sequences across entire genes or at a predefined subset of sites.
"""
import argparse
import Bio
import Bio.Phylo
from collections import defaultdict
import json
import numpy as np
import os
import sys

from .reconstruct_sequences import load_alignments


def read_masks(mask_file):
    ha_masks = {}
    with open(mask_file) as f:
        for line in f:
            (key, value) = line.strip().split()
            ha_masks[key] = np.frombuffer(value.encode(), 'S1').astype(int).astype(bool)

    return ha_masks


def mask_sites(aa, mask):
    return aa[mask[:len(aa)]]


def mask_distance(aaA, aaB, mask):
    """Return distance of sequences aaA and aaB by comparing sites in the given binary mask.

    >>> aaA = np.array(["A", "B", "C"], dtype="S1")
    >>> aaB = np.array(["A", "B", "D"], dtype="S1")
    >>> mask = np.array([0, 1, 1], dtype=np.bool)
    >>> mask_distance(aaA, aaB, mask)
    1
    >>> aaB = np.array(["A", "B", "X"], dtype="S1")
    >>> mask_distance(aaA, aaB, mask)
    0
    """
    sites_A = mask_sites(aaA, mask)
    sites_B = mask_sites(aaB, mask)

    # Count sites that differ between sequences excluding undetermined residues.
    distance = int(np.sum((sites_A != sites_B) & (sites_A != b"X") & (sites_B != b"X")))

    return distance


def register_arguments(parser):
    parser.add_argument("--tree", help="Newick tree", required=True)
    parser.add_argument("--alignment", nargs="+", help="sequence(s) to be used, supplied as FASTA files", required=True)
    parser.add_argument('--gene-names', nargs="+", type=str, help="names of the sequences in the alignment, same order assumed", required=True)
    parser.add_argument("--output", help="JSON file with calculated distances stored by node name and attribute name", required=True)
    parser.add_argument("--attribute-names", nargs="+", help="names to store distances associated with the corresponding masks", required=True)
    parser.add_argument("--masks", help="tab-delimited mask definitions with mask name in first column and binary mask in second column")
    parser.add_argument("--mask-names", nargs="*", help="name of each mask to use from the given masks file. Distances are calculated for each given mask and stored in each of the corresponding attribute names. If no mask is provided, all sites will be used.")


def run(args):
    # Load tree.
    tree = Bio.Phylo.read(args.tree, "newick")

    # Identify the name of the root node.
    root_node_name = tree.root.name

    # Load sequences.
    alignments = load_alignments(args.alignment, args.gene_names)

    # Index sequences by node name and gene.
    sequences_by_node_and_gene = defaultdict(dict)
    for gene, alignment in alignments.items():
        for record in alignment:
            sequences_by_node_and_gene[record.name][gene] = str(record.seq)

    # Create a single amino acid sequence for the root node based on the order of the genes.
    # Convert the string into an array for mask comparisons.
    root_node_sequence = np.fromstring(
        "".join([
            sequences_by_node_and_gene[root_node_name][gene]
            for gene in args.gene_names
        ]),
        dtype="S1"
    )

    # Determine which sites to include in distance calculations.
    # If any masks are specified, calculate the distance for each mask and store it in the corresponding attribute name.
    if args.masks:
        # Load masks.
        masks = read_masks(args.masks)

        # Map masks to attribute names where distances will be stored.
        attributes_by_mask = dict(zip(args.mask_names, args.attribute_names))
    else:
        # If no masks are specified, calculate the distance using all sites.
        mask_name = "all_sites"
        masks = {mask_name: np.ones(len(root_node_sequence)).astype(bool)}

        # Use the requested attribute name to store distances for this mask.
        attributes_by_mask = {mask_name: args.attribute_names[0]}

    # Calculate Hamming distance between the root node's sequence and the single
    # amino acid sequence for each node in the tree.
    distances_by_node = {}
    for node, node_translations in sequences_by_node_and_gene.items():
        node_sequence = np.fromstring(
            "".join([
                node_translations[gene]
                for gene in args.gene_names
            ]),
            dtype="S1"
        )

        # Calculate the distance for this node for each requested mask and store
        # the result in the corresponding attribute name.
        distances_by_node[node] = {}
        for mask_name, attribute in attributes_by_mask.items():
            distances_by_node[node][attribute] = mask_distance(
                root_node_sequence,
                node_sequence,
                masks[mask_name]
            )

    # Export distances to JSON.
    with open(args.output, "w") as oh:
        json.dump({"nodes": distances_by_node}, oh, indent=1, sort_keys=True)
