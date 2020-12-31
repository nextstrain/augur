#!/usr/bin/env python3
# coding: utf-8
"""Identify emerging clades from previously defined clades based on a minimum
number of new mutations that have reached a minimum frequency in a given region.

Example use cases:
# Find subclades based on nucleotide mutations with defaults.
python3 scripts/identify_emerging_clades.py \
  --tree-url http://data.nextstrain.org/ncov_global.json \
  --frequencies-url http://data.nextstrain.org/ncov_global_tip-frequencies.json \
  --nextstrain-url https://nextstrain.org/ncov/global \
  --output-table nuc_subclades.tsv \
  --output-html nuc_subclades.html
# Find region-specific subclades with nucleotide mutations.
python3 scripts/identify_emerging_clades.py \
  --tree-url http://data.nextstrain.org/ncov_europe.json \
  --frequencies-url http://data.nextstrain.org/ncov_europe_tip-frequencies.json \
  --nextstrain-url https://nextstrain.org/ncov/europe \
  --filter-attribute region \
  --filter-value Europe \
  --output-table europe_nuc_subclades.tsv \
  --output-html europe_nuc_subclades.html
# Find subclades based on spike amino acid mutations.
python3 scripts/identify_emerging_clades.py \
  --tree-url http://data.nextstrain.org/ncov_global.json \
  --frequencies-url http://data.nextstrain.org/ncov_global_tip-frequencies.json \
  --nextstrain-url https://nextstrain.org/ncov/global \
  --mutation-region S \
  --minimum-mutations 1 \
  --minimum-frequency 0.1 \
  --output-table spike_subclades.tsv \
  --output-html spike_subclades.html
"""
import argparse
from augur.utils import json_to_tree
import json
import numpy as np
import pandas as pd
import requests
import sys


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Identify emerging clades from previously defined clades.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tree-url", required=True, help="URL for a Nextstrain tree JSON")
    parser.add_argument("--frequencies-url", required=True, help="URL for a Nextstrain tip frequencies JSON")
    parser.add_argument("--nextstrain-url", help="URL for the corresponding Nextstrain build")

    parser.add_argument("--clade-membership-attribute", default="clade_membership", help="name of the attribute in the tree JSON to use for clade membership")
    parser.add_argument("--mutation-region", default="nuc", help="region of the genome to inspect for mutations")
    parser.add_argument("--minimum-mutations", default=2, type=int, help="minimum number of mutations to require for subclade since its annotated parent clade")
    parser.add_argument("--minimum-frequency", default=0.2, type=float, help="minimum frequency that a subclade must have been observed at")
    parser.add_argument("--minimum-timepoints-at-frequency", default=2, type=int, help="minimum number of timepoints a given subclade must have met the frequency threshold")
    parser.add_argument("--filter-attribute", help="name of a node attribute in the tree JSON to filter tips by and correspondingly renormalize frequencies to only those tips")
    parser.add_argument("--filter-value", help="value of the associated node attribute in the tree JSON to filter tips by")

    parser.add_argument("--output-table", required=True, help="tab-delimited data frame with mutations per putative subclade")
    parser.add_argument("--output-html", help="optional HTML page with links to the given Nextstrain build highlighting positions from each putative subclade")

    args = parser.parse_args()

    # Define inputs.
    tree_url = args.tree_url
    frequencies_url = args.frequencies_url
    nextstrain_url = args.nextstrain_url

    # Define parameters.
    clade_membership_attribute = args.clade_membership_attribute
    mutation_region = args.mutation_region
    minimum_mutations = args.minimum_mutations
    minimum_frequency = args.minimum_frequency
    minimum_timepoints_at_frequency = args.minimum_timepoints_at_frequency

    filter_attribute = args.filter_attribute
    filter_value = args.filter_value

    # Define outputs.
    subclades_table = args.output_table
    subclades_links = args.output_html

    if args.output_html is not None and args.nextstrain_url is None:
        print("WARNING: HTML output requested, but a Nextstrain URL was not provided. Skipping HTML output.", file=sys.stderr)

    # Load data
    tree_json = json.loads(requests.get(tree_url).content)
    tree = json_to_tree(tree_json)
    frequencies_json = json.loads(requests.get(frequencies_url).content)

    # Convert frequency lists into numpy arrays for easier summing of frequencies per clade.
    frequencies = {}
    for key, values in frequencies_json.items():
        if isinstance(values, dict) and "frequencies" in values:
            frequencies[key] = np.array(values["frequencies"])

    pivots = np.array(frequencies_json["pivots"])

    # If the user has defined a filter on the tips of the tree, include only those
    # tips in the frequencies and renormalize them to sum to 1.
    filtered_frequencies = {}
    if filter_attribute is not None and filter_value is not None:
        for tip in tree.find_clades(terminal=True):
            if filter_attribute in tip.node_attrs and tip.node_attrs[filter_attribute]["value"] == filter_value:
                filtered_frequencies[tip.name] = frequencies[tip.name]

        if len(filtered_frequencies) > 0:
            # Renormalized the remaining frequencies to sum to 1.
            total_per_timepoint = sum(filtered_frequencies.values())
            for strain, strain_frequencies in filtered_frequencies.items():
                filtered_frequencies[strain] = strain_frequencies / total_per_timepoint

            # Confirm the normalized frequencies sum to nearly 1.
            assert all(
                np.isclose(
                    np.ones_like(pivots),
                    sum(filtered_frequencies.values())
                )
            )

            # Reassign the global frequencies variable to these filtered frequencies.
            frequencies = filtered_frequencies
        else:
            print("ERROR: None of the tips in the tree passed the given filter.")

    # Annotate the cumulative number of mutations per node from the root.
    # Find the internal node representing each distinct clade.
    clade_node_by_name = {}
    for node in tree.find_clades(terminal=False):
        if clade_membership_attribute in node.node_attrs:
            clade_name = node.node_attrs[clade_membership_attribute]["value"]
            if clade_name not in clade_node_by_name:
                clade_node_by_name[clade_name] = node

    # Find mutations in descendants of clades.
    #
    # For each major previously identified clade, look for descendant clades that
    # have accumulated a minimum number of mutations from the original parent clade
    # and that have been observed at a minimum frequency for a minimum number of
    # timepoints.
    #
    # Track all positions with mutations on the path from the parent clade to the
    # putative clade's node in the tree.
    subclades = []

    for clade_name, clade in clade_node_by_name.items():
        # Look at all internal nodes descending from the current clade.
        for node in clade.find_clades(terminal=False):
            # Skip internal nodes with one tip. These tend to be "travel history" placeholder nodes.
            if len(node.clades) == 1:
                continue

            # Skip nodes that belong to a different clade. This handles the case where currently
            # annotated clades are nested within each other.
            if node.node_attrs[clade_membership_attribute] != clade.node_attrs[clade_membership_attribute]:
                continue

            # The first node in the loop will be the current clade, so initialize the state of the mutation
            # count and sets.
            if node == clade:
                node.mutation_count = 0
                node.mutation_set = set()
            else:
                # Each descendant starts with the mutations found on the path
                # from the annotated clade node to the descendant's parent.
                node.mutation_count = node.parent.mutation_count
                node.mutation_set = node.parent.mutation_set.copy()

                # Extract positions of mutations in the requested region
                # (e.g., "nuc", "S", "E", etc.). Each mutation has the form
                # of "C1059T" where the first and last character are the
                # ancestral and derived alleles and the remaining characters
                # are the integer position in the region.
                if hasattr(node, "branch_attrs") and mutation_region in node.branch_attrs["mutations"]:
                    positions = set([
                        int(mutation[1:-1])
                        for mutation in node.branch_attrs["mutations"][mutation_region]
                    ])
                    node.mutation_count += len(positions)
                    node.mutation_set.update(positions)

            if node.mutation_count >= minimum_mutations:
                # Sum frequencies per timepoint of tips descending from this node to get the node's frequencies.
                node_frequencies = sum(
                    frequencies.get(tip.name, np.zeros_like(pivots))
                    for tip in node.find_clades(terminal=True)
                )

                # If this node has been observed at the minimum frequency
                # for the minimum number of timepoints (not necessarily consecutive),
                # add it to the list of subclade candidates.
                timepoints_above_frequency = sum(node_frequencies >= minimum_frequency)
                if timepoints_above_frequency >= minimum_timepoints_at_frequency:
                    # Create a comma-delimited list of positions for copying/pasting
                    # into Nextstrain's color by genotype field.
                    mutation_positions = ",".join([
                        str(position)
                        for position in sorted(node.mutation_set)
                    ])

                    subclades.append({
                        "parent_clade": clade_name,
                        "node": node.name,
                        "node_date": np.round(node.node_attrs["num_date"]["value"], 2),
                        "timepoints_above_frequency": timepoints_above_frequency,
                        "mutations": mutation_positions
                    })

    if len(subclades) == 0:
        print("ERROR: no putative subclades were found for the given parameters.", file=sys.stderr)
        sys.exit(1)


    # Create a data frame of putative subclades and drop redundant collections of
    # mutations, keeping only the earliest node with those mutations.
    subclades = pd.DataFrame(subclades)
    distinct_subclades = subclades.sort_values(["parent_clade", "node_date"]).groupby("mutations").first().reset_index()
    distinct_subclades["mutation_region"] = mutation_region

    # Save distinct subclades.
    distinct_subclades.to_csv(subclades_table, sep="\t", index=False)
    print(f"Found {distinct_subclades.shape[0]} distinct subclades")

    if args.output_html and args.nextstrain_url:
        # Create an HTML page with links out to Nextstrain colored by genotype
        # at the positions associated with each putative subclade.
        with open(subclades_links, "w") as oh:
            print("<ul>", file=oh)
            for index, row_df in distinct_subclades.iterrows():
                parent_clade = row_df["parent_clade"]
                mutations = row_df["mutations"]
                print(
                    f"<li><a target='_new' href='{nextstrain_url}?c=gt-{mutation_region}_{mutations}&p=grid'>{parent_clade}: {mutations}</a></li>",
                    file=oh
                )

            print("</ul>", file=oh)
