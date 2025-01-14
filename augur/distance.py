"""Calculate the distance between sequences across entire genes or at a predefined subset of sites.

Distance calculations require selection of a comparison method (to determine
which sequences to compare) and a distance map (to determine the weight of a
mismatch between any two sequences).

**Comparison methods**

Comparison methods include:

#. root: the root and all nodes in the tree (the previous default for all distances)
#. ancestor: each tip from a current season and its immediate ancestor (optionally, from a previous season)
#. pairwise: all tips pairwise (optionally, all tips from a current season against all tips in previous seasons)

Ancestor and pairwise comparisons can be calculated with or without information
about the current season. When no dates are provided, the ancestor comparison
calculates the distance between each tip and its immediate ancestor in the given
tree. Similarly, the pairwise comparison calculates the distance between all
pairs of tips in the tree.

When the user provides a "latest date", all tips sampled after that date belong
to the current season and all tips sampled on that date or prior belong to
previous seasons. When this information is available, the ancestor comparison
calculates the distance between each tip in the current season and its last
ancestor from a previous season. The pairwise comparison only calculates the
distances between tips in the current season and those from previous
seasons.

When the user also provides an "earliest date", pairwise comparisons exclude
tips sampled from previous seasons prior to the given date. These two date
parameters allow users to specify a fixed time interval for pairwise
calculations, limiting the computationally complexity of the comparisons.

For all distance calculations, a consecutive series of gap characters (`-`)
counts as a single difference between any pair of sequences. This behavior
reflects the assumption that there was an underlying biological process that
produced the insertion or deletion as a single event as opposed to multiple
independent insertion/deletion events.

**Distance maps**

Distance maps are defined in JSON format with two required top-level keys.
The `default` key specifies the numeric (floating point) value to assign to all mismatches by default.
The `map` key specifies a dictionary of weights to use for distance calculations.
These weights are indexed hierarchically by gene name and one-based gene coordinate and are assigned in either a sequence-independent or sequence-dependent manner.
The simplest possible distance map calculates Hamming distance between sequences without any site-specific weights, as shown below:

.. code-block:: json

    {
        "name": "Hamming distance",
        "default": 1,
        "map": {}
    }

To ignore specific characters such as gaps or ambiguous nucleotides from the
distance calculation, define a top-level `ignored_characters` key with a list of
characters to ignore.

.. code-block:: json

    {
        "name": "Hamming distance",
        "default": 1,
        "ignored_characters": ["-", "N"],
        "map": {}
    }

By default, distances are floating point values whose precision can be controlled with the `precision` key that defines the number of decimal places to retain for each distance.
The following example shows how to specify a precision of two decimal places in the final output:

.. code-block:: json

    {
        "name": "Hamming distance",
        "default": 1,
        "map": {},
        "precision": 2
    }

Distances can be reported as integer values by specifying an `output_type` as `integer` or `int` as follows:

.. code-block:: json

    {
        "name": "Hamming distance",
        "default": 1,
        "map": {},
        "output_type": "integer"
    }

Sequence-independent distances are defined by gene and position using a numeric
value of the same type as the default value (integer or float). The following
example is a distance map for antigenic amino acid substitutions near influenza
A/H3N2 HA's receptor binding sites. This map calculates the Hamming distance
between amino acid sequences only at seven positions in the HA1 gene:

.. code-block:: json

    {
        "name": "Koel epitope sites",
        "default": 0,
        "map": {
            "HA1": {
                "145": 1,
                "155": 1,
                "156": 1,
                "158": 1,
                "159": 1,
                "189": 1,
                "193": 1
            }
        }
    }

Sequence-dependent distances are defined by gene, position, and sequence pairs
where the `from` sequence in each pair is interpreted as the ancestral state and
the `to` sequence as the derived state. The following example is a distance map
that assigns asymmetric weights to specific amino acid substitutions at a
specific position in the influenza gene HA1:

.. code-block:: json

    {
        "default": 0.0,
        "map": {
           "HA1": {
               "112": [
                   {
                       "from": "V",
                       "to": "I",
                       "weight": 1.192
                   },
                   {
                       "from": "I",
                       "to": "V",
                       "weight": 0.002
                   }
               ]
           }
       }
    }

The distance command produces a JSON output file in standard "node data" format
that can be passed to `augur export`. In addition to the standard `nodes` field,
the JSON includes a `params` field that describes the mapping of attribute names
to requested comparisons and distance maps and any date parameters specified by
the user. The following example JSON shows a sample output when the distance
command is run with multiple comparisons and distance maps:

.. code-block:: json

    {
        "params": {
            "attributes": ["ep", "ne", "ne_star", "ep_pairwise"],
            "compare_to": ["root", "root", "ancestor", "pairwise"],
            "map_name": [
                "wolf_epitope",
                "wolf_nonepitope",
                "wolf_nonepitope",
                "wolf_epitope"
            ],
            "latest_date": "2009-10-01"
        },
        "nodes": {
            "A/Afghanistan/AF1171/2008": {
                "ep": 7,
                "ne": 6,
                "ne_star": 1,
                "ep_pairwise": {
                    "A/Aichi/78/2007": 1,
                    "A/Argentina/3509/2006": 2
                }
            }
        }
    }

"""
import Bio
import Bio.Phylo
from collections import defaultdict
import copy
from itertools import chain
import json
import pandas as pd
import sys

from .argparse_ import ExtendOverwriteDefault
from .frequency_estimators import timestamp_to_float
from .io.file import open_file
from .reconstruct_sequences import load_alignments
from .utils import annotate_parents_for_tree, first_line, read_node_data, write_json


def read_distance_map(map_file):
    """Read a distance map JSON into a dictionary and assert that the JSON follows
    the correct format. Coordinates should be one-based in the JSON and are
    converted to zero-based coordinates on load.

    Parameters
    ----------
    map_file : str
        name of a JSON file containing a valid distance map

    Returns
    -------
    dict :
        Python representation of the distance map JSON

    Examples
    --------
    >>> sorted(read_distance_map("tests/data/distance_map_weight_per_site.json").items())
    [('default', 0), ('map', {'HA1': {144: 1}})]
    >>> sorted(read_distance_map("tests/data/distance_map_weight_per_site_and_sequence.json").items())
    [('default', 0.0), ('map', {'SigPep': {0: {('W', 'P'): -8.3}}})]
    """
    # Load the JSON.
    with open_file(map_file, "r") as fh:
        json_distance_map = json.load(fh)

    # Confirm that all required fields are present.
    assert "default" in json_distance_map, "the 'default' field is required in the distance map JSON"
    assert "map" in json_distance_map, "the 'map' field is required in the distance map JSON"

    # Reconstruct the distance map with zero-based integer coordinates and
    # dictionaries for sequence-specific weights.
    distance_map = copy.deepcopy(json_distance_map)
    distance_map["map"] = {}
    for gene, site_weights in json_distance_map["map"].items():
        distance_map["map"][gene] = {}

        for site, weights in site_weights.items():
            # Convert sequence-specific weights from the JSON-compatible list of
            # dictionaries to a dictionary of sequence-pair-to-weight mappings.
            try:
                weights = {(weight["from"], weight["to"]): weight["weight"]
                           for weight in weights}
            except TypeError:
                # If weights are not iterable, keep them as site-specific
                # values.
                pass

            # Convert each one-based site string to a zero-based integer.
            distance_map["map"][gene][int(site) - 1] = weights

    # Return the distance map.
    return distance_map


def get_distance_between_nodes(node_a_sequences, node_b_sequences, distance_map, aggregate_function=max):
    """Calculate distance between the two given nodes using the given distance map.

    In cases where the distance map between sequences is asymmetric, the first
    node is interpreted as the "ancestral" sequence and the second node is
    interpreted as the "derived" sequence.

    Parameters
    ----------
    node_a_sequences, node_b_sequences : dict
        sequences by gene name for two nodes (samples) in a tree

    distance_map : dict
        definition of site-specific and, optionally, sequence-specific distances
        per gene

    Returns
    -------
    float :
        distance between node sequences based on the given map

    Examples
    --------
    >>> node_a_sequences = {"gene": "ACTG"}
    >>> node_b_sequences = {"gene": "ACGG"}
    >>> distance_map = {"default": 0, "map": {}}
    >>> get_distance_between_nodes(node_a_sequences, node_b_sequences, distance_map)
    0.0
    >>> distance_map = {"default": 1, "map": {}}
    >>> get_distance_between_nodes(node_a_sequences, node_b_sequences, distance_map)
    1.0
    >>> distance_map = {"default": 0.0, "map": {"gene": {3: 1.0}}}
    >>> get_distance_between_nodes(node_a_sequences, node_b_sequences, distance_map)
    0.0
    >>> distance_map = {"default": 0.0, "map": {"gene": {2: 3.14159, 3: 1.0}}}
    >>> get_distance_between_nodes(node_a_sequences, node_b_sequences, distance_map)
    3.14159
    >>> distance_map = {"default": 0.0, "precision": 2, "map": {"gene": {2: 3.14159, 3: 1.0}}}
    >>> get_distance_between_nodes(node_a_sequences, node_b_sequences, distance_map)
    3.14
    >>> distance_map = {"default": 0.0, "output_type": "integer", "map": {"gene": {2: 3.14159, 3: 1.0}}}
    >>> get_distance_between_nodes(node_a_sequences, node_b_sequences, distance_map)
    3
    >>> distance_map = {"default": 0.0, "output_type": "int", "map": {"gene": {2: 3.14159, 3: 1.0}}}
    >>> get_distance_between_nodes(node_a_sequences, node_b_sequences, distance_map)
    3
    >>> distance_map = {"default": 0.0, "output_type": "unsupported", "map": {"gene": {2: 3.14159, 3: 1.0}}}
    >>> get_distance_between_nodes(node_a_sequences, node_b_sequences, distance_map)
    Traceback (most recent call last):
        ...
    ValueError: Unsupported output type of 'unsupported' provided in the distance map

    For site- and sequence-specific maps, the order of the input sequences
    matters; the first sequence is treated as the ancestral sequence while the
    second is treated as the derived.

    >>> distance_map = {"default": 0.0, "map": {"gene": {2: {('T', 'G'): 0.5}}}}
    >>> get_distance_between_nodes(node_a_sequences, node_b_sequences, distance_map)
    0.5
    >>> distance_map = {"default": 0.0, "map": {"gene": {2: {('T', 'G'): 0.5}}}}
    >>> get_distance_between_nodes(node_b_sequences, node_a_sequences, distance_map)
    0.0

    Treat a single indel as one event.

    >>> node_a_sequences = {"gene": "ACTG"}
    >>> node_b_sequences = {"gene": "A--G"}
    >>> distance_map = {"default": 1, "map": {}}
    >>> get_distance_between_nodes(node_a_sequences, node_b_sequences, distance_map)
    1.0

    Use the maximum weight of all sites affected by an indel with a site-specific distance map.

    >>> distance_map = {"default": 0, "map": {"gene": {1: 1, 2: 2}}}
    >>> get_distance_between_nodes(node_a_sequences, node_b_sequences, distance_map)
    2.0

    Use the maximum weight of all mutations at all sites affected by an indel with a mutation-specific distance map.

    >>> distance_map = {
    ...     "default": 0,
    ...     "map": {
    ...         "gene": {
    ...             1: {
    ...                 ('C', 'G'): 1,
    ...                 ('C', 'A'): 2
    ...             },
    ...             2: {
    ...                 ('T', 'G'): 3,
    ...                 ('T', 'A'): 2
    ...             }
    ...         }
    ...     }
    ... }
    >>> get_distance_between_nodes(node_a_sequences, node_b_sequences, distance_map)
    3.0

    Use the maximum weight of gaps at all sites affected by an indel with a mutation-specific distance map.

    >>> distance_map = {
    ...     "default": 0,
    ...     "map": {
    ...         "gene": {
    ...             1: {
    ...                 ('C', '-'): 1,
    ...                 ('C', 'A'): 2
    ...             },
    ...             2: {
    ...                 ('T', 'G'): 3,
    ...                 ('T', '-'): 2
    ...             }
    ...         }
    ...     }
    ... }
    >>> get_distance_between_nodes(node_a_sequences, node_b_sequences, distance_map)
    2.0

    If the default value is greater than any of the site-specific mismatches and
    the specific mismatch does not have a weighted defined, use the default
    weight.

    >>> distance_map = {
    ...     "default": 4,
    ...     "map": {
    ...         "gene": {
    ...             1: {
    ...                 ('C', 'G'): 1,
    ...                 ('C', 'A'): 2
    ...             },
    ...             2: {
    ...                 ('T', 'G'): 3,
    ...                 ('T', 'A'): 2
    ...             }
    ...         }
    ...     }
    ... }
    >>> get_distance_between_nodes(node_a_sequences, node_b_sequences, distance_map)
    4.0

    Count mismatches adjacent to indel events.

    >>> node_a_sequences = {"gene": "ACTGTA"}
    >>> node_b_sequences = {"gene": "A--CCA"}
    >>> distance_map = {"default": 1, "map": {}}
    >>> get_distance_between_nodes(node_a_sequences, node_b_sequences, distance_map)
    3.0

    Ignore specific characters defined in the distance map.

    >>> node_a_sequences = {"gene": "ACTGG"}
    >>> node_b_sequences = {"gene": "A--GN"}
    >>> distance_map = {"default": 1, "ignored_characters":["-"], "map": {}}
    >>> get_distance_between_nodes(node_a_sequences, node_b_sequences, distance_map)
    1.0
    >>> distance_map = {"default": 1, "ignored_characters":["-", "N"], "map": {}}
    >>> get_distance_between_nodes(node_a_sequences, node_b_sequences, distance_map)
    0.0

    """
    distance = 0.0
    ignored_characters = distance_map.get("ignored_characters", [])
    output_type = distance_map.get("output_type", "float")
    precision = distance_map.get("precision")

    for gene in node_a_sequences:
        gene_length = len(node_a_sequences[gene])

        # In a first pass, find all mismatches between the given sequences. Use
        # this pass to identify all sites in insertion/deletion (indel) events,
        # so we can calculate an aggregate weight per event. Track each event by
        # its start site.
        mismatches_by_site = defaultdict(set)
        in_gap = False
        for site in range(gene_length):
            if (node_a_sequences[gene][site] != node_b_sequences[gene][site]
                and node_a_sequences[gene][site] not in ignored_characters
                and node_b_sequences[gene][site] not in ignored_characters):
                if node_a_sequences[gene][site] == "-" or node_b_sequences[gene][site] == "-":
                    if not in_gap:
                        gap_start = site
                        in_gap = True

                    mismatches_by_site[gap_start].add(site)
                else:
                    in_gap = False
                    mismatches_by_site[site].add(site)
            else:
                in_gap = False

        # Sum distances across mismatched sites, aggregating indel events by a
        # summary function (e.g., max, mean, etc.).
        for sites in mismatches_by_site.values():
            mismatch_distances = []
            for site in sites:
                if gene in distance_map["map"] and site in distance_map["map"][gene]:
                    # Distances can be provided as either site- and
                    # sequence-specific dictionaries of sequence pairs to
                    # weights or as site-specific weights. Check for
                    # dictionaries first.
                    if isinstance(distance_map["map"][gene][site], dict):
                        seq_ancestral = node_a_sequences[gene][site]
                        seq_derived = node_b_sequences[gene][site]

                        # Check first for a user-defined weight for the
                        # mismatched bases. This supports mismatch weights
                        # between specific characters including gaps.
                        if (seq_ancestral, seq_derived) in distance_map["map"][gene][site]:
                            mismatch_distances.append(
                                distance_map["map"][gene][site][(seq_ancestral, seq_derived)]
                            )
                        # Next, check whether the mismatch is a gap. We want to
                        # take the aggregate of all weights at this site.
                        elif seq_ancestral == "-" or seq_derived == "-":
                            mismatch_distances.append(
                                aggregate_function(
                                    chain(
                                        (distance_map["default"],),
                                        distance_map["map"][gene][site].values()
                                    )
                                )
                            )
                        # Finally, use the default weight, if no
                        # sequence-specific weights are defined.
                        else:
                            mismatch_distances.append(distance_map["default"])
                    else:
                        mismatch_distances.append(distance_map["map"][gene][site])
                else:
                    mismatch_distances.append(distance_map["default"])

            # Aggregate the distances for all sites in the current mismatch.
            distance += aggregate_function(mismatch_distances)

    if output_type in ("integer", "int"):
        return int(distance)
    elif output_type != "float":
        raise ValueError(f"Unsupported output type of '{output_type}' provided in the distance map")

    if precision is not None:
        distance = round(distance, precision)

    return distance


def get_distances_to_root(tree, sequences_by_node_and_gene, distance_map):
    """Calculate distances between all samples in the given sequences and the node
    of the given tree using the given distance map.

    Parameters
    ----------
    tree : Bio.Phylo.BaseTree.Tree
        a rooted tree whose node names match the given dictionary of sequences
        by node and gene

    sequences_by_node_and_gene : dict
        nucleotide or amino acid sequences by node name and gene

    distance_map : dict
        site-specific and, optionally, sequence-specific distances between two
        sequences

    Returns
    -------
    dict :
        distances calculated between the root sequence and each sample in the
        tree and indexed by node name

    """
    distances_by_node = {}

    # Find the root node's sequences.
    root_node_sequences = sequences_by_node_and_gene[tree.root.name]

    # Calculate distance between root and all other nodes.
    for node_name, node_sequences in sequences_by_node_and_gene.items():
        distances_by_node[node_name] = get_distance_between_nodes(
            root_node_sequences,
            node_sequences,
            distance_map
        )

    return distances_by_node


def get_distances_to_last_ancestor(tree, sequences_by_node_and_gene, distance_map, latest_date):
    """Calculate distances between each sample in the given sequences and its last
    ancestor in a previous season using the given distance map.

    Parameters
    ----------
    tree : Bio.Phylo.BaseTree.Tree
        a rooted tree whose node names match the given dictionary of sequences
        by node and gene

    sequences_by_node_and_gene : dict
        nucleotide or amino acid sequences by node name and gene

    distance_map : dict
        site-specific and, optionally, sequence-specific distances between two
        sequences

    latest_date : pandas.Timestamp
        latest date to consider a node as a potential ancestor of a given
        sample; used to define a previous season relative to the most recent
        samples in the given tree.

    Returns
    -------
    dict :
        distances calculated between each sample in the tree and its last
        ancestor sequence with distances indexed by node name

    """
    if latest_date is not None:
        latest_date = timestamp_to_float(latest_date)

    distances_by_node = {}

    # Calculate distance between each tip and its closest ancestor in the last
    # season as defined by the given latest date threshold.
    for node in tree.find_clades(terminal=True):
        # If the given latest date is not None, skip nodes that were sampled
        # prior to this date.
        if latest_date is not None and node.attr["num_date"] < latest_date:
            continue

        # Find the closest ancestor of this node that was also sampled prior to
        # the given latest date. Stop searching once we reach the root. If the
        # latest date requested is None, the immediate parent of each node will
        # be used.
        parent = node.parent
        while parent != tree.root and latest_date is not None and parent.attr["num_date"] > latest_date:
            parent = parent.parent

        # Calculate distance between current node and its ancestor.
        distances_by_node[node.name] = get_distance_between_nodes(
            sequences_by_node_and_gene[parent.name],
            sequences_by_node_and_gene[node.name],
            distance_map
        )

    return distances_by_node


def get_distances_to_all_pairs(tree, sequences_by_node_and_gene, distance_map, earliest_date=None, latest_date=None):
    """Calculate distances between each sample in the given sequences and all other
    samples in previous seasons using the given distance map.

    Parameters
    ----------
    tree : Bio.Phylo.BaseTree.Tree
        a rooted tree whose node names match the given dictionary of sequences
        by node and gene

    sequences_by_node_and_gene : dict
        nucleotide or amino acid sequences by node name and gene

    distance_map : dict
        site-specific and, optionally, sequence-specific distances between two
        sequences

    earliest_date, latest_date : pandas.Timestamp
        earliest or latest date to consider a node for comparison to a given
        sample; used to define a range of previous seasons relative to the most
        recent samples in the given tree. Dates are open intervals (inclusive)
        for interval of previous seasons. The latest date is a closed lower
        bound on the interval of the current season.

    Returns
    -------
    dict :
        distances calculated between each sample in the tree and all samples
        from previous samples with distances indexed by primary sample name and
        then past sample name

    """
    if earliest_date is not None:
        earliest_date = timestamp_to_float(earliest_date)

    if latest_date is not None:
        latest_date = timestamp_to_float(latest_date)

    distances_by_node = {}

    # Calculate distance between each tip and all tips in previous seasons as
    # defined by the given latest date threshold.
    for node in tree.find_clades(terminal=True):
        # Skip nodes that were sampled on or prior to this date.
        if latest_date is not None and node.attr["num_date"] < latest_date:
            continue

        # Distances between this node and other nodes are indexed by the other
        # node name.
        distances_by_node[node.name] = {}

        # Find all nodes that were sampled prior to the given latest date.
        for past_node in tree.find_clades(terminal=True):
            # Calculate distance between current node and the past node. Allow
            # comparison between any two nodes if an earliest or latest date is
            # not given.
            if ((earliest_date is None or past_node.attr["num_date"] >= earliest_date) and
                (latest_date is None or past_node.attr["num_date"] <= latest_date)):
                distances_by_node[node.name][past_node.name] = get_distance_between_nodes(
                    sequences_by_node_and_gene[past_node.name],
                    sequences_by_node_and_gene[node.name],
                    distance_map
                )

    return distances_by_node


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("distance", help=first_line(__doc__))
    parser.add_argument("--tree", help="Newick tree", required=True)
    parser.add_argument("--alignment", nargs="+", action=ExtendOverwriteDefault, help="sequence(s) to be used, supplied as FASTA files", required=True)
    parser.add_argument('--gene-names', nargs="+", action=ExtendOverwriteDefault, type=str, help="names of the sequences in the alignment, same order assumed", required=True)
    parser.add_argument("--attribute-name", nargs="+", action=ExtendOverwriteDefault, help="name to store distances associated with the given distance map; multiple attribute names are linked to corresponding positional comparison method and distance map arguments", required=True)
    parser.add_argument("--compare-to", nargs="+", action=ExtendOverwriteDefault, choices=["root", "ancestor", "pairwise"], help="type of comparison between samples in the given tree including comparison of all nodes to the root (root), all tips to their last ancestor from a previous season (ancestor), or all tips from the current season to all tips in previous seasons (pairwise)", required=True)
    parser.add_argument("--map", nargs="+", action=ExtendOverwriteDefault, help="JSON providing the distance map between sites and, optionally, sequences present at those sites; the distance map JSON minimally requires a 'default' field defining a default numeric distance and a 'map' field defining a dictionary of genes and one-based coordinates", required=True)
    parser.add_argument("--date-annotations", help="JSON of branch lengths and date annotations from augur refine for samples in the given tree; required for comparisons to earliest or latest date")
    parser.add_argument("--earliest-date", help="earliest date at which samples are considered to be from previous seasons (e.g., 2019-01-01). This date is only used in pairwise comparisons. If omitted, all samples prior to the latest date will be considered.")
    parser.add_argument("--latest-date", help="latest date at which samples are considered to be from previous seasons (e.g., 2019-01-01); samples from any date after this are considered part of the current season")
    parser.add_argument("--output", help="JSON file with calculated distances stored by node name and attribute name", required=True)
    return parser


def run(args):
    # Load tree and annotate parents.
    tree = Bio.Phylo.read(args.tree, "newick")
    tree = annotate_parents_for_tree(tree)

    # Load sequences.
    alignments = load_alignments(args.alignment, args.gene_names)

    # Index sequences by node name and gene.
    sequences_by_node_and_gene = defaultdict(dict)
    for gene, alignment in alignments.items():
        for record in alignment:
            sequences_by_node_and_gene[record.name][gene] = str(record.seq)

    # Prepare the earliest date to compare nodes against.
    if args.earliest_date is None:
        earliest_date = None
    else:
        earliest_date = pd.Timestamp(args.earliest_date)

    # Prepare the latest date to compare nodes against.
    if args.latest_date is None:
        latest_date = None
    else:
        latest_date = pd.Timestamp(args.latest_date)

    # Load date annotations if any date ranges have been requested.
    if earliest_date or latest_date:
        if args.date_annotations is None:
            print(
                "ERROR: date annotations JSON from augur refine (e.g., branch_lengths.json) is required",
                file=sys.stderr
            )
            sys.exit(1)

        date_annotations = read_node_data(args.date_annotations)

        # Annotate tree with date annotations.
        for node in tree.find_clades():
            node.attr = date_annotations["nodes"][node.name]
            node.attr["num_date"] = node.attr["numdate"]

    final_distances_by_node = {}
    distance_map_names = []
    for compare_to, attribute, distance_map_file in zip(args.compare_to, args.attribute_name, args.map):
        # Load the given distance map.
        distance_map = read_distance_map(distance_map_file)
        distance_map_names.append(distance_map.get("name", distance_map_file))

        # Use the distance map to calculate distances between all samples in the
        # given tree and the desired target(s).
        if compare_to == "root":
            # Calculate distance between the root and all samples.
            distances_by_node = get_distances_to_root(
                tree,
                sequences_by_node_and_gene,
                distance_map
            )
        elif compare_to == "ancestor":
            # Calculate distance between the last ancestor for each sample in a
            # previous season.
            distances_by_node = get_distances_to_last_ancestor(
                tree,
                sequences_by_node_and_gene,
                distance_map,
                latest_date
            )
        elif compare_to == "pairwise":
            # Calculate distance between each sample and all other samples in a
            # previous season.
            distances_by_node = get_distances_to_all_pairs(
                tree,
                sequences_by_node_and_gene,
                distance_map,
                earliest_date,
                latest_date
            )
        else:
            # If the command line argument choices are defined properly above,
            # this block should never execute.
            print("ERROR: the comparison method '%s' is not supported" % compare_to, file=sys.stderr)
            sys.exit(1)

        # Map distances to the requested attribute name.
        # Convert data like:
        # {
        #   "A/AbuDhabi/24/2017": 1
        # }
        # to data like:
        #
        # {
        #   "A/AbuDhabi/24/2017": {
        #     "ep": 1
        #   }
        # }
        #
        for node_name, values in distances_by_node.items():
            if node_name not in final_distances_by_node:
                final_distances_by_node[node_name] = {}

            final_distances_by_node[node_name][attribute] = values

    # Prepare params for export.
    params = {
        "attribute": args.attribute_name,
        "compare_to": args.compare_to,
        "map_name": distance_map_names,
        "earliest_date": args.earliest_date,
        "latest_date": args.latest_date
    }

    # Export distances to JSON.
    write_json({"params": params, "nodes": final_distances_by_node}, args.output)
