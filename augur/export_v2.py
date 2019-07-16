"""
Export JSON files suitable for visualization with auspice.
"""

import os, sys
import re
import time
import numpy as np
from Bio import Phylo
from collections import defaultdict
from .utils import read_metadata, read_node_data, write_json, read_config, read_lat_longs, read_colors

def convert_tree_to_json_structure(node, metadata, div=0, strains=None):
    """
    converts the Biopython tree structure to a dictionary that can
    be written to file as a json. This is called recursively.
    Creates the strain property & divergence on each node

    input
        node -- node for which top level dict is produced.
        div  -- cumulative divergence (root = 0)

    returns
        tree in JSON structure
        list of strains
    """
    node_struct = {
        'strain': node.name,
        'div': div
    }

    if strains is None:
        strains = [node_struct["strain"]]
    else:
        strains.append(node_struct["strain"])

    if node.clades:
        node_struct["children"] = []
        for child in node.clades:
            if 'mutation_length' in metadata[child.name]:
                cdiv = div + metadata[child.name]['mutation_length']
            elif 'branch_length' in metadata[child.name]:
                cdiv = div + metadata[child.name]['branch_length']
            node_struct["children"].append(convert_tree_to_json_structure(child, metadata, div=cdiv, strains=strains)[0])

    return (node_struct, strains)


def get_colorings(config, traits, provided_colors, node_metadata, mutations_present):
    def _get_values(key):
        values_in_tree = set()
        for node_properties in node_metadata.values():
            # beware of special cases!
            if key == "citekey":
                values_in_tree.add(node_properties.get(key))
            elif key == "num_date":
                if key in node_properties:
                    values_in_tree.add(node_properties.get(key).get("value"))
            elif key in node_properties['traits']:
                values_in_tree.add(node_properties['traits'][key]['value'])
        return values_in_tree

    def _guess_type(key, values):
        # TODO detect other types (ordinal, boolean)
        if all([ isinstance(n, float) if isinstance(n, float) else isinstance(n, int) for n in values ]):
            return "continuous"
        return "categorical"

    def _guess_title(key):
        if key == "clade_membership":
            return "Clade"
        return key

    # NOTE: Emma -- currently if "country" is defined as a command-line trait
    # but there is a "color_options" in the config without it, then "country"
    # won't be a color-by. Not sure what behavior is desired here.
    if config.get("color_options"):
        color_config = config["color_options"]
    else:
        color_config = traits

    colorings = {}
    # handle special cases
    if mutations_present:
        colorings["gt"] = config.get("gt") if config.get("citekey") else {'title': 'Genotype', 'type': 'ordinal'}
    if _get_values("citekey"): # check if any citekeys (authors) are set on the tree
        colorings["citekey"] = config.get("citekey") if config.get("citekey") else {'title': 'Authors', 'type': 'categorical'}
    if _get_values("num_date"): # TODO: check if tree has temporal inference (possible to not have)
        colorings['num_date'] = config.get("num_date") if config.get("num_date") else {'title': 'Sampling date', 'type': 'continuous'}
    # remove keys which may exist but are special cases and/or have been handled above
    color_config.pop("gt", None)
    color_config.pop("citekey", None)
    color_config.pop("authors", None) # this was v1 syntax
    color_config.pop("num_date", None)


    for key, data in color_config.items():
        trait_values = _get_values(key) # e.g. list of countries, list of citekeys etc
        if not trait_values:
            print("WARNING - you asked for a color by for '{}' but it has no values on the tree => Ignoring.".format(key))
            continue
        colorings[key] = {
            "title": data.get("title") if data.get("title") else _guess_title(key),
            "type": data.get("type") if data.get("type") else _guess_type(key, trait_values)
        }
        # set color maps if provided in the config (FYI - color maps are in lower case)
        if key.lower() in provided_colors:
            values_lower = {str(val).lower(): val for val in trait_values}
            colorings[key]["scale"] = {values_lower[m[0]]: m[1] for m in provided_colors[key.lower()] if m[0] in values_lower}

    return colorings


def process_geographic_info(jsn, lat_long_mapping, node_metadata=None, nodes=None):
    if jsn is None or len(jsn)==0 :
        return {}
    geo = defaultdict(dict)

    traits = jsn #jsn["geographic_info"]

    for trait in traits:
        demes_in_tree = {data["traits"][trait]["value"] for name, data in node_metadata.items() if trait in data['traits']}
        for deme in demes_in_tree:
            try:
                geo[trait][deme] = lat_long_mapping[(trait.lower(),deme.lower())]
            except KeyError:
                print("Error. {}->{} did not have an associated lat/long value (matching performed in lower case)".format(trait, deme))
    return geo

def process_annotations(node_data):
    # treetime adds "annotations" to node_data
    if "annotations" not in node_data: # if haven't run tree through treetime
        return None
    return node_data["annotations"]

def process_panels(user_panels, meta_json):
    try:
        panels = meta_json["panels"]
    except KeyError:
        panels = ["tree", "map", "entropy"]

    if user_panels is not None and len(user_panels) != 0:
        panels = user_panels

    if "geographic_info" in meta_json:
        geoTraits = meta_json["geographic_info"].keys()
    else:
        geoTraits = []

    if "genome_annotations" in meta_json:
        annotations = meta_json["genome_annotations"].keys()
    else:
        annotations = []

    if "entropy" in panels and len(annotations) == 0:
        panels.remove("entropy")
    if "map" in panels and len(geoTraits) == 0:
        panels.remove("map")

    return panels


def collect_strain_info(node_data, tsv_path):
    """
    Integrate TSV metadata to the per-node metadata structure
    """
    strain_info = node_data["nodes"]
    meta_tsv, _ = read_metadata(tsv_path)

    for strain, node in strain_info.items():
        if strain in meta_tsv:
            for field in meta_tsv[strain]:
                node[field] = meta_tsv[strain][field]
    return strain_info


def construct_author_info_and_make_keys(node_metadata, raw_strain_info):
    """Gather the authors which appear in the metadata and assign them
    to nodes on the tree. Produce a mapping of seen authors and the
    metadata associated with them.

    :param node_metadata:
    :type node_metadata: dict
    :param raw_strain_info:
    :type raw_strain_info: dict
    :returns: author information. See JSON schema for format.
    :rtype: dict
    """
    author_info = {}
    for strain, node in node_metadata.items():
        try:
            authors = raw_strain_info[strain]["authors"]
        except KeyError:
            continue # internal node / terminal node without authors

        data = {
            "author": authors
        }
        if "title" in raw_strain_info[strain]:
            data["title"] = raw_strain_info[strain]["title"].strip()
        if "journal" in raw_strain_info[strain]:
            data["journal"] = raw_strain_info[strain]["journal"].strip()
        if "paper_url" in raw_strain_info[strain]:
            data["paper_url"] = raw_strain_info[strain]["paper_url"].strip()

        # make unique key...
        key = authors.split()[0].lower()
        if "journal" in data:
            # extract the year out of the journal
            matches = re.findall(r'\([0-9A-Z-]*(\d{4})\)', data["journal"])
            if matches:
                key += matches[-1]
        if "title" in data:
            key += raw_strain_info[strain]["title"].strip().split()[0].lower()

        node["citekey"] = key

        if key not in author_info:
            author_info[key] = data
        elif author_info[key] != data:
            print("WARNING: Contradictory author information: {} vs {}".format(author_info[key], data))

    return author_info


def transfer_metadata_to_strains(strains, raw_strain_info, traits):
    node_metadata = {}
    for strain_name in strains:
        node = {"traits":{}}

        try:
            raw_data = raw_strain_info[strain_name]
        except KeyError:
            raise Exception("ERROR: {} is not found in the node_data[\"nodes\"]".format(strain_name))

        # TRANSFER MUTATIONS #
        if "aa_muts" in raw_data or "muts" in raw_data:
            node["mutations"] = {}
            if "muts" in raw_data and len(raw_data["muts"]):
                node["mutations"]["nuc"] = raw_data["muts"]
            if "aa_muts" in raw_data:
                aa = {gene:data for gene, data in raw_data["aa_muts"].items() if len(data)}
                node["mutations"].update(aa)
                #convert mutations into a label
                if aa:
                    aa_lab = '; '.join("{!s}: {!s}".format(key,', '.join(val)) for (key,val) in aa.items())
                    node["labels"] = { "aa": aa_lab }


        # TRANSFER NODE DATES #
        if "numdate" in raw_data or "num_date" in raw_data: # it's ok not to have temporal information
            node["num_date"] = {
                "value": raw_data["num_date"] if "num_date" in raw_data else raw_data["numdate"]
            }
            if "num_date_confidence" in raw_data:
                node["num_date"]["confidence"] = raw_data["num_date_confidence"]

        # TRANSFER VACCINE INFO #

        # TRANSFER LABELS #
        if "clade_annotation" in raw_data:
            if 'labels' in node:
                node['labels']['clade'] = raw_data["clade_annotation"]
            else:
                node["labels"] = { "clade": raw_data["clade_annotation"] }

        # TRANSFER NODE_HIDDEN PROPS #

        # TRANSFER AUTHORS #

        # TRANSFER GENERIC PROPERTIES #
        for prop in ["url", "accession"]:
            if prop in raw_data:
                node[prop] = raw_data[prop]

        # TRANSFER TRAITS (INCLUDING CONFIDENCE & ENTROPY) #
        for trait in traits:
            if trait in raw_data and raw_data[trait]:
                node["traits"][trait] = {"value": raw_data[trait]}
                if trait+"_confidence" in raw_data:
                    node["traits"][trait]["confidence"] = raw_data[trait+"_confidence"]
                if trait+"_entropy" in raw_data:
                    node["traits"][trait]["entropy"] = raw_data[trait+"_entropy"]

        node_metadata[strain_name] = node
    return node_metadata

def add_metadata_to_tree(node, metadata):
    node.update(metadata[node["strain"]])
    if "children" in node:
        for child in node["children"]:
            add_metadata_to_tree(child, metadata)

def get_traits(node_data):
    exclude = ['branch_length', 'num_date', 'raw_date', 'numdate', 'clock_length',
               'mutation_length', 'date', 'muts', 'aa_muts', 'sequence', 'aa_sequences']
    traits = []
    for seq, val in node_data['nodes'].items():
        for t in val.keys():

            if '_confidence' in t or '_entropy' in t:
                continue

            if t not in traits and t not in exclude:
                traits.append(t)

    return traits


def get_root_sequence(root_node, ref=None, translations=None):
    '''
    create a json structure that contains the sequence of the root, both as
    nucleotide and as translations. This allows look-up of the sequence for
    all states, including those that are not variable.

    Parameters
    ----------
    root_node : dict
    	data associated with the node
    ref : str, optional
        filename of the root sequence
    translations : str, optional
        file name of translations

    Returns
    -------
    dict
        dict of nucleotide sequence and translations
    '''
    root_sequence = {}
    if ref and translations:
        from Bio import SeqIO
        refseq = SeqIO.read(ref, 'fasta')
        root_sequence['nuc']=str(refseq.seq)
        for gene in SeqIO.parse(translations, 'fasta'):
            root_sequence[gene.id] = str(gene.seq)
    else:
        root_sequence["nuc"] = root_node["sequence"]
        root_sequence.update(root_node["aa_sequences"])

    return root_sequence


def convert_camel_to_snake_case(string):
    """
    Converts a string from camel case to snake case.

    This is used to allow existing `auspice-config` files for various builds
    to work with the v2 export schema.
    """
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', string)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()


def add_core_args(parser):
    core = parser.add_argument_group("REQUIRED")
    core.add_argument('--tree','-t', required=True, help="tree to perform trait reconstruction on")
    core.add_argument('--metadata', required=True, help="tsv file with sequence meta data")
    core.add_argument('--node-data', required=True, nargs='+', help="JSON files with meta data for each node")
    return core


def add_config_args(parser):
    config = parser.add_argument_group("CONFIG OPTIONS")
    # XXX TODO: make clear either use auspice-config or additional config options
    config.add_argument('--auspice-config', help="file with auspice configuration")
    config.add_argument('--title', default="Analysis", help="Title to be displayed by auspice")
    config.add_argument('--maintainers', default=[""], nargs='+', help="Analysis maintained by")
    config.add_argument('--maintainer-urls', default=[""], nargs='+', help="URL of maintainers")
    config.add_argument('--geography-traits', nargs='+', help="What location traits are used to plot on map")
    config.add_argument('--extra-traits', nargs='+', help="Metadata columns not run through 'traits' to be added to tree")
    config.add_argument('--panels', default=['tree', 'map', 'entropy'], nargs='+', help="What panels to display in auspice. Options are : xxx")
    return config

def add_option_args(parser):
    options = parser.add_argument_group("OPTIONS")
    options.add_argument('--colors', help="file with color definitions")
    options.add_argument('--lat-longs', help="file latitudes and longitudes, overrides built in mappings")
    options.add_argument('--tree-name', default=False, help="Tree name (needed for tangle tree functionality)")
    options.add_argument('--minify-json', action="store_true", help="export JSONs without indentation or line returns")
    options.add_argument('--output-sequence', help="JSON file name that is passed on to auspice (e.g., zika_seq.json).")
    options.add_argument('--reference', required=False, help="reference sequence for export to browser, only vcf")
    options.add_argument('--reference-translations', required=False, help="reference translations for export to browser, only vcf")
    return options


def register_arguments_v2(subparsers):
    v2 = subparsers.add_parser("v2", help="Export version 2 JSON schema")
    core = add_core_args(v2)
    config = add_config_args(v2)
    options = add_option_args(v2)
    core.add_argument('--output-main', help="Main JSON file name that is passed on to auspice (e.g., zika.json).")
    return v2


def run_v2(args):
    T = Phylo.read(args.tree, 'newick')
    node_data = read_node_data(args.node_data) # args.node_data is an array of multiple files (or a single file)
    nodes = node_data["nodes"] # this is the per-node metadata produced by various augur modules

    if args.minify_json:
        json_indent = None
    else:
        json_indent = 2

    # export reference sequence data including translations. This is either the
    # inferred sequence of the root, or the reference sequence with respect to
    # which mutations are made on the tree (including possible mutations leading
    # to the root of the tree -- typical case for vcf input data).
    if args.output_sequence:
        if T.root.name in nodes:
            root_sequence = get_root_sequence(nodes[T.root.name], ref=args.reference,
                                                translations=args.reference_translations)
        else:
            root_sequence = {}

        write_json(root_sequence, args.output_sequence)


    ## SCHEMA v2.0 ##
    auspice_json = {}
    config = {}
    auspice_json["version"] = "2.0"

    if args.auspice_config:
        config = read_config(args.auspice_config)
        auspice_json["title"] = config["title"]
        auspice_json["maintainers"] = [{ "name": config["maintainer"][0], "url": config["maintainer"][1]}]
        # Set up display defaults
        if config.get("defaults"):

            for key in config["defaults"]:
                new_key = convert_camel_to_snake_case(key)
                config["defaults"][new_key] = config["defaults"].pop(key)

            auspice_json["display_defaults"] = config["defaults"]

    else:
        auspice_json['title'] = args.title
        auspice_json['maintainers'] = [{'name': name, 'url':url} for name, url in zip(args.maintainers, args.maintainer_urls)]

    # get traits to colour by etc - do here before node_data is modified below
    # this ensures we get traits even if they are not on every node
    traits = get_traits(node_data)
    if config.get('color_options'):
        traits.extend(config['color_options'].keys())
    if args.extra_traits:
        traits.extend(args.extra_traits)
    # Automatically add any specified geo traits - otherwise won't work!
    if args.geography_traits:
        traits.extend(args.geography_traits)
        traits = list(set(traits)) #ensure no duplicates
    # remove keys which may look like traits but are not
    excluded_traits = [
      "clade_annotation", # Clade annotation is label, not colorby!
      "authors" # authors are set as a node property, not a trait property
    ]
    traits = [t for t in traits if t not in excluded_traits]

    raw_strain_info = collect_strain_info(node_data, args.metadata)
    auspice_json["tree"], strains = convert_tree_to_json_structure(T.root, raw_strain_info)

    node_metadata = transfer_metadata_to_strains(strains, raw_strain_info, traits)
    auspice_json["author_info"] = construct_author_info_and_make_keys(node_metadata, raw_strain_info)

    # Set up filters
    if config.get('filters'):
        auspice_json['filters'] = config['filters']
    else:
        # This check allows validation to complete ok - but check auspice can handle having no author info! (it can in v1 schema)
        if len(auspice_json["author_info"]) == 0:    # if no author data supplied
            del auspice_json["author_info"]
            auspice_json['filters'] = traits
        else:
            auspice_json['filters'] = traits + ['authors']

    add_metadata_to_tree(auspice_json["tree"], node_metadata)


    auspice_json["colorings"] = get_colorings(
        config=config,
        traits=traits,
        provided_colors=read_colors(args.colors),
        node_metadata=node_metadata,
        mutations_present=True # TODO
    )

    # Set up geographic info
    if config.get("geo"):
        geo_config = config["geo"]
    else:
        geo_config = args.geography_traits

    lat_long_mapping = read_lat_longs(args.lat_longs)
    auspice_json["geographic_info"] = process_geographic_info(geo_config, lat_long_mapping, node_metadata=node_metadata)

    auspice_json["updated"] = time.strftime('%Y-%m-%d')
    auspice_json["genome_annotations"] = process_annotations(node_data)
    auspice_json["panels"] = process_panels(args.panels, auspice_json)

    if args.tree_name:
        if not re.search("(^|_|/){}(_|.json)".format(args.tree_name), str(args.output_main)):
            print("Error: tree name {} must be found as part of the output string".format(args.tree_name))
            sys.exit(2)
        auspice_json["tree_name"] = args.tree_name

    write_json(auspice_json, args.output_main, indent=json_indent)
