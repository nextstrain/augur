"""
Export JSON files suitable for visualization with auspice.
"""

import os
import re
import time
import numpy as np
from Bio import Phylo
from collections import defaultdict
from .utils import read_metadata, read_node_data, write_json, read_config, read_lat_longs, read_colors

def convert_tree_to_json_structure(node, metadata, div=0, nextflu_schema=False, strains=None):
    """
    converts the Biopython tree structure to a dictionary that can
    be written to file as a json. This is called recursively.
    Creates the strain property & divergence on each node
    input
        node -- node for which top level dict is produced.
        div  -- cumulative divergence (root = 0)
        nextflu_schema -- use nexflu schema (e.g. node["attr"]). This is deprecated.

    returns
        tree in JSON structure
        list of strains
    """
    if nextflu_schema:
        node_struct = {
            'attr': {"div": div},
            'strain': node.name,
            'clade': node.clade
        }
        for attr in ['branch_length', 'tvalue', 'yvalue', 'xvalue']:
            try:
                node_struct[attr] = node.__getattribute__(attr)
            except AttributeError:
                pass
    else:
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
            node_struct["children"].append(convert_tree_to_json_structure(child, metadata, div=cdiv, nextflu_schema=nextflu_schema, strains=strains)[0])

    return (node_struct, strains)


def recursively_decorate_tree_json_nextflu_schema(node, node_metadata, decorations):
    """
    This function is deprecated and is used to produce the nextflu-compatable JSON format

    For given decorations, add information from node_metadata to
    each node in the tree.
    * decorations must have property "key" which is the key used to insert
    into the node and the default key used to access node_metadata
    * if decorations has property "lookup_key", this is used to access
    node meta_data instead
    * if decorations has property "is_attr" (and it's value is True)
    then the result is inserted into node["attr"]

    returns Null
    """
    try:
        metadata = node_metadata[node["strain"]]
        metadata["strain"] = node["strain"]
    except KeyError:
        raise Exception("ERROR: node %s is not found in the node metadata."%n.name)

    for data in decorations:
        val = None
        insert_key = data["key"]
        try:
            if "lookup_key" in data:
                val = metadata[data["lookup_key"]]
            else:
                val = metadata[insert_key]
        except KeyError:
            pass

        if val is not None:
            if "is_attr" in data and data["is_attr"]:
                node["attr"][insert_key] = val
            else:
                if insert_key == 'aa_muts':
                    val = {k:v for k,v in val.items() if len(v) }
                node[insert_key] = val

    if "children" in node:
        for child in node["children"]:
            recursively_decorate_tree_json_nextflu_schema(child, node_metadata, decorations)


def tree_layout(T):
    """
    calculate tree layout.
    This function is deprecated, and only used for the nextflu JSON format
    """
    yval=T.count_terminals()
    clade = 0;
    for n in T.find_clades(order='postorder'):
        n.clade=clade; clade+=1;
        if n.is_terminal():
            n.yvalue=yval
            yval-=1
        else:
            child_yvalues = [c.yvalue for c in n]
            n.yvalue=0.5*(np.min(child_yvalues)+np.max(child_yvalues))


def process_colorings(jsn, color_mapping, nodes=None, node_metadata=None, nextflu=False):
    if nextflu:
        if "color_options" not in jsn:
            print("WARNING: no color options were defined")
            return
        data = jsn["color_options"]
    else:
        #if "colorings" not in jsn:
        if jsn is None or len(jsn)==0:
            print("WARNING: no colorings were defined")
            return
        data = {k: {'type': 'categorical'} for k in jsn}
        #Always add in genotype and date to colour by
        data['gt'] = {'title': 'Genotype', 'type': 'ordinal'}
        #figure out how to see if num_date is there? Is it possible to not have?
        data['num_date']  = {'title': 'Sampling date', 'type': 'continuous'}

    for trait, options in data.items():
        if nextflu:
            if "legendTitle" not in options: options["legendTitle"] = trait
            if "menuItem" not in options: options["menuItem"] = trait
            if "key" not in options: options["key"] = trait
        else:
            if "title" not in options: options["title"] = trait
            if "type" not in options:
                raise Exception("coloring {} missing type...".format(trait))

        if trait.lower() in color_mapping:
            # remember that the color maps (from the TSV) are in lower case, but this is not how they should be exported
            if nodes:
                values_in_tree = {node[trait] for node in nodes.values() if trait in node}
            else:
                values_in_tree = {data["traits"][trait]["value"] for name, data in node_metadata.items()}
            case_map = {val.lower(): val for val in values_in_tree}

            if nextflu:
                options["color_map"] = [(case_map[m[0]], m[1]) for m in color_mapping[trait.lower()] if m[0] in case_map]
            else:
                options["scale"] = {case_map[m[0]]: m[1] for m in color_mapping[trait.lower()] if m[0] in case_map}

    return data

def process_geographic_info(jsn, lat_long_mapping, nextflu=False, node_metadata=None, nodes=None):
    if (nextflu and "geo" not in jsn) or (not nextflu and (jsn is None or len(jsn)==0)):
        return {}
    geo = defaultdict(dict)

    traits = jsn["geo"] if nextflu else jsn #jsn["geographic_info"]

    for trait in traits:
        if nextflu:
            demes_in_tree = {node[trait] for node in nodes.values() if trait in node}
        else:
            demes_in_tree = {data["traits"][trait]["value"] for name, data in node_metadata.items()}
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

def process_panels(user_panels, meta_json, nextflu=False):
    try:
        panels = meta_json["panels"]
    except KeyError:
        panels = ["tree", "map", "entropy"]

    if user_panels is not None and len(user_panels) != 0:
        panels = user_panels

    if nextflu:
        if "geo" in meta_json:
            geoTraits = meta_json["geo"].keys()
        else:
            geoTraits = []
        if "annotations" in meta_json:
            annotations = meta_json["annotations"].keys()
        else:
            annotations = []
    else:
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

def add_tsv_metadata_to_nodes(nodes, meta_tsv, meta_json, extra_fields=['authors', 'url', 'accession']):
    """
    Only used for nextflu schema compatability

    Add the relevent fields from meta_tsv to the nodes
    (both are dictionaries keyed off of strain names)
    * the relevent fields are found by scanning the meta json
    together with the extra_fields param
    """
    fields = [x for x in meta_json["color_options"].keys() if x != "gt"] + extra_fields
    if "geo" in meta_json:
        fields += meta_json["geo"]

    for strain, node in nodes.items():
        if strain not in meta_tsv:
            continue
        for field in fields:
            if field not in node and field in meta_tsv[strain] and meta_tsv[strain][field]:
                node[field] = meta_tsv[strain][field]


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


def construct_author_info_nexflu(metadata, tree, nodes):
    """
    author info maps the "authors" property present on tree nodes
    to further information about the paper etc
    """

    authorsInTree = set()
    for node in tree.find_clades(order='postorder'):
        if node.is_terminal and node.name in nodes and "authors" in nodes[node.name]:
            authorsInTree.add(nodes[node.name]["authors"])

    author_info = defaultdict(lambda: {"n": 0})
    for strain, data in metadata.items():
        if "authors" not in data:
            print("Error - {} had no authors".format(strain))
            continue
        if data["authors"] not in authorsInTree:
            continue
        authors = data["authors"]
        author_info[authors]["n"] += 1
        # add in extra attributes if they're present in the meta TSV (for this strain...)
        for attr in ["title", "journal", "paper_url"]:
            if attr in data:
                if attr in author_info[authors] and data[attr].strip() != author_info[authors][attr].strip():
                    print("Error - {} had contradictory {}(s): {} vs {}".format(authors, attr, data[attr], author_info[authors][attr]))
                author_info[authors][attr] = data[attr].strip()

    return author_info

def construct_author_info_and_make_keys(node_metadata, raw_strain_info):
    """
    For each node, gather the relevant author information (returned)
    and create a unique key inside node_metadata
    """
    author_info = {}
    for strain, node in node_metadata.items():
        try:
            authors = raw_strain_info[strain]["authors"]
        except KeyError:
            continue # internal node / terminal node without authors

        data = {
            "authors": authors
        }
        if "title" in raw_strain_info[strain]:
            data["title"] = raw_strain_info[strain]["title"].strip()
        if "journal" in raw_strain_info[strain]:
            data["journal"] = raw_strain_info[strain]["journal"].strip()
        if "url" in raw_strain_info[strain]:
            data["url"] = raw_strain_info[strain]["paper_url"].strip()

        # make unique key...
        key = authors.split()[0].lower()
        if "journal" in data:
            # extract the year out of the journal
            matches = re.findall(r'\([0-9A-Z-]*(\d{4})\)', data["journal"])
            if matches:
                key += matches[-1]
        if "title" in data:
            key += raw_strain_info[strain]["title"].strip().split()[0].lower()

        node["authors"] = key

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
            raise Exception("ERROR: %s is not found in the node_data[\"nodes\"]"%strain_name)

        # TRANSFER MUTATIONS #
        if "aa_muts" in raw_data or "muts" in raw_data:
            node["mutations"] = {}
            if "muts" in raw_data and len(raw_data["muts"]):
                node["mutations"]["nuc"] = raw_data["muts"]
            if "aa_muts" in raw_data:
                aa = {gene:data for gene, data in raw_data["aa_muts"].items() if len(data)}
                node["mutations"].update(aa)

        # TRANSFER NODE DATES #
        if "numdate" in raw_data or "num_date" in raw_data: # it's ok not to have temporal information
            node["num_date"] = {
                "value": raw_data["num_date"] if "num_date" in raw_data else raw_data["numdate"]
            }
            if "num_date_confidence" in raw_data:
                node["num_date"]["confidence"] = raw_data["num_date_confidence"]

        # TRANSFER VACCINE INFO #

        # TRANSFER LABELS #

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
        newT = [t for t in list(val.keys()) if t not in traits and t not in exclude]
        traits.extend(newT)
    traits = [x for x in traits if '_confidence' not in x and '_entropy' not in x]

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


def register_arguments(parser):
    parser.add_argument('--tree', '-t', required=True, help="tree to perform trait reconstruction on")
    parser.add_argument('--metadata', required=True, help="tsv file with sequence meta data")
    parser.add_argument('--reference', required=False, help="reference sequence for export to browser, only vcf")
    parser.add_argument('--reference-translations', required=False, help="reference translations for export to browser, only vcf")
    parser.add_argument('--node-data', required=True, nargs='+', help="JSON files with meta data for each node")
    parser.add_argument('--auspice-config', help="file with auspice configuration")
    parser.add_argument('--colors', help="file with color definitions")
    parser.add_argument('--lat-longs', help="file latitudes and longitudes, overrides built in mappings")
    parser.add_argument('--new-schema', action="store_true", help="export JSONs using nexflu schema")
    parser.add_argument('--output-main', help="Main JSON file name that is passed on to auspice (e.g., zika.json).")
    parser.add_argument('--output-tree', help="JSON file name that is passed on to auspice (e.g., zika_tree.json). Only used with --nextflu-schema")
    parser.add_argument('--output-sequence', help="JSON file name that is passed on to auspice (e.g., zika_seq.json). Only used with --nextflu-schema")
    parser.add_argument('--output-meta', help="JSON file name that is passed on to auspice (e.g., zika_meta.json). Only used with --nextflu-schema")
    parser.add_argument('--title', default="Analysis", help="Title to be displayed by auspice")
    parser.add_argument('--maintainers', default=[""], nargs='+', help="Analysis maintained by")
    parser.add_argument('--maintainer-urls', default=[""], nargs='+', help="URL of maintainers")
    parser.add_argument('--geography-traits', nargs='+', help="What location traits are used to plot on map")
    parser.add_argument('--extra-traits', nargs='+', help="Metadata columns not run through 'traits' to be added to tree")
    parser.add_argument('--panels', default=['tree', 'map', 'entropy'], nargs='+', help="What panels to display in auspice. Options are : xxx")
    parser.add_argument('--minify-json', action="store_true", help="export JSONs without indentation or line returns")


def run(args):
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

    if not args.new_schema:
        # This schema is deprecated. It remains because:
        # (1) auspice can't use schema 2.0 yet, (2) nexflu doesn't use schema 2.0
        # export the tree JSON first
        meta_json = read_config(args.auspice_config)
        meta_tsv, _ = read_metadata(args.metadata)
        add_tsv_metadata_to_nodes(nodes, meta_tsv, meta_json)

        tree_layout(T)
        tree_json, _ = convert_tree_to_json_structure(T.root, nodes, nextflu_schema=True)

        # now the messy bit about what decorations (e.g. "country", "aa_muts") do we want to add to the tree?
        # see recursively_decorate_tree_json to understand the tree_decorations structure
        tree_decorations = [
            {"key": "num_date", "lookup_key": "numdate", "is_attr": True},
            {"key": "muts", "is_attr": False},
            {"key": "aa_muts", "is_attr": False}
        ]
        traits_via_node_metadata = {k for node in nodes.values() for k in node.keys()}
        traits_via_node_metadata -= {'sequence', 'mutation_length', 'branch_length', 'numdate',
                                     'mutations', 'muts', 'aa_muts', 'aa_sequences'}
        for trait in traits_via_node_metadata:
            tree_decorations.append({"key": trait, "is_attr": True})

        recursively_decorate_tree_json_nextflu_schema(tree_json, nodes, decorations=tree_decorations)
        write_json(tree_json, args.output_tree, indent=json_indent)

        # Export the metadata JSON
        lat_long_mapping = read_lat_longs(args.lat_longs)
        color_mapping = read_colors(args.colors)
        meta_json["updated"] = time.strftime("%d %b %Y")
        meta_json["virus_count"] = len(list(T.get_terminals()))
        meta_json["author_info"] = construct_author_info_nexflu(meta_tsv, T, nodes)
        meta_json["color_options"] = process_colorings(meta_json, color_mapping, nodes=nodes, nextflu=True)
        meta_json["geo"] = process_geographic_info(meta_json, lat_long_mapping, nodes=nodes, nextflu=True)
        annotations = process_annotations(node_data)
        if annotations:
            meta_json["annotations"] = annotations
        meta_json["panels"] = process_panels(None, meta_json, nextflu=True)

        write_json(meta_json, args.output_meta, indent=json_indent)
        return 0

    ## SCHEMA v2.0 ##
    unified = {}

    unified['title'] = args.title
    unified['maintainers'] = [{'name': name, 'href':url} for name, url in zip(args.maintainers, args.maintainer_urls)]
    unified["version"] = "2.0"

    # get traits to colour by etc - do here before node_data is modified below
    # this ensures we get traits even if they are not on every node
    traits = get_traits(node_data)
    if args.extra_traits:
        traits.extend(args.extra_traits)

    raw_strain_info = collect_strain_info(node_data, args.metadata)
    unified["tree"], strains = convert_tree_to_json_structure(T.root, raw_strain_info)

    node_metadata = transfer_metadata_to_strains(strains, raw_strain_info, traits)
    unified["author_info"] = construct_author_info_and_make_keys(node_metadata, raw_strain_info)

    # This check allows validation to complete ok - but check auspice can handle having no author info! (it can in v1 schema)
    if len(unified["author_info"]) == 0:    # if no author data supplied
        del unified["author_info"]
        unified['filters'] = traits
    else:
        unified['filters'] = traits + ['authors']

    add_metadata_to_tree(unified["tree"], node_metadata)

    color_mapping = read_colors(args.colors)
    unified["colorings"] = process_colorings(traits, color_mapping, node_metadata=node_metadata)

    lat_long_mapping = read_lat_longs(args.lat_longs)
    unified["geographic_info"] = process_geographic_info(args.geography_traits, lat_long_mapping, node_metadata=node_metadata)

    unified["updated"] = time.strftime('%Y-%m-%d')
    unified["genome_annotations"] = process_annotations(node_data)
    unified["panels"] = process_panels(args.panels, unified)

    write_json(unified, args.output_main, indent=json_indent)
