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
import warnings

# Set up warnings & exceptions
warn = warnings.warn
deprecationWarningsEmitted = False

def deprecated(message):
    warn(message, DeprecationWarning, stacklevel=2)
    global deprecationWarningsEmitted
    deprecationWarningsEmitted=True


def fatal(message):
    print("FATAL ERROR: {}".format(message))
    sys.exit(2)

class InvalidOption(Exception):
    pass

def customformatwarning(message, category, filename, lineno, line=None):
    if category.__name__ == "UserWarning":
        return "WARNING: {}\n\n".format(message)
    if category.__name__ == "DeprecationWarning":
        return "DEPRECATED: {}\n\n".format(message)
    return "{}\n".format(message)

warnings.formatwarning = customformatwarning
warnings.simplefilter("default") # show DeprecationWarnings by default


def convert_tree_to_json_structure(node, metadata, div=0, strains=None):
    """
    converts the Biopython tree structure to a dictionary that can
    be written to file as a json. This is called recursively.
    Creates the strain property & divergence on each node

    input
        node -- node for which top level dict is produced.
        div  -- cumulative divergence (root = 0). False → divergence won't be exported.

    returns
        tree in JSON structure
        list of strains
    """

    # Does the tree have divergence? (BEAST trees may not)
    # only calculate this for the root node!
    if div == 0 and 'mutation_length' not in metadata[node.name] and 'branch_length' not in metadata[node.name]:
        div = False

    node_struct = {'strain': node.name}
    if div is not False: # div=0 is ok
        node_struct["div"] = div

    if strains is None:
        strains = [node_struct["strain"]]
    else:
        strains.append(node_struct["strain"])

    if node.clades:
        node_struct["children"] = []
        for child in node.clades:
            if div is False:
                cdiv=False
            else:                
                if 'mutation_length' in metadata[child.name]:
                    cdiv = div + metadata[child.name]['mutation_length']
                elif 'branch_length' in metadata[child.name]:
                    cdiv = div + metadata[child.name]['branch_length']
            node_struct["children"].append(convert_tree_to_json_structure(child, metadata, div=cdiv, strains=strains)[0])

    return (node_struct, strains)


def get_colorings(config, traits, provided_colors, node_metadata, mutations_present):
    def _rename_authors_key(color_config):
        if not color_config.get("authors"):
            return
        if color_config.get("author"):
            warn("[config file] Both 'author' and 'authors' supplied as coloring options. Using 'author'.")
        deprecated("[config file] The 'authors' key is now called 'author'")
        color_config["author"] = color_config["authors"]
        del color_config["authors"]

    def _get_values(key):
        values_in_tree = set()
        for node_properties in node_metadata.values():
            # beware of special cases!
            if key in ["num_date", "author"]:
                values_in_tree.add(node_properties.get(key, {}).get("value"))
            elif key in node_properties['traits']:
                values_in_tree.add(node_properties['traits'][key]['value'])
        return values_in_tree

    def _get_type(key, config_data, trait_values):
        if config_data.get("type"):
            t = config_data.get("type")
            allowedTypes = ["continuous", "ordinal", "categorical", "boolean"]
            if t == "discrete":
                deprecated("[config file] Coloring types of 'discrete' are no longer valid. Please use one of {} instead. {} has been automatically set as 'categorical'".format(", ".join(allowedTypes), key))
                return "categorical"
            if t not in allowedTypes:
                warn("[config file] {}'s type: '{}' is not valid. Please choose from {}. This has been excluded!".format(key, t, ", ".join(allowedTypes)))
                raise InvalidOption()
            return t
        # no type supplied => try to guess
        if all([ isinstance(n, float) if isinstance(n, float) else isinstance(n, int) for n in trait_values ]):
            t = "continuous"
        else:
            t = "categorical"
        warn("[config file] {} was missing type information. We've guessed '{}'.".format(key, t))
        return t

    def _get_title(key, color_config):
        # preferentially get the title from the color_config if set
        if color_config.get(key):
            info = color_config.get(key)
            if "title" in info:
                return info["title"]
            if "menuItem" in info:
                deprecated("[config file] 'meunItem' has been replaced with 'title' (coloring '{}')".format(key))
                return info["menuItem"]
            if "legendTitle" in info:
                deprecated("[config file] 'legendTitle' has been replaced with 'title' (coloring '{}')".format(key))
                return info["legendTitle"]
        # hardcoded fallbacks:
        if key == "clade_membership":
            return "Clade"
        if key == "gt":
            return "Genotype"
        if key == "authors":
            return "Authors"
        if key == 'num_date':
            return 'Sampling date'
        # fallthrough
        return key

    # TODO: sort out how command line arguements play with colof_config, if defined
    if config.get("colorings"):
        color_config = config["colorings"]
    elif config.get("color_options"):
        color_config = config["color_options"]
        deprecated("[config file] 'color_options' has been replaced with 'colorings'")
    else:
        color_config = {t: {} for t in traits}

    # if the 'clade_membership' trait is defined then ensure it's also a coloring
    if 'clade_membership' in traits and 'clade_membership' not in color_config:
        color_config['clade_membership'] = {}

    # handle deprecated keys by updating to their new ones where possible
    _rename_authors_key(color_config)

    colorings = {}
    # handle special cases
    if mutations_present:
        colorings["gt"] = {'title': _get_title("gt", color_config), 'type': 'ordinal'}
    if _get_values("author"): # check if any nodes have author set
        colorings["author"] = {'title': _get_title("author", color_config), 'type': 'categorical'}
    if _get_values("num_date"): # TODO: check if tree has temporal inference (possible to not have)
        colorings['num_date'] = {'title': _get_title("num_date", color_config), 'type': 'continuous'}
    # remove these keys so they're not processed below
    color_config.pop("gt", None)
    color_config.pop("author", None)
    color_config.pop("num_date", None)
    color_config.pop("div", None) # makes no sense to have this

    # loop through provided colorings not handled above
    for key, config_data in color_config.items():
        trait_values = _get_values(key) # e.g. list of countries, regions etc
        if not trait_values:
            warn("you asked for a color by for '{}' but it has no values on the tree => Ignoring.".format(key))
            continue
        try:
            colorings[key] = {"title": _get_title(key, color_config), "type": _get_type(key, config_data, trait_values)}
        except InvalidOption:
            continue # a warning message will have been printed before this is thrown
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
            # deme may be numeric, or string
            try:
                deme_search_value = deme.lower()
            except AttributeError:
                deme_search_value = str(deme)

            try:
                geo[trait][deme] = lat_long_mapping[(trait.lower(), deme_search_value)]
            except KeyError:
                warn("{}->{} did not have an associated lat/long value (matching performed in lower case). Auspice won't be able to display this deme.".format(trait, deme))
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
    if meta_json.get("genome_annotations", None):
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
    if tsv_path:
        meta_tsv, _ = read_metadata(tsv_path) # meta_tsv={} if file isn't read / doeesn't exist

        for strain, node in strain_info.items():
            if strain in meta_tsv:
                for field in meta_tsv[strain]:
                    node[field] = meta_tsv[strain][field]
    return strain_info


def set_author_on_nodes(node_metadata, raw_strain_info):
    """Gather the authors which appear in the metadata and assign them
    to nodes on the tree.

    :param node_metadata:
    :type node_metadata: dict
    :param raw_strain_info:
    :type raw_strain_info: dict
    :returns: None
    :rtype: None
    """
    author_info = {}
    seen = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    for strain, node in node_metadata.items():
        author = raw_strain_info[strain].get("author")
        if not author:
            author = raw_strain_info[strain].get("authors")
        if not author:
            continue # internal node / terminal node without authors

        node["author"] = {"author": author}        
        if "title" in raw_strain_info[strain]:
            node["author"]["title"] = raw_strain_info[strain]["title"].strip()
        if "journal" in raw_strain_info[strain]:
            node["author"]["journal"] = raw_strain_info[strain]["journal"].strip()
        if "paper_url" in raw_strain_info[strain] and not raw_strain_info[strain]["paper_url"].strip("/").endswith("pubmed"):
            node["author"]["paper_url"] = raw_strain_info[strain]["paper_url"].strip()

        # add to `seen` which will later be used to create the unique value which auspice will display
        year_matches = re.findall(r'\([0-9A-Z-]*(\d{4})\)', node["author"].get("journal", ""))
        year = str(year_matches[-1]) if year_matches else "unknown"
        seen[author][year][node["author"].get("title", "unknown")].append(node)

    # turn "seen" into a unique "nice" string for auspice to display
    for author in seen.keys():
        for year in seen[author].keys():
            titles = sorted(seen[author][year].keys())
            for idx, title in enumerate(titles):
                value = author.split()[0].lower().capitalize()
                if year != "unknown":
                    value += " ({})".format(year)
                if len(titles) > 1:
                    value += " {}".format("abcdefghij"[idx])
                for node in seen[author][year][title]:
                    node["author"]["value"] = value

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

        # TRANSFER NODE.HIDDEN PROPS #
        hidden = raw_data.get("hidden", None)
        if hidden:
            if hidden in ["always", "divtree", "timetree"]:
                node["hidden"] = hidden
            elif hidden is True or str(hidden) == "1": # interpret this as hidden in both div + time tree
                node["hidden"] = "always"
            else:
                warn("Hidden node trait of {} is invalid. Ignoring.".format(hidden))

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
               'mutation_length', 'date', 'muts', 'aa_muts', 'sequence', 'aa_sequences',
               'hidden']
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


def register_arguments_v2(subparsers):
    v2 = subparsers.add_parser("v2", help="Export version 2 JSON schema")

    required = v2.add_argument_group(
        title="REQUIRED"
    )
    required.add_argument('--tree','-t', metavar="newick", required=True, help="Phylogenetic tree, usually output from `augur refine`")
    required.add_argument('--node-data', metavar="JSON", required=True, nargs='+', help="JSON files containing metadata for nodes in the tree")
    required.add_argument('--output-main', metavar="JSON", required=True, help="Ouput file for auspice")

    config = v2.add_argument_group(
        title="CONFIG OPTIONS",
        description="These control the display settings for auspice. \
            You can supply a config JSON (which has all available options) or command line arguments (which are more limited but great to get started). \
            Supplying both is fine too -- command line args will overrule what was set in the config file!"
    )
    config.add_argument('--auspice-config', metavar="JSON", help="Auspice configuration")
    config.add_argument('--title', nargs='+', help="Title to be displayed by auspice")
    config.add_argument('--maintainers', metavar="name", nargs='+', help="Analysis maintained by")
    config.add_argument('--maintainer-urls', metavar="url", nargs='+', help="URL of maintainers")
    config.add_argument('--geography-traits', metavar="trait", nargs='+', help="What location traits are used to plot on map")
    config.add_argument('--extra-traits', metavar="trait", nargs='+', help="Metadata columns not run through 'traits' to be added to tree")
    config.add_argument('--panels', default=['tree', 'map', 'entropy'], nargs='+', help="Restrict panel display in auspice. Options are %(default)s. Ignore this option to display all available panels.")

    optional_inputs = v2.add_argument_group(
        title="OPTIONAL INPUTS"
    )
    optional_inputs.add_argument('--metadata', metavar="TSV", help="Additional metadata for strains in the tree")
    optional_inputs.add_argument('--colors', metavar="TSV", help="Custom color definitions")
    optional_inputs.add_argument('--lat-longs', metavar="TSV", help="Latitudes and longitudes for geography traits (overrides built in mappings)")
    optional_inputs.add_argument('--reference', metavar="JSON", required=False, help="reference sequence for export to browser, only vcf")
    optional_inputs.add_argument('--reference-translations', metavar="???", required=False, help="reference translations for export to browser, only vcf")

    optional_settings = v2.add_argument_group(
        title="OPTIONAL SETTINGS"
    )
    optional_settings.add_argument('--tree-name', metavar="name", default=False, help="Tree name (needed for tangle tree functionality)")
    optional_settings.add_argument('--minify-json', action="store_true", help="export JSONs without indentation or line returns")

    optional_outputs = v2.add_argument_group(
        title="OPTIONAL OUTPUTS"
    )
    optional_outputs.add_argument('--output-sequence', metavar="JSON", help="(reconstructed) sequences for each node")
    
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

        # v2 schema uses 'display_defaults' not 'defaults' (v1).
        if config.get("defaults"):
            deprecated("[config file] 'defaults' has been replaced with 'display_defaults'")
            config["display_defaults"] = config["defaults"]

        # Set up display defaults
        if config.get("display_defaults"):

            for key in config["display_defaults"]:
                new_key = convert_camel_to_snake_case(key)
                config["display_defaults"][new_key] = config["display_defaults"].pop(key)

            auspice_json["display_defaults"] = config["display_defaults"]

    # Get title - config file overwrites command-line args.
    if config.get("title"):
        auspice_json['title'] = config['title']
    elif args.title:
        auspice_json['title'] = ' '.join(args.title)

    # Get maintainers. Config file overwrites command-line args.
    # Recognises and implements v1-style config spec without warnings.
    if config.get("maintainer"): #v1-type specification
        auspice_json["maintainers"] = [{ "name": config["maintainer"][0], "url": config["maintainer"][1]}]
    elif config.get("maintainers"): #v2-type specification (proposed by Emma)
        auspice_json['maintainers'] = [{'name': n[0], 'url': n[1]} for n in config['maintainers']]
    elif args.maintainers and args.maintainer_urls:
        auspice_json['maintainers'] = [{'name': name, 'url':url} for name, url in zip(args.maintainers, args.maintainer_urls)]

    # get traits to colour by etc - do here before node_data is modified below
    # this ensures we get traits even if they are not on every node
    traits = get_traits(node_data)
    if config.get('colorings'):
        traits.extend(config['colorings'].keys())
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
    set_author_on_nodes(node_metadata, raw_strain_info)

    # Set up filters
    if config.get('filters'):
        auspice_json['filters'] = config['filters']
        if "authors" in auspice_json['filters']:
            del auspice_json['filters'][auspice_json['filters'].index("authors")]
            auspice_json['filters'].append("author")

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

    # Set up panels for both config and command-line
    if config.get("panels"):
        panels = config["panels"]
    else:
        panels = args.panels
    auspice_json["panels"] = process_panels(panels, auspice_json)

    if args.tree_name:
        if not re.search("(^|_|/){}(_|.json)".format(args.tree_name), str(args.output_main)):
            fatal("tree name {} must be found as part of the output string".format(args.tree_name))
        auspice_json["tree_name"] = args.tree_name

    write_json(auspice_json, args.output_main, indent=json_indent)

    if deprecationWarningsEmitted:
        print("------------------------")
        print("There were deprecation warnings displayed. They have been fixed but these will likely become breaking errors in a future version of augur.")
        print("------------------------")
    print("")