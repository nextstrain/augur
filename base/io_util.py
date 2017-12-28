from __future__ import division, print_function
import Bio.Phylo

def myopen(fname, mode='r'):
    if fname[-2:] == 'gz':
        from gzip import open as gopen
        return gopen(fname, mode)
    else:
        return open(fname, mode)

def make_dir(dname):
    import os
    if not os.path.isdir(dname):
        try:
            os.makedirs(dname)
        except OSError as e:
            print("Cannot create run_dir",e)

def remove_dir(dname, max_attempts=5):
    import os, shutil, time
    if os.path.isdir(dname):
        # Try to remove the given directory repeatedly to compensate for NFS
        # latency that can result in OSError exceptions when a directory appears
        # not to be empty when shutil attempts to remove it.
        for i in xrange(max_attempts):
            try:
                shutil.rmtree(dname)
                return
            except OSError, e:
                time.sleep(i)

        # Try one last time and let OS exception propagate up.
        shutil.rmtree(dname)

def write_json(data, file_name, indent=1):
    import json
    try:
        handle = open(file_name, 'w')
    except IOError:
        pass
    else:
        json.dump(data, handle, indent=indent)
        handle.close()

def tree_to_json(node, extra_attr = []):
    tree_json = {}
    str_attr = ['strain']
    num_attr = ['xvalue', 'yvalue', 'tvalue', 'num_date']
    if hasattr(node, 'name'):
        tree_json['strain'] = node.name

    for prop in str_attr:
        if hasattr(node, prop):
            tree_json[prop] = node.__getattribute__(prop)
    for prop in num_attr:
        if hasattr(node, prop):
            try:
                tree_json[prop] = round(node.__getattribute__(prop),5)
            except:
                print("cannot round:", node.__getattribute__(prop), "assigned as is")
                tree_json[prop] = node.__getattribute__(prop)

    for prop in extra_attr:
        if len(prop)==2 and callable(prop[1]):
            if hasattr(node, prop[0]):
                tree_json[prop] = prop[1](node.__getattribute__(prop[0]))
        else:
            if hasattr(node, prop):
                tree_json[prop] = node.__getattribute__(prop)

    if node.clades:
        tree_json["children"] = []
        for ch in node.clades:
            tree_json["children"].append(tree_to_json(ch, extra_attr))
    return tree_json


def json_to_tree(json_dict):
    """Returns a Bio.Phylo tree corresponding to the given JSON dictionary exported
    by `tree_to_json`.

    >>> import json
    >>> json_fh = open("tests/json_tree_to_nexus/flu_h3n2_ha_3y_tree.json", "r")
    >>> json_dict = json.load(json_fh)
    >>> tree = json_to_tree(json_dict)
    >>> tree.name
    u'NODE_0002020'
    >>> len(tree.clades)
    2
    >>> tree.clades[0].name
    u'NODE_0001489'
    >>> hasattr(tree, "attr")
    True
    >>> "dTiter" in tree.attr
    True
    """
    node = Bio.Phylo.Newick.Clade()
    node.name = json_dict["strain"]

    if "children" in json_dict:
        # Recursively add children to the current node.
        node.clades = [json_to_tree(child) for child in json_dict["children"]]

    # Assign all non-children attributes.
    for attr, value in json_dict.iteritems():
        if attr != "children":
            setattr(node, attr, value)

    return node
