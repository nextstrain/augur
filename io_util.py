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
            print "Cannot create run_dir",e

def remove_dir(dname):
    import os, shutil
    if os.path.isdir(dname):
        import shutil
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
    str_attr = ['country','region','clade','strain', 'date', 'muts']
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
                print "cannot round:", node.__getattribute__(prop), "assigned as is"
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
