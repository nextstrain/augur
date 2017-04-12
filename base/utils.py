import matplotlib as mpl
from matplotlib import pyplot as plt

def generate_cmap(data, discrete):
    '''
    data is set or list (e.g. of countries)
    discrete is bool
    returns dict of (e.g.) country -> hex
    '''
    norm = mpl.colors.Normalize(0, len(data) - 1)
    cmap = mpl.cm.get_cmap("viridis")
    if discrete:
        if len(data) <= 10:
            cmap = mpl.cm.get_cmap("Vega10")
        elif len(data) <= 20:
            cmap = mpl.cm.get_cmap("Vega20")

    ret = {}
    for idx, val in enumerate(list(data)):
        ret[val] = mpl.colors.to_hex(cmap(norm(idx)))
    return ret
