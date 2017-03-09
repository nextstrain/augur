import json
from pprint import pprint
import pdb
import os

metaFiles = ['build/flu_H7N9_HA_meta.json', 'build/flu_H7N9_NA_meta.json'];
for metaFile in metaFiles:
    os.system("cp "+metaFile+" "+metaFile+".bak")

    with open(metaFile) as data_file:
        data = json.load(data_file)

    # print what we've got so far...
    for key in data["color_options"]:
        if "color_map" in data["color_options"][key]:
            print key
            print "\t" + "\n\t".join([x[0] for x in data["color_options"][key]["color_map"]])

    cmap = {
        "country": {
            "canada": "#ffffbf",
            "hong_kong": "#66c2a5",
            "china": "#d53e4f",
            "japan": "#3288bd",
            "taiwan": "#fdae61",
            "?": "#FFFFFF"
        },
        "host": {
            "avian": "#5e4fa2",
            "gallusgallus": "#3288bd",
            "passermontanus": "#e6f598",
            "goose": "#abdda4",
            "otheravian": "#66c2a5",
            "human": "#d53e4f",
            "duck": "#5e4fa2",
            "chicken": "#3288bd",
            "laboratoryderived": "#fee08b"
        },
        "division": {
            "?": "#FFFFFF",
            "beijing": "#d53e4f",
            "shandong": "#f46d43",
            "hunan": "#66c2a5",
            "guangxi": "#5e4fa2",
            "hebei": "#f46d43",
            "japan": "#4d4d4d",
            "taiwan": "#3288bd",
            "guangdong": "#5e4fa2",
            "jilin": "#c51b7d",
            "anhui": "#fee08b",
            "shanghai": "#ffffbf",
            "british_columbia": "#4d4d4d",
            "jiangsu": "#fdae61",
            "hubei": "#66c2a5",
            "zhejiang": "#e6f598",
            "xinjiang": "#d53e4f",
            "hong_kong": "#5e4fa2",
            "henan": "#f46d43",
            "guizhou": "#66c2a5",
            "fujian": "#abdda4",
            "jiangxi": "#66c2a5"
        }
    }



    for key in data["color_options"]:
        if "color_map" in data["color_options"][key] and key in cmap:
            print "changing color for " + key
            for idx in xrange(0, len(data["color_options"][key]["color_map"])):
                if (data["color_options"][key]["color_map"][idx][0] in cmap[key]):
                    data["color_options"][key]["color_map"][idx][1] = cmap[key][data["color_options"][key]["color_map"][idx][0]]

    with open(metaFile, 'w') as f:
        json.dump(data, f, indent=2, sort_keys=True)
