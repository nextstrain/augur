from augur.utils import load_features
# NOTE - this is from the installed augur, not the (potentially modified) source copy
# unless you are using an editable pip install
import argparse

parser = argparse.ArgumentParser(description="Test reading of a provided annotation (reference) file")
parser.add_argument('input', type=str, help="Reference file (GFF / GenBank)")
args = parser.parse_args()

features = load_features(args.input)
print("The following features were extracted from", args.input, "\n")
print(f"{'name':10}{'type':10}{'location':20}{'qualifiers found'}")
for name, data in features.items():
  print(f"{name:10}{data.type:10}{str(data.location):20}"+", ".join([k for k in data.qualifiers.keys()]))
  # print(data.qualifiers)



# Looks at feature qualifiers (in order) ["gene"], then "gene_name" then "locus_tag" FOR EACH entry. Last entry wins.

