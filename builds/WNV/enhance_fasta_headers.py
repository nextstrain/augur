from Bio import SeqIO
import os, sys, re
sys.path.append('../../../sacra') # sacra src directory
from src.entrez import query_genbank
from src.utils.genbank_parsers import choose_best_reference

fasta = {}
with open("data/full_dataset.singleline.aligned.pipeChar.fasta", "rU") as f:
    count = 0
    for record in SeqIO.parse(f, "fasta"):
        fasta[record.description.split('|')[0]] = record
        count += 1
        # if count > 30:
        #     break

## HOST & GENOTYPE INFO
hosts = {}
lineage = {}
with open("data/2018.04.08_WNV_headers_v2.csv", "rU") as f:
    for line in f:
        fields = line.strip().split(',')
        try:
            assert(len(fields) == 7)
        except AssertionError:
            print("Too many fields")
            print(line)
            sys.exit(2)
        hosts[fields[0]] = fields[5]
        lineage[fields[0]] = fields[6]

unknown = "Unknown"

for strain,record in fasta.iteritems():
    if strain not in hosts:
        print("{} not in hosts. setting to {}".format(strain, unknown))
        record.description += "|" + unknown
    else:
        record.description += "|" + hosts[strain]
    if strain not in lineage:
        print("{} not in lineage. setting to {}".format(strain, unknown))
        record.description += "|" + unknown
    else:
        record.description += "|" + lineage[strain]

## AUTHOR INFO. We want to add authors, journal, title, URL
accessions = fasta.keys()
gb = query_genbank([x for x in accessions if len(x) > 5]) # don't query obviously non-accession strings
refs = {accession: choose_best_reference(record) for accession, record in gb.iteritems()}
for accession,record in fasta.iteritems():
    if re.match(r'^W\d+$', accession):
        author = "Grubaugh et al"
        journal = "Unpublished"
        title = "West Nile virus genomic data from California"
        url = "https://andersen-lab.com/secrets/data/west-nile-genomics/"
    elif accession not in refs:
        print("{} has no genbank data. Setting as empty!".format(accession))
        record.description += "||||"
        continue
    else:
        reference = refs[accession]
        author = re.match(r'^([^,]*)', reference.authors).group(0) + " et al"
        journal = reference.journal
        title = reference.title

        if not reference.pubmed_id:
            print("no pubmed_id for ", author, title, "using accession")
            url = "https://www.ncbi.nlm.nih.gov/nuccore/" + accession
        else:
            url = "https://www.ncbi.nlm.nih.gov/pubmed/" + reference.pubmed_id


    record.description += "|{}|{}|{}|{}".format(author, journal, title, url)
    # print(record.description)


## clean for outp:
with open("data/WNV.fasta", "w") as fh:
    for strain,record in fasta.iteritems():
        fh.write(">{}\n".format(record.description))
        fh.write(str(record.seq))
        fh.write("\n")
