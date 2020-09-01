from Bio import Entrez
import re as re

# phrase = "photosystem II"
phrase = 'redox potential'
Entrez.email = "kenneth.mcguinness@gmail.com"  # Always tell NCBI who you are

handle = Entrez.esearch(db="pubmed", term=phrase, retmax=5000)
record = Entrez.read(handle)
handle.close()
idlist = record["IdList"]
# print(idlist)

from Bio import Medline
handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline",retmode="text")
records = Medline.parse(handle)

for record in records:
    # print(dir(record))
    # print(record.keys())
    r = re.compile('protein.+mV')
    rd = re.compile('redox')
    if r.search(record.get('AB','?')) and rd.search(record.get('AB','?')):
        # print("title:", record.get("TI", "?"))

        # print(record.get('AB','?'))
        print("source:", record.get("SO", "?"))
        print("")

    # print(found)

    # print("abstract",record.get("AB","?"))
    # print("authors:", record.get("AU", "?"))
    # exit()
    # exit()
handle.close()
