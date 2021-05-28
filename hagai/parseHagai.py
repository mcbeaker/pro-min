import pandas as pd
from Bio import PDB
from progress.bar import Bar
import os

f = "/home/kenneth/proj/proXtal/proteins/hagai/hagai.csv"
pdir = "/home/kenneth/proj/proXtal/proteins/hagai/pdb"
df = pd.read_csv(f)

feDF = df[df.Metal == "FE"]

pdbs = feDF["Microenvironment ID"].apply(lambda x: x[0:4]).sort_values().drop_duplicates()

# pdbs.to_csv(os.path.join(pdir,"pdbList.csv"),index=False)
bar = Bar("Processing",max=len(pdbs),fill='*',suffix='%(perccent).1f%% - %(eta)ds')

for i in pdbs:
    PDB.PDBList(verbose=False).retrieve_pdb_file(pdb_code=i,file_format="mmCif",pdir=pdir,obsolete=True)
    bar.next()
bar.finish()