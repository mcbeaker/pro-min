from bioservices.uniprot import UniProt #imports uniprot 
import pandas as pd 

u = UniProt(verbose=True)
for i, row in key_value.iterrows():
    # print(row.pdbID)
    map = u.mapping("PDB_ID","ACC", row.pdbID)
    # print(map)
    if(len(map[row.pdbID])>1):

        for ac in map[row.pdbID]:
            print(map)
            print(ac)
            fasta = u.get_fasta_sequence(ac)
            print(fasta)