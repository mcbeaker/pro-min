# import read_cif
import os
from CifFile import ReadCif
import pandas as pd
import numpy as np
from Bio.PDB import MMCIF2Dict
from Bio.PDB.PDBIO import Select
from Bio import PDB
from diffpy import structure 
import gemmi
from progress.bar import Bar
import subprocess
import inspect
import re
import multiprocessing as mp
# from elements import ELEMENTS
from collections import defaultdict
import json
from tqdm import tqdm

# strElements = [element.symbol.upper() for element in ELEMENTS]

# import urllib.request

import glob

#global variables
f = "/home/kenneth/proj/proXtal/proteins/hagai/hagai.csv"
pdir = "/home/kenneth/proj/proXtal/proteins/hagai/pdb"

global hagaiNames 
hagaiNames = ["Microenvironment_ID","Binding_motif","Cofactor","Metal","cofactor_group","EC","Head","Molecule","Organism_scientific","no_rank","superkingdom","phylum","class_","order","family","genus","organism_taxid","name","chains","Resolution","Structure_method","Keywords","Journal_reference","Release_date"]
df = pd.read_csv(f,names=hagaiNames,index_col=False)
logFile = [line.rstrip() for line in open('/home/kenneth/proj/proXtal/analyzeHagai.log').readlines()]
global feDF 

feDF = df[df['Metal'].str.contains("FE") & df['Structure_method'].str.contains("x-ray") & ~df['Microenvironment_ID'].isin(logFile)]

class SlowBar(Bar):
    suffix = '%(index)d/%(max)d - %(remaining_hours)d hours remaining'
    @property
    def remaining_hours(self):
        return self.eta // 3600

# global pbar
# pbar = tqdm(total=len(feDF.index))
# pbar = tqdm(total=10)

# bar = SlowBar("Processing",max=len(feDF.index),fill='*') #,suffix='%(percent).1f%% - %(eta)ds')


# # if os.path.exists(cifPath) == False:
# #     print(pdb + " needed")
# #     urllib.request.urlretrieve("https://files.rcsb.org/download/"+pdb.upper()+".cif",os.path.join(pdir,pdb+".cif"))

def get_cif_STR(cifPath):
        STR = PDB.MMCIFParser(QUIET=True).get_structure("pdb",cifPath)
        # DICT = PDB.MMCIF2Dict.MMCIF2Dict(cifPath)
        # print(DICT)
        return STR #,DICT

def chunkDict(inputDict):
    # print(inputDict)
    items = list(inputDict.items())
    chunksize = 24
    chunks = [items[i:1+chunksize] for i in range(0,len(items),chunksize)]
    return chunks

def list_tuples_to_dict(lTuple):
    # print(lTuple)
    output = dict(item for item in lTuple)
    # print(output)
    return output

def f_init(q):
    process_queue.q = q


#test one to three
def process_queue(id):

    # print(id) 
    # print(c)
    # print(i)
    # print("i:", i)
    # iDict = list_tuples_to_dict(chunks)
    # print(type(iDict))
    # print(iDict.keys())
    # iDict = list_to_dict(i)
    # print(id)
    # pbar.update(1)
    series = feDF.loc[feDF['Microenvironment_ID'] == id].to_dict()
    # print(list(series['Metal'].values())[0])
    # annotations = {attr:series[attr] for attr in hagaiNames}

    pdb,LIG,chain,resID,function = re.split('[._]',id)
    # # print(pdb,LIG,chain,resID,function,i.Metal)
    cifPath = os.path.join(pdir,pdb+".cif")
    STR = get_cif_STR(cifPath)

    try:

        RES = STR[0][chain][("H_"+LIG,int(resID)," ")]

    except Exception as e:
            print("  ",series['Microenvironment_ID'],("H_"+LIG,int(resID)," "))
            print(e)
            # continue
        
    resName = RES.get_resname()
    feATOMS = [atom for atom in RES.get_atoms() if "FE" in atom.get_name()] 

    atom_list = PDB.Selection.unfold_entities(STR[0],"A")
    ns = PDB.NeighborSearch(atom_list)
    
    #return dataframe
    output = []
    
    for atom in feATOMS:
        atomName = atom.get_name()
        feAtom = {"id_hagai":id,"ligName":resName,"ligChain":chain,"resID":resID,"ligAtom":atomName,"Metal":list(series['Metal'].values())[0],"ligX":atom.get_coord()[0],"ligY":atom.get_coord()[1],"ligZ":atom.get_coord()[2]}#,"micro_annotations":annotations}  
        neighbors = ns.search(atom.get_coord(),3,level='A')
        # output.append((id,len(neighbors)))
        c = 0 
        for neigh in neighbors: 
            if neigh-atom != 0:
                feAtom["N_"+str(c)+"_resName"]=neigh.get_parent().get_resname()
                feAtom["N_"+str(c)+"_atomName"] = neigh.get_name()
                feAtom["N_"+str(c)+"_resNum"] = neigh.get_parent().get_id()[1]
                feAtom["N_"+str(c)+"_x"] = neigh.get_coord()[0]
                feAtom["N_"+str(c)+"_y"] = neigh.get_coord()[1]
                feAtom["N_"+str(c)+"_z"] = neigh.get_coord()[2]
                feAtom["N_"+str(c)+"_dist"] = neigh-atom
                feAtom["N_"+str(c)+"_dx"] = neigh.get_coord()[0]-atom.get_coord()[0]
                feAtom["N_"+str(c)+"_dy"] = neigh.get_coord()[1]-atom.get_coord()[1]
                feAtom["N_"+str(c)+"_dz"] = neigh.get_coord()[2]-atom.get_coord()[2]
                count += 1
        output.append(feAtom)
        # neigh_dict = {id:{"feAtom":feAtom,"Microenvironment_annotations":annotations,"Neighbors":{"Neighbor_"+str(idx): {"N_"+str(idx)+"_resName":neigh.get_parent().get_resname(),"N_"+str(idx)+"_atomName":neigh.get_name(),"N_"+str(idx)+"_resNum":neigh.get_parent().get_id()[1],"N_"+str(idx)+"_atomCoord":neigh.get_coord(),"N_"+str(idx)+"_distance":neigh-atom,"N_"+str(idx)+"_dx":neigh.get_coord()[0]-atom.get_coord()[0],"N_"+str(idx)+"_dy":neigh.get_coord()[1]-atom.get_coord()[1],"N_"+str(idx)+"_dz":neigh.get_coord()[2]-atom.get_coord()[2]} for idx, neigh in enumerate(neighbors) if neigh-atom != 0 }}}
        # print(neigh_dict)
    process_queue.q.put(output)
    return output



def listener():
    '''listens for messages on a the q, writes to a file'''
    #https://stackoverflow.com/questions/13446445/python-multiprocessing-safely-writing-to-a-file
    
    pdir = "/home/kenneth/proj/proXtal/proteins/hagai/pdb"

    # if os.path.exists(os.path.join(pdir,"micro_multiproc.csv")):
    #     print("Output file exists already")
    #     exit()
    
    # print("Output jobs to file\n")

    with open(os.path.join(pdir,"micro_multiproc_2.csv"),"w") as outWriter:
        while 1:
            m = process_queue.q.get()
            # print(type(m))
            if m == "kill":
                outWriter.write('killed')
                break
            for i in m:
                outWriter.write(str(i) + '\n')
                outWriter.flush()

# print(logFile)
# print(feDF)
def main():
    
    manager = mp.Manager()
    q = manager.Queue()
    #init objects
    cores = mp.cpu_count()
    # print(cores)
    # https://stackoverflow.com/questions/3827065/can-i-use-a-multiprocessing-queue-in-a-function-called-by-pool-imap
    pool = mp.Pool(cores,f_init,[q])

    #put listeners to work
    watcher = pool.apply_async(listener)

    #fire off workers
    microIDs = feDF['Microenvironment_ID'].tolist()
    jobs = []
    output = pd.DataFrame()

    # bar = SlowBar("Processing",max=len(microIDs),fill='*') #,suffix='%(percent).1f%% - %(eta)ds')
    # pbar = tqdm(total=len(microIDs))
    # count = 0

    # for id in microIDs:
        # print(len(id))
    # count +=1
    # print(id)
    with tqdm(total=len(microIDs)) as pbar: 
        for i, job in tqdm(enumerate(pool.imap_unordered(process_queue,microIDs))):
            jobs.append(job)
            # print(job)
            pbar.update()
            # if i < 10:
                # jobs.append(job)
                # break
                # pbar.update()
            # print(type(job))

       
    q.put("kill")
    pool.close()
    pool.join()
    pbar.close()

main()