#!/usr/bin/env python
from Bio.PDB.PDBParser import PDBParser
import os
import pandas as pd
import numpy as np
from Bio.PDB.PDBIO import Select
from Bio import PDB
# from diffpy import structure 
# import gemmi
import inspect
import re
import multiprocessing as mp
from collections import defaultdict
from tqdm import tqdm
from collections import Counter
# Import SearchIO and suppress experimental warning
# from Bio import PDBConstructionWarning
# with warnings.catch_warnings():
#     warnings.simplefilter('ignore', PDBConstructionWarning)
#     from Bio import SearchIO

import glob
import sys 
sys.path.append(os.path.abspath("/home/kenneth/proj/proMin/code/tools"))
import tools
#global variables

global dfHag
dfHag = tools.load_findgeo_hagai_df()
print((dfHag.columns))
# exit()
def f_init(q):
    process_queue.q = q

def get_environment(atoms):

    metals = set(["CO","CU","FE","MN","MO","NI","V","W"]) #NO MG

    # elements = [get_element(atom.get_name().upper()) for atom in STR.get_atoms()]
    elements = [atom.element for atom in atoms]
    
    countElements = Counter(elements)
    # return countElements
    strCountElements = ""

    #     #sort based on element, returns a list of sorted keys 
    sortCountElement_keys = sorted(countElements.keys(),reverse=True)
    
    for m in set(sortCountElement_keys).intersection(metals):
        strCountElements += m+str(countElements[m])

    #   #string rest of elements
    for m in sorted(set(sortCountElement_keys).difference(metals)):
        strCountElements += m+str(countElements[m])

    return strCountElements

#test one to three
def process_queue(id):
    pdir = "/home/kenneth/proj/proMin/proteins/hagai/2018/pdbs"
    df = dfHag.iloc[id]
    # print(df)
    pdbFile = os.path.join(pdir,df.pdb+".pdb")
    #Create a PDBParser object
    parser = PDBParser() #PERMISSIVE = 0 will list all errors with PDB file
    STR = parser.get_structure(df.pdb, os.path.join(pdbFile))
    envComp = ""

    atoms = STR[0][df.chain].get_atoms()
    atomID = ""

    for atom in atoms:
        # print(atom.get_fullname())
        # print(atom.get_fullname().replace(' ',''))
        if atom.get_name() == df.atomName:
            resID = atom.get_parent().get_full_id()
            # print(resID,'blah')
            if (resID[3][0].find(df.resName) > -1) and (int(resID[3][1]) == df.resNum):
                atomID = resID[3]
                if atomID == '':
                    print(df.coord,resID,'blank')
                    continue

                chain = STR[0][df.chain]
                atom = chain[atomID][atom.get_name()]
                atom_list = PDB.Selection.unfold_entities(STR[0],'A')
                ns = PDB.NeighborSearch(atom_list)
                atoms = ns.search(atom.get_coord(),2.8,level='A')
                # print(atoms)
                envComp = df.coord + "," + df.resName + "," + ",".join(map(str,list(resID))) + "," + get_environment(atoms)

    # print(envComp)
    process_queue.q.put(envComp)
    return envComp



def listener():
    '''listens for messages on a the q, writes to a file'''
    #https://stackoverflow.com/questions/13446445/python-multiprocessing-safely-writing-to-a-file
    rdir = "/home/kenneth/proj/proMin/proteins/hagai/2018/pdbs/findgeo/envComp"

    # if os.path.exists(os.path.join(pdir,"micro_multiproc.csv")):
    #     print("Output file exists already")
    #     exit()
    
    # print("Output jobs to file\n")

    with open(os.path.join(rdir,"envComp.csv"),"w") as outWriter:
        while 1:
            m = process_queue.q.get()
            # print(type(m))
            if m == "kill":
                outWriter.write('killed')
                break
            else:
                outWriter.write(str(m) + '\n')
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
    jobs = []

    with tqdm(total=len(dfHag.index)) as pbar: 
        for i, job in tqdm(enumerate(pool.imap_unordered(process_queue,dfHag.index))):
            jobs.append(job)
            pbar.update()
    q.put("kill")
    pool.close()
    pool.join()
    pbar.close()

main()