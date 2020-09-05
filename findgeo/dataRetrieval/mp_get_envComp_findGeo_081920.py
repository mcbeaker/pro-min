#!/usr/bin/env python
import os
# from CifFile import ReadCif
import pandas as pd
import numpy as np
from Bio.PDB import MMCIF2Dict
from Bio.PDB.PDBIO import Select
from Bio import PDB
# from diffpy import structure 
# import gemmi
from progress.bar import Bar
import subprocess
import inspect
import re
import multiprocessing as mp
# from elements import ELEMENTS
from collections import defaultdict
import json
from tqdm import tqdm
from collections import Counter
import glob
import math
# from elements import ELEMENTS

#global variables
global pdbPaths
pdbPaths = {}

pdbPathsLoc='/home/kenneth/proj/proMin/proteins/rcsb/pdbs_0_1.5/findgeo/combineFindGeoResults'
# pdbPathsLoc='/home/kenneth/proj/proMin/minerals/database/amcsd_pdb/pdb_reSeq_res_atom/combineFindGeoResults'
count = 0 

# metals = ["FE_"]
for pdbFile in glob.iglob(pdbPathsLoc+'/*.pdb'):
    # if any(metal in pdbFile for metal in metals):
        # print(pdbFile)
    pdbPaths[count] = pdbFile
    count += 1
    # else:
        # continue
        

def get_pdb_STR(pdbPath):
    STR = PDB.MMCIFParser(QUIET=True).get_structure("pdb",pdbPath)
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

def get_metalRow(listAtoms,metalName):
    index = 0
    for i, atomRow in enumerate(listAtoms):
        if metalName in atomRow.get_name().upper():
            index = i
            return index
    return index

def get_element(atomName):
    elements = [e.symbol.upper() for e in ELEMENTS] 

    elem1 = ["C","N","O","F","P","S","W","V"]
    elem2 = ["FE","CO","MN","CU","NI","MO"]

    elem = ""

    if atomName[0:2] in elem2:
        elem = atomName[0:2]
    elif atomName[0:1] in elem1:
        elem = atomName[0:1]
    else:
        # elem = atomName[0:2]
        # print("Did not find element"+ str(atomName))
        elem=elem
    # print(elem)
    return elem

def get_atomDistances(structure,metalName):
    # print(metalName) 
    metals = ["FE", "CO", "MN", "CU", "NI", "MO","W", "V"]

    metalRow = get_metalRow(list(structure.get_atoms()),metalName)
    #print(atoms)
    distances = {} # create a dictionary to hold the dij values of this structure
    elemLigand = "" # this will be used to store the terminal elements name
    atoms = list(structure.get_atoms())
    numAtoms = len(atoms)

    # get dij distances for the structure and store in dictionary. leaving out metal atom: metal atom distance
    distance12 = {} # empty dictionary to hold information about dij's (the dij dictionaries will be made FROM this later)
    
    for idx in range(0,numAtoms):
        if idx != metalRow:
            distance = atoms[metalRow] - atoms[idx]
            atomNames = metalName+"_"+atoms[idx].get_name().upper()
            # distance12[atomNames] = (get_element(metalName)+":"+get_element(atoms[idx].get_name().upper()),distance)
            distance12[atomNames] = (atoms[metalRow].element+":"+atoms[idx].element,distance)
    return distance12
#
def get_valParms(metElem,ligElem):
    # need metal elements
    # metals = ["CO","CU","FE","MN","MO","NI","V","W"] #NO MG
    # ligands = ["C","N","O","F","P","S"]
    # parmFile = '/Users/ken/Box/proj/proXtal/nolan_valance/valenceParms/bvparm2006.csv'
	# print(metElem)
	# print(ligElem)
	parmFile='/home/kenneth/proj/proMin/valenceParms/bvparm2006.csv'
	# parmFile = '/home/kenneth/box/proj/proXtal/nolan_valance/valenceParms/bvparm2006.csv'
	dfParms = pd.read_csv(parmFile,sep=",",index_col=False,usecols=["valence_param_atom_1",
                                            "valence_param_atom_1_valence",
                                            "valence_param_atom_2",
                                            "valence_param_atom_2_valence",
                                            "valence_param_Ro",
                                            "valence_param_ref_id"])  
	dfParms['valence_param_atom_1'] = dfParms.valence_param_atom_1.str.upper()
	dfParms['valence_param_atom_2'] = dfParms.valence_param_atom_2.str.upper()
    # print(dfParms.valence_param_atom_1)
	dfParms = dfParms[(dfParms.valence_param_atom_1.isin(metElem)) & (dfParms.valence_param_atom_2.isin(ligElem))]
	dfParms = dfParms[dfParms.valence_param_atom_1_valence <= 6]
	# get distinct values of the dataframe based on column

	return dfParms

def get_valencyStates(dfParm):

    metalVal = {}

    for metal in list(dfParm.valence_param_atom_1):
        group = dfParm.groupby(['valence_param_atom_1','valence_param_atom_1_valence']).apply(lambda x:x['valence_param_atom_1_valence'].unique()) #get only unique valency for each metal
        metalVal[metal] = list(group[metal].valence_param_atom_1_valence)
    return metalVal

def calc_valence(r0,measured,ox_num):
	r0 = float(r0)
	# print(r0)
	measured = float(measured)
	ox_num = int(ox_num)
	# print(r0,measured,ox_num)
	b = 0.37    
	valence = math.exp((r0 - measured) / b)
	return valence

# defining the VALENCY calculation function --------------------------------------------------------------------------
# Author: Nolan Fehon/Kenneth
def get_valency(dists):  # inpu
	
	numDist = len(dists.keys())
	mElem = [value[0].split(":")[0] for value in dists.values()]
	lElem = [value[0].split(":")[1] for value in dists.values()] 
	
	# key, metElem, ligElem, dist
	minLigDist = [(key,value[0].split(":")[0],value[0].split(":")[1],value[1]) for key,value in zip(dists.keys(),dists.values())]
	
    # #keep track of metal and ligand name to get distance values
	dfParms = get_valParms(mElem,lElem)
	# print(dfParms)
	
	# get distinct values of the dataframe based on column
	if ('Fe' not in mElem) and ('N' not in lElem):
		dfParms = dfParms.drop_duplicates(subset = ["valence_param_atom_1","valence_param_atom_1_valence","valence_param_atom_2"], keep='first')
	
	#deal with Fe LS and Fe HS - Zheng Acta Cryst 2017 D73 316-325
	if ('Fe' in mElem) and ('N' in lElem): #check if Fe and N are in list
		#find index
		FeN = dfParms[(dfParms.valence_param_atom_1 == 'Fe') & (dfParms.valence_param_atom_2 == 'N')] #get rows with Fe and N
		dfParms = dfParms.drop_duplicates(subset = ["valence_param_atom_1","valence_param_atom_1_valence","valence_param_atom_2"], keep='first')
		dfParms.drop((dfParms.valence_param_atom_1 == 'Fe') & (dfParms.valence_param_atom_2_valence == 'N'), inplace = True) #drop all rows with Fe and N
		dfParms = dfParms.append(FeN,ignore_index = True)

	# print(dfParms)
	metVal = {}

	# loop through key, metElem, ligElem, dist
	#keep track of OxNum intersection
	oxInt = {}
	#count the keys
	count = 0
	for key,metElem,ligElem, dist in minLigDist:

		#get oxidative states
		dfOxParm = dfParms[(dfParms.valence_param_atom_1 == metElem) & (dfParms.valence_param_atom_2 == ligElem)]

		metVal[key] = {}
		metVal[key]['mElem'] = metElem
		metVal[key]['lElem'] = ligElem
		metVal[key]['Ox'] = []
		metVal[key]['Valence'] = []
		metVal[key]['r0'] = []
		# key,metElem,ligElem,dist,ox,r0,b
		# loop through standard oxidation values
		# print('test')
		for	i,oxParm in dfOxParm.iterrows():
			# print(oxParm)
			r0 = oxParm['valence_param_Ro']
			Ox = oxParm['valence_param_atom_1_valence']
			metVal[key]['r0'].append(r0)
			metVal[key]['Ox'].append(Ox)
			metVal[key]['Valence'].append(calc_valence(r0,dist,Ox))

		#first time through
		if count == 0:
			oxInt = set(metVal[key]['Ox'])
		#multiple times through
		else:
			oxInt = oxInt.intersection(set(metVal[key]['Ox']))
			# print('test2')
		count +=1
	#how close
	
	#loop through oxNum first for all atom_pairs
	feIndx = []
	for ox in oxInt:
		val = 0
		oxInd = 0
		
			if key == 'oxNum' or key == 'Valency':
				continue
			elif (key.find('FE') > -1) and (key.find('N') > -1):
				feIndx.append(i)
				continue
			else:
				# print(key)
				# print(type(metVal[key]))
				#need to deal with Fe2 and Fe3 low and high spin# not sure how to do that - only looking at one now
				oxInd = metVal[key]['Ox'].index(ox) 
				# print(oxInd)	
				val += float(metVal[key]['Valence'][oxInd])
				# print('test3')
	#deal with FE

		
		#update values 
		if abs(val-ox) < dif:
			valency = val
			oxNum = ox
			metVal['oxNum'] = oxNum
			# print('test4')
			# print(valency,ox)
			dif = abs(val-ox)
	
	for feIdx in feIndx:
	
	
	
	
	metVal['Valency'] = valency
	# print(metVal)

	return metVal

def calc_vecsum(structure,metalName,valenceDictionary):
	# print(valenceDictionary.keys())
	metals = ["FE", "CO", "MN", "CU", "NI", "MO","W", "V"]
	atoms = list(structure.get_atoms())
	metalRow = get_metalRow(list(structure.get_atoms()),metalName)
	metalAtom = atoms[metalRow]
	numAtoms = len(atoms)

	vecsum = 0
	fij = PDB.Vector(x=0,y=0,z=0)
	for idx in range(0,numAtoms):
		if idx != metalRow:
			# print('blah')
			atomNames = metalName+"_"+atoms[idx].get_name().upper()
			ligandAtom = atoms[idx]
			distance = abs(ligandAtom - metalAtom)
			vec = (ligandAtom.get_vector() - metalAtom.get_vector())
			rij = vec.__truediv__(distance)
			ligOcc = ligandAtom.get_occupancy()
			# print('ligOCC: ',ligOcc)
			# print('valence: ',valenceDictionary[atomNames]['Valence'])
			oxInd = valenceDictionary[atomNames]['Ox'].index(valenceDictionary['oxNum'])
			bondValence = float(valenceDictionary[atomNames]['Valence'][oxInd])
			# print('blha: ' + str(bondValence))
			sij = float(ligOcc) * bondValence
			# print('sij',sij)
			# raise TypeError('somethingHappend ' + str(ij))
			fij = fij.__add__(np.multiply(rij.get_array(),sij))
			# print('fij: ',fij)
	vecsum =  math.sqrt(fij.__mul__(fij)) / float(valenceDictionary['Valency'])
	# print('vecsum: ',vecsum)
	return vecsum


def get_rmsd(pdbFileName):

	pdb,metalID,geoShort,ext = re.split('[.]',pdbFileName)
	atomName,resNum,atomNum,chain = re.split('_',metalID)
	atomName = atomName.upper()

	pdbPathsLoc='/home/kenneth/proj/proMin/proteins/rcsb/pdbs_0_1.5/findgeo'
	# pdbPathsLoc='/home/kenneth/proj/proMin/proteins/hagai/pdbs'
	# pdbPathsLoc='/home/kenneth/proj/proMin/minerals/database/amcsd_pdb/pdb_reSeq_res_atom'

	rmsdFile = os.path.join(pdbPathsLoc,pdb,metalID,geoShort+".out")
	# print(rmsdFile)
	RMSD = 0
	# print(rmsdFile+"\n")
	if os.path.exists(rmsdFile) == True:
		with open(rmsdFile, "r") as a_file:
			for line in a_file:
				if "Total RMSD" in line: # only acting on line with RMSD present
					stripped_line = line.strip() # removing extra space characters from the line (extra)
					# print(stripped_line)
					found = re.search("(?!Total RMSD\w+)\d+.\d+", stripped_line) # pulling RMSD. won't take "." as a digit
					RMSD = found.group()
	# Findgeo RMSD (angstroms) - quantitative measure of the similarity between the configuration of the coordinating ligands and the various geometries possible for the metal coordination number, which can then be ranked to identify the best geometry assignment
	else:
		RMSD = "NA"
	# # convert radians to degrees
	# if RMSD != "NA":
	# 	# return RMSD
	# 	RMSD = str(np.rad2deg(float(RMSD)))
	
	return RMSD

def get_environment(STR,pdbFileName):

    pdb,metalID,geoShort,ext = re.split('[.]',pdbFileName)
    atomName,resNum,atomNum,chain = re.split('_',metalID)
    atomName = atomName.upper()
    metals = set(["CO","CU","FE","MN","MO","NI","V","W"]) #NO MG

    # elements = [get_element(atom.get_name().upper()) for atom in STR.get_atoms()]
    elements = [atom.element for atom in STR.get_atoms()]
    
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

	pdbPath = pdbPaths[id]
	# print(pdbPath)
	pdbFileName = os.path.basename(pdbPath)
	# pdb,amcsdID,metalName,resNum,atomNum,chain,geoShort,ext = re.split('[._]',pdbFileName)
	pdb,metalName,resNum,atomNum,chain,geoShort,ext = re.split('[._]',pdbFileName)
	metalName = metalName.upper()


	import warnings
	warnings.filterwarnings("error")
	# print(pdb,metalName)
	output=""
	try:

		STR = PDB.PDBParser(QUIET=False).get_structure(pdb,pdbPath)

		environment = get_environment(STR,pdbFileName)
		# print(environment)
		gRMSD = get_rmsd(pdbFileName)
		# print(gRMSD)
		atomDist = get_atomDistances(STR,metalName)
		# print(atomDist)
	# print(atomDist)
	# print(metalElements)
	# print(ligElements)
	# print(atomDist)
		vecsum =0 
		dict_metalValence = calc_valency(atomDist)
		if dict_metalValence['Valency'] != 0:
			vecsum = calc_vecsum(STR,metalName,dict_metalValence)

		# print(vecsum)
# print(valency)
# output = ""
# for i,dic in enumerate(atom_dist.items()):
#     if i == 0:
#         # print(dic)
#         output = pdbFileName + " " + dic[1][0]
#     else:
#         output += " " + dic[1][0]

	# print(output)
		# output = ""
		import json
		output = os.path.splitext(pdbFileName)[0] + " " + environment + " " + gRMSD + " " + str(dict_metalValence['Valency'])+ " " + str(vecsum)
		print('pdbID env gRMSD Valency vecsum')
		print(output+'\n')
		print(json.dumps(dict_metalValence,indent=4))

		# process_queue.q.put(output)
	except PDB.PDBExceptions.PDBConstructionWarning as e:
		print(e)
		print(pdbPath)
	# except TypeError as t:
	# 	print(t,'blak')
	# 	print(pdbPath)

	# print(output)
	return output



def listener():
	'''listens for messages on a the q, writes to a file'''
	#https://stackoverflow.com/questions/13446445/python-multiprocessing-safely-writing-to-a-file

	pdir = '/home/kenneth/proj/proMin/proteins/rcsb/pdbs_0_1.5/findgeo/combineFindGeoResults/analysis'
	# pdir = '/home/kenneth/proj/proMin/proteins/hagai/pdbs/combineFindGeoResults/analysis'
	# pdir = '/home/kenneth/proj/proMin/minerals/database/amcsd_pdb/pdb_reSeq_res_atom/combineFindGeoResults/analysis'
	outFile='pro_fg_dictionaryValues_vecsum.csv'
	outFile = os.path.join(pdir,outFile)
	if os.path.exists(pdir) == False:
		# print('exists ' + outFile)
		os.mkdir(pdir)
    #     print("Output file exists already")
    #     exit()
    
    # print("Output jobs to file\n")
	# print(outFile)
	with open(outFile,"w") as outWriter:
		while 1:
			m = process_queue.q.get()
			# print(m)
            # print(type(m))
			if m == "kill":
				outWriter.write('killed')
				break
			else:
				outWriter.write(m + '\n')
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

    with tqdm(total=len(pdbPaths.keys())) as pbar: 
        for i, job in tqdm(enumerate(pool.imap_unordered(process_queue,list(pdbPaths.keys())[0:10]))):
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
# print(get_valencyStates())