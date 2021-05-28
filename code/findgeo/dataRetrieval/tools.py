#!/usr/bin/env python
from Bio import PDB
import pandas as pd
import os
from collections import Counter
import math 
import re 
import glob
import numpy as np

def load_mineral_oxidation():
    transMetals = ['Fe','Cu','Mn','Ni','Mo','Co','V','W']
    transMetals.sort()
    dic = {}
    # f = '/Users/ken/Box/proj/proXtal/justin_redoxPotential/data/mineralname_mineralcharge.csv'
    f = '/home/kenneth/proj/proMin/minerals/mineralname_mineralcharge.csv'
    df = pd.read_csv(f,header=0,index_col=False)
    for i,series in df.iterrows():
        # print(series.Mineral)
        dic[series['Mineral']] = {}
        dic[series['Mineral']]['Formula'] = series['Formula']
        # if series['Formula'].find('Fe ') > -1:
            # print(series['Mineral'],series['Formula'])
        # print(dic)
        for metal in transMetals:
            # print(metal)
            if metal == 'Cu' and series[metal] == 'none':
                dic[series['Mineral']][metal.upper()] = '1+'
            else:
                dic[series['Mineral']][metal.upper()] = series[metal]
    return dic

def get_mineral_pdb_names():
    d = '/Users/ken/Box/proj/proXtal/data/proteins/code-pro-min/data/minerals'
    import glob as glob
    fList = [os.path.basename(pdbFile)for pdbFile in glob.iglob(d+'/*.pdb')]
    return fList


# mineral,Mn_30_30_A_irr,geoShort,ext = re.split('[.]',pdbFileName)
	# metalName,resNum,atomNum,chain = re.split('[._]',Mn_30_30_A_irr)
	
def get_pdbFilePaths(pdbPathsLoc,pdbFileNames):
	count = 0 
 
	# metals = ["FE_"]
	with open(os.path.join(pdbPathsLoc,pdbFileNames),"w") as output:
		
		for pdbFile in glob.iglob(pdbPathsLoc+'/*.pdb'):
			# if any(metal in pdbFile for metal in metals):
				# print(pdbFile)
			output.write(str(count)+","+pdbFile+"\n")
			count += 1
	# outFile = pd.read_csv(os.path.join(pdbPathsLoc,pdbFileNames),index_col=0)
	# return outFile

def get_pdb_STR(pdbPath):
    STR = PDB.MMCIFParser(QUIET=True).get_structure("pdb",pdbPath)
    # DICT = PDB.MMCIF2Dict.MMCIF2Dict(cifPath)
    # print(DICT)
    return STR #,DICT

def list_tuples_to_dict(lTuple):
    # print(lTuple)
    output = dict(item for item in lTuple)
    # print(output)
    return output

def get_metalRow(listAtoms,metalName):
    index = 0
    for i, atomRow in enumerate(listAtoms):
        if metalName in atomRow.get_name().upper():
            index = i
            return index
    return index

def get_angles(STR,metalName):
	atoms = list(STR.get_atoms())
	numAtoms = len(atoms)
	metalRow = get_metalRow(atoms,metalName)

	angles = []
	for i in range(0,numAtoms):
		for j in range(i+1,numAtoms):
			if (i != metalRow) & (j != metalRow):
				atm1 = atoms[i]
				# print("test\n")
				metal = atoms[metalRow]
				# print("test2\n")
				atm2 = atoms[j]
				# print("test3\n")
				angle = round(PDB.vectors.calc_angle(atm1.get_vector(),metal.get_vector(),atm2.get_vector()) * (180/np.pi),2)
				# print(str(angle)+"\n")
				atomNames = atm1.get_name().upper()+"_"+metalName+"_"+atm2.get_name().upper()
				# print("test5\n")
				residueNames = atm1.get_parent().get_resname() + "_" + metal.get_parent().get_resname() + "_" + atm2.get_parent().get_resname()
				# print("test6\n")
				residueNumbers = str(atm1.get_parent().get_full_id()[3][1]) + "_" +str(metal.get_parent().get_full_id()[3][1]) + "_" + str(atm2.get_parent().get_full_id()[3][1]) 
				# print("test7\n")
				residueChain = str(atm1.get_parent().get_full_id()[2]) + "_" + str(metal.get_parent().get_full_id()[2]) + "_" + str(atm2.get_parent().get_full_id()[2]) 
				# print("test8\n")
				occ = str(atm1.get_occupancy()) + "_" + str(metal.get_occupancy()) + "_" + str(atm2.get_occupancy()) 	
				elements = atm1.element + "_" + metal.element + "_" + atm2.element  
				angleInfo = [metal.element,elements,angle,atomNames,residueNames,residueNumbers,residueChain,occ]
				angleStr = "" #string to easily output to file
				for l in range(0,len (angleInfo)):
					if l == 0:
						angleStr = str(angleInfo[l])
					if l > 0:
						angleStr += "," + str(angleInfo[l])
				angles.append(angleStr)
	return angles

def get_atomDistances(structure,metalName):
    # print(metalName) 
    metals = ["FE", "CO", "MN", "CU", "NI", "MO","W", "V"]
    atoms = list(structure.get_atoms())     #print(atoms)

    metalRow = get_metalRow(atoms,metalName)

    distances = {} # create a dictionary to hold the dij values of this structure
    elemLigand = "" # this will be used to store the terminal elements name
    
    numAtoms = len(atoms)

    # get dij distances for the structure and store in dictionary. leaving out metal atom: metal atom distance
    distance12 = {} # empty dictionary to hold information about dij's (the dij dictionaries will be made FROM this later)
    
    for idx in range(0,numAtoms):
        if idx != metalRow:
            distance = atoms[metalRow] - atoms[idx]
            atomNames = metalName+"_"+atoms[idx].get_name().upper()
            residueNames = atoms[metalRow].get_parent().get_resname() + "_" + atoms[idx].get_parent().get_resname()
            residueNumbers = str(atoms[metalRow].get_parent().get_full_id()[3][1]) + "_" + str(atoms[idx].get_parent().get_full_id()[3][1]) 
            residueChain = str(atoms[metalRow].get_parent().get_full_id()[2]) + "_" + str(atoms[idx].get_parent().get_full_id()[2]) 
            metOcc = atoms[metalRow].get_occupancy()
            ligOcc = atoms[idx].get_occupancy() #use for valency calc
            distance12[atomNames] = (atoms[metalRow].element+":"+atoms[idx].element,distance,metOcc,ligOcc,residueNames,residueNumbers,residueChain)
    return distance12
#
def get_valParms(metElem,ligElem):
    # need metal elements
    # metals = ["CO","CU","FE","MN","MO","NI","V","W"] #NO MG
    # ligands = ["C","N","O","F","P","S"]
	parmFile='/home/kenneth/proj/proMin/valenceParms/bvparm2006.csv' #includes 2017 Zheng Fe - spin II/III
	
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

# defining the VALENCY calculation function --------------------------------------------------------------------------
# Author: Nolan Fehon/Kenneth
def calc_bond_valence(r0,measured,ox_num):
	r0 = float(r0)
	measured = float(measured)
	ox_num = int(ox_num)
	# print(r0,measured,ox_num)
	b = 0.37    
	valence = math.exp((r0 - measured) / b)
	return valence

def calc_bvs(dists,oxStates):  # inpu
	#calculates bond valence sum
	numDist = len(dists.keys())
	mElem = [value[0].split(":")[0] for value in dists.values()]
	lElem = [value[0].split(":")[1] for value in dists.values()] 
	
	# key, metElem, ligElem, dist,ligOcc
	minLigDist = [(key,value[0].split(":")[0],value[0].split(":")[1],value[1],value[3]) for key,value in zip(dists.keys(),dists.values())]
	
    # #keep track of metal and ligand name to get distance values
	dfParms = get_valParms(mElem,lElem)
	# print(dfParms)
	
	metVal = {}

	# loop through key, metElem, ligElem, dist
	#keep track of OxNum intersection
	oxInt = {}
	#count the keys
	count = 0
	for key,metElem,ligElem,dist,ligOcc in minLigDist:

		#get oxidative states
		dfOxParm = dfParms[(dfParms.valence_param_atom_1 == metElem) & (dfParms.valence_param_atom_2 == ligElem)]
		ox = oxStates[metElem]
		print(metElem,oxStates)
		# quit()
		metVal[key] = {}
		metVal[key]['mElem'] = metElem
		metVal[key]['lElem'] = ligElem
		metVal[key]['Ox'] = [] 
		metVal[key]['Valence'] = []
		# key,metElem,ligElem,dist,ox,r0,b
		# loop through standard oxidation values
		
		for	i,oxParm in dfOxParm.iterrows():
			# print(oxParm)
			r0 = oxParm['valence_param_Ro']
			Ox = oxParm['valence_param_atom_1_valence']
			metVal[key]['Ox'].append(Ox)
			# print(type(ligOcc))
			metVal[key]['Valence'].append(calc_bond_valence(r0,dist,Ox)*ligOcc)
		#first time through
		if count == 0:
			oxInt = set(metVal[key]['Ox'])
		#multiple times through
		else:
			oxInt = oxInt.intersection(set(metVal[key]['Ox']))
		count +=1
	#how close
	dif = 99
	valency = 0
	oxNum = 0

	#loop through oxNum first for all atom_pairs
	for ox in oxInt:
		val = 0
		for key in metVal.keys():
			oxInd = metVal[key]['Ox'].index(ox)	
			val += float(metVal[key]['Valence'][oxInd])
		
		#update values 
		if abs(val-ox) < dif:
			valency = val
			oxNum = ox
			# print(valency,ox)
			dif = abs(val-ox)
	
	metVal['valency'] = valency

	return metVal

def calc_vecsum(valency,STR,metalName):
 
	metalRow = get_metalRow(list(STR.get_atoms()),metalName)
	metals = ["FE", "CO", "MN", "CU", "NI", "MO","W", "V"]

    #print(atoms)
	distances = {} # create a dictionary to hold the dij values of this structure
	elemLigand = "" # this will be used to store the terminal elements name
	atoms = list(STR.get_atoms())
	numAtoms = len(atoms)

	# get dij distances for the structure and store in dictionary. leaving out metal atom: metal atom distance
	distance12 = {} # empty dictionary to hold information about dij's (the dij dictionaries will be made FROM this later)
	
	for idx in range(0,numAtoms):
		if idx != metalRow:
			distance = atoms[metalRow] - atoms[idx]
			atomNames = metalName+"_"+atoms[idx].get_name().upper()
			metOcc = atoms[metalRow].get_occupancy()
			ligOcc = atoms[idx].get_occupancy()
			# distance12[atomNames] = (get_element(metalName)+":"+get_element(atoms[idx].get_name().upper()),distance)
			distance12[atomNames] = (atoms[metalRow].element+":"+atoms[idx].element,distance,metOcc,ligOcc)
	return distance12

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
	else:
		RMSD = "NA"
	#RMSD is in Angstroms
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
