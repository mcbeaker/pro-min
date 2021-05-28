#!/usr/bin/env python
from Bio import PDB
import pandas as pd
import os
from collections import Counter
import math 
import re 
import glob
import numpy as np
import random

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
			# if metal == 'Cu' and series[metal] == 'none':
			#     dic[series['Mineral']][metal.upper()] = '1'
			# else:
			dic[series['Mineral']][metal.upper()] = series[metal]
	print(dic['Greigite'])
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
			# test_list = ['irr','4wes.FE_501_15503_B.oct.pdb','4tku.FE_605_16385_D.oct.pdb',
						# '4tku.FE_605_16316_B.oct.pdb','4wes.FE_502_15504_B.oct.pdb','6o7m.FE_602_31697_B.coc.pdb']
			# res = any(ele in pdbFile for ele in test_list) 
			# if  res == True : # any(metal in pdbFile for metal in metals):
				# continue
			output.write(str(count)+","+pdbFile+"\n")
			count += 1
	outFile = pd.read_csv(os.path.join(pdbPathsLoc,pdbFileNames),index_col=0,header=None)
	return outFile

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
			metVec = atoms[metalRow].get_vector()
			ligVec = atoms[idx].get_vector()
			atomNames = metalName+"_"+atoms[idx].get_name().upper()
			residueNames = atoms[metalRow].get_parent().get_resname() + "_" + atoms[idx].get_parent().get_resname()
			residueNumbers = str(atoms[metalRow].get_parent().get_full_id()[3][1]) + "_" + str(atoms[idx].get_parent().get_full_id()[3][1]) 
			residueChain = str(atoms[metalRow].get_parent().get_full_id()[2]) + "_" + str(atoms[idx].get_parent().get_full_id()[2])
			metOcc = atoms[metalRow].get_occupancy()
			ligOcc = atoms[idx].get_occupancy() #use for valency calc
			ligAtomNum = atoms[idx].get_serial_number() 
			if atomNames in distance12:
				n = random.randint(0,50)
				distance12[atomNames + "_" + str(n)] = (atoms[metalRow].element+":"+atoms[idx].element,distance,metOcc,ligOcc,residueNames,residueNumbers,residueChain,metVec,ligVec,ligAtomNum)
			else:
				distance12[atomNames] = (atoms[metalRow].element+":"+atoms[idx].element,distance,metOcc,ligOcc,residueNames,residueNumbers,residueChain,metVec,ligVec,ligAtomNum)
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

def calc_protein_bvs(ID,dists,env):  # inpu
	#calculates bond valence sum
	mElem = [value[0].split(":")[0] for value in dists.values()]
	lElem = [value[0].split(":")[1] for value in dists.values()] 
	
	# key, metElem, ligElem, dist,ligOcc
	minLigDist = [(key,value[0].split(":")[0],value[0].split(":")[1],value[1],value[3],value[-3],value[-2]) for key,value in zip(dists.keys(),dists.values())]
	
	# #keep track of metal and ligand name to get distance values
	dfParms = get_valParms(mElem,lElem)
	# print(dfParms)
	
	metVal = {}
	metVal[ID] = {}
	metVal[ID]['env']=env
	metVal[ID]['Valence'] = 0.0
	metVal[ID]['Ox'] = 0

	# loop through key, metElem, ligElem, dist
	#keep track of OxNum intersection
	oxInt = {}
	#count the keys
	count = 0
	print('test1')
	for key,metElem,ligElem,dist,ligOcc,metVec,ligVec in minLigDist:

		#get oxidative states
		dfOxParm = dfParms[(dfParms.valence_param_atom_1 == metElem) & (dfParms.valence_param_atom_2 == ligElem)].drop_duplicates(subset=['valence_param_atom_1_valence'])
		

		metVal[ID][key] = {}
		metVal[ID][key]['mElem'] = metElem
		metVal[ID][key]['lElem'] = ligElem
		metVal[ID][key]['Ox'] = []
		metVal[ID][key]['bond-valence'] = []
		metVal[ID]['coordNum'] = len(lElem)
		metVal[ID][key]['metVec'] = metVec
		metVal[ID][key]['ligVec'] = ligVec
		metVal[ID][key]['ligOcc'] = ligOcc
		metVal[ID][key]['dist'] = dist
		# key,metElem,ligElem,dist,ox,r0,b
		# loop through standard oxidation values
		
		for	i,oxParm in dfOxParm.iterrows():
			# print(oxParm)
			r0 = oxParm['valence_param_Ro']
			Ox = oxParm['valence_param_atom_1_valence']
			metVal[ID][key]['Ox'].append(Ox)
			# print(type(ligOcc))
			metVal[ID][key]['bond-valence'].append(calc_bond_valence(r0,dist,Ox)*ligOcc)
		#first time through
		if count == 0:
			oxInt = set(metVal[ID][key]['Ox'])
		#multiple times through
		else:
			oxInt = oxInt.intersection(set(metVal[ID][key]['Ox']))
		count +=1
	#how close
	dif = 99
	valency = 0
	oxNum = 0
	# print(dif)
	#loop through oxNum first for all atom_pairs
	for ox in oxInt:
		val = 0
		for key in metVal[ID].keys():
			if key not in ['Ox','Valence','coordNum','env']:
				oxInd = metVal[ID][key]['Ox'].index(ox)	
				val += float(metVal[ID][key]['bond-valence'][oxInd])
				# print(val)
		#update values 
		if abs(val-ox) < dif:
			valency = val
			oxNum = ox
			metVal[ID]['Ox'] = oxNum
			# print(valency,ox)
			dif = abs(val-ox)
	
	metVal[ID]['Valence'] = valency
	# print(metVal)
	return metVal
def calc_bvs(pdbFileName,dists,oxStates):  # inpu
	#calculates bond valence sum
	
	mElem = [value[0].split(":")[0] for value in dists.values()]
	lElem = [value[0].split(":")[1] for value in dists.values()] 
	oxState = oxStates[mElem[0]]
	print(oxState)
	# key, metElem, ligElem, dist,ligOcc
	minLigDist = [(key,value[0].split(":")[0],value[0].split(":")[1],value[1],value[3],value[-3],value[-2]) for key,value in zip(dists.keys(),dists.values())]
	
	# #keep track of metal and ligand name to get distance values
	dfParms = get_valParms(mElem,lElem)
	
	bond_valence = {}
	
	#loop over oxidation states 1 or 2, #_#
	for ox in oxState.split('_'):
		# print(ox)
		if ox == 'none':
			print('Blah',pdbFileName,ox)

		bond_valence[ox] = {}
		bond_valence[ox]['valence'] = 0.0

		for key,metElem,ligElem,dist,ligOcc,metVec,ligVec in minLigDist:

			bond_valence[ox][key] = {} 

			#get oxidative states
			dfOxParm = dfParms[(dfParms.valence_param_atom_1 == metElem) \
						& (dfParms.valence_param_atom_2 == ligElem) \
						& (dfParms.valence_param_atom_1_valence == int(ox))]

			r0 = dfOxParm.iloc[0]['valence_param_Ro'] #select the first in the list
			bond_val = calc_bond_valence(r0,dist,ox)*ligOcc
			bond_valence[ox][key]['bond_val'] = bond_val
			bond_valence[ox]['coordNum'] = len(lElem)
			bond_valence[ox][key]['metVec'] = metVec
			bond_valence[ox][key]['ligVec'] = ligVec
			bond_valence[ox][key]['ligOcc'] = ligOcc
			bond_valence[ox][key]['r0'] = r0
			bond_valence[ox][key]['dist'] = dist
			bond_valence[ox]['valence'] += bond_val

	return bond_valence

def calc_vecsum(metVal,ox):
	# print(valenceDictionary.keys())
	# The borderline and outlier thresholds are >0.10 and >0.23, respectively, for nVECSUM, 
	# >10% and >25%, respectively, for the vacancy parameter, which is the percentage of all expected coordination sites left vacant (Supplementary Fig. 2 and Supplementary Table 2). For example, ions with all coordination sites occupied by ligands (vacancy = 0) are classi- fied as acceptable. For geometry with an expected coordination number greater than four, metals with one vacant coordina- tion site (vacancy ≤ 25%) are borderline, and metals with two or more vacant coordination sites (vacancy > 25%)
	vecsum = 0
	fij = PDB.Vector(x=0,y=0,z=0)
	bonds = [key for key in metVal[ox].keys() if key not in ['coordNum','valence']]
	for bond in bonds:
		distance = metVal[ox][bond]['dist']
		metVec = metVal[ox][bond]['metVec']
		ligVec = metVal[ox][bond]['ligVec']
		# print('metVec',metVec)
		# print('ligVec',ligVec)
		vec = (ligVec - metVec)
		rij = vec.__truediv__(distance)
		ligOcc = metVal[ox][bond]['ligOcc']
		bondValence = metVal[ox][bond]['bond_val']
	# print('blha: ' + str(bondValence))
		sij = float(ligOcc) * bondValence
	# print('sij',sij)
	# raise TypeError('somethingHappend ' + str(ij))
		fij = fij.__add__(np.multiply(rij.get_array(),sij))
		# print('fij: ',fij)
	vecsum =  math.sqrt(fij.__mul__(fij)) / metVal[ox]['valence']
	# print('vecsum: ',vecsum)
	return vecsum
def calc_protein_vecsum(metVal,ID,ox):
	# print(metVal)
	# {'4RT5.NI_201_1655_A.hvp.pdb': {'env': 'NI1C1N2O4', 'Valence': 0, 'Ox': 0, 
	# 'NI_NE2': {'mElem': 'NI', 'lElem': 'N', 'Ox': [2, 3], 'bond-valence': [0.24315593231842803, 0.2644062425923594]}, 
	# 'NI_NE2_43': {'mElem': 'NI', 'lElem': 'N', 'Ox': [2, 3], 'bond-valence': [0.41219316064853856, 0.44821626924826263]}, 
	# 'NI_OE2': {'mElem': 'NI', 'lElem': 'O', 'Ox': [2, 3, 4], 'bond-valence': [0.22138086502773974, 0.28696034883432053, 0.3111966818850807]}, 
	# 'NI_C2': {'mElem': 'NI', 'lElem': 'C', 'Ox': [], 'bond-valence': []}, 
	# 'NI_O2': {'mElem': 'NI', 'lElem': 'O', 'Ox': [2, 3, 4], 'bond-valence': [0.310535384208486, 0.402525042833765, 0.43652183381558524]}, 
	# 'NI_O3': {'mElem': 'NI', 'lElem': 'O', 'Ox': [2, 3, 4], 'bond-valence': [0.10584098720469944, 0.1371941816444871, 0.14878144062458418]}, 
	# 'NI_O': {'mElem': 'NI', 'lElem': 'O', 'Ox': [2, 3, 4], 'bond-valence': [0.08400752496077941, 0.10889300963040541, 0.1180899849582612]}}}
	# print(valenceDictionary.keys())
	# The borderline and outlier thresholds are >0.10 and >0.23, respectively, for nVECSUM, 
	# >10% and >25%, respectively, for the vacancy parameter, which is the percentage of all expected coordination sites left vacant (Supplementary Fig. 2 and Supplementary Table 2). For example, ions with all coordination sites occupied by ligands (vacancy = 0) are classi- fied as acceptable. For geometry with an expected coordination number greater than four, metals with one vacant coordina- tion site (vacancy ≤ 25%) are borderline, and metals with two or more vacant coordination sites (vacancy > 25%)
	vecsum = 0
	fij = PDB.Vector(x=0,y=0,z=0)
	# print(metVal[ID].keys())
	bonds = [key for key in metVal[ID].keys() if key not in ['env','Valence','Ox','coordNum']]
	# print(bonds)
	# oxInd = metVal[ID]['Ox'].index(ox)	
	for bond in bonds:
		# print(ox)
		oxInd = metVal[ID][bond]['Ox'].index(ox)
		# print(oxInd)
		distance = metVal[ID][bond]['dist']
		# print(distance)
		metVec = metVal[ID][bond]['metVec']
		# print(metVec)
		ligVec = metVal[ID][bond]['ligVec']
		# print(ligVec)
		vec = (ligVec - metVec)
		# print(vec) 
		rij = vec.__truediv__(distance)
		# print(rij)
		ligOcc = metVal[ID][bond]['ligOcc']
		# print(ligOcc)
		bondValence = metVal[ID][bond]['bond-valence'][oxInd]
		# print(bondValence)

		# print(bondValence)
	# print('blha: ' + str(bondValence))
		sij = float(ligOcc) * bondValence
		# print(sij)
	# print('sij',sij)
	# raise TypeError('somethingHappend ' + str(ij))
		fij = fij.__add__(np.multiply(rij.get_array(),sij))

		# print('fij: ',fij)
	vecsum =  math.sqrt(fij.__mul__(fij)) / metVal[ID]['Valence']
	# print('vecsum: ',vecsum)
	return vecsum

def get_rmsd(pdbFileName,pdbPathsLoc):

	pdb,metalID,geoShort,ext = re.split('[.]',pdbFileName)
	atomName,resNum,atomNum,chain = re.split('_',metalID)
	atomName = atomName.upper()

	# pdbPathsLoc='/home/kenneth/proj/proMin/proteins/feS/clusters/findGeo'
	# pdbPathsLoc='/home/kenneth/proj/proMin/proteins/rcsb/pdbs_0_1.5/findgeo'
	# pdbPathsLoc='/home/kenneth/proj/proMin/proteins/hagai/pdbs'
	# pdbPathsLoc='/home/kenneth/proj/proMin/minerals/database/amcsd_pdb/pdb_reSeq_res_atom'

	rmsdFile = os.path.join(pdbPathsLoc,pdb,metalID,geoShort+".out")
	print(rmsdFile)
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

	# pdb,rNum,metalID,geoShort,ext = re.split('[.]',pdbFileName)
	# atomName,resNum,atomNum,chain = re.split('_',metalID)
	# atomName = atomName.upper()
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
