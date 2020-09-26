#!/usr/bin/env python
import os
# from CifFile import ReadCif
import pandas as pd
import numpy as np
from Bio import PDB
import subprocess
import re
import multiprocessing as mp
from collections import defaultdict
from tqdm import tqdm
from collections import Counter
import sys
import tools

#global variables
global pdbPaths
pdbPathsLoc='/home/kenneth/proj/proMin/proteins/rcsb/pdbs_0_1.5/findgeo/combineFindGeoResults'
# pdbPathsLoc = '/home/kenneth/proj/proMin/minerals/database/amcsd_pdb/pdb_reSeq_res_atom/combineFindGeoResults'
# tools.get_pdbFilePaths(pdbPathsLoc,'pdbFileNames.csv')
pdbPaths = pd.read_csv(os.path.join(pdbPathsLoc,'pdbFileNames.csv'),index_col=0)
# print(pdbPaths)
# exit()

def f_init(q):
	process_queue.q = q

#test one to three
def process_queue(id):

	pdbPath = pdbPaths.iloc[id]['pdbPath']
	# print(pdbPath)
	pdbFileName = os.path.basename(pdbPath)
	pdb,metalName,resNum,atomNum,chain,geoShort,ext = re.split('[._]',pdbFileName)
	# print(pdbPath)
	# print(re.split('[.]',pdbFileName))
	
	# mineral,Mn_30_30_A_irr,geoShort,ext = re.split('[.]',pdbFileName)
	# metalName,resNum,atomNum,chain = re.split('[._]',Mn_30_30_A_irr)
	
	metalName = metalName.upper()

	import warnings
	warnings.filterwarnings("error")
	# print(pdb,metalName)
	output=""
	try:
		STR = PDB.PDBParser(QUIET=False).get_structure('pdb',pdbPath)

		environment = tools.get_environment(STR,pdbFileName)

		gRMSD = tools.get_rmsd(pdbFileName)

		angles = tools.get_angles(STR,metalName)

		for i in angles:
			# print('Protein,'+i)
			# process_queue.q.put('Mineral,'+pdbFileName + ',' + i + ',' + environment)
			process_queue.q.put('Protein,'+pdbFileName + ',' + i + ',' + environment)

		# atomDist = tools.get_atomDistances(STR,metalName)

		# for key in atomDist.keys():
		# 	metName,ligName = key.split("_")
		# 	metElm, ligElm = atomDist[key][0].split(':')
		# 	dist = atomDist[key][1]
		# 	metOcc = atomDist[key][2]
		# 	ligOcc = atomDist[key][3]
		# 	metResName, ligResName = atomDist[key][4].split("_")
		# 	metResNum, ligResNum = atomDist[key][5].split("_") 
		# 	metChain, ligChain = atomDist[key][6].split("_") 
		# 	output = pdbFileName + ',' + metElm + "," + ligElm + "," + str(round(dist,2)) + "," + metName + "," + ligName + "," +  str(metOcc) + "," + \
		# 				str(ligOcc) + "," + metResName + "," + ligResName + "," + str(metResNum) + "," + ligResNum + "," + \
		# 				metChain + "," + ligChain 
		# 	process_queue.q.put(output)


					
		# print(atomDist)
	# print(atomDist)
	# print(metalElements)
	# print(ligElements)
	# print(atomDist)
		# metVal = tools.calc_bvs(atomDist)

		# vecsum = tools.calc_vecsum(metVal,STR,metalName)

		# output = os.path.splitext(pdbFileName)[0] + " " + environment + " " + str(gRMSD) + " " + str(metVal['valency'])
		
		# process_queue.q.put(output)
		
	except PDB.PDBExceptions.PDBConstructionWarning as e:
		print(e)
		print(pdbPath)
	except TypeError as t:
		print(t)
		print(pdbPath)
		exc_type, exc_obj, exc_tb = sys.exc_info()
		fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
		print(exc_type, fname, exc_tb.tb_lineno)

	# print(output)
	return output

def listener():
	'''listens for messages on a the q, writes to a file'''
	#https://stackoverflow.com/questions/13446445/python-multiprocessing-safely-writing-to-a-file

	pdir = '/home/kenneth/proj/proMin/proteins/rcsb/pdbs_0_1.5/findgeo/combineFindGeoResults/analysis'
	# pdir = '/home/kenneth/proj/proMin/proteins/hagai/pdbs/combineFindGeoResults/analysis'
	# pdir = '/home/kenneth/proj/proMin/minerals/database/amcsd_pdb/pdb_reSeq_res_atom/combineFindGeoResults/analysis'
	# outFile='min_fg_dictionaryValues_ligOcc_angles.csv'
	outFile='pro_fg_dictionaryValues_ligOcc_angles.csv'
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
			print(m)
            # print(type(m))
			if m == "kill":
				outWriter.write('killed')
				break
			else:
				outWriter.write(m + '\n')
				outWriter.flush()

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
        for i, job in tqdm(enumerate(pool.imap_unordered(process_queue,list(pdbPaths.index)))):
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