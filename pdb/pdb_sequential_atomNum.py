from p3d.protein import Protein
import glob
from shutil import copyfile
import os
import pandas as pd
#needs to be on the level of findgeo.input


basePath='/home/kenneth/proj/proMin/minerals/database/amcsd_pdb'

pdbNames=os.path.join(basePath,'pdbName.csv')
dfPDBs=pd.read_csv(pdbNames,header=None,index_col=False)

for i,row in dfPDBs.iterrows():

    pdbFile=os.path.join(basePath,row[0]+'.pdb')
    outPDB = os.path.join(basePath,'pdb_reSeq_res_atom',row[0]+'.pdb')
    copyfile(pdbFile,pdbFile+'.bak')
    pdb = Protein(pdbFile)

    for n in range(0,len(pdb.atoms)):
        pdb.atoms[n].idx = n
        pdb.atoms[n].resid = n
        # print(pdb.atoms[n].aa)
        # print(pdb.atoms[n].atype)
        # print(dir(pdb.atoms[n].info))

        # print(dir(pdb.atoms[n]))
        # exit()
        if pdb.atoms[n].atype.strip().find('Wat') > -1:
            pdb.atoms[n].elementType = 'O'
        hAt = ['Hl','HhR','HhL','Hh1R','HhR,''Hl','Hwl','Hw2','HOh1','HOh8','Hw11','Hw12','Hw21','Hw22','HwL','HwR','Hw1L','Hw1R','Hh1L','HW1a','HW1b','HW2a','HW2b','HW3']
        # print(type(pdb.atoms[n].atype))
        # print(pdb.atoms[n].atype)
        # t = ['Pb1']
        # if str(pdb.atoms[n].atype.strip()) in t :
            # print("Blah")
        # exit()
        if pdb.atoms[n].atype.strip() in hAt:
            pdb.atoms[n].elementType = 'H'
        if pdb.atoms[n].atype.strip().find('Op') > -1:
            pdb.atoms[n].elementType = 'P'
        if pdb.atoms[n].atype.strip().find('Zn') > -1:
            pdb.atoms[n].elementType = 'ZN'
        if pdb.atoms[n].atype.strip().find('Pb') > -1:
            pdb.atoms[n].elementType = 'PB'
        if pdb.atoms[n].atype.strip().find('Bi') > -1:
            pdb.atoms[n].elementType = 'BI'
        if pdb.atoms[n].atype.strip().find('Ba') > -1:
            pdb.atoms[n].elementType = 'BA'
        if pdb.atoms[n].atype.strip().find('Al') > -1:
            pdb.atoms[n].elementType = 'AL'
        if pdb.atoms[n].atype.strip().find('Sb') > -1:
            pdb.atoms[n].elementType = 'SB'
        if pdb.atoms[n].atype.strip().find('As') > -1:
            pdb.atoms[n].elementType = 'AS'
        if pdb.atoms[n].atype.strip().find('Ti') > -1:
            pdb.atoms[n].elementType = 'TI'
        

        
    pdb.writeToFile(outPDB) 

# pdbPathsLoc='/home/kenneth/proj/proMin/minerals/database/amcsd_pdb/FindGeoSummative_noHeader.csv'
# df=pd.read_csv(pdbPathsLoc,header=None,index_col=False)

# for i,row in df.iterrows():
#     pdb = row[0].split('.')[0]
#     metal = row[4]
#     geo = row[2]
#     pdbPath=os.path.join(basePath,pdb,metal,"findgeo.input")
#     pdbOUT = os.path.join(basePath,'combineFindGeoResults',pdb+'.'+metal+'.'+geo+'.pdb')
#     print(pdbPath)
#     print(pdbOUT)
    

#         # print(dir(pdb.atoms[0]))
#         exit()
    # exit()

# count = 0 

# # metals = ["FE_"]
# for pdbFile in glob.iglob(pdbPathsLoc+'/*.pdb'):
#     #101m.A_FE_155_FE_FE,1270,spy,square pyramid (regular),FE_155_1270_A
#     #echo $line
#     line=`echo $line | tr -d '\n'`
#     # echo $line
#     pdb=`echo $line | cut -d'.' -f 1`
#     geo=`echo $line | cut -d',' -f 3`
#     metal=`echo $line| cut -d',' -f 5`
#     pdb,metalID,geoShort,ext = re.split('[.]',pdbFileName)
# 	atomName,resNum,atomNum,chain = re.split('_',metalID)
	
#     pdbPath=$basePath/$pdb/$metal/"findgeo.input"

#     outPDB=$combineFindGeoFolder/$pdb.$metal.$geo.pdb
#     atomName = atomName.upper()
#     backup = pdbFile+'.bak'
#     copyfile(pdbFile,backup)
#     print(pdbFile)    
#     pdb = Protein(pdbFile)

#     for n in range(0,len(pdb.atoms)):
#     # print(pdbFile)
#     # print(dir(pdb))
#     # print(type(pdb))
#     # print(dir(at))
#         pdb.atoms[n].idx = n
#     # print(at.pos_in_list)
#         pdb.atoms[n].resid = n
#     # exit()
#     pdb.writeToFile(pdbFile)