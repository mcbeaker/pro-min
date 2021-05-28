from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.MMCIFParser import MMCIFParser

import os

def get_STR(filePath,fileType='pdb'):
    if fileType == "pdb":
        STR = MMCIFParser(QUIET=True).get_structure("pdb",filePath)
        return STR

    if fileType == "cif":
        DICT = MMCIF2Dict(filePath)
        # print(DICT)
        return DICT 
    else:
        raise TypeError("%s is not a valid fileType" %fileType)
import glob

pdir = "/home/kenneth/proj/proXtal/proteins/hagai/cif"



for pdb in os.listdir(pdir):
    if pdb.endswith(".cif"):
        STR = get_STR(os.path.join(pdir,pdb),fileType='cif')

        # print(dir(STR))
        # print(STR['_pdbx_database_PDB_obs_spr.pdb_id'],STR['_entity.pdbx_ec'],sep=',')
        # print(STR.keys())
        print(str(STR['_entry.id']).lower(),STR['_entity.pdbx_ec'][0])

# _entity.pdbx_ec