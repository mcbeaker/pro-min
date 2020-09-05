#!/usr/bin/env python
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBList import PDBList
from Bio.PDB import vectors
import numpy as np
import os

# Reading a PDB file

#get PDB from rcsb.org
basePath = '/home/kenneth/proj/proMin/code/biopython'
filename ='pdb1fdn.ent'
pdb=PDBList().retrieve_pdb_file("1FDN",file_format='pdb')

#Create a PDBParser object
parser = PDBParser() #PERMISSIVE = 0 will list all errors with PDB file
structure = parser.get_structure("1FDN", os.path.join(basePath,filename))

# print(type(structure)) #what type of object did the parser return
# <class 'Bio.PDB.Structure.Structure'>

# print(dir(structure)) #check what attributes exist
# ['__class__', '__contains__', '__delattr__', '__delitem__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__getitem__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__iter__', '__le__', '__len__', '__lt__', '__module__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', '_generate_full_id', '_id', '_reset_full_id', 'add', 'child_dict', 'child_list', 'copy', 'detach_child', 'detach_parent', 'full_id', 'get_atoms', 'get_chains', 'get_full_id', 'get_id', 'get_iterator', 'get_level', 'get_list', 'get_models', 'get_parent', 'get_residues', 'has_id', 'header', 'id', 'insert', 'level', 'parent', 'set_parent', 'transform', 'xtra']

# # PDB Structure object, layers 1)model which contains 2)chains which contains 3)residues which contains 4)atoms
model = structure[0]
chain = model["A"]
# print(list(chain.get_residues()))

#dealing with hetero atom - http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec201
residue = chain[("H_SF4",56," ")] #H_SF4 signifies that this is a non-standard residue in the PDB file, this stands for Hetero_SF4 
# PDB information directly from the file
# HETATM  385 FE1  SF4 A  56       2.260  14.394  74.719  1.00 10.77          FE
# HETATM  386 FE2  SF4 A  56       0.728  16.628  74.564  1.00 10.52          FE
# HETATM  387 FE3  SF4 A  56      -0.061  14.381  73.264  1.00 10.25          FE
# HETATM  388 FE4  SF4 A  56      -0.098  14.443  76.081  1.00 11.02          FE
# HETATM  389  S1  SF4 A  56      -1.348  15.643  74.660  1.00 10.48           S
# HETATM  390  S2  SF4 A  56       0.583  12.793  74.702  1.00 10.90           S
# HETATM  391  S3  SF4 A  56       1.661  15.767  76.430  1.00  9.67           S
# HETATM  392  S4  SF4 A  56       1.808  15.575  72.849  1.00  9.54           S

# print(residue)
# <Residue SF4 het=H_SF4 resseq=56 icode= >
# print(dir(residue))
# ['__class__', '__contains__', '__delattr__', '__delitem__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__getitem__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__iter__', '__le__', '__len__', '__lt__', '__module__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', '_generate_full_id', '_id', '_reset_full_id', 'add', 'child_dict', 'child_list', 'copy', 'detach_child', 'detach_parent', 'disordered', 'flag_disordered', 'full_id', 'get_atom', 'get_atoms', 'get_full_id', 'get_id', 'get_iterator', 'get_level', 'get_list', 'get_parent', 'get_resname', 'get_segid', 'get_unpacked_list', 'has_id', 'id', 'insert', 'is_disordered', 'level', 'parent', 'resname', 'segid', 'set_parent', 'sort', 'transform', 'xtra']

fe1_atom = residue['FE1']
# print(fe1_atom)
# <Atom FE1>
# print(dir(fe1_atom))
# ['__class__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__le__', '__lt__', '__module__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__sub__', '__subclasshook__', '__weakref__', '_assign_atom_mass', '_assign_element', '_sorting_keys', 'altloc', 'anisou_array', 'bfactor', 'coord', 'copy', 'detach_parent', 'disordered_flag', 'element', 'flag_disorder', 'full_id', 'fullname', 'get_altloc', 'get_anisou', 'get_bfactor', 'get_coord', 'get_full_id', 'get_fullname', 'get_id', 'get_level', 'get_name', 'get_occupancy', 'get_parent', 'get_serial_number', 'get_sigatm', 'get_siguij', 'get_vector', 'id', 'is_disordered', 'level', 'mass', 'name', 'occupancy', 'parent', 'serial_number', 'set_altloc', 'set_anisou', 'set_bfactor', 'set_coord', 'set_occupancy', 'set_parent', 'set_serial_number','set_sigatm', 'set_siguij', 'sigatm_array', 'siguij_array', 'transform', 'xtra']

# print(fe1_atom.get_vector)
# <bound method Atom.get_vector of <Atom FE1>>

fe1Vec = fe1_atom.get_vector()
# print(fe1_atom.get_vector())
# <Vector 2.26, 14.39, 74.72>

# fe2_atom = residue['FE2']
# fe2Vec = fe2.get_vector()

#calculate S3-FE4-S2  angle
S3_atom = residue['S3'] #need to know which atoms (atom name i.e. "FE4") are connected and which atom is the center atom in angle
S3Vec = S3_atom.get_vector()

FE4_atom = residue['FE4']
FE4Vec = FE4_atom.get_vector()

S2_atom = residue['S2']
S2Vec = S2_atom.get_vector()


angleS3_SF4_S2 = vectors.calc_angle(S3Vec,FE4Vec,S2Vec)
# 1.867073476442546 #in radians
angleDegrees = (angleS3_SF4_S2 * 180/np.pi) #convert to degrees
print(angleDegrees)


# residue = chain
# print(residue)
# atom = residue["SF4"]

# print(atom.get_list())