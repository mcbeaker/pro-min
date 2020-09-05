import read_cif
import os
from CifFile import ReadCif
import pandas as pd
import numpy as np
from Bio.PDB import MMCIF2Dict
from Bio.PDB.PDBIO import Select
from Bio import PDB
from diffpy import structure 
import gemmi
from progress.bar import Bar, IncrementalBar
import subprocess
import inspect
import re

from elements import ELEMENTS

strElements = [element.symbol.upper() for element in ELEMENTS]

# import urllib.request

import glob

# # if os.path.exists(cifPath) == False:
# #     print(pdb + " needed")
# #     urllib.request.urlretrieve("https://files.rcsb.org/download/"+pdb.upper()+".cif",os.path.join(pdir,pdb+".cif"))

def get_cif_STR(cifPath):
        STR = PDB.MMCIFParser(QUIET=True).get_structure("pdb",cifPath)
        # DICT = PDB.MMCIF2Dict.MMCIF2Dict(cifPath)
        # print(DICT)
        return STR #,DICT

f = "/home/kenneth/proj/proXtal/proteins/hagai/hagai.csv"
pdir = "/home/kenneth/proj/proXtal/proteins/hagai/pdb"

hagaiNames = ["Microenvironment_ID","Binding_motif","Cofactor","Metal","cofactor_group","EC","Head","Molecule","Organism_scientific","no_rank","superkingdom","phylum","class_","order","family","genus","organism_taxid","name","chains","Resolution","Structure_method","Keywords","Journal_reference","Release_date"]
df = pd.read_csv(f,names=hagaiNames)
feDF = df[df['Metal'].str.contains("FE") & df['Structure_method'].str.contains("x-ray")]


# print(len(df.index))
bar = IncrementalBar("Processing",max=len(feDF.index))#,fill='*')#),suffix='%(percent).1f%% - %(eta)ds')

with open(os.path.join(pdir,"micro.csv"),"w+") as outWriter:
    for i in feDF.itertuples():
        pdb,LIG,chain,resID,function = re.split('[._]',i.Microenvironment_ID)
        # print(pdb,LIG,chain,resID,function,i.Metal)
        cifPath = os.path.join(pdir,pdb+".cif")
        STR = get_cif_STR(cifPath)

        try:

            RES = STR[0][chain][("H_"+LIG,int(resID)," ")]

        except Exception as e:
                print("  ",i.Microenvironment_ID,("H_"+LIG,int(resID)," "))
                print(e)
                continue
            
        resName = RES.get_resname()
        feATOMS = [atom for atom in RES.get_atoms() if "FE" in atom.get_name()] 

        atom_list = PDB.Selection.unfold_entities(STR,"A")
        ns = PDB.NeighborSearch(atom_list)

        for atom in feATOMS:
            atomName = atom.get_name()
            feAtom = {"ligName":resName,"ligChain":chain,"resID":resID,"ligAtom":atomName,"Metal":i.Metal,"ligCoords":atom.get_coord()}  
            annotations = {attr:getattr(i,attr) for attr in hagaiNames}
            neighbors = ns.search(atom.get_coord(),3,level='A')
            neigh_dict = {i.Microenvironment_ID:{"feAtom":feAtom,"Microenvironment_annotations":annotations,"Neighbors":{"Neighbor_"+str(idx): {"neigh_"+str(idx)+"_resName":neigh.get_parent().get_resname(),"neigh_"+str(idx)+"_atomName":neigh.get_name(),"neigh_"+str(idx)+"_resNum":neigh.get_parent().get_id()[1],"neigh_"+str(idx)+"_atomCoord":neigh.get_coord(),"neigh_"+str(idx)+"_distance":neigh-atom} for idx, neigh in enumerate(neighbors) if neigh-atom != 0 }}}
            outWriter.write(str(neigh_dict)+"\n")
            # print(neigh_dict)
                    # print(neigh_dict['feAtom'].keys())
                    # for feKey,feVal,nKey,nVal in zip(neigh_dict['feAtom'].keys(),neigh_dict['feAtom'].values(), \ 
                    #                                 neigh_dict["Neighbors"].keys(),neigh_dict["Neighbors"].values()):
                        
                        # outWriter.writerow([feKey,feVal,nKey,nVal])

        bar.next()
    bar.finish()



# for f in feCIFFiles:




#     #loop through lig atoms, get fe coords, get nearest neighbors of fe's
#     for lig in ligands:

#         for atom in lig.get_atoms():
#             # print(atom.get_name()) 
#             if "FE" in atom.get_name():

#                  

#                 


#         minName = os.path.splitext(os.path.basename(f))[0]
        
#         minPath = f

#         STR = PDB.MMCIFParser(QUIET=True).get_structure("STR",minPath)
#         print([atom for atom in grPDB.get_atoms()])
#         # gr_atom_list = PDB.Selection.unfold_entities(grPDB)
#         # print(gr_atom_list)
#         # gr_Fe = [atom for atom in grPDB.get_atoms() if "FE" in atom.get_name().upper()]
#         # gr_ns = PDB.NeighborSearch(gr_atom_list)
#         # gr_feNeighbors = [gr_ns.search(atom.get_coord(),5,level='A') for atom in gr_Fe]
#         # # print(gr_feNeighbors[0])
#         # # print(gr_feNeighbors[1])

#         # for id,a1 in enumerate(gr_Fe):

#         #     for a2 in gr_feNeighbors[id]:
#         #         if a1-a2 < 3 and a1-a2 != 0:
#         #             print("%s %s %s %s" %(id,a1,a2,a1-a2))
#         #     break

#         # if os.path.exists(outPDB) == True:
#         # # if os.path.exists(outXYZ) == False:
#         #     try:
                
#         #         # print(minPath)
#         #         Str = structure.loadStructure(f)
#         #         # Str.assignUniqueLabels()
#         #         Str.write(outPDB,"pdb")
#         # #         Str.write(outXYZ,"xyz")
#         # #         # print(outCIF)
#         # #         # cif = ReadCif(minPath)
#         # #         # print(dir(cif))
#         # #         # break

#         #     except Exception as e:
#         #         print(minPath)
#         #         print(e)
#         #         exit()
#         bar.next()
# bar.finish()




# print(grStr.element)
# print([id for id,elem in enumerate(grStr.element) if "FE" in elem.upper()])
# grSubCells = gemmi.SubCells(grStr[0],grStr.cell,5)
# print(dir(grStr[0]))

# gr = "/home/kenneth/proj/proXtal/amcsd/cif/mineral/Fe/pdb/greigite_2.pdb"
# grPDB = PDB.PDBParser().get_structure("gr","/home/kenneth/proj/proXtal/amcsd/cif/mineral/Fe/pdb/greigite_2.pdb")
# gr = "/home/kenneth/Downloads/greigite.pdb"
# grStr = structure.loadStructure(gr)
# grStr.assignUniqueLabels()
# print(grStr.label)
# grStr.write("/home/kenneth/proj/proXtal/amcsd/cif/mineral/Fe/pdb/greigite_2.pdb","pdb")

# gr = "/home/kenneth/proj/proXtal/amcsd/cif/mineral/Fe/pdb/greigite_2.pdb"
# grPDB = PDB.PDBParser().get_structure("gr",gr)
# gr_atom_list = PDB.Selection.unfold_entities(grPDB,'A')
# gr_Fe = [atom for atom in grPDB.get_atoms() if "FE" in atom.get_name().upper()]
# gr_ns = PDB.NeighborSearch(gr_atom_list)
# gr_feNeighbors = [gr_ns.search(atom.get_coord(),5,level='A') for atom in gr_Fe]
# # print(gr_feNeighbors[0])
# # print(gr_feNeighbors[1])

# for id,a1 in enumerate(gr_Fe):
#     for a2 in gr_feNeighbors[id]:
#         if a1-a2 < 2.5 and a1-a2 != 0:
#             print("%s %s %s %s" %(id,a1,a2,a1-a2))
#     break


# print(dir(fdnPDB))
# print([chain for chain in fdnPDB.get_chains()])

# gr = ReadCif(gr)
# fdn = PDB.MMCIFParser().get_structure("fdn","/home/kenneth/proj/proXtal/proteins/1fdn.cif")
# print()
# print([chain for chain in fdn.get_chains()])
# print(fdn[0])

# resi = ['SF4']

# FE_resi = [resi for resi in fdn.get_residues() if resi.get_id()[0] == "H_SF4"]

# for feRes in FE_resi:
#     for a1 in feRes.get_atoms():
#         for a2 in feRes.get_atoms():

#             if a1-a2 < 3.0 and a1-a2 != 0: #and ("S" not in a1.get_name() and "S" not in a2.get_name()):

#                 print("%s %s %s %s %s" %(feRes.get_id()[1],feRes.get_resname(),a1.get_name(),a2.get_name(),a1-a2))


# print(FE_resi)
# S = [atom for atom in fdn.get_atoms() if "S" in atom.get_name() and "SG" not in atom.get_name()]

# for f,s in zip(FE,S):
        # print("%s %s %s" %(f.get_name(),s.get_name(),f-s))

# metals = ["Mg", "Ca", "Fe", "Zn", "Cu", "Mn", "Co"]





# for model in fdn:
#     for chain in model:
#         for residue in chain:
#             for atom in residue:
#                 print(atom)

# print(inspect.getmembers(fdn_doc.get_structure))
# fdnStruct = fdn_doc.get_structure
# for i in fdnStruct.get_atoms():
#     print(i)
# gr_block = gr_doc.sole_block()
# fdn_block = fdn_doc.sole_block()

# gr_elems = set(gr_block.find_loop("_atom_site.type_symbol"))
# print("gr" + ' ' + ' '.join(gr_elems))

# fdn_elem = set(fdn_block.find_loop("_atom_site.type_symbol"))
# print("fdn" + ' ' + ' '.join(fdn_elem))


# for f in os.listdir(d):

#     filePath = d+f
#     print(filePath)
#     df_cif = read_cif.diff_read_cif(filePath)

#     # from diffpy import lattice
#     from scipy import spatial

    
#     # print(spatial.distance.cdist(df_cif.xyz_cartn,df_cif.xyz_cartn,'euclidean'))
#     # print(df_cif.element)
#     # print(df_cif.label)

#     df = pd.DataFrame(data=np.triu(spatial.distance.cdist(df_cif.xyz_cartn,df_cif.xyz_cartn,'euclidean')),columns=df_cif.element,index=df_cif.label)

#     df.to_csv(d+f+'.csv')

#     # print(df_cif.lattice.distance('K','K_2'))

#     # print(dir(df_cif))
#     # print(df_cif)
#     # print(df_cif.xyz_cartn)
#     # print(df_cif.label)
#     # print(df_cif.element)

#     break
    # cf = ReadCif(d+f)
    
    # if '_atom_site_occupancy' in cf['global']:

    #     mineral = cf['global']['_chemical_name_mineral']
    #     print(mineral)
    #     a = float(cf['global']['_cell_length_a'])
    #     b = float(cf['global']['_cell_length_b'])
    #     c = float(cf['global']['_cell_length_c'])

    #     alpha = float(cf['global']['_cell_angle_alpha'])
    #     beta = float(cf['global']['_cell_angle_beta'])
    #     gamma = float(cf['global']['_cell_angle_beta'])


    #     # print(a,b,c,alpha,beta,gamma)
    #     labels = cf['global']['_atom_site_label']
    #     us = cf['global']['_atom_site_fract_x']
    #     vs = cf['global']['_atom_site_fract_y']
    #     ws = cf['global']['_atom_site_fract_z']
    #     occ = cf['global']['_atom_site_occupancy']

    #     processed = zip(labels,us,vs,ws,occ)

    #     lxyz = []

    #     #only keep the highest occupancy if sharing the same xyz coords

    #     setXYZO = set()
    #     dirXYZ = {}
    #     for label,u,v,w,occ in processed:
    #         r = read_cif.fract_to_cart(a,b,c,float(u),float(v),float(w))
    #         if (r[0],r[1],r[2]) not in setXYZO:
    #             dirXYZr[0],r[1],r[2])] = float(occ)
    #             setXYZO.add((r[0],r[1],r[2]))
    #         else:
    #             if float(occ) > dirXYZ[(r[0],r[1],r[2])]:
    #                 dirXYZ[(label,r[0],r[1],r[2])] = float(occ)
    #     print(dirXYZ)
            # lxyz.append([label,r[0],r[1],r[2],occ])
            
        # print(fe_lxyz)
        # print(other_lxyz)
            # print(xyz)
            # print("%s %s %s %s %s %s %s %s %s %s %s %s %s" % (label,a,b,c,alpha,beta,gamma,u,v,w,r[0],r[1],r[2]))
            # print(xyz)


    #     for i in range(0,len(lxyz)):
    #         for j in range(i+1,len(lxyz)):
    #             l1 = lxyz[i][0]
    #             l2 = lxyz[j][0]
    #             x1 = lxyz[i][1]
    #             y1 = lxyz[i][2]
    #             z1 = lxyz[i][3]
    #             x2 = lxyz[j][1]
    #             y2 = lxyz[j][2]
    #             z2 = lxyz[j][3]
    #             print("%s %s %s %s %s %s %s %s %s\n" %(mineral,l1,l2,x1,x2,y1,y2,z1,z2))
    #             xyz = read_cif.xyz_dist(x1,y1,z1,x2,y2,z2)
    #             if xyz != 0:
    #                 outFile.write("%s %s %s %s %s %s %s\n" %(mineral,l1,l2,x2-x1,y2-y1,z2-z1,xyz))
    # outFile.close()        