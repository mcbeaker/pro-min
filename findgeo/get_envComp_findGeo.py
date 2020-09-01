from progress.bar import Bar
import shutil
import pandas as pd
from collections import Counter
from Bio.PDB import MMCIF2Dict
from Bio.PDB.PDBIO import Select
from Bio import PDB
import os
import re

'''
#Format
#pdbID,atomNum,geoShort,geoLong,metalID (FE_155_1270_A) - metalAtomName_resNum_atomNum_chain
#Format pdbID "${pdb}.${chain}_${atomName}_${resNum}_${atomName}_${atomName},
# ${atomName},${geoAbbr},${geoFull},
# ${metalID}" Format- metalAtomName_resNum_atomNum_chain

#ouput
#pdbID,element,envComp,geoShort,geoLong,resNum,atomNum,atomName,resName
'''
pwd = os.join('home','kenneth','Box','proj','proXtal','data','proteins','findGeoResults')
proFile = os.join('home','kenneth','Box','proj','proXtal','data','proteins','FindGeoSummative_wo_CA_K_NA_ZN.csv')
df = pd.read_csv(proFile,header=0,index_col=False)

#metal information
metals = ["CO","CU","FE","MG","MN","MO","NI","V","W"] #added MG 06/18/2020

outPDBTypeError = open(os.path.join(pwd,"outPDBTypeError.txt"),'w') #these PDBs do not exist

#new columns
atomEnvironmentComposition = []
atomNum = []
atomNames = []
elements = []
resNames = []
resNums = []
pdbIDs = []
chains = []

#Get atom environment composition (element, residueName)

# count = 0
with Bar('Processing',max=len(df.index)) as bar:

    for row in df.itertuples():
        
        pdb,chain_metal_resNum_atomName_element = re.split('[.]',row.pdbID)
        chain,metal,resNum,atomName,element = re.split('_',chain_metal_resNum_atomName_element)

        #update columns
        pdbIDs.append(pdb)
        chains.append(chain)
        atomNames.append(atomName.upper())


    #     if pdbPath not in sortNotFound and pdbPath not in outPDBTypeError:
    # #             #print(pdbPath)

    #             # try:
        print(pdbPath)
    #     STR = PDB.PDBParser(QUIET=False).get_structure("pdb",pdbPath)
        
    #     # df_STR = pd.DataFrame(columns=['resName','resNum','atomName','atomNum','atomElement'])
    #     df_STR = pd.DataFrame(columns=['atomName','atomNum','atomElement'])

    #     #get residueName,residueNum,atomName,atomNum,atomElement
    #     # print(pdb,[a.element for a in STR.get_atoms()],'\n')
    #     print(pdb,[a.get_serial_number() for a in STR.get_atoms()],'\n')
    #     for atom in STR.get_atoms():
            
    #         atomDict = {#'resName':atom.get_parent().get_resname(),
    #                     #'resNum': atom.get_full_id()[3][1],
    #                     'atomName':atom.get_id(),
    #                     'atomNum': atom.get_serial_number(),
    #                     'atomElement': atom.element}
    #         # print(atomDict)
    #         df_STR = df_STR.append(atomDict,ignore_index=True)

    #     # #position of metal
    #     element = df_STR.loc[df_STR['atomName'] == atomName]['atomElement']
    #     # # print(element)
    #     elements.append(element.item())
    #     # # print(element)
    #     countElements = Counter(df_STR.atomElement)
    #     # # print(df_STR.atomNum,'\n')
    #     # # print(df_STR.atomName,'\n')
    #     # # print(df_STR.atomElement,'\n')
    #     # print(pdbPath,'\n', countElements)
        
    #     # resName = df_STR.loc[df_STR['atomName'] == atomName]['resName']
    #     # resNames.append(resName)
    #     # resNum = df_STR.loc[df_STR['atomName'] == atomName]['resNum']
    #     # resNums.append(resNum)
        
    #     # print(df_STR)
            

    #     strCountElements = ""

    # # #     #sort based on element, returns a list of sorted keys 
    #     sortCountElement_keys = sorted(countElements.keys(),reverse=True)
        
    # #     #string metals first
    # #     print(str(set(countElements.keys()).intersection(metals)))
    #     for m in set(countElements.keys()).intersection(metals):
    #         strCountElements += m+str(countElements[m])

    #     # if len(set(countElements.keys()).intersection(metals)) > 1:
    #         # print(pdbID,set(countElements.keys()).intersection(metals))


    # #     #string rest of elements
    #     for m in set(countElements.keys()).difference(metals):
    #         strCountElements += m+str(countElements[m])
            
    # # #     # if pdb == '1a6l':
    # # #             # print(countElements,strCountElements)
    # # #     # print(countElements)
    # # #     # print(sortCountElement)            

    # # #     # print(strCountElements)
    #     atomEnvironmentComposition.append(strCountElements)
    # #     # else:
    # #     #     atomEnvironmentComposition.append("NAN")
    # #     #     elements.append("NAN")
    # #     #     resNames.append("NAN")
    # #     #     atomNames.append("NAN")
    # #     #     pdbIDs.append("NAN")

        bar.next()
        # print('\n')
bar.finish()

# add new column to data fralt
# me with env. composition
# df['envComposition'] = atomEnvironmentComposition
# df['element'] = elements
# # df['resName'] = resNames
# df['pdbID'] = pdbIDs
# df['atomName'] = atomNames 

# df.sort_values(by=['element','geoShort','envComposition','pdbID'],ascending=True,inplace=True)

#did not execute first time
# df = df[['atomID','resNum','atomNum','chain']] = df.metalID.str.split('_',expand=True)
# df = df.drop_duplicates(subset=['pdbID','envComposition','element','geoShort'],keep='first')

# df = df.drop(['metalID'],axis=1)
# df = df.dropna(subset=['envComposition','element'])
df.to_csv(os.path.join(pwd,'FindGeoSummative_wEnvComp.csv'),index=False)

# df.to_csv('/home/kenneth/software/findgeo/FindGeoSummative_wEnvComp.csv',index=False)

#             outListDict.append({"ID":row.pdbID+"_"+str(row.atomNum),"pdb":pdb,"chain":chain,"metal":metal,"resNum":str(resNum),
#             "resID":resID,"atomNum":str(row.atomNum),"atomID":atomName,
#             "geoShort":row.geoShort,"geoLong":row.geoLong,"strCountElements":strCountElements})

# with open(os.path.join(pwd,"FindGeoAtomAnalysis.csv"),"w") as file:
#     file.write(json.dump(outListDict))

#         # except TypeError:
#             # outPDBTypeError.write(pdbPath + '\n')
            
#         #    print(atom.element)
#             #exit(0)

# #"101m.A_FE_155_FE_FE,1270,spy,square pyramid (regular),FE_155_1270_A"
# #print(df.columns)
# #print(df[3].value_counts())
# #print(df[df[3] == 'Irregular geometry'])

# metals = ["FE","CU","NI","MO"]