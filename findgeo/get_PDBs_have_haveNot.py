import os

wd='/Users/ken/Box/proj/proXtal/data/proteins/findGeoResults'
proFile='/Users/ken/Box/proj/proXtal/FindGeoSummative_wo_CA_K_NA_ZN.csv'

pdbHave = set()
pdbNotHave = set()

with open(proFile,"r") as findGeo:
    
    for i in findGeo.readlines()[1:]:
        i = i.rstrip()
        pdb = i.split(".")[0]
        geo = i.split(",")[2]
        metal = i.split(",")[4]
        # print(repr(metal))
        pdbPath = os.path.join(wd,pdb,metal,"findgeo.input")
        # outPDB= os.path.join(wd,"combFindGeoPDB",pdb+"."+metal+"."+geo+".pdb")
        # print(pdbPath)    
        if os.path.exists(pdbPath):
            print(pdb)
            pdbHave.add(pdb)

        else:
            print(pdbPath)
            pdbNotHave.add(pdb)

with open(wd+"/"+"havePDBs.txt","w+") as outPDB:
    for i in pdbHave:
        outPDB.write(i+"\n")

with open(wd+"/"+"pdbNotHave.txt","w+") as outPDB:
    for i in pdbNotHave:
        outPDB.write(i+"\n")

