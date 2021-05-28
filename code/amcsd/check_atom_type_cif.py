from diffpy import structure 
from progress.bar import Bar
import os
from collections import Counter

finalDir = '/home/kenneth/proj/proXtal/amcsd/xCif/singleOccupancy/natural/containsTransitionMetals/finalDataset'

finalData = '/home/kenneth/proj/proXtal/amcsd/xCif/singleOccupancy/natural/containsTransitionMetals/finalDataset/amcsd_final_dataset.txt'
atom_labels = set()

atomTypeSub = {"ALL":"AL","ASL":"AS","BAL":"BL","BIA":"BA",\
                "BIB":"BB","BIM":"BM","HHL":"BL","HHR":"HR",\
                "HWL":"HW","HWR":"HR","PBB":"PB","PBL":"PL",\
                "PBM":"PM","SBL":"SL","TIL":"TL","ZNI":"ZN"}

metals = ["CO","CU","FE","MN","MO","NI","V","W","ZN"]

with Bar('Processing',max=369) as bar:        
    with open(finalData,"r") as f:
        for line in f.readlines():
            # print(line[:-1])
            fileName = line[:-1]
            if 'Blyth' not in fileName:
                cf = structure.loadStructure(os.path.join(finalDir,fileName))
                cf.assignUniqueLabels

                # print(cf[1].__sizeof__())
                # print(cf.label)
                for i in cf:
                    # print(i.label)
                    i.__setattr__('label',\
                        i.label.upper().replace('_',"").replace("-","").replace("OH","H").replace("(","").replace(")","").replace("WAT","OW").replace(",","").replace(".",""))
                    i.__setattr__('element',i.element.upper())
                    if i.label[0:4] == 'ZNII':
                        i.__setattr__('label',i.label.replace('ZNII','ZN'))
                    elif i.label[0:3] in atomTypeSub.keys():
                        i.__setattr__('label',i.label.replace(i.label[0:3],atomTypeSub[i.label[0:3]]))

                #-------------------------------------------------------------------------#
                #generate supercell
                elements = cf.element

                countElements = Counter(elements)
                    # print(countElements)

                strCountElements = ""

                    #     #sort based on element, returns a list of sorted keys 
                # if 'Blyth' in fileName:
                    # print(fileName)
                    # print(cf) 
                    # print(countElements)
                sortCountElement_keys = sorted(countElements.keys(),reverse=True)
                metalsMineral =  [(el,countElements[el]) for el in set(countElements.keys()).intersection(metals)]
                countMetals = [c[1] for c in metalsMineral]
                size = max(max(countMetals),len(set(countElements.keys()).intersection(metals)))
                # print(countMetals)
                # print(line[:-1])
                print(fileName)
                print(metalsMineral)
                print(size)

            
            # atom_labels.update(cf.label)

# with open('/home/kenneth/proj/proXtal/amcsd/xCif/singleOccupancy/natural/containsTransitionMetals/finalDataset/atom_labels.txt','w') as f:
    # for line in atom_labels:
        # f.write(line+'\n')


