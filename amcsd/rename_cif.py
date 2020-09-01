from CifFile import ReadCif
import pandas as pd

import os

finData = {}
    
with open('/home/kenneth/proj/proXtal/amcsd/xCif/min_amcsd.txt') as f:
    next(f)
    for line in f:
        amcsd = str(line.split(',')[1].rstrip())
        mineral = line.split(',')[0]
        # print(amcsd)
        finData.update({amcsd:mineral})
# print(finData.keys())

# print(finData['amcsd'])

# print(finData['0002004'])

# exit()
test = '/home/kenneth/proj/proXtal/amcsd/xCif/singleOccupancy/natural/containsTransitionMetals/finalDataset/test'

d = '/home/kenneth/proj/proXtal/amcsd/database/'
noGo = ['xx14764.cif','xx12015.cif','xx14771.cif','xx14774.cif','xx14760.cif',\
'xx18483.cif','xx15076.cif','xx14768.cif','xx14762.cif','xx15880.cif','xx15075.cif','xx7038.cif','xx18398.cif',\
    'xx14759.cif','xx14770.cif','xx18396.cif','xx14761.cif','xx14765.cif','xx18487.cif','xx7290.cif','xx18484.cif',\
        'xx18485.cif','xx14767.cif','xx14763.cif','xx14766.cif','xx18486.cif','xx20064.cif','xx14769.cif','xx7289.cif']


amcsdCodes = []
from progress.bar import Bar

with Bar('Processing',max=len(os.listdir(d))) as bar:
    for f in os.listdir(d):

        if 'xx' in os.path.splitext(f)[0] and f not in noGo: 
            try:
                # print(os.path.join(d,f))
                fileName = os.path.join(d,f)
                cf = ReadCif(fileName)              

                # print(cf['global']['_chemical_name_mineral'])ds
                if '_database_code_amcsd' in cf['global']:
                    # print(cf['global']['_chemical_name_mineral'])
                    code = str(cf['global']['_database_code_amcsd'])

                    if code in finData.keys():
                        mineral = finData[code]
                    
                        outName = os.path.join(d,'formatName',mineral+"_"+code+".cif")

                        with open(outName,"w") as outfile:
                            outfile.write(cf.WriteOut())

                    # print(code + '\n')
                    # print(repr(code)) 
                    # amcsdCodes.append(code)
                    # if code in finData.keys():
                    #     print(code + ' yes')
                    #     exit()
                    # if code not in finData.keys():
                    #     print(code +' no')
                    # print(code)
                    
                    # print(type(finData['amcsd'][1]))
                    # print('blah')
                    # print(code)
                    # print(finData.loc[finData['amcsd'] == 11781 ])
                    # mineral =  finData[code]#['min']                                                                                                                                            
                    # print(mineral['min'])
                    # exit()
                    # outName = os.path.join(test,mineral+'.cif')
                    # # print(outName)
                    # with open(outName,"w") as outfile:
                    #     outfile.write(cf.WriteOut())
            except Exception as e:
            #     # print(f)
                # print(e)
                continue
        bar.next()
bar.finish()

# with open('/home/kenneth/proj/proXtal/amcsd/xCif/amcsdCode.txt','w') as f:
#     for i in amcsdCodes:

#         f.write(i + '\n')

# # 14579 6001 0
# # print("%s %s %s" % (countMineralName,countAMCSDName,countNone))
# # print(amcsdCodes)
# # print(set(finData.keys()).intersection(set(amcsdCodes)))
# # # noGoLen = len(noGo)
# # # # print(noGoLen)
# # # countMineralName = 0
# # countAMCSDName = 0
# # countNone = 0
# # countExists = 0
