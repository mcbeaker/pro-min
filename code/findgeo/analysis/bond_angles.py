#!/usr/bin/env python
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import seaborn as sb
import numpy as np
import scipy

rdir = '/home/kenneth/proj/proMin/results'

#combine min and pro
dfPro = pd.read_csv(os.path.join(rdir,'pro_fg_dictionaryValues_ligOcc_angles.csv'),header=0,index_col=False)
dfPro['Type'] = 'Protein'

dfMin = pd.read_csv(os.path.join(rdir,'min_fg_dictionaryValues_ligOcc_angles.csv'),header=0,index_col=False)
dfMin['Type'] = 'Mineral'

dfPro_Min = dfPro.append(dfMin)
dfPro_Min.replace(' ', '', regex=True)

dfMinElem = set(dfPro_Min[dfPro_Min.Type == "Mineral"].Angle.unique())
dfProElem = set(dfPro_Min[dfPro_Min.Type == "Protein"].Angle.unique())
overlap = dfMinElem.intersection(dfProElem)
transElems = sorted(set([bond.split("_")[1] for bond in overlap]))
# print(transElems)
# print(overlap)
# exit()
elem = 0
# notBonds = ['CU_CL','FE_C','FE_CL','MN_CL','NI_SE','MO_S','MN_C','V_N']
notAngles = ['CL_CU_O','CL_FE_O', 'CL_MN_CL','O_FE_CL','S_W_S',\
    'C_MN_O','C_W_C','O_MN_C','O_FE_C','S_FE_O','O_FE_S','S_V_S','O_MN_S','S_MN_O','O_CU_S']
fig, ax = plt.subplots(figsize=(30, 10),num=None, dpi=300, facecolor='w', edgecolor='k')
# fig.subplots_adjust(hspace=0, wspace=0)
all_angles = ['O_CO_O','O_CU_O','S_CU_S','O_FE_O','S_FE_S','O_MN_O','O_MO_O','S_MO_S','O_NI_O','S_NI_S','O_W_O','O_V_O']
O_angles = ['O_CO_O','O_CU_O','O_FE_O','O_MN_O','O_MO_O','O_NI_O','O_W_O','O_V_O']
df = dfPro_Min[dfPro_Min['Angle'].isin(O_angles)] #& (dfPro_Min.Degrees.between(dfPro_Min.Degrees.quantile(0.05), dfPro_Min.Degrees.quantile(.98)))]

g = sb.violinplot(x='Angle', y="Degrees",hue='Type',bw=0.2,data=df,split=True,order=O_angles)
g.set_xticklabels(g.get_xticklabels(), rotation=30)
g.set_xlabel(xlabel='Angle',fontsize=20)
g.set_ylabel(ylabel='Degrees',fontsize=20)
g.tick_params(labelsize=20)
# g.set(ylim=(,160))
outBL = os.path.join(rdir,'O_angles_pro_min_split_bondAngles.png')
plt.tight_layout()
sb.set_style('whitegrid')
plt.savefig(outBL)# 

for i in transElems:
    elemComb = [angle for angle in overlap if (angle.split("_")[1] == transElems[elem]) & (angle not in notAngles)]
    # print(elemComb)
    df = dfPro_Min[dfPro_Min['Angle'].isin(elemComb)].sort_values('Angle',ascending=False) #& (dfPro_Min.Degrees.between(dfPro_Min.Degrees.quantile(0.05), dfPro_Min.Degrees.quantile(.98)))]
    
    # print(df)
    # exit()
    
    # print(elemComb)
    # plt.text(0.5, 0.5, str((4, 4, i)), fontsize=18, ha='center')
    # print(df)
    
    for eComb in elemComb:
            minVal = df[(df.Type == 'Mineral') & (df.Angle == eComb)]['Degrees']#.dropna()
            print(eComb + ' minVal descriptors\n')
            print(minVal.describe())
            print('minVal Median ' + str(np.median(minVal)) + '\n')
            proVal = df[(df.Type == 'Protein') & (df.Angle == eComb)]['Degrees']#.dropna()
            print(eComb + ' proVal descriptors\n')
            print(proVal.describe())
            print('proVal Median ' + str(np.median(proVal)) + '\n')
            # print(proVal)
            print(eComb + ' val t-test\n')
            print(scipy.stats.ttest_ind(minVal,proVal))
            print(eComb + ' ks test\n')
            print(scipy.stats.ks_2samp(minVal,proVal))
    
    elem += 1