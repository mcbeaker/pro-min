#!/usr/bin/env python
# %%
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import seaborn as sb
import numpy as np
import scipy
# %%
rdir = '/home/kenneth/proj/proMin/results'

#combine min and pro
cols = ['pdbFileName','metElm','ligElm','dist','metName','ligName', \
        'metOcc' ,'ligOcc' ,'metResName' ,'ligResName','metResNum' ,'ligResNum', 'metChain','ligChain'] 
dfPro = pd.read_csv(os.path.join(rdir,'pro_fg_dictionaryValues_ligOcc_distances.csv'),header=None,index_col=False,names=cols)
dfPro['Type'] = 'Protein'

dfMin = pd.read_csv(os.path.join(rdir,'min_fg_dictionaryValues_ligOcc_distances.csv'),header=None,index_col=False,names=cols)
dfMin['Type'] = 'Mineral'

dfPro_Min = dfPro.append(dfMin)
dfPro_Min['elemComb'] = dfPro_Min.metElm + "_" + dfPro_Min.ligElm  

# %% 
# print(dfPro_Min.elemComb.unique())
# ['FE_S', 'FE_N', 'FE_O', 'CU_N', 'CU_S', 'CU_O', 'NI_O', 'FE_C',
#        'NI_N', 'CO_O', 'CO_N', 'NI_C', 'MN_N', 'MN_O', 'MO_O', 'W_O',
#        'NI_S', 'V_O', 'MN_CL', 'CU_C', 'CO_C', 'MN_C', 'NI_ZN', 'NI_SE',
#        'CO_CL', 'CU_ZN', 'NI_P', 'MN_P', 'V_S', 'CO_S', 'MO_N', 'MO_S',
#        'MN_S', 'MN_F', 'FE_CL', 'CU_CL', 'MN_MG', nan, 'V_N', 'FE_SI',
#        'W_S', 'CU_AU', 'FE_SB', 'NI_SN', 'NI_PB', 'MN_BE', 'CU_PD',
#        'FE_PT', 'NI_TE', 'NI_SB', 'NI_AS', 'W_C', 'CU_AS', 'CO_AS',
#        'CU_SB', 'CO_SB', 'CU_PT', 'MN_SI', 'FE_SE', 'CU_SE', 'FE_AS',
#        'CU_SN', 'CO_SE', 'FE_TE', 'CO_TE', 'CU_I', 'CU_U', 'FE_AL'
dfMinElem = set(dfPro_Min[dfPro_Min.Type == "Mineral"].elemComb.unique())
dfProElem = set(dfPro_Min[dfPro_Min.Type == "Protein"].elemComb.unique())
overlap = dfMinElem.intersection(dfProElem)
transElems = sorted(set([bond.split("_")[0] for bond in overlap]))
# {'CO_O','CO_S','CU_CL','CU_O','CU_S','FE_C','FE_CL','FE_O','FE_S','MN_C','MN_CL',
#  'MN_O','MN_S','MO_O','MO_S','NI_O','NI_S','NI_SE','V_N','V_O','V_S','W_O'}
# %%
fig,axs = plt.subplots(2,4, sharex=False,sharey= False)
fig.subplots_adjust(hspace=0.4, wspace=0.4)
elem = 0
notBonds = ['CU_CL','FE_C','FE_CL','MN_CL','NI_SE','MO_S','MN_C','V_N']
for row in range(0,2):
    for col in range(0,4):
        elemComb = [bond for bond in overlap if (bond.split("_")[0] == transElems[elem]) & (bond not in notBonds)]
        df = dfPro_Min[dfPro_Min['elemComb'].isin(elemComb)]
        
        print(elemComb)
    #     plt.text(0.5, 0.5, str((4, 4, i)), fontsize=18, ha='center')
        df = df[(df.metElm.str.find(transElems[elem]) != -1) & (df.dist.between(df.dist.quantile(0.05), df.dist.quantile(.98)))]
        g = sb.catplot(x="elemComb", y="dist",hue='Type',scale_hue=False,scale='count',bw=0.2,data=df,split=True,kind='violin',ax = axs[row,col])
        for eComb in elemComb:
                minVal = df[(df.Type == 'Mineral') & (df.elemComb == eComb)]['dist']#.dropna()
                print(eComb + ' minVal descriptors\n')
                print(minVal.describe())
                print('minVal Median ' + str(np.median(minVal)) + '\n')
                proVal = df[(df.Type == 'Protein') & (df.elemComb == eComb)]['dist']#.dropna()
                print(eComb + ' proVal descriptors\n')
                print(proVal.describe())
                print('proVal Median ' + str(np.median(proVal)) + '\n')
                # print(proVal)
                print(eComb + ' val t-test\n')
                print(scipy.stats.ttest_ind(minVal,proVal))
                print(eComb + ' ks test\n')
                print(scipy.stats.ks_2samp(minVal,proVal))
        outBL = os.path.join(rdir,transElems[elem]+"_"+'split_VP_noScaleHue_count_0.2bw_pro_min_split_bondLengths.png')
        plt.savefig(outBL)# 
        elem += 1

# dfPro['elemComb'] = dfPro.metElm + "_" + dfPro.ligElm
# dfPro = dfPro[(dfPro.metElm == 'CU') & (dfPro.dist.between(dfPro.dist.quantile(0.05), dfPro.dist.quantile(.98)))] # & ((dfPro.elemComb == 'FE_S') | (dfPro.elemComb == 'FE_O')) \

# %%
dfPro['elemComb'] = dfPro.metElm + "_" + dfPro.ligElm
dfPro = dfPro[(dfPro.metElm == 'FE') & ((dfPro.elemComb == 'FE_S') | (dfPro.elemComb == 'FE_O')) \
        & (dfPro.dist.between(dfPro.dist.quantile(0.05), dfPro.dist.quantile(.98)))]
# dfPro = dfPro.groupby(by=['metElm'])
# g = sb.catplot(x="elemComb", y="dist",col="metElm", data=dfPro)

# axValType.legend(('Protein','Mineral'))
# outBL = os.path.join(rdir,'BP_FE_S_O_pro_bondLengths.png')
# plt.savefig(outBL)
# %%
dfMin['elemComb'] = dfMin.metElm + "_" + dfMin.ligElm
dfMin = dfMin[(dfMin.metElm == 'FE') & ((dfMin.elemComb == 'FE_S') | (dfMin.elemComb == 'FE_O'))\
        & (dfMin.dist.between(dfMin.dist.quantile(0.05), dfMin.dist.quantile(.98)))]
dfMin['Type'] = 'Mineral'
# dfPro = dfPro.groupby(by=['metElm'])
# g = sb.catplot(x="elemComb", y="dist",col="metElm", data=dfMin)

# axValType.legend(('Protein','Mineral'))
# outBL = os.path.join(rdir,'BP_FE_S_O_min_bondLengths.png')
# plt.savefig(outBL)
# %%
dfPro_Min = dfPro.append(dfMin)
g = sb.catplot(x="elemComb", y="dist",col="metElm",hue='Type',split=True,data=dfPro_Min,kind='violin')
outBL = os.path.join(rdir,'BP_FE_S_O_pro_min_split_bondLengths.png')
plt.savefig(outBL)

# %%
minVal = dfMin[(dfMin.Type == 'Mineral') & (dfMin.elemComb == 'FE_S')]['dist']#.dropna()
print('minVal descriptors\n')
print(minVal.describe())
proVal = dfPro[(dfPro.Type == 'Protein') & (dfPro.elemComb == 'FE_S')]['dist']#.dropna()
print('proVal descriptors\n')
print(proVal.describe())
# print(proVal)
print('val t-test\n')
print(scipy.stats.ttest_ind(minVal,proVal))
print('ks test\n')
print(scipy.stats.ks_2samp(minVal,proVal))
# %%
minVal = dfMin[(dfMin.Type == 'Mineral') & (dfMin.elemComb == 'FE_O')]['dist']#.dropna()
print('minVal descriptors\n')
print(minVal.describe())
proVal = dfPro[(dfPro.Type == 'Protein') & (dfPro.elemComb == 'FE_O')]['dist']#.dropna()
print('proVal descriptors\n')
print(proVal.describe())
# print(proVal)
print('val t-test\n')
print(scipy.stats.ttest_ind(minVal,proVal))
print('ks test\n')
print(scipy.stats.ks_2samp(minVal,proVal))# %%

# %%
