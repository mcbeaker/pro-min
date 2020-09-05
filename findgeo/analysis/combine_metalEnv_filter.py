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


#1a00.FE_142_4387_A.spy FE1N5 0.593 2.1477735644671845
# pdir = '/home/kenneth/proj/proMin/proteins/hagai/findgeo/total_findgeo'
#all results go here
rdir = '/home/kenneth/proj/proMin/results'

mdir = '/home/kenneth/proj/proMin/minerals/database/amcsd_pdb/pdb_reSeq_res_atom/combineFindGeoResults/analysis'
mdir = '/home/kenneth/proj/proMin/minerals/database/amcsd_pdb/pdb_reSeq_res_atom/combineFindGeoResults/analysis'
minFile='pro_fg_dictionaryValues.csv'
minFile = os.path.join(mdir,minFile)
columns = ['metalID','envComp','gRMSD','valence']
dfMin = pd.read_csv(minFile,sep=' ',header=None,index_col=False,names=columns)
# mineral format
# Naquite_0008632.Fe_4_4_A.irr FE1SI3 NA 0
dfMin['type'] = "Mineral"
splitMinColumns = ['mineral','amcsdID','atomName','resNum','atomNum','chain','geo']
dfMin[splitMinColumns] = dfMin.metalID.str.split('\.|_',expand=True)
dfMin.drop(columns='metalID',inplace=True)
dfMin['valence'] = dfMin.valence.astype(float)
dfMin['gRMSD'] = dfMin.gRMSD.astype(float)
# print('val is null\n') #1841
# print(len(dfMin[dfMin.valence.isnull()]))
# print(len(dfMin.valence.isnull()))

# print('geo is irr')
# print(len(dfMin[dfMin.geo  == 'irr'])) #788
# exit()
def get_description(df,name):
    print(name + " description\n")
    print(df.describe())
    print(name + " val == 0\n")
    print(len(df[df.valence.isnull()])) 
    print(name + " gRMSD == 'irr'\n")
    print(len(df[df.geo == 'irr']))
# print("Describe initial metal sites\n")
# print(dfMin.describe())
# count  1053.000000  1841.000000
# mean      0.302689     0.838647
# std       0.217892     0.855241
# min       0.000000     0.000000
# 25%       0.137000     0.337388
# 50%       0.267000     0.668119
# 75%       0.442000     1.013470
# max       0.966000     6.052019
# exit()

#remove rows that have both val==0 and have irr geo
# dfMin = dfMin[(dfMin.valence != 0) & (dfMin.geo != 'irr')]
dfMin['atomName'] = dfMin.atomName.str.replace('\d+|\+','').str.upper()
dfMin['pdb'] = dfMin[['mineral','amcsdID']].agg('-'.join, axis=1)
# get_description(dfMin,"filter geo != 'irr' & val != null" )
# gRMSD     valence
# count  972.000000  972.000000
# mean     0.302587    1.144416
# std      0.209076    0.953424
# min      0.000000    0.043816
# 25%      0.145750    0.520008
# 50%      0.267000    0.877834
# 75%      0.435750    1.365295
# max      0.966000    5.926891
# exit()
# dfMin.sort_values(['pdb','atomName','resNum','atomNum','chain'],axis=0,inplace=True)

# exit()
#--------------------------------------
# pdir='/home/kenneth/proj/proMin/proteins/hagai/pdbs/combineFindGeoResults/analysis'
pdir='/home/kenneth/proj/proMin/proteins/rcsb/pdbs_0_1.5/findgeo/combineFindGeoResults/analysis'
proFile='pro_fg_dictionaryValues.csv'
proFile = os.path.join(pdir,proFile)
columns = ['metalID','envComp','gRMSD','valence']
dfPro = pd.read_csv(proFile,sep=' ',header=None,index_col=False,names=columns)
dfPro['type'] = 'Protein'
#protein format
splitProColumns = ['pdb','atomName','resNum','atomNum','chain','geo']
dfPro[splitProColumns] = dfPro.metalID.str.split('\.|_',expand=True)
dfPro.drop(columns='metalID',inplace=True)
dfPro['valence'] = dfPro.valence.astype(float)
dfPro['gRMSD'] = dfPro.gRMSD.astype(float)
# print(dfPro.columns)
# print(dfPro[dfPro.pdb == '1gn8'])
# dfPro = dfPro[(dfPro.valence != 0) & (dfPro.geo != 'irr')]
# print(dfPro[dfPro.pdb == '1gn8'])
# exit()
dfPro['atomName'] = dfPro.atomName.str.replace('\d+|\+','').str.upper()
# dfPro['pdb'] = dfPro[['mineral','amcsdID']].agg('-'.join, axis=1) 
# dfPro.sort_values(['pdb','atomName','resNum','atomNum','chain'],axis=0,inplace=True)
dfPro['mineral'] = np.nan
dfPro['amcsdID'] = np.nan
print('Starting pro metal sites\n')
# print(len(dfPro.index)) 
# 31431

#combine min and pro
dfAll = pd.concat([dfMin,dfPro],ignore_index=True)
dfAll.sort_values(['pdb','chain','resNum','atomName','atomNum'],axis=0,inplace=True)
dfAll.to_csv(os.path.join(rdir,'min_pro_0_1.5_dictionaryValues.csv'),header=True,index=False,na_rep=np.nan,float_format='%.2f')

# exit()
# print("Combined min/pro sites\n")

# print(len(dfAll.index))
# 32403

order = ["CO","CU","FE","MN","MO","NI","V","W"]

#distributions of mineral and protein bulk valence
# bulkVal = dfAll.groupby(by=['pdb','atomName','resNum','envComp'],sort=False,as_index=False).first()
# bulkValData= bulkVal[bulkVal.valence.between(bulkVal.valence.quantile(.05), bulkVal.valence.quantile(.98))]
# axBulkVal = plt.axes(label="Bulk Valence")
# fig_bulkVal = sb.violinplot(x="type",y="valence",data=bulkValData,ax=axBulkVal).get_figure()
# outVal = os.path.join(rdir,'min_pro_val_percentile.png')
# fig_bulkVal.savefig(outVal)

#descriptive statistics
# minVal = bulkValData.where(bulkValData.type == 'Mineral')['valence'].dropna()
# print('minVal descriptors\n')
# print(minVal.describe())
# proVal = bulkValData.where(bulkValData.type == 'Protein')['valence'].dropna()
# print('proVal descriptors\n')
# print(proVal.describe())

# print('val t-test\n')
# print(scipy.stats.ttest_ind(minVal,proVal))
# exit()

#distribution of bulk mineral and protein gRMSD
# bulkRMSD = dfAll.groupby(by=['pdb','atomName','resNum','envComp'],sort=False,as_index=False).first()
# bulkRMSDData= bulkRMSD[bulkRMSD.gRMSD.between(bulkRMSD.gRMSD.quantile(.05), bulkRMSD.gRMSD.quantile(.98))]
# axBulkRMSD = plt.axes(label="Bulk gRMSD")
# fig_bulkRMSD = sb.violinplot(x="type",y="gRMSD",data=bulkRMSDData,ax=axBulkRMSD).get_figure()
# outRMSD = os.path.join(rdir,'min_pro_gRMSD_percentile.png')
# fig_bulkRMSD.savefig(outRMSD)

#gRMSD data for t-test
# minRMSD = bulkRMSDData.where(bulkRMSDData.type == 'Mineral')['gRMSD'].dropna()
# print('minRMSD descriptors\n')
# print(minRMSD.describe())
# proRMSD = bulkRMSDData.where(bulkRMSDData.type == 'Protein')['gRMSD'].dropna()
# print('proRMSD descriptors\n')
# print(proRMSD.describe())

# print('gRMSD t-test\n')
# print(scipy.stats.ttest_ind(minRMSD,proRMSD))
# exit()

gLAll = dfAll.groupby(by=['type','pdb','atomName','resNum','envComp'],sort=False,as_index=False).first()
print('Starting group by type/pdb/atomName/resNum/envComp metal sites\n')
print(len(gLAll.index)) 
# 16321
exit()
# gLAll = grAll.agg({'chain':lambda x: ','.join(set(x)),'valence':lambda x:np.round(np.average(x),decimals=2),'gRMSD':lambda x:np.round(np.average(x),decimals=2),'geo':lambda x: ",".join(set(x)) })#,'envComp':lambda x: ','.join(set(x)),'geo':lambda x: ",".join(set(x)) })
# with pd.option_context('display.max_rows', None, 'display.max_columns', 10):  # more options can be specified also
    # print(gLAll)

# axValType = plt.axes(label='Valence Type')
valData= gLAll[gLAll.valence.between(gLAll.valence.quantile(.05), gLAll.valence.quantile(.98))]
# print(valData.groupby('type').valence.describe())
# fig_valType = sb.violinplot(x="atomName",y="valence",order=order,hue="type",data=valData,split=True,ax=axValType).get_figure()
# outValType = os.path.join(rdir,'byType_min_pro_val_percentile.png')
# fig_valType.savefig(outValType)

#-------------------------------
# axRMSDType = plt.axes(label='gRMSD Type')
rmsdData= gLAll[gLAll.gRMSD.between(gLAll.gRMSD.quantile(0.1), gLAll.gRMSD.quantile(.98))]
# print(rmsdData.groupby('type').gRMSD.describe())
# fig_RMSDType = sb.violinplot(x="atomName",y="gRMSD",order=order,hue="type",data=rmsdData,split=True,ax=axRMSDType).get_figure()
# outRMSDType = os.path.join(rdir,'byType_min_pro_gRMSD_percentile.png')
# fig_RMSDType.savefig(outRMSDType)
# exit()

#valency data for t-test
minVal = valData[(valData.type == 'Mineral') & (valData.valence != np.nan)].groupby(by='atomName')['valence']#.dropna()#.groupby(by='atomName')
print('minVal descriptors\n')
print(minVal.describe())

proVal = valData[(valData.type == 'Protein') & (valData.valence != np.nan)].groupby(by='atomName')['valence']#.dropna()#.groupby(by='atomName')
print('proVal descriptors\n')
print(proVal.describe())

for atom in order:

    minVal = valData.where((valData.type == 'Mineral') & (valData.atomName == atom))['valence'].dropna()
    proVal = valData.where((valData.type == 'Protein') & (valData.atomName == atom))['valence'].dropna()
    print(atom+' val t-test\n')
    print(scipy.stats.ttest_ind(minVal,proVal))

    minRMSD = rmsdData.where((rmsdData.type == 'Mineral') & (rmsdData.atomName == atom))['gRMSD'].dropna()
    proRMSD = rmsdData.where((rmsdData.type == 'Protein') & (rmsdData.atomName == atom))['gRMSD'].dropna()
    print(atom+' gRMSD t-test\n')
    print(scipy.stats.ttest_ind(minRMSD,proRMSD))
exit()
# proVal = rmsdData.where(rmsdData.type == 'Protein')['valence'].dropna().groupby(by='atomName')
# print('proVal descriptors\n')
# print(proVal.describe())

# print('val t-test\n')
# print(scipy.stats.ttest_ind(minVal,proVal))

# #gRMSD data for t-test
# minRMSD = rmsdData.where(rmsdData.type == 'Mineral')['gRMSD'].dropna().groupby(by='atomName')
# print('minRMSD descriptors\n')
# print(minRMSD.describe())
# proRMSD = rmsdData.where(rmsdData.type == 'Protein')['gRMSD'].dropna().groupby(by='atomName')
# print('proRMSD descriptors\n')
# print(proRMSD.describe())

# print('gRMSD t-test\n')
# print(scipy.stats.ttest_ind(minRMSD,proRMSD))
exit()
# ax2.set(ylim=(-0.2,1.1))
# ax2 = sb.violinplot(x="atomName",y="gRMSD",order=order,hue="type",data=rmsdData,split=True,ax=ax2).get_figure()

# rmsd = ax2.get_figure()
# outRMSD = os.path.join(pdir,'min_pro_gRMSD_percentile.png')
# ax2.savefig(outRMSD)


#-----------------------

# val = sb.violinplot(x=gLVal.atomName,y=gLVal.valence).get_figure()
# val = sb.violinplot(x=gLMin.atomName,y=gLMin.valence).get_figure()
# # outVal = os.path.join(pdir,'min_val_percentile_05_98.png')

# # rmsd = sb.violinplot(x=gLgRMSD.atomName,y=gLgRMSD.gRMSD).get_figure()
# rmsd = sb.violinplot(x=gLMin.atomName,y=gLMin.gRMSD).get_figure()
# # outRMSD = os.path.join(pdir,'min_RMSD_percentile_05_98')
# outRMSD = os.path.join(pdir,'min_RMSD')
# rmsd.savefig(outRMSD)







# print(gL)
# gLVal = gL[gL.valence.between(gL.valence.quantile(.05), gL.valence.quantile(.98))]
# gLVal = gLVal.groupby(by='atomName',as_index=False)
# print(gLVal.groups)
# exit()

# gLgRMSD = gL[gL.gRMSD.between(gL.gRMSD.quantile(.05), gL.gRMSD.quantile(.98))]
# gLgRMSD = gLgRMSD.groupby(by='atomName',as_index=False)


# print(list(gL.atomName.drop_duplicates()))


# plt.savefig('/home/kenneth/testgRMSD_05_98.png')
# gLgRMSD.boxplot(by='atomName',column='gRMSD',return_type=None,figsize=(20,20))
# gLVal.boxplot(by='atomName',column='valence',return_type=None,figsize=(20,20))
# gL.plot.box(by=gL.atomName,columns=gL.valence,return_type=None)
# # plt.show()
# # print(gL[gL.atomName == 'W10,W1,W11,W12,W13,W14,W15,W16,W17,W18,W2,W3,W4,W5,W6,W7,W8,W9'])

# # atomName = grouped_lists.groupby(by='atomName',as_index=False,sort=False)
# # print(atomName.head())
# # atomName.groups['FE1,FE2,FE3,FE4,FE5,FE6,FE7,FE8']
# # print(atomName.get_group('FE1,FE2,FE3,FE4,FE5,FE6,FE7,FE8').describe())
# # print(atomName.get_group('FE').describe())



# for i in grouped_lists:
#     print(i)


# unq = df_filt.drop_duplicates(subset=splitColumns)

# print(unq)

