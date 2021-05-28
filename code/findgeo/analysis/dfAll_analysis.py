#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib
matplotlib.use('Agg')
import seaborn as sb
import numpy as np
import scipy

# rdir = '/home/kenneth/proj/proMin/results'
rdir = '/Volumes/proMin/results'

#combine min and pro
# dfAll = pd.read_csv(os.path.join(rdir,'ec_min_pro_fg_dictionaryValues.csv'),header=0,index_col=False)
# dfAll = pd.read_csv(os.path.join(rdir,'nr_40','nr_40_pro_min.csv'),header=0,index_col=0)
dfAll = pd.read_csv(os.path.join(rdir,'pdbs_0_1.5','min_pro_0_1.5_dictionaryValues.csv'),header=0,index_col=0)

exit()

# mdir = '/home/kenneth/proj/proMin/minerals/database/data'
# dfMin = pd.read_csv(os.path.join(mdir,'RRUFF_AMCSD_Fe_Mineral_IDs_3-shaunna-retrieved.csv'),header=0,index_col=False)

#merge 
# dfAgeMinPro = pd.merge(dfAll,dfMin,how='outer',left_on='mineral',right_on='Mineral_Name_plain').Oldest_Known_Age_Ma.dropna()
# print(dfAgeMinPro)
# print(dfAgeMinPro[~dfAgeMinPro.mineral.isnull()].Oldest_Known_Age_Ma.isnull())

# exit()

# print(dfAll.Resolution)
#plot histogram of resolution
'''resData = dfAll[(dfAll.type == 'Protein')  & (dfAll.Resolution != 'None')]
resData = resData.astype({'Resolution':float})
#distributions of mineral and protein bulk valence
resGroup = resData.groupby(by=['pdb','atomName','resNum','envComp'],sort=False,as_index=False).first()

axRes = plt.axes(label="Res")
axRes.tick_params(axis='both', which='major', labelsize=18)
# plt.rcParams['ytick.labelsize'] = 'large'
sb.violinplot(x="atomName",y="Resolution",data=resData,ax=axRes)
axRes.set_ylabel('Resolution'+ u' \u212B',fontsize=18)
axRes.set_xlabel('Metal',fontsize=18)
outRes = os.path.join(rdir,'VP_pdb_resolution.png')
plt.savefig(outRes)'''
#resolution
# exit()

# print(dfAll.head().gRMSD)
dfAll.gRMSD = dfAll.gRMSD*(180/np.pi)
# print(dfAll.head().gRMSD)

order = ["CO","CU","FE","MN","MO","NI","V","W"]

#distributions of mineral and protein bulk valence
bulkVal = dfAll.groupby(by=['pdb','atomName','resNum','envComp'],sort=False,as_index=False).first()
# print(bulkVal.sort_values(['valence'],axis=0,inplace=False,ascending=False))
bulkValData= bulkVal[bulkVal.valence.between(bulkVal.valence.quantile(.05), bulkVal.valence.quantile(.98))]
# bulkValData= bulkVal[bulkVal.valence.between(bulkVal.valence.quantile(0), bulkVal.valence.quantile(1))]
bulkValData=bulkValData[((bulkValData.type == 'Mineral') | (bulkValData.type == 'Protein')) & (bulkValData.valence != 0)]
bulkValData['val_ENV'] = 'Metal Coordination Environment'
# print(bulkValData.sort_values(['valence'],axis=0,inplace=False,ascending=False))

# bulkValData=bulkValData[(bulkValData.type == 'Protein') & (bulkValData.valence != 0)]
# print(bulkValData)
axBulkVal = plt.axes(label="Bulk Valence")
axBulkVal.tick_params(axis='both', which='major', labelsize=18)
# plt.rcParams['ytick.labelsize'] = 'large'
sb.violinplot(x="val_ENV",y="valence",hue='type',split=True,data=bulkValData,ax=axBulkVal,legend=['Protein','Mineral'])
axBulkVal.set_ylabel('Valence',fontsize=18)
axBulkVal.set_xlabel('',fontsize=18)

# axBulkVal.yaxis.get_label().set_fontsize(18)

# axBulkVal.set_xlabel('Valence',fontsize=24)

# sb.stripplot(x="type",y="valence",data=bulkValData,ax=axBulkVal,zorder=1,color='grey')
# print('test')
# fig_bulkVal = sb.swarmplot(x="type",y="valence",data=bulkValData,ax=axBulkVal).get_figure()
outVal = os.path.join(rdir,'pdbs_0_1.5','bulkVal_VP_per_min_pro.png')
plt.savefig(outVal)
# exit()

#descriptive statistics
minVal = bulkValData[(bulkValData.type == 'Mineral')]['valence']#.dropna()
print('minVal descriptors\n')
print(minVal.describe())
proVal = bulkValData[(bulkValData.type == 'Protein')]['valence']#.dropna()
print('proVal descriptors\n')
print(proVal.describe())
# print(proVal)
print('val t-test\n')
print(scipy.stats.ttest_ind(minVal,proVal))
# exit()

#distributions of mineral and protein bulk valence
bulkRMSD = dfAll.groupby(by=['pdb','atomName','resNum','envComp'],sort=False,as_index=False).first()
bulkRMSDData= bulkRMSD[bulkRMSD.gRMSD.between(bulkRMSD.gRMSD.quantile(.05), bulkRMSD.gRMSD.quantile(.98))]
# bulkRMSDData= bulkRMSD[bulkRMSD.gRMSD.between(bulkRMSD.gRMSD.quantile(0), bulkRMSD.gRMSD.quantile(1))]
bulkRMSDData=bulkRMSDData[((bulkRMSDData.type == 'Mineral') | (bulkRMSDData.type == 'Protein')) & ~(bulkRMSDData.gRMSD.isnull())]
bulkRMSDData['gRMSD_ENV'] = 'Metal Coordination Environment'
# bulkRMSDData=bulkRMSDData[(bulkRMSDData.type == 'Protein') & (bulkRMSDData.valence != 0)]
# print(bulkRMSDData)
axBulkRMSD = plt.axes(label="Bulk RMSDence")
axBulkRMSD.tick_params(axis='both', which='major', labelsize=18)
axBulkRMSD.legend(('Protein','Mineral'))
sb.violinplot(x='gRMSD_ENV',y="gRMSD",hue='type',split=True,data=bulkRMSDData,ax=axBulkRMSD)
axBulkRMSD.set_ylabel('gRMSD' +u' \N{DEGREE SIGN}',fontsize=18)
axBulkRMSD.set_xlabel(xlabel='')

# sb.stripplot(x="type",y="valence",data=bulkRMSDData,ax=axBulkRMSD,zorder=1,color='grey')
# print('test')
# fig_bulkRMSD = sb.swarmplot(x="type",y="valence",data=bulkRMSDData,ax=axBulkRMSD).get_figure()
outRMSD = os.path.join(rdir,'pdbs_0_1.5','bulkRMSD_VP_per_min_pro.png')

plt.savefig(outRMSD)
# exit()
#descriptive statistics
minRMSD = bulkRMSDData[(bulkRMSDData.type == 'Mineral')]['gRMSD']#.dropna()
print('minRMSD descriptors\n')
print(minRMSD.describe())
proRMSD = bulkRMSDData[(bulkRMSDData.type == 'Protein')]['gRMSD']#.dropna()
print('proRMSD descriptors\n')
print(proRMSD.describe())
print(proRMSD)
# print('val t-test\n')
print(scipy.stats.ttest_ind(minRMSD,proRMSD))

# exit()

#group statistics

gLAll = dfAll.groupby(by=['type','pdb','atomName','resNum','envComp'],sort=False,as_index=False).first()
# print('Starting group by type/pdb/atomName/resNum/envComp metal sites\n')
# print(gLAll.Cofactor.unique())
# exit() 





#------VAL-bytype-violin
axValType = plt.axes(label='Valence Type')
valData= gLAll[gLAll.valence.between(gLAll.valence.quantile(0.05), gLAll.valence.quantile(.98))]
# valData= gLAll[gLAll.valence.between(gLAll.valence.quantile(0), gLAll.valence.quantile(1))]
print(valData.groupby('atomName').valence.describe())#.count)

sb.violinplot(x="atomName",y="valence",order=order,hue="type",data=valData,split=True,ax=axValType,legend=['Protein','Mineral'])
# axValType.legend(('Protein','Mineral'))
axValType.set_ylabel('Valence (vu)',fontsize=18)
axValType.set_xlabel(xlabel='Metal')
outValType = os.path.join(rdir,'pdbs_0_1.5','val_VP_type_per_min_pro.png')
plt.savefig(outValType)
# exit()
#-------Val-bytype-violin

#-------Val-byCofactor-violin




#-------Val-byCofactor-violin


#-------Val-bytype-barplot
# axValBarType = plt.axes(label='Valence Bar Type')

labels = order
minBarData = valData[valData.type == 'Mineral'].groupby('atomName').count()
# print(minBarData)
proBarData = valData[valData.type == 'Protein'].groupby('atomName').count()
x = np.arange(len(labels[0:7]))  # the label locations - proteins do not have W
width = 0.35  # the width of the bars

fig, axBPVal = plt.subplots()
print(np.shape(proBarData[0:7]),np.shape(minBarData[0:7])) #proBarData[0:7] does not have W - [0:8] includes W
# exit()
rects1 = axBPVal.bar(x + width/2, proBarData.valence[0:7], width, label='Proteins')
rects2 = axBPVal.bar(x - width/2, minBarData.valence[0:7], width, label='Minerals')
print(rects1,rects2)
# exit()
# Add some text for labels, title and custom x-axis tick labels, etc.
axBPVal.set_ylabel('Count of Valency calculations per metal type')
axBPVal.set_xlabel("Metal")
# axBPVal.set_title('Number of Metal Coordination sites for Proteins and Minerals ')
axBPVal.set_xticks(x)
axBPVal.set_xticklabels(labels[0:7])
axBPVal.legend(('Protein','Mineral'))


def autolabel(rects,ax):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
autolabel(rects1,axBPVal)
autolabel(rects2,axBPVal)

fig.tight_layout()
outValBarType = os.path.join(rdir,'pdbs_0_1.5','val_BP_type_per_min_pro.png')
plt.savefig(outValBarType)
#-------Val-bytype-barplot

#-------RMSD-byType-violinPlot------------------------
axRMSDType = plt.axes(label='gRMSD Type')
rmsdData= gLAll[gLAll.gRMSD.between(gLAll.gRMSD.quantile(0.05), gLAll.gRMSD.quantile(.98))]
# rmsdData= gLAll[gLAll.gRMSD.between(gLAll.gRMSD.quantile(0), gLAll.gRMSD.quantile(1))]
# print(rmsdData.groupby('atomName').gRMSD.describe())
fig_RMSDType = sb.violinplot(x="atomName",y="gRMSD",order=order[0:7],hue="type",data=rmsdData,split=True,ax=axRMSDType,legend=['Protein','Mineral'])
# axRMSDType.legend(('Protein','Mineral'))
axRMSDType.set_ylabel('gRMSD' +u' \N{DEGREE SIGN}',fontsize=18)
axRMSDType.set_xlabel(xlabel='Metal')
# outRMSDType = os.path.join(rdir,'nr40_gRMSD_VP_type_per_min_pro.png')
outRMSDType = os.path.join(rdir,'pdbs_0_1.5','gRMSD_VP_type_per_min_pro.png')
plt.savefig(outRMSDType)
#-------RMSD-byType-violinPlot------------------------

#-------RMSD-bytype-barplot
# axRMSDBarType = plt.axes(label='RMSD Bar Type')
minBarData = rmsdData[rmsdData.type == 'Mineral'].groupby('atomName').count()
# print(minBarData:)
proBarData = rmsdData[rmsdData.type == 'Protein'].groupby('atomName').count()
x = np.arange(len(labels[0:7]))  # the label locations
width = 0.35  # the width of the bars

figRMSD, axRMSD = plt.subplots()
rects1 = axRMSD.bar(x + width/2, proBarData.gRMSD[0:7], width, label='Proteins')
rects2 = axRMSD.bar(x - width/2, minBarData.gRMSD[0:7], width, label='Minerals')


# Add some text for labels, title and custom x-axis tick labels, etc.
axRMSD.set_ylabel('Count of gRMSD calculations per metal type')
axRMSD.set_xlabel("Metal")
axRMSD.xaxis.set_major_locator(ticker.FixedLocator(x))
axRMSD.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
# axRMSD.set_title('Number of Val:
# axRMSD.set_xticklabels(labels)
axRMSD.legend(('Protein','Mineral'))

autolabel(rects1,axRMSD)
autolabel(rects2,axRMSD)

figRMSD.tight_layout()
# outRMSDBarType = os.path.join(rdir,'nr40_gRMSD_BP_type_per_min_pro.png')
outRMSDBarType = os.path.join(rdir,'pdbs_0_1.5','gRMSD_BP_type_per_min_pro.png')
plt.savefig(outRMSDBarType)
#-------RMSD-bytype-barplot


for atom in order:

    minRMSD = rmsdData.where((rmsdData.type == 'Mineral') & (rmsdData.atomName == atom))['gRMSD'].dropna()
    proRMSD = rmsdData.where((rmsdData.type == 'Protein') & (rmsdData.atomName == atom))['gRMSD'].dropna()
    print(atom+' gRMSD t-test')
    print(scipy.stats.ttest_ind(minRMSD,proRMSD))
    # print("\n")

for atom in order:

   minVal = valData.where((valData.type == 'Mineral') & (valData.atomName == atom))['valence'].dropna()
   proVal = valData.where((valData.type == 'Protein') & (valData.atomName == atom))['valence'].dropna()
   print(atom+' val t-test')
   print(scipy.stats.ttest_ind(minVal,proVal))
   print("\n")
# exit()
# ec number analysis
