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
dfAll = pd.read_csv(os.path.join(rdir,'ec_min_pro_fg_dictionaryValues.csv'),header=0,index_col=False)
print(dfAll.head().gRMSD)
dfAll.gRMSD = dfAll.gRMSD*(180/np.pi)
print(dfAll.head().gRMSD)
exit()
order = ["CO","CU","FE","MN","MO","NI","V","W"]

#distributions of mineral and protein bulk valence
bulkVal = dfAll.groupby(by=['pdb','atomName','resNum','envComp'],sort=False,as_index=False).first()
# print(bulkVal.sort_values(['valence'],axis=0,inplace=False,ascending=False))
bulkValData= bulkVal[bulkVal.valence.between(bulkVal.valence.quantile(.05), bulkVal.valence.quantile(.98))]
bulkValData=bulkValData[((bulkValData.type == 'Mineral') | (bulkValData.type == 'Protein')) & (bulkValData.valence != 0)]
# print(bulkValData.sort_values(['valence'],axis=0,inplace=False,ascending=False))

# bulkValData=bulkValData[(bulkValData.type == 'Protein') & (bulkValData.valence != 0)]
# print(bulkValData)
axBulkVal = plt.axes(label="Bulk Valence")
sb.violinplot(x="type",y="valence",data=bulkValData,ax=axBulkVal,legend=['Protein','Mineral'])
# sb.stripplot(x="type",y="valence",data=bulkValData,ax=axBulkVal,zorder=1,color='grey')
# print('test')
# fig_bulkVal = sb.swarmplot(x="type",y="valence",data=bulkValData,ax=axBulkVal).get_figure()
outVal = os.path.join(rdir,'VP_min_pro_bulkVal_percentile.png')
plt.savefig(outVal)

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
bulkRMSDData=bulkRMSDData[((bulkRMSDData.type == 'Mineral') | (bulkRMSDData.type == 'Protein')) & ~(bulkRMSDData.gRMSD.isnull())]
# bulkRMSDData=bulkRMSDData[(bulkRMSDData.type == 'Protein') & (bulkRMSDData.valence != 0)]
# print(bulkRMSDData)
axBulkRMSD = plt.axes(label="Bulk RMSDence")
axBulkRMSD.legend(('Protein','Mineral'))
sb.violinplot(x="type",y="valence",data=bulkRMSDData,ax=axBulkRMSD)
# sb.stripplot(x="type",y="valence",data=bulkRMSDData,ax=axBulkRMSD,zorder=1,color='grey')
# print('test')
# fig_bulkRMSD = sb.swarmplot(x="type",y="valence",data=bulkRMSDData,ax=axBulkRMSD).get_figure()
outRMSD = os.path.join(rdir,'VP_min_pro_bulkRMSD_percentile.png')
plt.savefig(outRMSD)

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
print('Starting group by type/pdb/atomName/resNum/envComp metal sites\n')

#------VAL-bytype-violin
axValType = plt.axes(label='Valence Type')
valData= gLAll[gLAll.valence.between(gLAll.valence.quantile(0.05), gLAll.valence.quantile(.98))]
print(valData.groupby('atomName').valence.describe())#.count)

sb.violinplot(x="atomName",y="valence",order=order,hue="type",data=valData,split=True,ax=axValType,legend=['Protein','Mineral'])
# axValType.legend(('Protein','Mineral'))
outValType = os.path.join(rdir,'VP_byType_min_pro_val_percentile.png')
plt.savefig(outValType)
#-------Val-bytype-violin

#-------Val-bytype-barplot
# axValBarType = plt.axes(label='Valence Bar Type')

labels = order
minBarData = valData[valData.type == 'Mineral'].groupby('atomName').count()
# print(minBarData)
proBarData = valData[valData.type == 'Protein'].groupby('atomName').count()
x = np.arange(len(labels))  # the label locations
width = 0.35  # the width of the bars

fig, axBPVal = plt.subplots()
rects1 = axBPVal.bar(x + width/2, proBarData.valence, width, label='Proteins')
rects2 = axBPVal.bar(x - width/2, minBarData.valence, width, label='Minerals')


# Add some text for labels, title and custom x-axis tick labels, etc.
axBPVal.set_ylabel('Count of Valency calculations per metal type')
axBPVal.set_xlabel("Metal")
# axBPVal.set_title('Number of Metal Coordination sites for Proteins and Minerals ')
axBPVal.set_xticks(x)
axBPVal.set_xticklabels(labels)
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
outValBarType = os.path.join(rdir,'BP_byType_min_pro_val_percentile.png')
plt.savefig(outValBarType)
#-------Val-bytype-barplot

#-------RMSD-byType-violinPlot------------------------
axRMSDType = plt.axes(label='gRMSD Type')
rmsdData= gLAll[gLAll.gRMSD.between(gLAll.gRMSD.quantile(0.05), gLAll.gRMSD.quantile(.98))]
# print(rmsdData.groupby('atomName').gRMSD.describe())
fig_RMSDType = sb.violinplot(x="atomName",y="gRMSD",order=order,hue="type",data=rmsdData,split=True,ax=axRMSDType,legend=['Protein','Mineral'])
# axRMSDType.legend(('Protein','Mineral'))
outRMSDType = os.path.join(rdir,'VP_byType_min_pro_gRMSD_percentile.png')
plt.savefig(outRMSDType)
#-------RMSD-byType-violinPlot------------------------

#-------RMSD-bytype-barplot
# axRMSDBarType = plt.axes(label='RMSD Bar Type')
minBarData = rmsdData[rmsdData.type == 'Mineral'].groupby('atomName').count()
# print(minBarData:)
proBarData = rmsdData[rmsdData.type == 'Protein'].groupby('atomName').count()
x = np.arange(len(labels))  # the label locations
width = 0.35  # the width of the bars

figRMSD, axRMSD = plt.subplots()
rects1 = axRMSD.bar(x + width/2, proBarData.gRMSD, width, label='Proteins')
rects2 = axRMSD.bar(x - width/2, minBarData.gRMSD, width, label='Minerals')


# Add some text for labels, title and custom x-axis tick labels, etc.
axRMSD.set_ylabel('Count of gRMSD calculations per metal type')
axRMSD.set_xlabel("Metal")
# axRMSD.set_title('Number of Val:
axRMSD.set_xticklabels(labels)
axRMSD.legend(('Protein','Mineral'))

autolabel(rects1,axRMSD)
autolabel(rects2,axRMSD)

figRMSD.tight_layout()
outRMSDBarType = os.path.join(rdir,'BP_byType_min_pro_RMSD_percentile.png')
plt.savefig(outRMSDBarType)
#-------RMSD-bytype-barplot


for atom in order:

    minRMSD = rmsdData.where((rmsdData.type == 'Mineral') & (rmsdData.atomName == atom))['gRMSD'].dropna()
    proRMSD = rmsdData.where((rmsdData.type == 'Protein') & (rmsdData.atomName == atom))['gRMSD'].dropna()
    print(atom+' gRMSD t-test')
    print(scipy.stats.ttest_ind(minRMSD,proRMSD))
    print("\n")

for atom in order:

   minVal = valData.where((valData.type == 'Mineral') & (valData.atomName == atom))['valence'].dropna()
   proVal = valData.where((valData.type == 'Protein') & (valData.atomName == atom))['valence'].dropna()
   print(atom+' val t-test')
   print(scipy.stats.ttest_ind(minVal,proVal))
   print("\n")

# ec number analysis
axEC1Type = plt.axes(label='EC1')
# proECData= gLAll[(gLAll.EC1==1)&(gLAll.type=='Protein') & ~(gLAll.EC1.isnull()) & (gLAll.valence.between(gLAll.valence.quantile(0.05), gLAll.valence.quantile(.98)))]
proECData= gLAll[(gLAll.atomName=='FE')&(gLAll.type=='Protein') & ~(gLAll.EC1.isnull()) & (gLAll.valence.between(gLAll.valence.quantile(0.05), gLAll.valence.quantile(.98)))]
print(proECData.groupby('EC1').valence.describe())#.count)

sb.violinplot(x="atomName",y="valence",hue="EC1",data=proECData,ax=axEC1Type,legend=['Protein'])
# axValType.legend(('Protein','Mineral'))
outEC1Type = os.path.join(rdir,'VP_byFe_EC_pro_val_percentile.png')
plt.savefig(outEC1Type)
# exit()
ec = ['EC1','EC2','EC3','EC4','EC5','EC6','EC7']
for e1 in range(0,len(ec)):
    for e2 in range(e1+1,len(ec)):
        e1Data = proECData.where((proECData.EC1 == e1))['valence'].dropna()
        e2Data = proECData.where((proECData.EC1 == e2))['valence'].dropna()
        print(str(e1)+' ' + str(e2) + ' val t-test')
        print(scipy.stats.ttest_ind(e1Data,e2Data))
        print("\n")

#facet grid for EC1 and atomType
axFacet = plt.axes(label="facet")
proEC_Metal_Data = gLAll[(gLAll.type=='Protein') & ~(gLAll.EC1.isnull()) & (gLAll.valence.between(gLAll.valence.quantile(0.05), gLAll.valence.quantile(.98)))]
g = sb.FacetGrid(proEC_Metal_Data, col="atomName", col_wrap=4, height=2, col_order=order,)
g.map(sb.violinplot, "EC1", "valency", order=order,color=".3", ci=None,ax=axFacet)
outFacetECType = os.path.join(rdir,'VP_Facet_EC_pro_val_percentile.png')
plt.savefig(outFacetECType)