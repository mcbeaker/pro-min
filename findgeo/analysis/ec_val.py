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
gLAll = dfAll.groupby(by=['type','pdb','atomName','resNum','envComp'],sort=False,as_index=False).first()
print('Starting group by type/pdb/atomName/resNum/envComp metal sites\n')
print("EC Analysis\n")
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
for e1 in range(1,8):
    for e2 in range(e1+1,8):
        e1Data = proECData.where((proECData.EC1 == e1))['valence'].dropna()
        e2Data = proECData.where((proECData.EC1 == e2))['valence'].dropna()
        # print(str(e1)+' ' + str(e2) + ' val t-test')
        # print(scipy.stats.ttest_ind(e1Data,e2Data))
        # print("\n")
# exit()
#EC on x, and valence on Y
# ec number analysis
axEC17Type = plt.axes(label='EC17')
# proECData= gLAll[(gLAll.EC1==1)&(gLAll.type=='Protein') & ~(gLAll.EC1.isnull()) & (gLAll.valence.between(gLAll.valence.quantile(0.05), gLAll.valence.quantile(.98)))]
proEC17Data= gLAll[(gLAll.type=='Protein') & ~(gLAll.EC1.isnull()) & (gLAll.valence.between(gLAll.valence.quantile(0.05), gLAll.valence.quantile(.98)))]
proEC17Data = proEC17Data.astype({"EC1":int})
proEC17Data = proEC17Data.astype({"EC1":str})
# print(proECData.groupby('EC1').valence.describe())#.count)

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 28}

SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 30

# plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
# plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
# plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
# plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
# plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# matplotlib.rc('figure', titlesize=BIGGER_SIZE)

sb.violinplot(x="EC1",y="valence",data=proEC17Data,ax=axEC17Type)

# axValType.legend(('Protein','Mineral'))
outEC17Type = os.path.join(rdir,'VP_EC17_pro_val_percentile.png')
plt.savefig(outEC17Type)

#facet grid for EC1 and atomType
proEC_Metal_Data = gLAll[(gLAll.EC1 != 7) & (gLAll.type=='Protein') & ~(gLAll.EC1.isnull()) & (gLAll.valence.between(gLAll.valence.quantile(0.05), gLAll.valence.quantile(.98)))]
proEC_Metal_Data = proEC_Metal_Data.astype({"EC1":int})
# g = sb.catplot(x="EC1", y="valence",col="atomName",data=proEC_Metal_Data, kind="violin",height=4, aspect=.7,col_order=order,col_wrap=4)
proEC_Metal_Data = proEC_Metal_Data.astype({"EC1":str})
proEC_Metal_Data = proEC_Metal_Data.astype({"atomName":str})

# print(proEC_Metal_Data.atomName.dtype)
sb.set(font_scale = 2.5)
g = sb.catplot(x="atomName", y="valence",orient='v',data=proEC_Metal_Data, kind="violin",height=10, aspect=1,order=order, col_wrap =3, col='EC1',sharex=False)

# g = sb.FacetGrid(proEC_Metal_Data, col="atomName", col_wrap=4, height=2, col_order=order)
# g.map(sb.violinplot, x="EC1", y="valence",color=".3", order = ec, ci=None,ax=axFacet)
outFacetECType = os.path.join(rdir,'VP_colAtomName_EC_pro_val_percentile.png')
plt.savefig(outFacetECType)

#compare 1 and 3 EC
axEC13Type = plt.axes(label='EC13')
# proECData= gLAll[(gLAll.EC1==1)&(gLAll.type=='Protein') & ~(gLAll.EC1.isnull()) & (gLAll.valence.between(gLAll.valence.quantile(0.05), gLAll.valence.quantile(.98)))]
pro13ECData= gLAll[((gLAll.EC1==1)|(gLAll.EC1==3))&(gLAll.type=='Protein') & 
    ~(gLAll.EC1.isnull()) & (gLAll.valence.between(gLAll.valence.quantile(0.05), gLAll.valence.quantile(.98)))]
print(pro13ECData.groupby('EC1').valence.describe())#.count)

sb.violinplot(x="atomName",y="valence",hue="EC1",split=True,data=pro13ECData,ax=axEC13Type,legend=['EC1','EC3'])
# axValType.legend(('Protein','Mineral'))
axEC13Type.set_ylabel('Valence (vu)',fontsize=18)
axEC13Type.set_xlabel(xlabel='Metal')
outEC13Type = os.path.join(rdir,'VP_13EC_pro_val_percentile.png')
plt.savefig(outEC13Type)

for o in range(0,len(order)):
    e1Data = pro13ECData.where((pro13ECData.EC1 == 1) & (pro13ECData.atomName == order[o]))['valence'].dropna()
    e2Data = pro13ECData.where((pro13ECData.EC1 == 3) & (pro13ECData.atomName == order[o]))['valence'].dropna()
    print(str(order[o]) + ' ' + str(order[o]))
    print(str(1)+' ' + str(3) + ' val t-test')
    print(scipy.stats.ttest_ind(e1Data,e2Data))
    print("\n")