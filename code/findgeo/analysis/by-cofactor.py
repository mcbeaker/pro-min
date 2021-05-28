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
import matplotlib.ticker as ticker

rdir = '/home/kenneth/proj/proMin/results'
HEME = ['HEA','SRM','HAS','HEM','DHE','HEC']
SF = ['ICS','FS3','CFM','CLF','SF4','FES']
manganese = ['MN3']
# MN3; Manganese (III) Ion

copper = ['CUA','CU1','CU']
# CUA;DINUCLEAR COPPER ION
# CU1
# CU

# FEO; MU-OXO-DIIRON

nickel = ['F43','3NI']
# F43; factor 430
# 3NI; Nickel (III) Ion

tungsten = ['W02','W']
# W02; OCTADECATUNGSTENYL DIPHOSPHATE
# W; Tungsten Ion

molybdenum = ['MOS','MO','MOO']
# MOS; DIOXOTHIOMOLYBDENUM(VI) ION
# MO; Molybdenum atom
# MOO; molybdenum ion

cobolt = ['CNC','B12','NCO','3CO','COH']
# COH; PROTOPORPHYRIN IX CONTAINING CO
# CNC; CO-CYANOCOBALAMIN
# B12; COBALAMIN
# NCO; COBALT HEXAMMINE(III)
# 3CO; Cobalt (III) ION

vanadium = ['V70','VN4']
# V70; Meta Vanadate
# VN4; oxido(dioxo)vandium

# FC6; HEXACYANOFERRATE(3-)

# HEA; Heme-A
# SRM; SIROHEME
# HAS; HEME-AS
# HEM; Heme
# DHE; Heme D
# HEC; Heme C

# ICS; iron-sulfur-molybdenum cluster with interstitial carbon
# FS3; FE3-S4
# CFM; FE-MO-S CLUSTER
# CLF; FE(8)-S(7) CLUSTER

#combine min and pro
dfAll = pd.read_csv(os.path.join(rdir,'nr_40','nr_40_pro_min.csv'),header=0,index_col=0)
gLAll = dfAll.groupby(by=['type','pdb','atomName','resNum','envComp'],sort=False,as_index=False).first()
gLAll = gLAll[~gLAll.valence.le(0)]
dfFE = gLAll[((gLAll.Cofactor == 'FE') | (gLAll.Cofactor == 'FEO')) & (gLAll.valence < 8)].copy(deep=True)
dfFE.Cofactor = 'FE'
# print(dfFe.valence > 4)
# exit()
# exit()
print('Starting group by type/pdb/atomName/resNum/envComp metal sites\n')
# print("EC Analysis\n")
axCoType = plt.axes(label='cofactor')
proCoData= gLAll[(gLAll.atomName == 'FE') & (gLAll.Cofactor.isin(HEME) | gLAll.Cofactor.isin(SF)) & (gLAll.type=='Protein') & (gLAll.valence.between(gLAll.valence.quantile(0.05), gLAll.valence.quantile(.98)))]
print(len(proCoData.atomName))
ifHeme = proCoData.Cofactor.isin(HEME)
# exit()
ifFeS = proCoData.Cofactor.isin(SF)
ifMN = gLAll.Cofactor.isin(manganese)
ifNI = gLAll.Cofactor.isin(nickel)
ifV = gLAll.Cofactor.isin(vanadium)
ifW = gLAll.Cofactor.isin(tungsten)
ifCO = gLAll.Cofactor.isin(cobolt)
ifMO = gLAll.Cofactor.isin(molybdenum)
ifCU = gLAll.Cofactor.isin(copper)

MN = gLAll[ifMN].copy(deep=True)
MN['Cofactor'] = 'Mn'

NI = gLAll[ifNI].copy(deep=True)
NI['Cofactor'] = 'Ni'
V = gLAll[ifV].copy(deep=True)
V['Cofactor'] = 'V'
W = gLAll[ifW].copy(deep=True)
W['Cofactor'] = 'W'
CO = gLAll[ifCO].copy(deep=True) 
CO['Cofactor'] = 'CO'
MO = gLAll[(ifMO)].copy(deep=True)
MO['Cofactor'] = 'MO'
CU = gLAll[(ifCU)].copy(deep=True)
CU['Cofactor'] = 'CU'
# print(CU[CU.valence > 2][['valence','Microenvironment_ID','envComp']])
# valence               Microenvironment_ID      envComp
# 76       2.06  1cyx.CUA_A_316_electrontransport  CU1C1N1O1S2
# 939      2.04     2j5w.CU_A_3046_oxidoreductase      CU1N2S1
# 1320     2.42       3cdz.CU_A_757_bloodclotting    CU1C3N1S1
# 1822     2.88      3x1e.CU_A_407_oxidoreductase        CU1O2
# 2133     2.08           4oak.CU_A_201_hydrolase    CU1C2N3O3

heme=proCoData[ifHeme].copy(deep=True)
# print(len(heme.index))
heme['Cofactor'] = 'Heme'
FeS=proCoData[ifFeS].copy(deep=True)
# print(len(FeS.index))
FeS['Cofactor'] = 'FeS'

# print('Heme descriptors\n')
# print(heme.describe())
# print()

# print('FeS descriptors\n')

# print(FeS.describe())
# print('heme_fes t-test\n')
# print(scipy.stats.ttest_ind(heme.valence,FeS.valence))

print('FE descriptors\n')
print(dfFE.describe())

# print(MN)
order = ['Heme','FeS','Fe','Mn','Mo','Cu','Co','V','W','Ni']

heme_FeS = pd.concat([heme,FeS,dfFE,MN,MO,CU,CO,V,W,NI])
heme_FeS = heme_FeS[(heme_FeS.valence.between(heme_FeS.valence.quantile(.05), heme_FeS.valence.quantile(.98)))]

# for o in range(0,len(order)):
#     for p in range(o+1,len(order)):
#         c1 = heme_FeS.where(heme_FeS.Cofactor.str.upper() == order[o].upper())['valence'].dropna()
#         c2 = heme_FeS.where(heme_FeS.Cofactor.str.upper() == order[p].upper())['valence'].dropna()
#         # print(c1,c2)
#         print(str(order[o]) + ' ' + str(order[p]) + ' val t-test')
#         print(scipy.stats.ttest_ind(c1,c2))
#         print("\n")

# print(heme_FeS[heme_FeS.valence < 0.2][['valence','Microenvironment_ID','envComp']])
# 1302     0.14  3bvf.FE_A_1006_oxidoreductase    FE1O2 # bound/close to FE1007
# 1646     0.17     3pqi.FE_A_246_viralprotein    FE1N2 #probably 
# 2084     0.17   4mhu.FE_A_401_oxidoreductase  FE1N2O1

sb.violinplot(x="Cofactor",y="valence",data=heme_FeS,scale='area',ax=axCoType)

# sb.stripplot(x="atomName", y="valence",hue="Cofactor",dodge=False, data=heme_FeS[(heme_FeS.valence != 0)], jitter=1,ax=axCoType,zorder=1,palette=[(0,0,0),(0.5,0.5,0.5)],alpha=0.5)
# axCoType.legend(coOrder)
axCoType.set_ylabel('Valence (vu)')
axCoType.set_xlabel("Cofactor",labelpad=10)
# axBPVal.set_title('Number of Metal Coordination sites for Proteins and Minerals ')
x = np.arange(len(order))  # the label locations
axCoType.xaxis.set_major_locator(ticker.FixedLocator(x))
axCoType.xaxis.set_major_formatter(ticker.FixedFormatter(order))
plt.tight_layout()
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=True,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=True) # labels along the bottom edge are off
outCoType = os.path.join(rdir,'nr_cof_per_VP_pro.png')
plt.savefig(outCoType)

# HEM = ['DHE','HEC',]
# ['HEM' 'CU1' 'MN' 'SF4' nan 'DHE' 'CU' 'WO4' 'FE' 'CO' 'NI' 'FE2' 'VO4'
#  'C2O' 'HEC' 'CNC' 'CUA' 'FES' 'B12' 'FEO' 'F43' 'F3S' 'MOO' 'CFM' 'CLF'
#  'WO2' 'MN3' 'MOS' 'HEA' 'MO' '3NI' 'NCO' 'V7O' 'FC6' 'SRM' '3CO' 'VN4'
#  'HAS' 'COH' 'W' 'ICS']



# exit()
# ec = ['EC1','EC2','EC3','EC4','EC5','EC6','EC7']
# for e1 in range(1,8):
#     for e2 in range(e1+1,8):
#         e1Data = proECData.where((proECData.EC1 == e1))['valence'].dropna()
#         e2Data = proECData.where((proECData.EC1 == e2))['valence'].dropna()
        # print(str(e1)+' ' + str(e2) + ' val t-test')
        # print(scipy.stats.ttest_ind(e1Data,e2Data))
        # print("\n")
# exit()
#EC on x, and valence on Y
# ec number analysis
# axEC17Type = plt.axes(label='EC17')
# # proECData= gLAll[(gLAll.EC1==1)&(gLAll.type=='Protein') & ~(gLAll.EC1.isnull()) & (gLAll.valence.between(gLAll.valence.quantile(0.05), gLAll.valence.quantile(.98)))]
# proEC17Data= gLAll[(gLAll.type=='Protein') & ~(gLAll.EC1.isnull()) & (gLAll.valence.between(gLAll.valence.quantile(0.05), gLAll.valence.quantile(.98)))]
# proEC17Data = proEC17Data.astype({"EC1":int})
# proEC17Data = proEC17Data.astype({"EC1":str})
# print(proECData.groupby('EC1').valence.describe())#.count)

# font = {'family' : 'normal',
#         'weight' : 'normal',
#         'size'   : 28}

# SMALL_SIZE = 8
# MEDIUM_SIZE = 10
# BIGGER_SIZE = 30

# plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
# plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
# plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
# plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
# plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# matplotlib.rc('figure', titlesize=BIGGER_SIZE)

# sb.violinplot(x="EC1",y="valence",data=proEC17Data,ax=axEC17Type)

# # axValType.legend(('Protein','Mineral'))
# outEC17Type = os.path.join(rdir,'VP_EC17_pro_val_percentile.png')
# plt.savefig(outEC17Type)

# #facet grid for EC1 and atomType
# proEC_Metal_Data = gLAll[(gLAll.EC1 != 7) & (gLAll.type=='Protein') & ~(gLAll.EC1.isnull()) & (gLAll.valence.between(gLAll.valence.quantile(0.05), gLAll.valence.quantile(.98)))]
# proEC_Metal_Data = proEC_Metal_Data.astype({"EC1":int})
# # g = sb.catplot(x="EC1", y="valence",col="atomName",data=proEC_Metal_Data, kind="violin",height=4, aspect=.7,col_order=order,col_wrap=4)
# proEC_Metal_Data = proEC_Metal_Data.astype({"EC1":str})
# proEC_Metal_Data = proEC_Metal_Data.astype({"atomName":str})

# # print(proEC_Metal_Data.atomName.dtype)
# sb.set(font_scale = 2.5)
# g = sb.catplot(x="atomName", y="valence",orient='v',data=proEC_Metal_Data, kind="violin",height=10, aspect=1,order=order, col_wrap =3, col='EC1',sharex=False)

# # g = sb.FacetGrid(proEC_Metal_Data, col="atomName", col_wrap=4, height=2, col_order=order)
# # g.map(sb.violinplot, x="EC1", y="valence",color=".3", order = ec, ci=None,ax=axFacet)
# outFacetECType = os.path.join(rdir,'VP_colAtomName_EC_pro_val_percentile.png')
# plt.savefig(outFacetECType)

# #compare 1 and 3 EC
# axEC13Type = plt.axes(label='EC13')
# # proECData= gLAll[(gLAll.EC1==1)&(gLAll.type=='Protein') & ~(gLAll.EC1.isnull()) & (gLAll.valence.between(gLAll.valence.quantile(0.05), gLAll.valence.quantile(.98)))]
# pro13ECData= gLAll[((gLAll.EC1==1)|(gLAll.EC1==3))&(gLAll.type=='Protein') & 
#     ~(gLAll.EC1.isnull()) & (gLAll.valence.between(gLAll.valence.quantile(0.05), gLAll.valence.quantile(.98)))]
# print(pro13ECData.groupby('EC1').valence.describe())#.count)

# sb.violinplot(x="atomName",y="valence",hue="EC1",split=True,data=pro13ECData,ax=axEC13Type,legend=['EC1','EC3'])
# # axValType.legend(('Protein','Mineral'))
# axEC13Type.set_ylabel('Valence (vu)',fontsize=18)
# axEC13Type.set_xlabel(xlabel='Metal')
# outEC13Type = os.path.join(rdir,'VP_13EC_pro_val_percentile.png')
# plt.savefig(outEC13Type)

# for o in range(0,len(order)):
#     e1Data = pro13ECData.where((pro13ECData.EC1 == 1) & (pro13ECData.atomName == order[o]))['valence'].dropna()
#     e2Data = pro13ECData.where((pro13ECData.EC1 == 3) & (pro13ECData.atomName == order[o]))['valence'].dropna()
#     print(str(order[o]) + ' ' + str(order[o]))
#     print(str(1)+' ' + str(3) + ' val t-test')
#     print(scipy.stats.ttest_ind(e1Data,e2Data))
#     print("\n")