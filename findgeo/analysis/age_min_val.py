
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
import statsmodels

rdir = '/home/kenneth/proj/proMin/results'

#combine min and pro
dfAll = pd.read_csv(os.path.join(rdir,'ec_min_pro_fg_dictionaryValues.csv'),header=0,index_col=False)

mdir = '/home/kenneth/proj/proMin/minerals/database/data'
dfMin = pd.read_csv(os.path.join(mdir,'RRUFF_AMCSD_Fe_Mineral_IDs_3-shaunna-retrieved.csv'),header=0,index_col=False)

#merge 
dfAgeMinPro = pd.merge(dfAll,dfMin,how='outer',left_on='mineral',right_on='Mineral_Name_plain')
# print(dfAgeMinPro)
#distributions of mineral and protein bulk valence
bulkVal = dfAgeMinPro.groupby(by=['type','pdb','envComp'],sort=False,as_index=False).first()
print(bulkVal)
# print(bulkVal.sort_values(['valence'],axis=0,inplace=False,ascending=False))
bulkValData= bulkVal[bulkVal.valence.between(bulkVal.valence.quantile(.05), bulkVal.valence.quantile(.95))]
bulkValData=bulkValData[(bulkValData.type == 'Mineral') & (bulkValData.valence != 0) & (~bulkValData.Oldest_Known_Age_Ma.isnull())]
bulkValData['val_ENV'] = 'Metal Coordination Environment'
print(np.unique(list(bulkValData.pdb)))
# print(bulkValData.sort_values(['valence'],axis=0,inplace=False,ascending=False))

# bulkValData=bulkValData[(bulkValData.type == 'Protein') & (bulkValData.valence != 0)]
# print(bulkValData)
plt.figure(num=None, figsize=(8, 8), dpi=300, facecolor='w', edgecolor='k')
axBulkVal = plt.axes(label="Bulk Valence")
axBulkVal.tick_params(axis='both', which='major', labelsize=18)
# plt.rcParams['ytick.labelsize'] = 'large'
from scipy import stats

X = bulkValData["Oldest_Known_Age_Ma"]
Y = bulkValData["valence"]

slope, intercept, r_value, p_value, slope_std_error = stats.linregress(X,Y)
predict_y = slope * X + intercept
print(slope,intercept,r_value, p_value, slope_std_error)
p = sb.regplot(x="Oldest_Known_Age_Ma",y="valence",ci=100,ax=axBulkVal,fit_reg=True,robust=True,
            data=bulkValData,
            line_kws={'label':'$y=%1.2e*x+%.2f$'%(-6.53e-05, 1.00298638524996497345e+00)})

slope, intercept, r_value, p_value, slope_std_error = stats.linregress(p.get_lines()[0].get_xdata(),p.get_lines()[0].get_ydata(),)
# predict_y = slope * X + intercept
# print(slope,intercept,r_value, p_value, slope_std_error)
print("slope: %0.2e    intercept: %0.20e" % (slope, intercept))

axBulkVal.legend()
axBulkVal.set_ylabel('Valence',fontsize=18)
axBulkVal.set_xlabel('Oldest Known Age of Mineral (Mya)',fontsize=18)

# axBulkVal.yaxis.get_label().set_fontsize(18)

# axBulkVal.set_xlabel('Valence',fontsize=24)
# cmap = sb.cubehelix_palette(as_cmap=True)

# sb.stripplot(x="type",y="Oldest_Known_Age_Ma",data=bulkValData,ax=axBulkVal,zorder=1,color='grey',cmap=cmap)
# print('test')
# fig_bulkVal = sb.swarmplot(x="type",y="valence",data=bulkValData,ax=axBulkVal).get_figure()
outMinAgeVal = os.path.join(rdir,'SP_min_Age_bulkVal_95percentile.png')
plt.savefig(outMinAgeVal)

#descriptive statistics
# minVal = bulkValData[(bulkValData.type == 'Mineral')]['valence']#.dropna()
# print('minVal descriptors\n')
# print(minVal.describe())
# proVal = bulkValData[(bulkValData.type == 'Protein')]['valence']#.dropna()
# print('proVal descriptors\n')
# print(proVal.describe())
# # print(proVal)
# print('val t-test\n')
# print(scipy.stats.ttest_ind(minVal,proVal))
# exit()