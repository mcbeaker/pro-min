#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import pandas as pd

import os

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')


rdir = '/home/kenneth/proj/proMin/results'

#combine min and pro
dfAll = pd.read_csv(os.path.join(rdir,'ec_min_pro_fg_dictionaryValues.csv'),header=0,index_col=False)

mdir = '/home/kenneth/proj/proMin/minerals/database/data'
dfMin = pd.read_csv(os.path.join(mdir,'RRUFF_AMCSD_Fe_Mineral_IDs_3-shaunna-retrieved.csv'),header=0,index_col=False)

#merge 
dfAgeMinPro = pd.merge(dfAll,dfMin,how='outer',left_on='mineral',right_on='Mineral_Name_plain')
print(dfAgeMinPro)
#distributions of mineral and protein bulk valence
bulkVal = dfAgeMinPro.groupby(by=['pdb','atomName','resNum','envComp'],sort=False,as_index=False).first()
# print(bulkVal.sort_values(['valence'],axis=0,inplace=False,ascending=False))
bulkValData= bulkVal[bulkVal.valence.between(bulkVal.valence.quantile(.05), bulkVal.valence.quantile(.98))]
bulkValData=bulkValData[(bulkValData.type == 'Mineral') & (bulkValData.valence != 0) & (~bulkValData.Oldest_Known_Age_Ma.isnull())]
bulkValData['val_ENV'] = 'Metal Coordination Environment'

x = bulkValData.valence

fig, ax = plt.subplots()
violins = ax.violinplot(x)

# ymin, ymax = ax.get_ylim()
# xmin, xmax = ax.get_xlim()

# # create a numpy image to use as a gradient
# Nx,Ny=1,1000
# imgArr = np.tile(np.linspace(0,1,Ny), (Nx,1)).T
# cmap = 'hsv'

# for violin in violins['bodies']:
#     path = Path(violin.get_paths()[0].vertices)
#     patch = PathPatch(path, facecolor='none', edgecolor='none')
#     ax.add_patch(patch)
#     img = ax.imshow(imgArr, origin="lower", extent=[xmin,xmax,ymin,ymax], aspect="auto",
#                     cmap=cmap,
#                     clip_path=patch)

# # colorbar
# ax_divider = make_axes_locatable(ax)
# cax = ax_divider.append_axes("right", size="5%", pad="2%")
# norm = matplotlib.colors.Normalize(vmin=ymin, vmax=ymax)
# cb = matplotlib.colorbar.ColorbarBase(cax, cmap=matplotlib.cm.get_cmap(cmap),
#                                 norm=norm,
#                                 orientation='vertical')
plt.savefig(os.path.join(rdir,'gradient.png'))