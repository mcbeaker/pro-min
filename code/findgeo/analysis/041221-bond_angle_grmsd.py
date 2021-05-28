
from matplotlib.markers import MarkerStyle
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import pandas as pd
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sn
import importlib 
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.size'] = 12

matplotlib.use('Agg')
rdir='/home/kenneth/proj/proMin/results'

getData = importlib.import_module("get_data")

dfBL = getData.get_bondLengthInfo()
# print(dfBL)
dfAngle = getData.getAngle()
# print(dfAngle)
dfGrmsd = getData.grmsd()
# print(dfGrmsd)
# exit()
fig, ax = plt.subplots(1, 3,figsize=(12, 5), sharey=False)

# bond len
sn.swarmplot(data=dfBL,y='Bond_Length',x='geo',hue='Environment',size=3.5,ax=ax[0],marker=MarkerStyle(marker='o', fillstyle='none'))
ax[0].set_title('A.',y=1.02, fontsize = 16,loc='left',x=-0.15)
ax[0].set_yticks([2.0,2.2,2.4,2.6,2.8])
ax[0].yaxis.set_minor_locator(AutoMinorLocator(2))
ax[0].set_xlabel("Fe Coordination Geometry")
ax[0].set_ylabel("Fe-S Bond Length (Å)")

# bond angle
sn.swarmplot(data=dfAngle,y='Degrees',x='geoShort',hue='Environment',size=3.5,ax=ax[1],marker=MarkerStyle(marker='o', fillstyle='none'))
ax[1].set_title('B.',y=1.02, fontsize = 16,loc='left',x=-0.15)
ax[1].set_yticks(list(range(80,140,10)))
ax[1].yaxis.set_minor_locator(AutoMinorLocator(2))
ax[1].set_ylabel("S-Fe-[S/N]* Bond Angle (°)")
ax[1].set_xlabel("Fe Coordination Geometry")

# grmsd
sn.swarmplot(data=dfGrmsd,y='gRMSD',x='geoShort',hue='Environment',size=3.5,ax=ax[2],marker=MarkerStyle(marker='o', fillstyle='none'))
ax[2].set_title('C.',y=1.02, fontsize = 16,loc='left',x=-0.15)
ax[2].set_yticks([0,0.1,0.2,0.3,0.4,0.5])
ax[2].yaxis.set_minor_locator(AutoMinorLocator(2))
ax[2].set_ylabel("gRMSD (Å)")
ax[2].set_xlabel("Fe Coordination Geometry")
# for tick in ax[2].get_xticklabels():
#     tick.set_fontname("Arial")
# for tick in ax[2].get_yticklabels():
#     tick.set_fontname("Arial")
plt.tight_layout()
plt.savefig(os.path.join(rdir,'041321-bl_angle_grmsd.png'),dpi=400)
