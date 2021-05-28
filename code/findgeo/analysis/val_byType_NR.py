import sys 
import os
sys.path.append(os.path.abspath("/home/kenneth/proj/proMin/code/tools"))
import tools
import matplotlib.pyplot as plt
import seaborn as sb

rdir = '/home/kenneth/proj/proMin/results'

df = tools.load_dfAll_NR()
order = ["CO","CU","FE","MN","MO","NI","V","W"]

df = df.groupby(by=['type','pdb','atomName','resNum','envComp'],sort=False,as_index=False).first()
print(df[df.valence > 8][['envComp','valence','pdb','Microenvironment_ID']])
print(df[df.valence < 0][['envComp','valence','pdb','Microenvironment_ID']])
# print(gLAll.mean().valence)
# print(gLAll.first().valence)
# exit()

#------VAL-bytype-violin
ax = plt.axes(label='ax')
df = df[df.valence.between(df.valence.quantile(0.0), df.valence.quantile(1))]
# valData= gLAll[gLAll.valence.between(gLAll.valence.quantile(0), gLAll.valence.quantile(1))]
print(df.groupby('atomName').valence.describe())#.count)

sb.violinplot(x="atomName",y="valence",order=order,hue="type",data=df,split=True,ax=ax,legend=['Protein','Mineral'])

ax.set_ylabel('Valence (vu)',fontsize=18)
ax.set_xlabel(xlabel='Metal')
outValType = os.path.join(rdir,'nr40_valFirst_VP_type_per_min_pro.png')
plt.savefig(outValType)
# -------Val-bytype-violin

#negative valence for MO and W
# Metals within the bounds of metals