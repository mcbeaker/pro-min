#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.ticker as ticker
matplotlib.use('Agg')
import seaborn as sb
import numpy as np
import scipy
import re as re
from matplotlib.lines import Line2D

transform = {'CO1O6':'CoO','CU1O2':'CuO','CU1O3':'CuO','CU1O4':'CuO','CU1O5':'CuO','CU1O6':'CuO',
'CU1S2':'CuS','CU1S3':'CuS','CU1S4':'CuS',
'FE1O2':'FeO','FE1O3':'FeO','FE1O4':'FeO','FE1O5':'FeO','FE1O6':'FeO',
'FE1S2':'FeS','FE1S3':'FeS','FE1S4':'FeS',
'MN1O2':'MnO','MN1O3':'MnO','MN1O4':'MnO','MN1O5':'MnO','MN1O6':'MnO',
'MO1O2':'MoO','MO1O4':'MoO','MO1O5':'MoO','MO1O6':'MoO','MO1S6':'MoS',
'NI1O6':'NiO','NI1S2':'NiS','NI1S4':'NiS','V1O4':'VO','V1O5':'VO','W1O4':'WO','W1O6':'WO'}

CoO = (19/255,0/255,255/255,1)
CuO =  (16/255,204/255,205/255,1)
CuS = (184/255,115/255,51/255,1)
NiO = (0/255,141/255,124/255,1)
NiS = (114/256,116/256,114/256,1)
FeO = (253/255,230/255,56/255,1)
FeS =  (255/255,0/255,0/255,1)
MnO = (0/255,255/255,0/255,1) #green
MoO = (0/255,0/255,255/255,1)
MoS = (192/255,192/255,192/255,1)
VO = (255/255,128/255,0/255,1)
WO = (255/255, 214/255, 170/255,1) #tungsten 100w

labels = {'CoO':CoO,'CuO':CuO,'CuS':CuS,'NiO':NiO, 
'NiS':NiS,'FeO':FeO,'FeS':FeS,'MnO':MnO,'MoO':MoO,'MoS':MoS,'VO':VO,'WO':WO}


rdir = '/home/kenneth/proj/proMin/results'

#combine min and pro
# dfAll = pd.read_csv(os.path.join(rdir,'min_pro_fg_dictionaryValues.csv'),header=0,index_col=False)
dfNR = pd.read_csv(os.path.join(rdir,'nr_40','nr_40_pro_min.csv'),header=0,index_col=0)
dfNR['findgeo'] = 'yes'
dfNR['envSimple']= dfNR['envComp'].map(transform)
dfNR = dfNR[(~dfNR.envSimple.isnull()) & (dfNR.valence.between(dfNR.valence.quantile(0.05), dfNR.valence.quantile(.98)))]
dfNR['color']= dfNR['envSimple'].map(labels)

dfNR = dfNR.groupby(by=['type','pdb','atomName','resNum','envComp'],sort=False,as_index=False).first()
print(dfNR.type)
# exit()

axEnv = plt.axes(label='envComp')
sb.violinplot(x="envSimple",y="valence",split=True,hue='type',data=dfNR,scale='area',ax=axEnv)#,palette=labels)

# sb.stripplot(x="atomName", y="valence",hue="Cofactor",dodge=False, data=heme_FeS[(heme_FeS.valence != 0)], jitter=1,ax=axCoType,zorder=1,palette=[(0,0,0),(0.5,0.5,0.5)],alpha=0.5)
# axCoType.legend(coOrder)
axEnv.set_ylabel('Valence (vu)')
axEnv.set_xlabel("Environment",labelpad=10)
# axBPVal.set_title('Number of Metal Coordination sites for Proteins and Minerals ')
x = np.arange(len(labels.keys()))  # the label locations
axEnv.xaxis.set_major_locator(ticker.FixedLocator(x))
axEnv.xaxis.set_major_formatter(ticker.FixedFormatter(list(labels.keys())))
plt.tight_layout()
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=True,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=True) # labels along the bottom edge are off
outEnv = os.path.join(rdir,'nr_first_envSimple_per_VP.png')
plt.savefig(outEnv)

exit()


rdir = '/home/kenneth/proj/proMin/results'
#combine min and pro
# dfMerge.to_csv(os.path.join(rdir,'ec_min_pro_fg_dictionaryValues.csv'))

# print(dfMerge)

# exit()




# print(len(dfMerge.index))
#extra 10k rows
# notRetrievedDF = dfMerge[(dfMerge.hagai == 'yes') & (dfMerge.envComp.isnull())]

#connect proteins and minerals if they have the same envComp
#source (Mineral), Target (Protein ID,module)
#filter geo irr
gr = dfMerge[(dfMerge.geo != 'irr') & (((dfMerge.type == 'Protein') & (~dfMerge.Binding_motif.isnull())) | (dfMerge.type == "Mineral"))]
# print(gr.geo)
# exit()

grGroup = gr.groupby(['envComp','geo']).agg({"gRMSD":lambda x: x.mean(),"valence": lambda x: x.mean(),'type':lambda x:",".join(np.unique(x)),'pdb':lambda x:",".join(np.unique(x)),'Binding_motif':lambda x:list(map(int,np.unique(x.dropna())))}).reset_index()
# print(grGroup)

# print(grGroup.Binding_motif)
# exit()
grGroup['pdb'] = [x.split(',') for x in grGroup.pdb]

dfMP = grGroup[grGroup.type == 'Mineral,Protein'].reset_index()

# print(dfMP.pdb.apply(lambda col: ",".join(col[i] for i, s in enumerate(col) if '-' in col)))
# print(dfMP.pdb.apply(lambda col: )
# print(dfMP.apply(lambda x: 1 if '-' in x['pdb'] else 0, axis=1))#.apply(lambda x : x.str.contains('-'))]
dfMP['minCoLoc'] = dfMP.pdb.apply(lambda x: [i for i in x if "-" in i])

#combinations of moduleID and minerals
from itertools import product
# dfModMin = df(list(product(list(dfMP.Binding_motif),list(dfMP.minCoLoc))))
def f(x):
    return list(product(*[x.Binding_motif,x.minCoLoc]))

from collections import defaultdict, Counter

d = defaultdict(list)
for k, v in list(dfMP.apply(f,axis=1).explode()):
    # print(k,v)
    d[k].append(v)
# print(Counter(d))
# dfMP['Binding_motif'] = [x.split(',') for x in dfMP.Binding_motif]
# dfMP['Binding_motif'] = dfMP.Binding_motif.astype(int)

# print(dfMP)
# exit()

# print(dfMP.pdb.value_counts(normalize=True))

dfMP.to_csv(os.path.join(rdir,'nr_groupedMinPro.csv'))



# print(dfMP)
# exit()
# print(dfMP.iloc[43]['pdb'])
# print(np.unique(dfMP.envComp))
# exit()
# colors = # 
CoO = (19/255,0/255,255/255,1)
CuO =  (16/255,204/255,205/255,1)
CuS = (184/255,115/255,51/255,1)
NiO = (0/255,141/255,124/255,1)
NiS = (114/256,116/256,114/256,1)
FeO = (253/255,230/255,56/255,1)
FeS =  (255/255,0/255,0/255,1)
MnO = (0/255,255/255,0/255,1) #green
MoO = (0/255,0/255,255/255,1)
MoS = (192/255,192/255,192/255,1)
VO = (255/255,128/255,0/255,1)
WO = (255/255, 214/255, 170/255,1) #tungsten 100w

labels = {'CoO':CoO,'CuO':CuO,'CuS':CuS,'NiO':NiO, 
'NiS':NiS,'FeO':FeO,'FeS':FeS,'MnO':MnO,'MoO':MoO,'MoS':MoS,'VO':VO,'WO':WO}

# order = ['FES','HEME','FE','NI','MN','CU','CO','MO','W','V']
# colors = [FES,HEME,FE,NI,MN,CU,CO,MO,W,V]
# print(dfMP.envComp)
# format
# envComp             type                                                pdb                                      Binding_motif
# 129     CO1O2  Mineral,Protein  [1hv9, 2oi6, 2v0h, 2vor, 2vos, 2wke, 3ka4, 3ka...  141.0,194.0,451.0,533.0,1190.0,1330.0,1441.0,1...

# edges = set()
nodesMin = set()
nodesPro = set()
# graph
import networkx as nx
from itertools import product
G = nx.Graph()   
# edges format (min,mod,env=env)
# get all minerals from list by retrieving strings with "-"
count = 0
for i,row in dfMP.iterrows():
    # print(row.Binding_motif)
    minList = [x for x in row.pdb if x.find('-') > -1]
    modList = row.Binding_motif#.split(",") 
    # print(minList)
    # print(modList)
    edges = product(*[minList,modList])
    if row.envComp.find('MN') > -1:
        count += 1
    # print(list(edges))
    G.add_edges_from(list(edges),geo=row.geo,envComp=row.envComp,label=transform[row.envComp],color=labels[transform[row.envComp]])
    # print(G.edges(data=True))
    G.add_nodes_from(minList,type='Mineral',bipartite=0,envComp=row.envComp)
    G.add_nodes_from(modList,type='Protein',bipartite=1,envComp=row.envComp)
    # exit()

print(len(G.nodes))
print(len(G.edges))
# print(G.number_of_nodes())
# print(G.number_of_edges())
minNodes = {n for n, d in G.nodes(data=True) if d['bipartite']==0}
proNodes = set(G) - minNodes
# print(G.edges())
pos = nx.bipartite_layout(G,minNodes,scale=1,aspect_ratio=4/3)
# print(pos)
#filter edges

# for env in colors.keys():
#     print(env)
#     selected_edges = [(u,v) for u,v,e in G.edges(data=True) if e['envComp'] == env]
#     # print(selected_edges)
#     d_edges = nx.draw_networkx_edges(G,pos=pos,edgelist=selected_edges,width=0.01,edge_color=labels[colors[env]])


    # print(d_edges)
# print (selected_edges)
edges = G.edges(data=True)
# print([(u,v,e) for u,v,e in edges if e['label'] == 'FeS'])
# print(edges)
labels = {e['label'] for u,v,e in edges}
colors = {e['color'] for u,v,e in edges}
d_edges = nx.draw_networkx_edges(G, pos, edges=edges, edge_color=colors)

from collections import Counter

# print(Counter([transform[e['envComp']] for u,v,e in edges]))

d_nodes = nx.draw_networkx_nodes(G,pos=pos,node_shape='o',align='horizontal',node_size=0.1)
# plt.rcParams['figure.dpi'] = 300
# plt.figure()
def make_proxy(clr,**kwargs):
    return Line2D([0, 1], [0, 1], color=clr, **kwargs)

# generate proxies with the above function
proxies = [make_proxy(eval(clr), lw=5) for clr in labels]
#https://stackoverflow.com/questions/19877666/add-legends-to-linecollection-plot - uses plotted data to define the color but here we already have colors defined, so just need a Line2D object.
plt.legend(proxies, labels,loc='right',ncol=1)

plt.tight_layout()
plt.savefig(os.path.join(rdir,'nr40_nx_min_pro_geo_mod_envComp.png'),dpi=300)

nx.write_edgelist(G,os.path.join(rdir,'nr_40_edgeList_protein_geo_modules_envComp.csv'),data=['env','geo','size'])
exit()
with open(os.path.join(rdir,'nr_40_nodes_list.csv'),'w') as fileWriter:
    fileWriter.write('id, type\n')
    for l in G.nodes:
        print(type(l))
        if '-' in l:
            fileWriter.write(l + ',0\n')
        else:
            fileWriter.write(l + ',1\n')

exit()

# exit()


 



#,'type','pdb'])#.Binding_motif.unstack()
# print(gr['pdb','Binding_motif'].agg('envComp',lambda x: ",".join(set(x)),'type',lambda x: ",".join(set(x)),'pdb',lambda x: ",".join(set(x)),'Binding_motif',lambda x: ",".join(set(x))))
grAgg = type(gr['pdb','Binding_motif'].agg(lambda x: list(set(x))))
print(grAgg.groupby('envComp').agg({'type':','.join(),'pdb':','.join(),'Binding_motif':",".join(list(map(str,np.unique(x))))}))
exit()
# gLAll = grAll.agg({'chain':lambda x: ','.join(set(x)),'valence':lambda x:np.round(np.average(x),decimals=2),'gRMSD':lambda x:np.round(np.average(x),decimals=2),'geo':lambda x: ",".join(set(x)) })#,'envComp':lambda x: ','.join(set(x)),'geo':lambda x: ",".join(set(x)) }





# print(gr.agg({'Binding_motif':lambda x: ",".join((map(str,x))),'envComp':lambda x: ",".join(x),'type':lambda x:','.join(x),'pdb':lambda x:','.join(x)}))
# print(grAgg)

#list) #,'mineral']].apply(lambda x: (x['pdb']))    #list)#['Binding_motif'].apply(list) 
# print(pd.DataFrame(gr))
#keys
# envComp
#values
# [['Mineral', 'Paganoite-0006855'], ['Protein', '1nzc'], ['Protein', '1uiv'], ['Protein', '2ixl'], ['Protein', '2j5s'], ['Protein', '2jih'], ['Protein', '2y1v'], ['Protein', '3foo'], ['Protein', '3mmu']]
# print(type(gr))
# print(gr['W1S3','Mineral','Catamarcaite-0006132' ])
# print(gr['W1S3','Protein','1xdy' ])
exit()
from collections import defaultdict
d = defaultdict(list)

listTuples = []
seen = 0
# print(gr.keys())
for k, v in gr.keys():
    print(v)
    d[k].append(v)
    # if len(k[1][0]) == 2:
        # print(k)
print(d.values())
# print(d.keys())
exit()
from itertools import product,chain
envComp_min_prot = [(i[0],gr[(i[0],'Mineral')],gr[(i[0],'Protein')]) for i in list(d.items()) if len(i[1])>1]

uniqueEnvProMin = [(envComp,list(np.unique(pro)),list(np.unique(min))) for envComp,min,pro in envComp_min_prot]

edgeList = []
import networkx as nx
G = nx.Graph()
nodes = set()

from progress.bar import Bar
bar = Bar("Processing",max=len(uniqueEnvProMin),fill='*',suffix='%(percent).1f%% - %(eta)ds')
for envComp,min,pro in uniqueEnvProMin:
    # print(type(envComp))
    # print(pro)
    # print(min)
    # print([x for a in permutations(pro,len(min))])
    # print([envComp]*10)
    #add nodes
    edges = list(product(*[pro,min]))
    nodes.update(list(chain(*edges)))
    G.add_edges_from(edges,envComp=envComp)
    # print(G.edges())
    # exit()
    bar.next()
G.add_nodes_from(nodes)
bar.finish()


nx.write_edgelist(G,os.path.join(rdir,'edgeList_protein_modules_envComp.csv'),data=['envComp'])

with open(os.path.join(rdir,'nodes_list.csv'),'w') as fileWriter:
    fileWriter.write('id, type\n')
    for l in nodes:
        if '-' in l:
            fileWriter.write(l + ',0\n')
        else:
            fileWriter.write(l + ',1\n')


exit()
# print(edgeList)
# print([ ( envComp,[ list(zip(x,list(np.unique(pro)))) for x in permutations(list(np.unique(min)),len(list(np.unique(pro))))]) for envComp,min,pro in envComp_min_prot])
# for e,m,p in envComp_min_prot:
#     print(type(e))
    # print(','.join(m)+'\n')
    # print(','.join(p)+'\n')

# print(envComp_min_prot)



# print(groupedMinProEnvComp[('FE1S3SB3', 'Mineral')])
# for key in groupedMinProEnvComp.keys():
    # print(key )
    

# print(grouped)
exit()
# print(dfMerge.columns)
# print(dfMerge[dfMerge.type == 'Protein']['Binding_motif'])

# merge if same chain,resnum
# print(dfMerge.groupby(by=['envComp','Binding_motif']))






grouped = dfMerge[(dfMerge.type == 'Protein') & (~dfMerge.Binding_motif.isnull())].groupby('envComp')['Binding_motif'].apply(list) 
print(grouped)



# for key in grouped.keys():
#     grouped.at[key] = np.unique(grouped.get(key))
#     combTuple = combinations(grouped.get(key),2)
#     G.add_edges_from(combTuple,envComp=key)
#     # print(G.edges())
# nx.write_edgelist(G,os.path.join(rdir,'edgeList_protein_modules_envComp.txt'),data=['envComp'])
# print(dir(grouped))

# for series in grouped: 
#     print(series.keys())
#     exit()
    # G.add_edge(l) 