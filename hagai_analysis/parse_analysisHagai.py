import multiprocessing as mp
import ast
import os
from progress.bar import Bar
import pandas as pd
import string
from tqdm import tqdm
import glob
import re

pdir = "/home/kenneth/proj/proXtal/proteins/hagai/pdb"
os.chdir(pdir)

class SlowBar(Bar):
    suffix = '%(index)d/%(max)d - %(remaining_hours)d hours remaining'
    @property
    def remaining_hours(self):
        return self.eta // 3600

fileName = os.path.join(pdir,"micro_multiproc_wo_alt.csv")
df = pd.DataFrame()

with open(fileName) as f:
    l = list(map(ast.literal_eval,f.readlines()))
    df = pd.DataFrame(l)
    df.to_csv(os.path.join(pdir,"df_micro.csv"),index=True)


# outFile = os.path.join(pdir,"df_micro.csv") 
# df = pd.read_csv(outFile,index_col=0)
# df = df.iloc[df.isnull().sum(axis=1).argsort()]
# df_notHeme = df[~df.ligName.str.contains("HE")].dropna(axis='columns',how='all')
# df_FE = df[df.ligName.str.match("FE")].dropna(axis='columns',how='all') 
# # print(dir(df.ligName.str)) #  ("FE")].dropna(axis='columns',how='all') 
# # print(df[df.id_hagai.str.contains("1eb7")])
# print(df_FE[:1])
# print(df_FE.iloc[1].to_string())


# df = df[~df['id_hagai'].str.contains("3p9p",na=False) & ~df['id_hagai'].str.contains('2vxi',na=False)]
# df = df.dropna(axis='columns',how='all')
# print(df)
# df = df.dropna(axis='columns',how='all')
# print(df.isnull().sum(axis=1).head())
# print(df.loc[df.loc[17201].isnull())
# df.to_csv(os.path.join(pdir,"df_micro_wrangled.csv"))
# .sort_values(ascending=True).index])
# spl = lambda x : pd.Series(re.split('\._',x))
# df[['pdb','LIG','chain','resID','function']] = df.id_hagai.apply(spl)
# print(df)

# print(df.dropna())
# d = df.drop(index=[17395,17063,17140,17232,17310,17376,11763,17035,17045,17102])
# d = df.dropna(axis=1,how='all')
# dx = d[d.columns[d.columns.to_series().str.contains('dx')]]
# dy = d[d.columns[d.columns.to_series().str.contains('dy')]]
# dz = d[d.columns[d.columns.to_series().str.contains('dz')]]
# print(d)
# print(dx.dropna())
# print(d.loc[11763])
# print(d.loc[17035])
# print(d.loc[17045])
# print(d.loc[17102])

    # print(df)

# files = glob.glob("df_x*.csv")
# df = pd.DataFrame()

# for xF in files:
#     print(xF)
#     # if os.path.exists(os.path.join(pdir,"df_"+xF+".csv")) == False:

#     with open(xF,'r') as f:
#         lines = pd.read_csv(xF,index_col=False)
#         with tqdm(total=len(files)) as pbar: 
#             df = pd.concat([df,lines],sort=True)
#             pbar.update(1)

    #             # print(len(df.index))
    #             # if idx == 1000:
    #                 # break
    #     # print(s)
    #     # table = s.maketrans("","",)
    #     # l = s.translate(table,bad_chars)
    #     # print(l)
    #     # df = pd.concat(ast.literal_eval(l)[0],sort=True)
    #     # print(df)
    # # print(len(df.index))

# df.to_csv(os.path.join(pdir,"df_micro_combine.csv"),index=False)




    