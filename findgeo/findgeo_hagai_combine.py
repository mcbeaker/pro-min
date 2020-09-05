import pandas as pd
import numpy as np

hagaiID = '/home/kenneth/proj/proXtal/proteins/hagai/id_bindingMotif.csv'
hagaiID_df = pd.read_csv(hagaiID,header=0,index_col=False)
#hagaiID,bindingMotif

#------format hagai data---------
hagaiID_df[['pdbID','res_chain_resNum_function']] = hagaiID_df.hagaiID.str.split('.',expand=True)
hagaiID_df[['resID','chain','resNum','function']] = hagaiID_df.res_chain_resNum_function.str.split('_',expand=True)
hagaiID_df = hagaiID_df.drop(columns=['res_chain_resNum_function'])
hagaiID_df = hagaiID_df.drop(columns=['hagaiID'])
hagaiID_df['resNum'] = pd.to_numeric(hagaiID_df['resNum'])
# print(hagaiID_df)

df = pd.read_csv('/home/kenneth/software/findgeo/FindGeoSummative_wEnvComp_reducedDuplicates.csv',header=0,index_col=False)

# print(hagaiID_df.columns)
# Index(['bindingMotif', 'pdbID', 'resID', 'chain', 'resNum', 'function'], dtype='object')
# print(df.columns)
#"pdbID,atomNum,geoShort,geoLong,envComposition,element,resNum,chain"

# df.index_col = pd.Int64Index(range(1,len(df)+1))
print(df.head())
print(hagaiID_df.head())
# df['bindingMotiff'] = ""
# df['function'] = ""

print(type(df.resNum[1]))
print(type(hagaiID_df.resNum[1]))
new_df = pd.merge(df,hagaiID_df,on=['pdbID','resNum','chain'],how='inner')
print(new_df)

new_df.to_csv('/home/kenneth/software/findgeo/findgeo_hagai_combined.csv',index=False)
