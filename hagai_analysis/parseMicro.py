import ast
import pandas as pd

f = "/home/kenneth/proj/proXtal/proteins/hagai/hagai.csv"
pdir = "/home/kenneth/proj/proXtal/proteins/hagai/pdb"

df = pd.read_csv(f)
# print(df.keys)
feXrayDF = df[df['Metal'].str.contains("FE") & df['Structure_method'].str.contains("x-ray")]
pdbs = feXrayDF["Microenvironment ID"].apply(lambda x: x[0:4]).sort_values().drop_duplicates()


# 4rvy.HEM_f_101_oxidoreductase
# 3wr9.B_502_FE2_FE

print(feXrayDF)

with open(os.path.join(pdir,"micro.csv"),"w+") as micro:
    dicts = micro.readlines()

    for line in dicts: 

        dct = ast.literal_eval(line)
        print(dct.keys())


for i in feDF.itertuples():
    pdb,LIG,chain,resID,function = re.split('[._]',i.Microenvironment_ID)
