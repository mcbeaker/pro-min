import pandas as pd

unqList = []
with open("/Users/ken/Box/proj/proXtal/data/proteins/hagai_analysis/2018/unique_modules_SPAN.csv",'r') as unq:
    unqList = [i.rstrip() for i in unq.readlines()]


df = pd.read_csv("/Users/ken/Box/proj/proXtal/data/proteins/hagai_analysis/2018/network.csv")

SPAN_mods = df[df['Binding_motif'].isin(unqList)]

outFile = "/Users/ken/Box/proj/proXtal/data/proteins/hagai_analysis/2018/microspheresIn_SPAN.csv"
SPAN_mods.to_csv(outFile)