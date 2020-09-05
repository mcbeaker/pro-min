#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib
matplotlib.use('Agg')
import seaborn as sb
import numpy as np
import scipy


def load_findgeo_hagai2018():
    f = "/home/kenneth/proj/proMin/proteins/hagai/2018/pdbs/findgeo/envComp"
    fileName = os.path.join(f,os.path.join(f,"findgeoID.csv"))
    df = pd.read_csv(fileName,header=0,index_col=0)
    # print(df.columns)
    return df 

def load_hagai2018_findgeo_metal_coord():
    f = "/home/kenneth/proj/proMin/proteins/hagai/2018/pdbs/findgeo/envComp"
    fileName = os.path.join(f,"pdb_atomID_resNum_atomNum_chain_geo.txt")
    df = pd.read_csv(fileName,header=0,index_col=None)
    # print(df.columns)
    df[['pdb','atomName','resNum','atomNum','chain','geo','ext']] = df.coord.str.split('\.|_',expand=True)
    df.to_csv(os.path.join(f,"findgeoID.csv"))
    return df

def load_hagai2018():
    f = "/home/kenneth/proj/proMin/proteins/hagai/2018/data/hagai.csv"
    df = pd.read_csv(f,index_col=None)
    # print(df.columns)
    df[['pdb','resName','chain','resNum','function']] = df['Microenvironment_ID'].str.split('\.|_',expand=True)
    df = df.astype({'resNum':'int64'})
    df['hagai'] = 'yes'
    #columns
    # ,Microenvironment_ID,Binding_motif,Cofactor,Metal,cofactor_group,Head,Molecule,Organism_scientific,
    # no_rank,superkingdom,phylum,class_,order,family,genus,organism_taxid,name,chains,Resolution,Structure_method,Keywords,
    # Journal_reference,Release_date,pdb,resName,chain,resNum,function,hagai
    return df

def load_dfAll_NR():
    rdir = '/home/kenneth/proj/proMin/results'
    dfAll = pd.read_csv(os.path.join(rdir,'nr_40','nr_40_pro_min.csv'),header=0,index_col=0)
    return dfAll
    
def load_dfAll():
    rdir = '/home/kenneth/proj/proMin/results'
    dfAll = pd.read_csv(os.path.join(rdir,'ec_min_pro_fg_dictionaryValues.csv'),header=0,index_col=False)
    return dfAll

def autolabel(rects,ax):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

def merge_findgeo_hagai_df():
    f = "/home/kenneth/proj/proMin/proteins/hagai/2018/pdbs/findgeo/envComp"
    dfHag = load_hagai2018()
    # print(dfHag.resNum.dtype)
    dfFindGeo = load_findgeo_hagai2018()
    # print(dfFindGeo.resNum.dtype)
    df3 = pd.merge(dfHag,dfFindGeo,how="inner",
        right_on=['pdb','resNum','chain'],
        left_on= ['pdb','resNum','chain'])
    
    # df3['coord'] = df3['resName'] + "_" + df3['coord']
    
    df3.to_csv(os.path.join(f,"hagai_findgeoID.csv"))

    # return df3
def load_findgeo_hagai_df():
    f = "/home/kenneth/proj/proMin/proteins/hagai/2018/pdbs/findgeo/envComp"
    fileName = os.path.join(f,os.path.join(f,"hagai_findgeoID.csv"))
    df = pd.read_csv(fileName,header=0,index_col=0)
    # print(df.columns)
    return df  

def main():
    # load_hagai2018_findgeo_metal_coord()
    merge_findgeo_hagai_df()

main()
