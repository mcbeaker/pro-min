import pandas as pd
import os

def read_hagai_2018():
    f = "/home/kenneth/proj/proMin/proteins/hagai/data/hagai.csv"
    df = pd.read_csv(f,index_col=None)
    print(df.columns)
    df[['pdb','resName','chain','resNum','function']] = df['Microenvironment_ID'].str.split('\.|_',expand=True)
    df['hagai'] = 'yes'
    #columns
    # ,Microenvironment_ID,Binding_motif,Cofactor,Metal,cofactor_group,Head,Molecule,Organism_scientific,
    # no_rank,superkingdom,phylum,class_,order,family,genus,organism_taxid,name,chains,Resolution,Structure_method,Keywords,
    # Journal_reference,Release_date,pdb,resName,chain,resNum,function,hagai

    return df

def read_ec_findgeo_2018():
    rdir="/home/kenneth/proj/proMin/results"
    f = "/home/kenneth/proj/proMin/results/ec_min_pro_fg_dictionaryValues.csv"
    dfHag = pd.read_csv(f,index_col=0,header=0)
    cols_to_order = ['pdb', 'chain']
    new_columns = cols_to_order + (dfHag.columns.drop(cols_to_order).tolist())
    dfHag = dfHag[new_columns]
    dfHag.drop(columns=['PDB','CHAIN','LIG','resID'],inplace=True)
    dfHag.to_csv(os.path.join(rdir,"pdb_chain_data_081020.csv"))
    print(dfHag.columns)
    # 'pdb', 'chain', 'amcsdID', 'atomName', 'atomNum', 'envComp', 'gRMSD',
    #    'geo', 'mineral', 'resNum', 'type', 'valence', 'findgeo',
    #    'Microenvironment_ID', 'Binding_motif', 'Cofactor', 'Metal',
    #    'cofactor_group', 'EC', 'Head', 'Molecule', 'Organism_scientific',
    #    'no_rank', 'superkingdom', 'phylum', 'class_', 'order', 'family',
    #    'genus', 'organism_taxid', 'name', 'chains', 'Resolution',
    #    'Structure_method', 'Keywords', 'Journal_reference', 'Release_date',
    #    'function', 'hagai', 'ACCESSION', 'EC_NUMBER', 'EC1', 'EC2', 'EC3',
    #    'EC4']
    return dfHag


def merge_uniprot_pdb(df):

   f = "/home/kenneth/proj/proMin/proteins/hagai/data/pdb_chain_uniprot.csv"
   dfUni = pd.read_csv(f,comment="#",index_col=None) #ignore first line
   dfUni.rename(inplace=True,columns={"PDB": "pdb", "CHAIN": "chain"})

    #PDB,CHAIN,SP_PRIMARY,RES_BEG,RES_END,PDB_BEG,PDB_END,SP_BEG,SP_END
    #combine pdb and chain


   dfHag_Uni = pd.merge(left = dfUni, right = df, right_on = ["pdb",'chain'], left_on = ['pdb','chain'], how = 'right')

   print(df.columns)
#    print(len(dfHag_Uni.pdb.unique())) #9530
#    print(len(dfHag_Uni.SP_PRIMARY.unique())) #2636
   print(dfHag_Uni[~dfHag_Uni.SP_PRIMARY.isnull()].SP_PRIMARY) #null 965, #not null 32865
   return dfHag_Uni

def get_uniprot_df(accList):
    from bioservices.uniprot import UniProt #imports uniprot 
 
    u = UniProt(verbose=True)
    df = u.get_df(accList)
    return df

def test():
    # print(merge_uniprot_pdb_return_unique_seqs().columns)
    # print(merge_uniprot_pdb_return_unique_seqs().head(25))
    # print(merge_uniprot_pdb().SP_PRIMARY)
    ddir = '/home/kenneth/proj/proMin/proteins/hagai/data'
    dfHag = read_ec_findgeo_2018()

    # dfHag = read_hagai_2018()

    cols_to_order = ['pdb', 'chain']
    new_columns = cols_to_order + (dfHag.columns.drop(cols_to_order).tolist())
    dfHag = dfHag[new_columns]
    dfHag.to_csv(os.path.join(ddir,'2018_hag_pdb_chain.csv'),index=False)
    exit()
    print(len(dfHag.index))
    dfHag_Uni = merge_uniprot_pdb(dfHag)
    print(len(dfHag_Uni.index))
    print(dfHag_Uni.duplicated())
    exit()
    dfHag_Uni.to_csv(os.path.join(ddir,'2018_hag_uni_pdb.csv'))
    # exit()
    print(dfHag_Uni.SP_PRIMARY.isnull())
    
    #get unique SP_PRIMARY
    uniqSP_PRIMARY = dfHag_Uni.SP_PRIMARY.unique()
    dfHag_Uni_nonRed = dfHag_Uni[(dfHag_Uni['SP_PRIMARY'].isnull())]
    
    # dfHag_Uni.drop_duplicates(,inplace=True)

    print(dfHag_Uni_nonRed[dfHag_Uni_nonRed.SP_PRIMARY.isnull()])

    exit()
    print(len(dfHag[~dfHag.SP_PRIMARY.isnull()].index))
    # # print(dfHag)
    dfHag.drop_duplicates(subset=['SP_PRIMARY'],inplace=True)
    print(len(dfHag[~dfHag.SP_PRIMARY.isnull()].index))
    # dfHag.dropna(subset=['SP_PRIMARY'],inplace=True)
    #  print(dfHag.iloc[0])

    # dfUni = get_uniprot_df(list(dfHag.SP_PRIMARY))
    dfUni = pd.read_csv(os.path.join(ddir,'2018_hag_uni_pdb.csv'),index_col=0)
    print(len(dfUni[~dfUni.Sequence.isnull()].index))

    dfUni_Hag = pd.merge(left=dfHag,right=dfUni,left_on=['SP_PRIMARY'],right_on=['Entry'],how='outer')
    print(len(dfUni_Hag.index)) 
    # dfUni_Hag.to_csv(os.path.join(ddir,'2018_hag_uni_pdb_sequence.csv'))
    # dfUni_Hag = pd.read_csv(os.path.join(ddir,'2018_hag_uni_pdb_sequence.csv'))

    # from Bio import SeqIO
    # from Bio.Seq import Seq
    # from Bio.SeqRecord import SeqRecord
    # for i,row in dfUni_Hag.iterrows():
        
    #     simple_seq = Seq(str(row.Sequence))
    #     simple_seq_r = SeqRecord(simple_seq,id=str(row.Entry))

    # # print(dfUni_Hag['Sequence similarities'])
    #     with open(os.path.join(ddir,"fasta",str(row.pdb)+"_"+str(row.chain)+"_"+str(row.Entry)+".fasta"), "w") as output_handle:
    #         SeqIO.write(simple_seq_r, output_handle, "fasta")
    # print(dfUni.sequence)
    # print(dfUni.columns)
    # dfHag_Uni = pd.merge(right = dfUni, left = dfHag, right_on = ['','chain'], left_on = ['pdb','chain'], how = 'outer')

    #remove duplicates and merge NaN 
    # df[(~df.duplicated()) | (df['col'].isnull())] # keep duplicates
    # df.drop_duplicates(subset=['pdb','SP_PRIMARY'],keep="first",inplace=True)
    # dfUni = get_uniprot_df(list(df.SP_PRIMARY))
    # print(dfUni)

def main():
    test()

main()
