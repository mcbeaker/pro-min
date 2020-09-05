#!/usr/bin/env python
from bioservices.uniprot import UniProt #imports uniprot 
import pandas as pd 
import os
from progress.bar import Bar
hdir = '/home/kenneth/proj/proMin/proteins/hagai/data'

hagaiNames = ["Microenvironment_ID","Binding_motif","Cofactor","Metal","cofactor_group","EC","Head","Molecule","Organism_scientific","no_rank","superkingdom","phylum","class_","order","family","genus","organism_taxid","name","chains","Resolution","Structure_method","Keywords","Journal_reference","Release_date"]
dfHag = pd.read_csv(os.path.join(hdir,'hagai.csv'),header=0,index_col=False)
dfHag[['pdb','LIG','chain','resID','function']] = dfHag.Microenvironment_ID.str.split('[._]',expand=True)
dfHag['resID'] = dfHag.resID.astype(int)
dfHag['Binding_motif'] = dfHag.Binding_motif.astype(int)
dfHag['hagai'] = 'yes'
import numpy as np
dfHag['ec'] = np.nan
u = UniProt(verbose=True)

bar = Bar("Processing",max=len(dfHag.index),fill='*',suffix='%(percent).1f%% - %(eta)ds')

map = u.mapping('PDB_ID','ACC',query=dfHag.pdb)
# print(map)
# exit()
df = u.get_df([id[0] for id in map.values()]) 
df.to_csv(os.path.join(hdir,'df_microPDBs.csv')) #returns dataframe with 
# Unnamed: 0
# Entry
# Entry name Gene names Gene names  (primary ) Gene names  (synonym ) Gene names  (ordered locus )
# Gene names  (ORF ) Organism Organism ID
# Protein names Proteomes Taxonomic lineage (ALL) Taxonomic lineage IDs Virus hosts Sequence Length Mass Gene encoded by Alternative products (isoforms) 
# Erroneous gene model prediction Erroneous initiation Erroneous termination Erroneous translation Frameshift Mass spectrometry Polymorphism RNA editing Sequence caution 
# Alternative sequence Natural variant Non-adjacent residues Non-standard residue Non-terminal residue Sequence conflict Sequence uncertainty 
# Version (sequence) Domains Domain count Domain [CC] Sequence similarities Coiled coil Compositional bias Domain [FT] Motif Region Repeat Zinc finger 
# EC number Absorption Catalytic activity Cofactor General annotation (ENZYME REGULATION) Function [CC] Kinetics Pathway Redox potential 
# Temperature dependence pH dependence Active site Binding site DNA binding Metal binding Nucleotide binding Site Gene ontology (GO) Gene ontology (biological process) Gene ontology (molecular function) 
# Gene ontology (cellular component) Gene ontology IDs InterPro Interacts with Subunit structure [CC] 
# PubMed ID Mapped PubMed ID Date of creation Date of last modification Date of last sequence modification Version (entry) 
# 3D Beta strand Helix Turn Subcellular location [CC] Intramembrane Topological domain Transmembrane Annotation Features Caution Tissue specificity 
# Miscellaneous [CC] Keywords Protein existence Status Sequence annotation (Features) Protein families Version Comments Cross-reference (null) 
# Keyword ID Pathway.1 Allergenic properties Biotechnological use Disruption phenotype Involvement in disease Pharmaceutical use Toxic dose 
# Post-translational modification Chain Cross-link Disulfide bond Glycosylation Initiator methionine Lipidation Modified residue Peptide 
# Propeptide Signal peptide Transit peptide Taxonomic lineage (all) Taxonomic lineage (SUPERKINGDOM) Taxonomic lineage (KINGDOM) 
# Taxonomic lineage (SUBKINGDOM) Taxonomic lineage (SUPERPHYLUM) Taxonomic lineage (PHYLUM) Taxonomic lineage (SUBPHYLUM) 
# Taxonomic lineage (SUPERCLASS) Taxonomic lineage (CLASS) Taxonomic lineage (SUBCLASS) Taxonomic lineage (INFRACLASS) 
# Taxonomic lineage (SUPERORDER) Taxonomic lineage (ORDER) Taxonomic lineage (SUBORDER) Taxonomic lineage (INFRAORDER) Taxonomic lineage (PARVORDER) 
# Taxonomic lineage (SUPERFAMILY) Taxonomic lineage (FAMILY) Taxonomic lineage (SUBFAMILY) Taxonomic lineage (TRIBE) Taxonomic lineage (SUBTRIBE) 
# Taxonomic lineage (GENUS) Taxonomic lineage (SUBGENUS) Taxonomic lineage (SPECIES GROUP) Taxonomic lineage (SPECIES SUBGROUP) Taxonomic lineage (SPECIES) 
# Taxonomic lineage (SUBSPECIES) Taxonomic lineage (VARIETAS) Taxonomic lineage (FORMA) Taxonomic lineage IDs (all) Taxonomic lineage IDs (SUPERKINGDOM) 
# Taxonomic lineage IDs (KINGDOM) Taxonomic lineage IDs (SUBKINGDOM) Taxonomic lineage IDs (SUPERPHYLUM) Taxonomic lineage IDs (PHYLUM) 
# Taxonomic lineage IDs (SUBPHYLUM) Taxonomic lineage IDs (SUPERCLASS) Taxonomic lineage IDs (CLASS) Taxonomic lineage IDs (SUBCLASS) 
# Taxonomic lineage IDs (INFRACLASS) Taxonomic lineage IDs (SUPERORDER) Taxonomic lineage IDs (ORDER) Taxonomic lineage IDs (SUBORDER) Taxonomic lineage IDs (INFRAORDER) 
# Taxonomic lineage IDs (PARVORDER) Taxonomic lineage IDs (SUPERFAMILY) Taxonomic lineage IDs (FAMILY) Taxonomic lineage IDs (SUBFAMILY) Taxonomic lineage IDs (TRIBE) 
# Taxonomic lineage IDs (SUBTRIBE) Taxonomic lineage IDs (GENUS) Taxonomic lineage IDs (SUBGENUS) Taxonomic lineage IDs (SPECIES GROUP) Taxonomic lineage IDs (SPECIES SUBGROUP) 
# Taxonomic lineage IDs (SPECIES) 
# Taxonomic lineage IDs (SUBSPECIES) Taxonomic lineage IDs (VARIETAS) Taxonomic lineage IDs (FORMA) Cross-reference (db_abbrev) Cross-reference (EMBL)
print(df)
exit()

for i, row in dfHag.iterrows():
    map = u.mapping('PDB_ID','ACC',query=row.pdb)
    # print(map[row.pdb][0])
    
    ec = u.search(map[row.pdb][0], columns="ec")
    # print(dir(ec))
    # print(type(ec))
    if ec.replace('\n','') != 'EC number':
        df = u.get_df(map[row.pdb][0])
        print(df['EC number'])
        # for c in df.columns:
        #     print(c+'\n')
        exit()
        ec = ec.replace('\n\n','\n').replace('EC number\n','')

        print(repr(ec))
        dfHag.iloc[i,dfHag.columns.get_loc('ec')] = [ec.replace(' ','').split(';')]
        xit()
    bar.next()
bar.finish()
print(dfHag.ec)
    # exit()
#     # print(row.pdbID)
 # Returns sequence on the ZAP70_HUMAN accession Id
# sequence = u.search("ZAP70_HUMAN", columns="sequence")
# print(u.search(map.values()))
    # if(len(map[row.pdbID])>1):

    #     for ac in map[row.pdbID]:
    #         print(map)
    #         print(ac)
    #         fasta = u.get_fasta_sequence(ac)
    #         print(fasta)

exit()