#!/bin/bash
pdbDir='/home/kenneth/proj/proMin/proteins/rcsb/pdbs_0_1.5/pdbs'
pdbs=$1
# echo $pdb
findDir='/home/kenneth/proj/proMin/proteins/rcsb/pdbs_0_1.5/findgeo_4ang'

# line=$1

IFS=''
while read pdb; do
    pdb=`echo $pdb | cut -f1 -d' '`
    echo $pdb
    rm -rf $findDir/$pdb
    mkdir $findDir/$pdb
    cp $pdbDir/$pdb.pdb $findDir/$pdb
    cd $findDir/$pdb
    findgeo.py -w $findDir/$pdb -p $pdb.pdb -e H -o -t 4 #threshold 4 
    # cd ..

done < $pdbs

# awk -F',' {'printf ("%s.%s\n", $1, $2)'} ../2018_hag_pdb_chain.csv | xargs -L 1 -P 12 -I pdbChain get_fasta_chain.sh pdbChain