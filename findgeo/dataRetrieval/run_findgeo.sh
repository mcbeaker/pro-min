#!/bin/bash
pdbDir='/home/kenneth/proj/proMin/proteins/rcsb/pdbs_0_1.5/pdbs'
# pdbs=$2
findDir='/home/kenneth/proj/proMin/proteins/rcsb/pdbs_0_1.5/findgeo'

line=$1

IFS=''
# while read line; do
    # echo $lined
    rm -rf $findDir/$line
    mkdir $findDir/$line
    cp $pdbDir/$line.pdb $findDir/$line
    cd $findDir/$line
    findgeo.py -w $findDir/$line -p $line.pdb -e H -o
    # cd ..

# done < $pdbs

# awk -F',' {'printf ("%s.%s\n", $1, $2)'} ../2018_hag_pdb_chain.csv | xargs -L 1 -P 12 -I pdbChain get_fasta_chain.sh pdbChain