#!/bin/bash



#!/bin/bash

basePath='/home/kenneth/proj/proMin/proteins/hagai/data'
fastaFolder=$basePath/fasta

IFS=''
while read line; do
    pdb=`echo $line | cut -d, -f1`
    chain=`echo $line | cut -d, -f2`
    echo $pdb.$chain

    wget -O $fastaFolder/${pdb}.${chain}.fasta www.rcsb.org/fasta/chain/${pdb}.${chain}/download
    sleep 0.5
    
done < $basePath/2018_hag_pdb_chain.csv