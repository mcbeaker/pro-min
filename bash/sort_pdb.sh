#!/bin/bash

wd='home/kenneth/Box/proj/proXtal/data/proteins/findGeoResults'
proFile='/Users/ken/Box/proj/proXtal/FindGeoSummative_wo_CA_K_NA_ZN.csv'

IFS=''
while read line; do

    pdb=`cut $line -d"." -f1`
    
done < $proFile
