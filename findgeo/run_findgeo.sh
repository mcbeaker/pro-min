#!/bin/bash
cd '/Users/ken/Box/proj/proXtal/data/minerals/amcsd_final/pdb'
pdbs='/Users/ken/Box/proj/proXtal/data/minerals/amcsd_final/pdb/pdbName.csv'


IFS=''
while read line; do
    # echo $line
    mkdir $PWD/${line%.pdb}
    cp $PWD/$line $PWD/${line%.pdb}
    cd $PWD/${line%.pdb}
    findgeo.py -w $PWD -p $line -e H -o
    cd .. 
done < $pdbs