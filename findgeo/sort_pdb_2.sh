#!/bin/bash


wd='/Users/ken/Box/proj/proXtal/data/proteins/findGeoResults'
proFile='/Users/ken/Box/proj/proXtal/FindGeoSummative_wo_CA_K_NA_ZN.csv'

IFS=''
while read line; do
    #format of line
    #101m.A_FE_155_FE_FE,1270,spy,square pyramid (regular),FE_155_1270_A
    #echo $line
    line=${line%$'\n'}
    pdb=`echo $line | cut -d'.' -f 1`
    geo=`echo $line | cut -d',' -f 3`
    metal=`echo $line| cut -d',' -f 5`
    
    if [[ ! "$metal" == *"ZN_"* ]] && [[ ! "$metal" == *"K_"* ]] && [[ ! "$metal" == *"CA_"* ]] && [[ ! "$metal" == *"NA_"* ]]; then
        pdbPath=$wd/$pdb/$metal/"findgeo.input"

        outPDB=$wd/"combFindGeoPDB"/$pdb.$metal.$geo.pdb
        # echo $outPDB

        if [[ ! -f $pdbPath ]]; then
             echo $line >> $wd/combFindGeoPDB/notHavePDB.txt
        fi

        if [[ -f $pdbPath ]]; then
            echo $pdb >> $wd/combFindGeoPDB/havePDB.txt
            # echo $outPDB
        fi

        # if [[ ! -f $outPDB ]] && [[ -f $pdbPath ]]; then
        #     pdb_sort $pdbPath > $outPDB
        # fi
    fi

    

done < $proFile
