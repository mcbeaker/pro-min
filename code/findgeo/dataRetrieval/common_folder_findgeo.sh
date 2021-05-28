#!/bin/bash


basePath=$1
findGeoSummativeFile=$2
# line=$1
# echo $line
combineFindGeoFolder=$basePath/combineFindGeoResults

#make output folder if it doesnt exist
# if [[ ! -d $combineFindGeoFolder ]]; then
    # mkdir $combineFindGeoFolder
# fi

IFS=''
while read line; do
    #format of line
    #101m.A_FE_155_FE_FE,1270,spy,square pyramid (regular),FE_155_1270_A
    # echo $line
    line=`echo $line | tr -d '\n'`
    # echo $line
    pdb=`echo $line | cut -d'.' -f 1`
    pdbLower=`echo $pdb | perl -ne 'print lc'`
    geo=`echo $line | cut -d',' -f 3`
    metalID=`echo $line | cut -d',' -f 5`
    chain=`echo $metalID | cut -d'_' -f4`
    metalCheck=`echo $line | cut -d',' -f 5| cut -d'_' -f 1|sed 's/[0-9]\+$//'`
    # echo $metalCheck
    metals=("FE" "CO" "MN" "CU" "NI" "MO" "W" "V")
    # echo "${pdbLower}_${chain}"
    # echo $pdb
    # echo $geo
    # echo $metalID
    # echo $metalID
    # if [[ "$metal" == "FE" ]] || [[ ! "$metal" == *"CO"* ]] && [[ ! "$metal" == *"CA_"* ]] && [[ ! "$metal" == *"NA_"* ]]; then
    if [[ " ${metals[@]} " =~ " ${metalCheck} " ]];then
        # echo $metals
        pdbPath=$basePath/$pdb/$metalID/"findgeo.input"
        # echo $pdbPath
        # exit
        outPDB=$combineFindGeoFolder/$pdb.$metalID.$geo.pdb
        # echo $outPDB

        # if [[ ! -f $pdbPath ]]; then
        #     #  echo $pdbPath
        #      echo $line >> $combineFindGeoFolder/combinnotHavePDB.txt
        # fi

        # if [[ -f $pdbPath ]]; then

        #     echo $line >> $combineFindGeoFolder/havePDB.txt
        #     # echo $outPDB
        # fi

        if [[ ! -f $outPDB ]] && [[ -f $pdbPath ]]; then
            # echo 'blah'
            #  pdb_sort $pdbPath > $outPDB
            cp $pdbPath $outPDB
        fi
        # exit
    fi

done < $findGeoSummativeFile
