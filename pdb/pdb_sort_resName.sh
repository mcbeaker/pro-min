#!/bin/bash


# findGeoSummativeFile=$1
# line=$1
# echo $line
basePath='/home/kenneth/proj/proMin/proteins/rcsb/pdbs_0_1.5/findgeo'
combineFindGeoFolder=$basePath/combineFindGeoResults

#make output folder if it doesnt exist
# if [[ ! -d $combineFindGeoFolder ]]; then
    # mkdir $combineFindGeoFolder
# fi

# IFS=''
for i in $(ls $combineFindGeoFolder/*.pdb); do
    cp $i $i.bak
     pdb_sort $i > $i.pdb
     rm $i.bak
     mv $i.pdb $i
done
    #format of line
    #101m.A_FE_155_FE_FE,1270,spy,square pyramid (regular),FE_155_1270_A
    #echo $line
#     line=`echo $line | tr -d '\n'`
#     # echo $line
#     pdb=`echo $line | cut -d'.' -f 1`
#     geo=`echo $line | cut -d',' -f 3`
#     metalID=`echo $line | cut -d',' -f 5`
#     metalCheck=`echo $line | cut -d',' -f 5| cut -d'_' -f 1|sed 's/[0-9]\+$//'`
#     # echo $metalCheck
#     metals=("FE" "CO" "MN" "CU" "NI" "MO" "W" "V")
#     # echo $pdb
#     # echo $geo
#     # echo $metalID
#     # echo $metalID
#     # if [[ "$metal" == "FE" ]] || [[ ! "$metal" == *"CO"* ]] && [[ ! "$metal" == *"CA_"* ]] && [[ ! "$metal" == *"NA_"* ]]; then
#     if [[ " ${metals[@]} " =~ " ${metalCheck} " ]];then
#         # echo $metals
#         pdbPath=$basePath/$pdb/$metalID/"findgeo.input"

#         outPDB=$combineFindGeoFolder/$pdb.$metalID.$geo.pdb
#         # echo $outPDB

#         # if [[ ! -f $pdbPath ]]; then
#         #     #  echo $pdbPath
#         #      echo $line >> $combineFindGeoFolder/combinnotHavePDB.txt
#         # fi

#         # if [[ -f $pdbPath ]]; then

#         #     echo $line >> $combineFindGeoFolder/havePDB.txt
#         #     # echo $outPDB
#         # fi

#         if [[ ! -f $outPDB ]] && [[ -f $pdbPath ]]; then
#             #  pdb_sort $pdbPath > $outPDB
#             cp $pdbPath $outPDB
#         fi
#     fi

# done < $findGeoSummativeFile
