#!/bin/bash

f='/home/kenneth/software/findgeo/FindGeoSummative.csv'
geoDir='/home/kenneth/software/findgeo/proGeo_pdb'

IFS=''
while read line; do
    #format of line
    #101m.A_FE_155_FE_FE,1270,spy,square pyramid (regular),FE_155_1270_A
    #echo $line
    pdb=`echo $line | cut -d'.' -f 1`
    geo=`echo $line | cut -d',' -f 3`
    metal=`echo $line| cut -d',' -f 5`
    
    if [ ! -d "/home/kenneth/software/findgeo/getPro_pdb/${pdb}" ]; then
        mkdir "/home/kenneth/software/findgeo/getPro_pdb/${pdb}"
    fi

    if [ ! -d "/home/kenneth/software/findgeo/getPro_pdb/${pdb}/${metal}" ]; then
        mkdir "/home/kenneth/software/findgeo/getPro_pdb/${pdb}/${metal}"
    fi
     
    if [ $geo == "irr" ]; then
        echo $line
    else
    #echo $geo
        cp "/home/kenneth/software/findgeo/${pdb}/${metal}/${geo}.pdb" ${geoDir} 
    fi

done < '/home/kenneth/software/findgeo/FindGeoSummative.csv'
