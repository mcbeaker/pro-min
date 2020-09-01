#!/bin/bash
pdb=$1

formatted=""

while read i;
do 
	id=`echo "$i" | cut -d':' -f1`; 
	metal=`echo "$id" | cut -d'_' -f1`
	residue=`echo "$id" | cut -d'_' -f2`
	atom=`echo "$id" | cut -d'_' -f3`
	chain=`echo "$id" | cut -d'_' -f4`
	
	etry=`echo "$i" | cut -d':' -f2`; 
	geoAbbr=`echo "$etry" | cut -d'-' -f1 | tr -d '[:space:]'`
	geoFull=`echo "$etry" | cut -d'-' -f2 | xargs`


	echo $i
	echo "${pdb}.${chain}_${metal}_${residue}_${metal}_${metal},${atom},${geoAbbr},${geoFull},${id}" >> "${pdb}_FindGeo.csv"
done < "$pdb/FindGeo.summary"

#MG_3102_81622_X: irr - Irregular etry
