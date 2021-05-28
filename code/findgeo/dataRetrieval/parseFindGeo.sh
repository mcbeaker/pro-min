#!/bin/bash
pdb=`echo "$1" | cut -d'.' -f1`

while read i;
do 
	id=`echo "$i" | cut -d':' -f1` 
	metal=`echo "$id" | cut -d'_' -f1`
	residue=`echo "$id" | cut -d'_' -f2`
	atom=`echo "$id" | cut -d'_' -f3`
	chain=`echo "$id" | cut -d'_' -f4`
	
	geometry=`echo "$i" | cut -d':' -f2`; 
	geoAbbr=`echo "$geometry" | cut -d'-' -f1 | tr -d '[:space:]'`
	geoFull=`echo "$geometry" | cut -d'-' -f2 | xargs`


	# echo $i
	echo "${pdb}.${chain}_${metal}_${residue}_${metal}_${metal},${atom},${geoAbbr},${geoFull},${id}"
done < "$pdb/FindGeo.summary"

#MG_3102_81622_X: irr - Irregular geometry
