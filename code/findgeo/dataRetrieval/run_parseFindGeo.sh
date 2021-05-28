#!/bin/bash
pdbsFile=$1
output='FindGeoSummative_noHeader.csv'
formatted=""

while read i;
do 
	parseFindGeo.sh $i >> $output
	# parseFindGeo.sh $pdbsFile >> $output
done < $pdbsFile

#MG_3102_81622_X: irr - Irregular geometry
