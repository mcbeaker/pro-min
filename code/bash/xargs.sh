cat FindGeoSummative_noHeader.csv | while read -r i; do printf "%q\n";done | xargs -L 1 -P 12 -I line common_folder_findgeo.sh line