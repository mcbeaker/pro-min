unzip cif_archive.zip 
csplit cifdata.txt /END/ {*}
sed -i.bak '/END/d' xx*
rm *.bak
sed -i.bak "/_chemical_formula_sum ''/d" xx*
rm *.bak
mkdir reformat
mv xx* reformat/
