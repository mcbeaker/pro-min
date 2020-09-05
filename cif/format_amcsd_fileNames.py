from progress.bar import Bar
from CifFile import ReadCif
import os
import shutil
import glob

d = "/home/kenneth/proj/proMin/minerals/database/data/xx"

files = glob.glob(os.path.join(d,"xx*"))

shaunna_list = []

#get shaunna's list of minerals
with open(os.path.join(d,'mineralNames.csv')) as f:
    next(f)
    for line in f:
        shaunna_list.append(line.rstrip())

with Bar('Processing',max=len(files)) as bar:
    for f in files:
        try:
            cf = ReadCif(f)

            #print(f)
            if '_chemical_name_mineral' in cf['global']: 

                if cf['global']['_chemical_name_mineral'] in shaunna_list:

                    code = str(cf['global']['_database_code_amcsd']).lstrip("0")
                    mineral = cf['global']['_chemical_name_mineral']
                    
                    if os.path.exists(os.path.join(d,mineral)) == False:
                        os.mkdir(os.path.join(d,mineral))

                    outName = os.path.join(d,mineral,mineral+"_"+code+".cif")
                    
                    shutil.copy(f,outName)
                    shutil.move(f,os.path.join(d,'csplit'))
                # names.append(cf['global']['_chemical_name_mineral'])
                # print(f)
        except Exception as e:
            print(f)
            print(e)
            continue
        bar.next()
    bar.finish()