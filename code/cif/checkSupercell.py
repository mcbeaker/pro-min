from progress.bar import Bar
from CifFile import ReadCif
import PyCif
import os
from shutil import copy2 #dest is directory
import glob
import subprocess
import sys

d = "/home/kenneth/proj/proMin/minerals/database/reformat"
outDir = "/home/kenneth/proj/proMin/minerals/supercell" 
error = os.path.join(outDir,'analysis','err.txt')
out = os.path.join(outDir,'analysis','out.txt')

shaunna_list = []


# log.write('Feasible,SuperCellSize,File,Combinations\n')
# log.flush()  # <-- here's something not to forget!
#get shaunna's list of minerals
# with open(error) as err:

with open(os.path.join(d,'dir.txt')) as f:
    direcs = [line.rstrip() for line in f.readlines()]
    with Bar('Processing',max=len(direcs)) as bar:
        count = 0
        for drc in direcs:
            #cd into new folder
            minDir = os.path.join(outDir,drc)
            os.chdir(minDir)
            cifFiles = glob.glob(os.path.join(minDir,"*.cif"))
            for cifFile in cifFiles:
                #get basename
                # print(cifFile)
                cifName=os.path.splitext(os.path.basename(cifFile))[0]
                dirMinCif = os.path.join(minDir,cifName)
                if os.path.isdir(dirMinCif) == False:
                    os.mkdir(dirMinCif)
                
                copy2(cifFile, dirMinCif)
                superCellCif = os.path.join(dirMinCif,cifName+'.cif')
                #change into individual mineral cif folders
                os.chdir(dirMinCif)
                #run supercell.sh
                try:
                    cf = ReadCif(superCellCif)
                    # c = subprocess.Popen(['supercellCall.sh',superCellCif,'2x2x2',out],stderr=error)
                except Exception as e:
                    print(superCellCif)
                    print(e)
                #     sys.stdout.write(e)
                #     sys.stdout.write(pdbPath) 
                # for line in c.stderr:
                #     # print(line)
                #     # print(str(line))
                #     sys.stdout.write(str(line))
                # c.wait()

                # c = subprocess.Popen(['supercellCall.sh',superCellCif,'2x2x2'], stdout=log, stderr=subprocess.PIPE)
            #change back to previous directory
                
                os.chdir(minDir)
            #change back into original directory
            os.chdir(d)
            bar.next()
        bar.finish()
                # count += 1
                # # print(dirMinCif+"\n")
                # if count == 10:
                #     exit()
