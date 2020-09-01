
import pandas as pd
import re
from pymol import cmd
from progress.bar import Bar
import __main__
# __main__.pymol_argv = ['pymol','-qci'] # Pymol: quiet and no GUI
import pymol
from time import sleep
import os

class SlowBar(Bar):
    suffix = '%(index)d/%(max)d - %(remaining_hours)d hours remaining'
    @property
    def remaining_hours(self):
        return self.eta // 3600




uniq_envFile = '/home/kenneth/proj/proXtal/proteins/hagai/data/uniq_environments.txt'
envFile = '/home/kenneth/proj/proXtal/proteins/hagai/data/df_micro.csv'

pdbs = '/home/kenneth/proj/proXtal/proteins/hagai/pdb/'

uniq_env = pd.read_csv(uniq_envFile,index_col=None,header=None)

env = pd.read_csv(envFile,index_col=0,header=0,na_filter=False)
cmd.reinitialize()

bindingMotif = ''

bar = SlowBar("Processing",max=len(uniq_env.index),fill='*') #,suffix='%(percent).1f%% - %(eta)ds')

for row,series in uniq_env.iterrows():
    pdb,ligChain,ligName,resNum,ligAtom,element = re.split('[._]',series[0])
    cmd.load(pdbs+pdb+".cif")
    
    neighs = env.loc[env.metalID == series[0],]
    metal = pdb+'//'+ligChain+'/'+resNum+'/'+ligAtom
    selection = '('+pdb+'//'+ligChain+'/'+resNum+'/'+ligAtom
    
    fileName = re.sub("//","_",metal)
    fileName = re.sub("/","_",fileName)
    
    if os.path.exists(fileName):
        continue
    
    else:
        for row, envi in neighs.iterrows():
            # print(type(series['NresNum']))
            # try:
            # print(str(envi['NresNum']))
            neigh = pdb+'//'+ligChain+'/'+str(envi['NresNum'])+'/'+envi['NatomName']
            selection += ' + ' + neigh

            cmd.distance(ligAtom+'_'+envi['NatomName'],metal,neigh)
            
                # print(neigh)
            
        # print(selection)
        selection += ")"
                # print(str(envi['NresNum']))
        cmd.select("envi",selection)
        cmd.zoom("envi")
        cmd.set("sphere_scale",0.1,"envi")
        cmd.show(representation="sticks",selection="envi")
        cmd.zoom("envi")
        
        cmd.save(pdbs+'pse/'+fileName+'.pse')
        cmd.reinitialize()
        sleep(0.1)
    bar.next()
bar.finish()
