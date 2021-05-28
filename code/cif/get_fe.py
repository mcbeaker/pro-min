from CifFile import ReadCif


import os
from shutil import copy

d = '/home/kenneth/proj/proXtal/amcsd/cif/mineral/'

for f in os.listdir(d):
    cf = ReadCif(d+f)

    labels = cf['global']['_atom_site_label']

    subs = "Fe"

    if len(list(filter(lambda x: subs in x, labels))) > 0:

        copy(f,d+'Fe/')

