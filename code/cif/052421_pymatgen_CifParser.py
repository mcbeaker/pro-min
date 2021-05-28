from itertools import groupby
from pymatgen.io.cif import CifParser 
import os
import pymatgen.symmetry  as pms
from pymatgen.core import Molecule

import pymatgen.analysis.bond_valence as pma
from pymatgen.core.periodic_table import Element, Species, get_el_sp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from monty.serialization import loadfn
import numpy as np
import pandas as pd
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt


# List of electronegative elements specified in M. O'Keefe, & N. Brese,
# JACS, 1991, 113(9), 3226-3229. doi:10.1021/ja00009a002.
ELECTRONEG = [
    Element(sym)
    for sym in [
        "H",
        "B",
        "C",
        "Si",
        "N",
        "P",
        "As",
        "Sb",
        "O",
        "S",
        "Se",
        "Te",
        "F",
        "Cl",
        "Br",
        "I",
    ]
]

module_dir = "/home/kenneth/software/pymatgen/pymatgen/analysis" 

# Read in BV parameters.
BV_PARAMS = {}
for k, v in loadfn(os.path.join(module_dir, "bvparam_1991.yaml")).items():
    BV_PARAMS[Element(k)] = v

# Read in yaml containing data-mined ICSD BV data.
all_data = loadfn(os.path.join(module_dir, "icsd_bv.yaml"))
ICSD_BV_DATAR= {Species.from_string(sp): data for sp, data in all_data["bvsum"].items()}
PRIOR_PROB = {Species.from_string(sp): data for sp, data in all_data["occurrence"].items()}


mins = ['Greigite_127','Pyrite_12728','Pyrrhotite_288','Mackinawite_14518', 'Marcasite_12726', 'Smythite_80', 'Tochilinite_15564', 'Troilite_4160']
d = '/home/kenneth/proj/proMin/minerals/database/reformat'
bv_elem = []
bv_sum = []
bv_type = []
bvs = []

for m in mins:
    print(m)
    mn,num = m.split('_')
    fi = os.path.join(d,mn,mn+'_'+str(num)+'.cif')
    # print(fi)
    # exit()
    parser = CifParser(fi)
    structure = parser.get_structures(primitive=True)[0]
    els = [Element(el.symbol) for el in structure.composition.elements]

    if not set(els).issubset(set(BV_PARAMS.keys())):
        raise ValueError("Structure contains elements not in set of BV parameters!")

    bv = pma.BVAnalyzer()
    
    # Perform symmetry determination and get sites grouped by symmetry.
    if bv.symm_tol:
        finder = SpacegroupAnalyzer(structure, bv.symm_tol)
        symm_structure = finder.get_symmetrized_structure()
        equi_sites = symm_structure.equivalent_sites
    else:
        equi_sites = [[site] for site in structure]

    # Sort the equivalent sites by decreasing electronegativity.
    equi_sites = sorted(equi_sites, key=lambda sites: -sites[0].species.average_electroneg)

    # Get a list of valences and probabilities for each symmetrically
    # distinct site.
    all_prob = []
    if structure.is_ordered:
        for sites in equi_sites:
            test_site = sites[0]
            # print(test_site)
            nn = structure.get_neighbors(test_site, 3.0)
        
            bv_sum.append(pma.calculate_bv_sum(test_site, nn))
            bv_elem.append(test_site.specie.symbol)
            bv_type.append('Mineral')

#protein
d = '/home/kenneth/proj/proMin/proteins/feS'
with open(os.path.join(d,'pdbNames.txt')) as f:
    for i in f.readlines():
        fi  = os.path.join(d,i.strip()+'.pdb')
        mol = Molecule.from_file(fi)
        # print(dir(parser))
        print('test1')
        print(list(BV_PARAMS.keys()))
        els = [Element(el.symbol) for el in mol.composition.elements if Element(el.symbol) in list(BV_PARAMS.keys())]
        print(els)
        exit()
        if not set(els).issubset(set(BV_PARAMS.keys())):
            raise ValueError("Structure contains elements not in set of BV parameters!")
        sites = mol.sites
        # Sort the equivalent sites by decreasing electronegativity.
        sites = sorted(sites, key=lambda sites: -sites[0].species.average_electroneg)
        print('test')
        for site in sites:
            test_site = site[0]
            # print(test_site)
            nn = mol.get_neighbors(test_site, 3.0)
        
            bv_sum.append(pma.calculate_bv_sum(test_site, nn))
            bv_elem.append(test_site.specie.symbol)
            bv_type.append('Protein')

        # bv = pma.BVAnalyzer()
        
        # # Perform symmetry determination and get sites grouped by symmetry.
        # if bv.symm_tol:
        #     finder = SpacegroupAnalyzer(mol, bv.symm_tol)
        #     symm_structure = finder.get_symmetrized_structure()
        #     equi_sites = symm_structure.equivalent_sites
        # else:
        #     equi_sites = [[site] for site in structure]
        # print(equi_sites)
        # exit()
# exit()
df = pd.DataFrame(np.random.randn(len(bv_sum),2),columns=['Element','Bond_Valence_Sum'])
df['Element'] = bv_elem
df['Bond_Valence_Sum'] = bv_sum
print(df.groupby('Element').describe())
# print(df['Element'])
df.boxplot(by='Element')
plt.show()


