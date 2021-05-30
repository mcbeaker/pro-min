from itertools import groupby
from pymatgen.io.cif import CifParser 
from pymatgen.io.xyz import XYZ
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
from glob import glob


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

module_dir = "/home/kenneth/.pyenv/versions/3.9.5/lib/python3.9/site-packages/pymatgen/analysis" 

# Read in BV parameters.
BV_PARAMS = {}
for k, v in loadfn(os.path.join(module_dir, "bvparam_1991.yaml")).items():
    BV_PARAMS[Element(k)] = v

# Read in yaml containing data-mined ICSD BV data.
all_data = loadfn(os.path.join(module_dir, "icsd_bv.yaml"))
ICSD_BV_DATAR= {Species.from_string(sp): data for sp, data in all_data["bvsum"].items()}
PRIOR_PROB = {Species.from_string(sp): data for sp, data in all_data["occurrence"].items()}


mins = ['Greigite_127','Pyrite_12728','Pyrrhotite_288','Mackinawite_14518', 'Marcasite_12726', 'Smythite_80', 'Tochilinite_15564', 'Troilite_4160']
d = '/home/kenneth/proj/pro-min/minerals'
bv_elem = []
bv_sum = []
bv_type = []
bvs = []
bv_name = []
bv_geo = []

for m in mins:
    print(m)
    mn,num = m.split('_')
    fi = os.path.join(d,mn+'_'+str(num)+'.cif')
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
            if test_site.specie.symbol in ['Fe']:
            # print(test_site)
                nn = structure.get_neighbors(test_site, 3.0)
            
                bv_sum.append(round(pma.calculate_bv_sum(test_site, nn),2))
                bv_elem.append(test_site.specie.symbol)
                bv_type.append('Mineral')
                bv_name.append(m)
# print(bv_sum)
# exit()
#protein

# d = '/home/kenneth/proj/pro-min/proteins/feS/clusters/findGeo/combineFindGeoResults'
d='/home/kenneth/proj/pro-min/proteins/feS/sites'
xyzs = glob(os.path.join(d,'*.xyz'))
for xyz in xyzs:
    print(xyz)
    if xyz != os.path.join(d,'pdb.xyz'):
        
        xyzObj = XYZ.from_file(xyz)
        # print(dir(xyzObj))
        mol = xyzObj.molecule
        
        for site in mol.sites:
            print(site)
            import re
            f = os.path.splitext(os.path.split(xyz)[1])[0]
            # print(re.split('[_\.]',f))
            # exit()
            pdb,resNum,elName,resID,atomID,chain,geo = re.split('[_\.]',f)
            
            # print(pdb,resNum,elName,resID,atomID,chain,geo)
            if (geo != 'irr' and site.specie.symbol == 'Fe'):
                nn = mol.get_neighbors(site, 3.0)
                print(nn)
                bv_sum.append(round(pma.calculate_bv_sum(site, nn),2))
                bv_elem.append(site.specie.symbol)
                bv_type.append('Protein')
                bv_name.append(f.upper())
  

df = pd.DataFrame(np.random.randn(len(bv_sum),4),columns=['Element','Bond_Valence_Sum','Type','Name'])
df['Element'] = bv_elem
df['Bond_Valence_Sum'] = bv_sum
df['Type'] = bv_type
df['Name'] = bv_name

out = '/home/kenneth/proj/pro-min/results'
df.to_csv(os.path.join(out,'053021_vp_min_pro_val.csv'),index=False,columns=['Name','Type','Bond_Valence_Sum','Element'])
# exit()
# print(df.Element)
# exit()
import seaborn as sb

# print(df.groupby('Element').describe())
# print(df['Element'])
print(df.groupby(['Type','Element']).describe()) #.reset_index().pivot(index='Element', values='Bond', columns='level_1')
# exit()
# el = ['Fe','S']
# df = df.loc[(df['Element'].isin(el)) & (df['Bond_Valence_Sum'] > -3)]
# print(df.groupby(['Type','Element']).describe()) #.reset_index().pivot(index='Element', values='Bond', columns='level_1')
# print(df.describe(include='all'))
# exit()
# print(df)
# exit()
print(df[df.Type=='Protein'])
vp = sb.violinplot(data=df, x='Element', y = 'Bond_Valence_Sum', hue='Type',split=True)
vp.set(ylabel='Bond Valence Sum')
plt.tight_layout()
plt.savefig(os.path.join(out,'053021-vp_min_pro_val.png'),dpi=400)
# plt.show()


