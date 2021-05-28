import os
import pandas as pd
from CifFile import ReadCif
import numpy as np
import sys

#paths and files
base = '/home/kenneth/proj/proMin/minerals/database/data/cif/feS'

minNameFile = 'feS_mineralNames.txt'
minerals = os.path.join(base,minNameFile)

# shaunnaFolder = 'data/shaunna'
# amcsdFolder = os.path.join(base,shaunnaFolder)
#---------------------

dfMin = pd.read_csv(minerals,usecols = ['MineralName'],index_col=False)
# print(list(dfMin.mineralName))

# cf['global'].keys()
dfOUTPUT = pd.DataFrame(columns=['_chemical_name_mineral',
'_publ_author_name','_journal_name_full','_journal_volume','_journal_year','_journal_page_first','_journal_page_last','_publ_section_title',
'_database_code_amcsd','_chemical_formula_sum','_cell_length_a','_cell_length_b',
'_cell_length_c','_cell_angle_alpha','_cell_angle_beta','_cell_angle_gamma','_cell_volume','_exptl_crystal_density_diffrn','_symmetry_space_group_name_h-m',
'dupLabOcc','Multiplicity','Pressure','Temperature'])
# if atoms contain same fractxyz, then list subs with occupancy

cifKeys = ['_chemical_name_mineral','_publ_author_name','_journal_name_full','_journal_volume','_journal_year','_journal_page_first','_journal_page_last',
'_publ_section_title','_database_code_amcsd','_chemical_formula_sum','_cell_length_a','_cell_length_b','_cell_length_c','_cell_angle_alpha',
'_cell_angle_beta','_cell_angle_gamma','_cell_volume','_exptl_crystal_density_diffrn','_symmetry_space_group_name_h-m']

for i,row in dfMin.iterrows():
		
	cifFolder = os.path.join(base,row.MineralName)
	if os.path.exists(cifFolder): #row.mineralName not in ['Angelaite','Arsenopyrite']:

		for cifFile in [f for f in os.listdir(cifFolder) if f.endswith(".cif")]:
			try:
				cf = ReadCif(os.path.join(cifFolder,cifFile))

				# '_atom_site_label','_atom_site_fract_x','_atom_site_fract_y','_atom_site_fract_z','_atom_site_occupancy','_atom_site_u_iso_or_equiv'])
				cifDic = {key:cf['global'][key] for key in cifKeys}
				cifDic['Temperature'] = 'False'
				cifDic['Pressure'] = 'False'

				#remove new line from '_publ_section_title'
				cifDic['_publ_section_title'] = cifDic['_publ_section_title'].replace("\n"," ").replace('\r'," ").strip()
				# print(cifDic)
				#with multiplicity
				if '_atom_site_occupancy' in cf['global']: #keep if it doesnt have multiple occupancies
					
					#find atom labels that have the same fractxyz
					metals = ['Co','Cu','Fe','Mn','Mo','Ni','V','W','Zn']
					data = zip(cf['global']['_atom_site_label'],cf['global']['_atom_site_fract_x'],cf['global']['_atom_site_fract_y'],
							cf['global']['_atom_site_fract_z'],cf['global']['_atom_site_occupancy'])
					dfMult = pd.DataFrame(data,columns=['label','x','y','z','occ'])

					#pull out first of each duplicate
					dfSameFirst = dfMult[dfMult.duplicated(subset=['x','y','z'],keep='first')] # not all multiplicities will have same metal - & (dfMult.label.str.contains('|'.join(metals)))]

					lab_occ = []

					#pull out rest of duplicates
					for i,row in dfSameFirst.iterrows():
						rest = dfMult[(dfMult['x'] == row['x']) & (dfMult['y'] == row['y']) & (dfMult['z'] == row['z'])]
						lab_occ.append((list(rest['label']),list(rest['occ'])))

					cifDic['dupLabOcc'] = lab_occ
				else:
					cifDic['dupLabOcc'] = np.nan
				# print(type(cifDic))
				if 'T =' in cifDic['_publ_section_title']:
					cifDic['Temperature'] = 'True'
					print(cifDic)
				if 'Pressure' in cifDic['_publ_section_title']:
					cifDic['Pressure'] = 'True'
				if 'P = ' in cifDic['_publ_section_title']:
					cifDic['Pressure'] = 'True'

				if 'Gpa' in cifDic['_publ_section_title']:
					cifDic['Pressure'] = 'True'
				# 	print(cifDic['_chemical_name_mineral'])
				# 	print(cifDic['_database_code_amcsd'])
				# 	print(cifDic['_publ_section_title']+'\n')
				if cifDic['Temperature'] != 'True' and cifDic['Pressure'] != 'True':
					dfOUTPUT = dfOUTPUT.append(cifDic,ignore_index=True)
				# print(type(dfOUTPUT))
				# exit()
			except Exception as e:
				print(e)
				print("Error on line {}".format(sys.exc_info()[-1].tb_lineno))
	else:
		print(cifFolder, " does not exist")
dfOUTPUT.to_csv(os.path.join(base,'feS_amcsd_stats_boolean_Press_Temp.csv'))