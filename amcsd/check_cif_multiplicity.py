from scipy import ndimage
import numpy as np
from CifFile import ReadCif
from collections import Counter
from diffpy import structure 
from diffpy.structure.expansion import supercell
import pandas as pd
import os
import shutil
import multiprocessing as mp
from tqdm import tqdm
import string

def closest_node(node, nodes):
    nodes = np.asarray(nodes)
    # print(np.asarray(nodes))
    deltas = nodes - node
    # print(len(deltas))
    # print(deltas[0])
    dist_2 = 0
    try: 
        dist_2 = np.einsum('ij,ij->i', deltas, deltas)
    except Exception as e:
        print(str(e) + "this is in the function")
    
    # print(len(dist_2))
    return np.argmin(dist_2)

    # print(help(Structure))
global data_out
data_out = []

global finalDir
finalDir = '/home/kenneth/proj/proXtal/amcsd/database/formatName/findgeo'

global finalData
homeDir = '/home/kenneth/proj/proXtal/amcsd/database/formatName'
finalData = [os.path.join(homeDir,line[:-1]) for line in open(os.path.join(homeDir,'formatName.txt'),'r').readlines()]
# print(finalData)

# global df_out
# df_out = [] #pd.DataFrame(columns=["mineral","element","x","y","z","id"])

#metal information
global metals
global atomTypeSub
metals = ["CO","CU","FE","MN","MO","NI","V","W","ZN"]
atomTypeSub = {"ALL":"AL","ASL":"AS","BAL":"BL","BIA":"BA",\
                "BIB":"BB","BIM":"BM","HHL":"BL","HHR":"HR",\
                "HWL":"HW","HWR":"HR","PBB":"PB","PBL":"PL",\
                "PBM":"PM","SBL":"SL","TIL":"TL","ZNI":"ZN"}

def f_init(q):
    process_queue.q = q

def process_queue(file_location):
    # print(fileName)
    data = ""
    try:
        file_name = os.path.basename(file_location)
        # print(file_name)
        cf = structure.PDFFitStructure().read(file_location,"cif")
        print(cf)
        # cf = structure.loadStructure(file_location,fmt='cif')
        # print(dir(stru))#.filename+str(cf.stru.label))
        # cf.stru.assignUniqueLabels()

        # # print(cf.stru.label)

        # # print(cf[1].__sizeof__())
        # # print(cf.label)
        # for i in cf.stru:
        #     # print(i.label)
        #     i.__setattr__('label',\
        #         i.label.upper().replace('_',"").replace("-","").replace("(","").replace(")","").replace("WAT","OW").replace(",","").replace(".",""))
        #     i.__setattr__('element',i.element.upper())
        #     if i.label[0:4] == 'ZNII':
        #         i.__setattr__('label',i.label.replace('ZNII','ZN'))
        #     elif i.label[0:3] in atomTypeSub.keys():
        #         i.__setattr__('label',i.label.replace(i.label[0:3],atomTypeSub[i.label[0:3]]))
        # # print("2"+cf.filename+str(cf.stru.label))

        # countElements = Counter(cf.stru.element)
        # # # print(countElements)
        # # print(cf.stru.label)

        # strCountElements = ""

        # # sortCountElement = sorted(countElements.keys(),reverse=True)
        # metalsMineral =  [(el,countElements[el]) for el in set(countElements.keys()).intersection(metals)]
        # countMetals = [c[1] for c in metalsMineral]
        # size = max(max(countMetals),len(set(countElements.keys()).intersection(metals)))
        # # # print(countMetals)
        # # # print(line[:-1])
        # # # print(metalsMineral)
        # # print(size)
        # size = 2

        # # for m in set(countElements.keys()).intersection(metals):
        # #     strCountElements += m+str(countElements[m])

        # #     # if len(set(countElements.keys()).intersection(metals)) > 1:
        # #         # print(pdbID,set(countElements.keys()).intersection(metals))
        
        # # #     #string rest of elements
        # # for m in set(countElements.keys()).difference(metals):
        # #     strCountElements += m+str(countElements[m])

        # # # print(strCountElements)

        # # # print(cfLabels)
        # transitionMetalsMineral = set()
        # for s in cf.stru.label:
        #     if any(substring in s for substring in metals):
        #         if s == 'W':
        #             transitionMetalsMineral.add(s)
        #         elif 'W' not in s:
        #             transitionMetalsMineral.add(s)
        
        # # # print(transitionMetalsMineral)
        # # #-------------------------------------------------------------------------#
        # # #generate supercell
        # # cf.stru = supercell(cf.stru,(size,size,size))
        # cf.stru = supercell(cf.stru,(size,size,size))
        # # print("3"+cf.filename+str(cf.stru.label))
        # # # print(scf)
        # # cf.stru.assignUniqueLabels()
        # # print("4"+cf.filename+str(cf.stru.label))
        # outDir = os.path.join(finalDir,file_name.strip('.cif'))
        # if os.path.exists(outDir == False):
        #     os.mkdir(outDir)

        # outName = os.path.join(outDir,file_name.strip('.cif')+'.pdb')
        # if os.path.exists(outName) == False:
        #     # print(finalDir)
        #     # print(outName) 
        #     cf.stru.write(outName,'pdb')
        # # # exit() 
        # # # scf.write(outName,'cif')
        # # print(cf.stru.xyz_cartn[:,:])
        # # rows = cf.stru.xyz_cartn.shape[0]
        # # cols = cf.stru.xyz_cartn.shape[1]
        # # cf.stru.xyz_cartn = np.reshape(cf.stru.xyz_cartn,(rows,cols))
        # # print(cf.stru.xyz_cartn)
        # COM = [0,0,0]

        # try:

        #     COM=np.mean(cf.stru.xyz_cartn[:,:],axis=0)
        # except Exception as d:
        #     print(str(d) + "this is mean")
        
        # # print(COM)
        # # # countIdx = 0
        # # #for each transition metal find the closest
        # data = ""
        # for m in transitionMetalsMineral:
        #     data = ""
        #     cf_idx_metal = []
        #     cf_coords = []
        #     cf_trans = []
        #     count = 0
        #     for i in range(0,len(cf.stru.label)):
        #         if m == cf.stru[i].label:
        #             # print(m + " "  + cf.stru[i].label)
        #             cf_idx_metal.append((cf.stru[i].label,i))
        #             cf_coords.append(cf.stru[i].xyz_cartn)#,i,m))
        #             cf_trans.append(cf.stru[i])
        #             # continue
        #     # print(cf_coords)

        #     closest = closest_node(COM,np.asarray(cf_coords))
        #     newData = file_name.strip('.cif')+','+ cf_trans[closest].element +','+ str(cf_idx_metal[closest][0])+"_"+"1"+"_"+str(cf_idx_metal[closest][1])+"_"+"A"
            
        #     if count == 0:
        #         data = newData
        #     else:
        #         data += "," + newData
        #     # print(closest)
        #     # print(scf_idx_metal[closest])
        #     # print(type(fileName.strip('.cif')))
        #     # print(type(scf_trans[closest].element))
        #     # print(type(scf_trans[closest].xyz_cartn[0]))
        #     # print(type(scf_trans[closest].xyz_cartn[1]))
        #     # print(type(scf_trans[closest].xyz_cartn[2]))
        #     # print(type(scf_idx_metal[closest][0]+"_"+"1"+"_"+str(scf_idx_metal[closest][1])+"_"+"A"))
            
        #     # print(data)
        #     #print(data)
        #     # print(list(data.items()))
        #     # global data_out
        #     # data_out = []
        #     # data_out = data_out.append(data)
        #     # print(data_out)
        #     process_queue.q.put(data)
    
    except Exception as e:
        print("Error: " + file_name + '_' + str(e))
        # print(data)
    return(data)
        
        

def listener():
    '''listens for messages on a the q, writes to a file'''
    #https://stackoverflow.com/questions/13446445/python-multiprocessing-safely-writing-to-a-file

    fileName = '/home/kenneth/proj/proXtal/amcsd/database/formatName/findgeo/amcsd_findgeo_id_mult.csv'
    # keys = ["mineral","element","x","y","z","id"]
        # print(data)
    with open(fileName,"w") as outWriter:
        while 1:
            m = process_queue.q.get()
            # print(m)
            # print(type(m))
            if m == "kill":
                outWriter.write('killed')
                break
            else:
                mSplit = m.split(',')
                for s in range(0,len(mSplit),3):
                    outStr = mSplit[s] + ',' + mSplit[s+1] + ',' + mSplit[s+2]
                    outWriter.write(outStr + '\n')

            # for i in m:
                # string = ''
                # for idx,key in enumerate(keys):
                    # if idx == 0:
                        # string == str(i[key][0])
                    # else:
                        # string += " " + str(i[key][0])
                # print(string)
                # outWriter.write(m + '\n')
                # outWriter.flush()

def main():

    manager = mp.Manager()
    q = manager.Queue()
    #init objects
    cores = mp.cpu_count()

    # print(cores)
    # https://stackoverflow.com/questions/3827065/can-i-use-a-multiprocessing-queue-in-a-function-called-by-pool-imap
    pool = mp.Pool(cores,f_init,[q])
    watcher = pool.apply_async(listener)

    jobs = []
    with tqdm(total=len(finalData)) as pbar: 
        for i, job in tqdm(enumerate(pool.imap_unordered(process_queue,finalData))):
            # print(i)
            # print(job)
            jobs.append(job)
            # print(job)
            pbar.update()
        # collect results from the workers through the pool result queue
        # for job in jobs: 
            # job.get()
    q.put("kill")
    pool.close()
    pool.join()
    pbar.close()

main()



        # exit()
# pd.to_csv('/home/kenneth/proj/proXtal/amcsd/xCif/singleOccupancy/natural/containsTransitionMetals/finalDataset/amcsd_findgeo_id.csv')

            # print(scf_coords[0])
            # print(scf[768].xyz_cartn,scf[768].label)

            # scf_idx = \
            #     [(scfLabels[s],s) for s in range(0,len(scfLabels)) if any(substring in scfLabels[s] for substring in transitionMetalsMineral)]

            # for i in scf_idx:
            #     print(i[1])
            # scf_coords_trans = [scf[x[1]].xyz_cartn for x in scf_idx_transitionMetalsMineral] 
            # print(scf[768].xyz_cartn) 

            # # print(scf_idx_transitionMetalsMineral)
            # # print(scf_coords_trans)
            
                    # for i in scf.xyz_cartn:        #     print(type(i))
            # # print(ndimage.measurements.center_of_mass(np.array(scf.xyz_cartn)))

            # print(dir(cf))
            # ['addNewAtom', 'angle', 'anisotropy', 'append', 'assignUniqueLabels', 'clear', 'composition', 'copy', 'count', 
            # 'distance', 'element', 'extend', 'getLastAtom', 'index', 'insert', 'label', 'lattice', 'occupancy', 'pdffit', 
            # 'placeInLattice', 'pop', 'read', 'readStr', 'remove', 'reverse', 'sort', 'title', 'tolist', 'write', 'writeStr', 
            # 'x', 'xyz', 'xyz_cartn', 'y', 'z']
            # "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(...)
            # # print(dir(cf))
            # for i in cf:
            #     # print(i.label)
            #     i.__setattr__('label',i.label.upper().replace('_',"").replace("-","").replace("OH","H"))
            #     i.__setattr__('element',i.element.upper())
            #     print(i)
            #     # print(i.label)
            #     # exit()
            # print(fileName.split('_')[0])
            # elements = cf.element
            # labels = cf.label.replace("_","").replace("-","")
            # print(cf.xyz_cartn)

# from progress.bar import Bar

# with Bar('Processing',max=len(finalSet)) as bar:
#     for d in naturalMinerals:
#         if d != 'minerals.txt':

#             for f in os.listdir(os.path.join(natDir,d)):
#                 cf = ReadCif(os.path.join(natDir,d,f))

#                 labels = cf['global']['_atom_site_label']

#                 if any(substring in labels for substring in metals):
#                     # print(d)
#                     # print(labels)
#                     newFile = os.path.join(containsMetalsMinerals,d+"_"+cf['global']['_database_code_amcsd']+".cif")
#                     oldFile = os.path.join(natDir,d,f)
#                     shutil.copyfile(oldFile,newFile)
#     bar.next()
# bar.finish()
                # print(newFile)
                # print(oldFile)
                # print('vi '+ d +'/'+f,cf['global']['_database_code_amcsd'])



# from progress.bar import Bar

# with Bar('Processing',max=len(os.listdir(d))) as bar:
#     for f in os.listdir(d):

#         if os.path.splitext(f)[1] == '.cif' and f not in noGo: 
#             try:
#                 cf = ReadCif(os.path.join(d,f))
#                 if '_atom_site_occupancy' not in cf['global']: #keep if it doesnt have multiple occupancies
#                     count += 1
#                     if os.path.exists(os.path.join(d,cf['global']['_chemical_name_mineral'])) == False:
#                         os.mkdir(os.path.join(d,cf['global']['_chemical_name_mineral']))
                    
#                     shutil.move(os.path.join(d,f),os.path.join(d,cf['global']['_chemical_name_mineral']),f)
#                     # names.append(cf['global']['_chemical_name_mineral'])
#                     # print(f)
#             except Exception as e:
#                 # print(f)
#                 # print(e)
#                 continue
#         bar.next()
#     bar.finish()
    # print(count)
# print(Counter(names))
            # print("test")
            # print(dir(cf))
            # print(dir(cf['global']))
            # print(cf['global'].keys())

            # print(cf['global'][''])
            # exit()
        #         # print(cf['global']['_chemical_name_mineral'])
        #         outName = d+"cif/mineral/test/"+cf['global']['_chemical_name_mineral']+'.cif'
        #         if os.path.exists(outName):
        #             print(outName + "exists")
        #             countExists +=1
        #         else:
        #             outfile = open(outName,"w")
        #             outfile.write(cf.WriteOut())
        #             outfile.close()
        #         # countMineralName += 1

        #     elif '_amcsd_formula_title' in cf['global']: 
        #         # print(cf['global']['_amcsd_formula_title']) 
        #         # countAMCSDName += 1
        #         outName = d+"cif/formTitle/test/"+cf['global']['_amcsd_formula_title']+'.cif'

        #         if os.path.exists(outName):
        #             print(outName + "exists")
        #             countExists +=1
        #         else:
        #             outfile = open(outName,"w")
        #             outfile.write(cf.WriteOut())
        #             outfile.close()
        #     else:
        #         countNone += 1
        #         print(f)
        # except Exception as e:
        #     # print(f)
        #     # print(e)
        #     continue
# print(countExists)
# 14579 6001 0
# print("%s %s %s" % (countMineralName,countAMCSDName,countNone))
