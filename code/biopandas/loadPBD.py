# from biopandas.pdb import PandasPdb
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from glob import glob
from tqdm import tqdm
import multiprocessing as mp
import os
global pdbs
global home
home='/media/kenneth/potions/pdb/050121'

pdbs = glob(os.path.join(home,'*.cif'))
# print(pdbs[0])

def f_init(q):
    process_queue.q = q

def process_queue(id):
    pdbPath = pdbs[id]
    # print(pdbPath)
    output = MMCIF2Dict(pdbPath)
    strct = ''
    if '_reflns.d_resolution_high' in output:
        if output['_reflns.d_resolution_high'][0] != '?': 
            if float(output['_reflns.d_resolution_high'][0]) <= 1.0:
                strct=output['_entry.id']
    # struct = parser.get_structure()
    # output = struct.get_header()
    # ppdb = PandasPdb().read_pdb(pdbPath)
    # print('blah')
    # output= ppdb.df['Others']
    return strct

def listener():
    outFile = os.path.join(home,'le1Ang','pdbsLE1Ang.txt')
    with open(outFile,"w") as outWriter:
        m = process_queue.q.get()
        # print(m+'\n')
        # print(type(m))
        if m == "kill":
            # outWriter.write('killed')
            exit()
        else:
            outWriter.write(m+'\n')
            outWriter.flush()

def main():
    

    manager = mp.Manager()
    q = manager.Queue()
    #init objects
    cores = mp.cpu_count()
    pool = mp.Pool(cores,f_init,[q])
    # watcher = pool.apply_async(listener)

    jobs = []
    with tqdm(total=len(pdbs)) as pbar:
        for i, job in tqdm(enumerate(pool.imap_unordered(process_queue,range(0,len(pdbs))))): #len(pdbs))))):
            jobs.append(job)
            pbar.update()


    q.put("kill")
    pool.close()
    pool.join()
    pbar.close()
main()