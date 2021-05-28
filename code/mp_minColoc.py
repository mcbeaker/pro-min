from tqdm import tqdm
import multiprocessing as mp
from os import path
import pandas as pd
import re as re

global kpath
kpath ='/Users/ken/googleDrive/proj/ENIGMA/RISE2020'
fLoc = path.join(kpath,"teams/Falade_Kenneth_RISE2020/data/mineralLocations.csv")
fInt = path.join(kpath,"teams/Falade_Kenneth_RISE2020/data/mineralNames.csv")

global dfLoc
dfLoc = pd.read_csv(fLoc,sep="\t",usecols=['MED_minerals_at_loc'])
dfLoc.dropna(inplace=True)

global dfInt
dfInt = pd.read_csv(fInt,usecols=['mineralName'])
# dfIntList = list(dfInt.mineralName)

global seen 
seen = []

# dfLoc = dfLoc[dfLoc.MED_minerals_at_loc.str.contains("|".join(dfIntList))]
# print(dfLoc)
# exit()lt


def f_init(q):
    process_queue.q = q

def process_queue(id):

    minInt = set([dfInt.iloc[id].mineralName])
    # print(dfInt.iloc[id].mineralName)
    coMinLo = pd.DataFrame()

    #filtered for rows that contain mineral of interest within this for loop
    # if dfLoc.MED_minerals_at_loc.str.contains(dfInt.iloc[id].mineralName).any():
    coMinLo = dfLoc[dfLoc.MED_minerals_at_loc.str.contains(dfInt.iloc[id].mineralName)]

    out_coLoc = set()
    # count = 0
    coLocMin = []
    for index,row in coMinLo.iterrows():
        # print(dfInt.iloc[id].mineralName)
        rowCoMinLoc = set(row.MED_minerals_at_loc.split(','))

        diff = rowCoMinLoc.difference(minInt)

        for minLoc in diff:
            pair = (dfInt.iloc[id].mineralName,minLoc)
            pair1 = (minLoc,dfInt.iloc[id].mineralName)

            # if pair not in seen and pair1 not in seen:
            #     seen.append(pair)
            #     out_coLoc.update([pair])

            # if pair in seen and pair1 not in seen:
            #     continue

    # print(out_coLoc)
    out_coLoc = [" ".join(tup) for tup in out_coLoc]
    # print(out_coLoc)
    # print(out_coLoc[0:1])
    process_queue.q.put(out_coLoc)
    return out_coLoc

def listener():
    '''listens for messages on a the q, writes to a file'''
    #https://stackoverflow.com/questions/13446445/python-multiprocessing-safely-writing-to-a-file
    
    pdir = path.join(kpath,'teams/Falade_Kenneth_RISE2020/data')

    # if os.path.exists(os.path.join(pdir,"micro_multiproc.csv")):
    #     print("Output file exists already")
    #     exit()
    
    # print("Output jobs to file\n")

    with open(path.join(pdir,"coLocMin.csv"),"w") as outWriter:
        while 1:
            proc = process_queue.q.get()
            # print(type(m))
            # print("this: ",proc,'\n')

            if proc == "kill":
                outWriter.write('killed')
                break
            else:
                # for i in re.findall("[^ ]+ [^ ]+",proc):
                # print(proc)
                # outWriter.write(str(proc)+'\n')
                # outWriter.flush()
                for i in proc:
            #         # print(i)
                    outWriter.write(str(i) + "\n")
                    outWriter.flush()

def main():
    
    manager = mp.Manager()
    q = manager.Queue()
    #init objects
    cores = mp.cpu_count()
    # print(cores)
    # https://stackoverflow.com/questions/3827065/can-i-use-a-multiprocessing-queue-in-a-function-called-by-pool-imap
    pool = mp.Pool(cores,f_init,[q])

    #put listeners to work
    watcher = pool.apply_async(listener)

    #fire off workers
    jobs = []

    # print(id)
    ind = dfInt.index.values.tolist()
    with tqdm(total=len(ind)) as pbar: 
        for i, job in tqdm(enumerate(pool.imap_unordered(process_queue,ind))):
            jobs.append(job)
            # print(job)
            pbar.update()
            # if i < 10:
                # jobs.append(job)
                # break
                # pbar.update()
            # print(type(job))

       
    q.put("kill")
    pool.close()
    pool.join()
    pbar.close()

main()