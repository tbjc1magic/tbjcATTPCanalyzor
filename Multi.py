
from Analyzor.DataFactory import DataFactory
from Analyzor import VertexAnalyzor
import time
from subprocess import call
from multiprocessing import Pool
import operator
import json
import numpy as np

def runProcess(exe):
    p = subprocess.Popen(exe,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()

    err = p.stderr.readlines()

    if len(err) != 0:
        for line in err:
            print line
        exit(1)

    out = p.stdout.readlines()
    return out

SQLpath ="data/SQL/10C/"

def ProcessFile(fname):

    r = fname.split('/')[-1].split('.')[0].split('_')
    runID,fID = map(int,r)

    dp = DataFactory(fname,SQLpath+'ProtoMap.db')
    dist = []
    for i in  sorted(dp.t3['EventID'].unique()):

        try:
            #print i
            image = dp.ConstructImage(i)
            image = VertexAnalyzor.FilterBackground(image)
            points,(xc,yc) = VertexAnalyzor.GetEventPositions(image,0)
            r = VertexAnalyzor.GetEventInfo(points,(xc,yc))
            dist.append({'runID':runID,'fileID':fID,
                'eventID':i,'data':r})
        except:
            pass

    return dist

if __name__ == "__main__":

    pool = Pool(processes=4)

    param = []
    runs = [85,88]

    for run in runs:
        for i in range(16):
            path = SQLpath+'{:04d}_{:04d}.db'.format(run,i)
            param.append(path)

    start_time = time.time()
    res = pool.map(ProcessFile,param)

    end_time = time.time()

    print "total process cost "+str(end_time - start_time)

    data = reduce(operator.add,res)

    with open('data.dat','w') as f:
        json.dump(data,f)
