from Analyzor.DataFactory import DataFactory
from Analyzor.VertexAnalyzor import GetRange,FilterBackground
import time
from subprocess import call
from multiprocessing import Pool
import operator

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
    print fname
    dp = DataFactory(fname,SQLpath+'ProtoMap.db')

    dist = []
    for i in  sorted(dp.t3['EventID'].unique()):

        try:
            #print i
            image = dp.ConstructImage(i)
            image = FilterBackground(image)
            r = GetRange(image.astype(np.uint8),0)
            dist.append(r)
        except:
            pass

    return dist

if __name__ == "__main__":

    pool = Pool(processes=8)

    param = []
    runs = [85]

    for run in runs:
        for i in range(2):
            path = SQLpath+'{:04d}_{:04d}.db'.format(run,i)
            param.append(path)

    start_time = time.time()
    res = pool.map(ProcessFile,param)

    end_time = time.time()

    print "total process cost "+str(end_time - start_time)

    data = [_ for _ in reduce(operator.add,res) if _<1e6]

    with open('data.dat','w') as f:
        json.dump(data,f)
