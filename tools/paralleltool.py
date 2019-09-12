import numpy as np

def chunks(l, n):
    n = max(1, n)
    lenl = len(l)
    chucklist = [l[i:i+n] for i in xrange(0, lenl, n)]
    return chucklist

    
def joindata(results):
    argsortlist=np.argsort([results[i][0] for i in range(len(results))])
    totdict = np.array([results[i][1] for i in range(len(results))])
    totdict=totdict[argsortlist]
    comdict={}
    for key in totdict[0]:
        totarray=np.array([])
        for indict in totdict:
            totarray = np.append(totarray,indict[key])
        comdict[key] = totarray.flatten()
    return comdict

def collect_results(result):
    results.extend(result)