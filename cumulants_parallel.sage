from time import time
import datetime
import multiprocessing as mp

global cumulants

def partitions(points):
    if len(points) == 1:
        yield [ points ]
        return
    first = points[0]
    for smaller in partitions(points[1:]):
        for m, subset in enumerate(smaller):
            yield smaller[:m] + [[ first ] + subset]  + smaller[m+1:]
        yield [ [ first ] ] + smaller

def nonflat(partition,r):
    p = []
    for j in partition:    
        seq = list(map(lambda x: (x-1)//r,j))
        p.append(len(seq) == len(set(seq)))
    return all(p)

def connected(partition,n,r):
    q = []; c = 0
    if n  == 1: return all([len(j)==1 for j in partition])
    for j in partition:
        jk = list(set(map(lambda x: (x-1)//r,j)))
        if(len(jk)>1):            
            if c == 0:
                q = jk; c += 1
            elif(set(q) & set(jk)):
                d=[y for y in (q+jk) if y not in q]
                q = q + d
    return n == len(set(q))

def connectednonflat(n,r):
    points = list(range(1,n*r+1))
    randd = []
    for m, p in enumerate(partitions(points), 1):
        randd.append(sorted(p))
    cnfp = [e for e in randd if (connected(e,n,r) and nonflat(e,r))]
    for rou in range(r,(r-1)*n+2): 
        rs = [d for d in cnfp if len(d)==rou]
        print("Connected non-flat partitions with",rou,"blocks:",len(rs))
    print("Connected non-flat set partitions:",len(cnfp))
    return cnfp

def graphs(G,E,setpartition,n):
    r=len(set(flatten(G)));rhoG = []
    for j in range(n):
        for hop in G: rhoG.append([r*j+hop[0],r*j+hop[1]])
        for l in range(len(E)):
            F=E[l]
            for i in F: rhoG.append([j*r+i,n*r+l+1]);
    for i in setpartition:
        if(len(i)>1):
            b = []
            for j in rhoG:
                b.append([i[0] if ele in i else ele for ele in j])
            rhoG = b
    for i in rhoG: i.sort()
    return rhoG

def inner(n,d,G,E,mu,H,setpartition,z,r):
    rhoG=graphs(G,E,setpartition,n)
    for ll in range(len(E)+1):
        for l in range(1,d+1): z[d*(n*r+ll)+l] = var(str(y)+str(ll)+str('_')+str(l))
    for key in range(1,n*r+1): 
        for l in range(1,d+1): z[key*d+l] = var(str(x)+str(key)+str(x)+str(l))
    edgesrhoG = [i for n, i in enumerate(rhoG) if i not in rhoG[:n]]
    vertrhoG = set(flatten(edgesrhoG));
    for ll in range(len(E)): vertrhoG.remove(n*r+ll+1);
    strr = '*λ'*len(vertrhoG)
    for i in vertrhoG:
        for l in range(1,d+1): strr = '*mu({},{},{})'.format(z[i*d+l],λ,β) + strr
        for l in range(1,d+1): strr = strr + ').integrate({},-infinity,+infinity)'.format(z[i*d+l])
    for i in edgesrhoG:
        for l in range(1,d+1): strr = '*H({},{},{})'.format(z[i[0]*d+l],z[i[1]*d+l],β) + strr
    strr = '('*len(vertrhoG)*d+strr[1:]
    return eval(preparse(strr))

def collect_result(result):
    global cumulants
    global iii
    global tim
    iii=iii+1;
    if (mod(iii,100)==0):
        tim=(time()-t_start2)*(lencnfp-iii)/iii/60
        print('[%d]\r'%(iii),'Est remaining time (minutes):%d'%(tim),end="")
    cumulants+=result
	
def c(n,d,G,E,mu,H):
    global cumulants
    global iii
    global t_start2
    t_start2 = time()
    d_start2 = datetime.datetime.now()
    r=len(set(flatten(G)));
    x,y=var("x,y")
    cumulants = 0; iii = 0
    z = dict(enumerate([str(x)+str(key)+str(x)+str(l) for key in range(0,n*r+1) for l in range(1,d+1)], start=1))
    global lencnfp
    cnfp=connectednonflat(n,r)
    lencnfp=len(cnfp)
    pool = mp.Pool(4) # pool = mp.Pool(mp.cpu_count())
    for setpartition in cnfp: 
        pool.apply_async(func = inner, args=(n,d,G,E,mu,H,setpartition,z,r), callback=collect_result)
    pool.close()
    pool.join()
    print("\n");
    d_end2 = datetime.datetime.now()
    print("Runtime is",(d_end2-d_start2))
    return cumulants._sympy_() 

