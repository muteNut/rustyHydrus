import sys,csv,re
def load_f(path):
    rows=[]
    for l in open(path):
        p=l.split()
        if len(p)>=20:
            try:
                rows.append([float(x) for x in p[:20]])
            except: pass
    return rows
f=load_f(sys.argv[1]+"/T_LEVEL.OUT")
r=list(csv.DictReader(open(sys.argv[2]+"/T_LEVEL.csv")))
print("fortran rows",len(f),"rust rows",len(r))
names=["time","rTop","rRoot","vTop","vRoot","vBot","sum_rTop","sum_rRoot","sum_vTop","sum_vRoot","sum_vBot","hTop","hRoot","hBot","runoff","sum_runoff","volume","sum_infil","sum_evap","tlevel"]
keys=["time","rTop","rRoot","vTop","vRoot","vBot","sum_rTop","sum_rRoot","sum_vTop","sum_vRoot","sum_vBot","hTop","hRoot","hBot","runoff","sum_runoff","volume","sum_infil","sum_evap","tlevel"]
# match by tlevel
rd={int(float(x["tlevel"])):x for x in r}
maxd={k:0 for k in keys}
n=0
for row in f:
    tl=int(row[19])
    if tl not in rd: continue
    n+=1
    x=rd[tl]
    for i,k in enumerate(keys[:19]):
        a=row[i]; b=float(x[k])
        d=abs(a-b)/max(abs(a),abs(b),1e-3)
        if d>maxd[k]: maxd[k]=d
print("matched",n)
for k in keys[:19]: print("%-10s max rel diff %.2e"%(k,maxd[k]))
# time mismatch
tf=[row[0] for row in f]; tr=[float(x["time"]) for x in r]
print("last times",tf[-1],tr[-1])
