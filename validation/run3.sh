#!/bin/bash
# usage: run3.sh name hyst top bot hinit
export PATH=/opt/r91/bin:$PATH
n=$1; rm -rf $n ${n}_f ${n}_r; python3 gen3.py $n $2 $3 $4 $5
cp -r $n ${n}_f; (cd ${n}_f && timeout 200 /home/claude/fortran/hydrus_ref . < /dev/null 2>&1 | grep -i -E "converged|error" | head -2)
/home/claude/hydrus-rs/target/debug/hydrus1d-cli $n -o ${n}_r -q 2>&1 | tail -2
python3 - <<PY
import csv
f=[]
for l in open("${n}_f/T_LEVEL.OUT"):
    p=l.split()
    if len(p)>=20:
        try: f.append([float(x) for x in p[:20]])
        except: pass
r=list(csv.DictReader(open("${n}_r/T_LEVEL.csv")))
def near(t):
    a=min(f,key=lambda x:abs(x[0]-t)); b=min(r,key=lambda x:abs(float(x['time'])-t)); return a,b
for t in [10,20,30]:
    a,b=near(t)
    print("t=%g F sumvTop=%.4f sumvRoot=%.4f sumvBot=%.4f vol=%.4f hTop=%.4g hBot=%.4g | R %.4f %.4f %.4f %.4f %.4g %.4g"%(t,a[8],a[9],a[10],a[16],a[11],a[13],float(b['sum_vTop']),float(b['sum_vRoot']),float(b['sum_vBot']),float(b['volume']),float(b['hTop']),float(b['hBot'])))
print("rows F/R",len(f),len(r))
PY
