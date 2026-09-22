import sys,csv
case=sys.argv[1]; ns=int(sys.argv[2]); neq=(len(sys.argv)>3 and sys.argv[3]=="neq")
blocks={};cur=None
for l in open(case+"_f/NOD_INF.OUT"):
    if l.strip().startswith("Time:"):
        cur=float(l.split()[1]); blocks[cur]=[]
    else:
        p=l.split()
        if cur is not None and len(p)>=11+ns:
            try: blocks[cur].append([float(x) for x in p])
            except: pass
r=list(csv.DictReader(open(case+"_r/NOD_INF.csv")))
for t in [10.0,20.0,30.0]:
    fb=blocks[t]; rb=[x for x in r if abs(float(x['time'])-t)<1e-6]
    out=[]
    md={"head":0,"temp":0}
    for j in range(ns): md["conc%d"%(j+1)]=0
    if neq:
        for j in range(ns): md["sorb%d"%(j+1)]=0
    scale={}
    for a,b in zip(fb,rb):
        md["head"]=max(md["head"],abs(a[2]-float(b["head"]))/max(abs(a[2]),1))
        md["temp"]=max(md["temp"],abs(a[10]-float(b["temp"])))
        for j in range(ns):
            k="conc%d"%(j+1)
            mx=max(x[11+j] for x in fb)
            md[k]=max(md[k],abs(a[11+j]-float(b[k]))/max(mx,1e-9))
            if neq:
                k2="sorb%d"%(j+1)
                mx=max(abs(x[11+ns+j]) for x in fb)+1e-12
                md[k2]=max(md[k2],abs(a[11+ns+j]-float(b[k2]))/mx)
    print(t,{k:"%.1e"%v for k,v in md.items()})
# solute cumulative
def load(path):
    rows=[]
    for l in open(path):
        p=l.split()
        if len(p)>=14:
            try: rows.append([float(x) for x in p])
            except: pass
    return rows
f=load(case+"_f/solute1.out"); 
print("solute1.out last row F:",f[-1][:14])
tl=list(csv.DictReader(open(case+"_r/T_LEVEL.csv")))[-1]
print("R last: cvTop=%s cvBot=%s sum_cvTop1=%s sum_cvBot1=%s ch0=%s ch1=%s cTop=%s cBot=%s"%(tl['cvTop1'],tl['cvBot1'],tl['sum_cvTop1'],tl['sum_cvBot1'],tl['sum_cvCh0_1'],tl['sum_cvCh1_1'],tl['cTop1'],tl['cBot1']))
