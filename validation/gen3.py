import os,sys
def gen(d, hyst=0, top="atm", bot="free", hinit="dry", dtfixed=False):
    os.makedirs(d,exist_ok=True)
    topinf = "t" if top=="atm" else "f"
    kodtop = -1 if top=="atm" else 1
    atmbc = "t"  if top=="atm" else "f"
    botline = {"free":"f f t f -1 f 0","seep":"f f f t -1 f 0","gwl":"f f f f -1 f 0".replace("f f f f","f t f f"),"drain":"f f f f -1 t 0","consth":"f f f f 1 f 0","constq":"f f f f -1 f 0"}[bot]
    if hyst==0:
        mats="""   thr     ths    Alfa      n         Ks       l
   0.078   0.43   0.036   1.56     24.96     0.5
   0.02    0.35   0.145   2.68    712.8      0.5
"""
        hy="Model   Hysteresis\n  0          0\n"
    else:
        hy=f"Model   Hysteresis\n  0          {hyst}\nIKappa\n -1\n"
        mats="""   thr     ths    Alfa      n         Ks       l  thm thsw alfaw ksw
   0.078   0.43   0.036   1.56     24.96     0.5  0.43 0.40 0.072 15.0
   0.02    0.35   0.145   2.68    712.8      0.5  0.35 0.33 0.29 500.0
"""
    extra=""
    if top=="const" or False: pass
    # constant flux read: when !TopInf and KodTop==-1  or bottom const q
    needconst = (top!="atm" and kodtop==-1) or (bot=="constq")
    if needconst: extra += "rTop rBot rRoot\n 0 -0.05 0\n"
    if bot=="gwl": extra += "GWL0L Aqh Bqh\n 0 -0.5 -0.02\n"
    if bot=="drain": extra += "iPosDr\n 1\nzBotDr spacing entres\n 100 800 0\nKhTop\n 30\n"
    dt = "0.01      0.01       0.01" if dtfixed else "0.001      1e-05       5"
    sel = f"""Pcp_File_Version=4
*** BLOCK A: BASIC INFORMATION *****************************************
Heading
Test3
LUnit  TUnit  MUnit
cm
days
mmol
lWat   lChem lTemp  lSink lRoot lShort lWDep lScreen lVariabBC lEquil lInverse
 t     f     f      {'t' if top=='atm' else 'f'}     f     f      f     f       {atmbc}         t      f
lSnow  lHP1   lMeteo  lVapor lActiveU lFluxes lIrrig  lDummy  lDummy  lDummy
 f       f       f       f       f       f       f       f       f       f
NMat    NLay  CosAlpha
  2       2       1
*** BLOCK B: WATER FLOW INFORMATION ************************************
MaxIt   TolTh   TolH
  20    0.001     1
TopInf WLayer KodTop InitCond
 {topinf}     f      {kodtop}       f
BotInf qGWLF FreeD SeepF KodBot DrainF  hSeep
 {botline}
{extra}hTab1   hTabN
    0.001   200
{hy}{mats}*** BLOCK C: TIME INFORMATION ******************************************
        dt       dtMin       dtMax     DMul    DMul2  ItMin ItMax  MPL
       {dt}        1.3     0.7     3     7     3
      tInit        tMax
        0          30
  lPrintD  nPrStep tPrintInterval lEnter
     f       1         1            f
 TPrint
       10   20   30
"""
    if top=="atm":
        sel += """*** BLOCK G
     Model   OmegaC
       0                                   1
P0       P2H      P2L      P3       r2H      r2L
-10     -200     -800     -8000     0.5      0.1
POptm
-25 -25
"""
    open(d+"/Selector.in","w").write(sel)
    lines=["Pcp_File_Version=4","    1","    x h","  101 0 0 0.0"]
    for i in range(101):
        x=-i*1.0
        mat=1 if x>-50 else 2
        beta=1.0 if x>-30 else 0.0
        if hinit=="dry": h=-200.0 if mat==1 else -80.0
        elif hinit=="hydro": h=-70.0-x
        else: h=-5.0 if x<-50 else -100
        if top=="const" and i==0: h=-10.0
        lines.append("%d %.5f %.3f %d %d %.4f 1 1 1"%(i+1,x,h,mat,mat,beta))
    lines += ["    1","    51"]
    open(d+"/Profile.dat","w").write("\n".join(lines)+"\n")
    if top=="atm":
        T=[(2,3.0,0.3,0.4),(4,0.0,0.5,0.5),(8,0.0,0.5,0.5),(12,8.0,0.3,0.4),(20,0.0,0.6,0.6),(30,0.0,0.6,0.6)]
        if hyst>0: T=[(3,3.0,0.0,0.0),(6,0.0,0.6,0.5),(9,4.0,0.0,0.0),(12,0.0,0.7,0.5),(16,5.0,0.0,0.0),(20,0.0,0.7,0.5),(25,3,0,0),(30,0.0,0.5,0.4)]
        if bot=="seep" or hinit=="wet": T=[(2,15.0,0.3,0.4),(6,0.0,0.3,0.4),(10,12.0,0.3,0.4),(30,0.0,0.6,0.6)]
        rows="".join(f" {t} {p} {e} {r} 10000 0 0 0\n" for t,p,e,r in T)
        open(d+"/ATMOSPH.IN","w").write("Pcp_File_Version=4\n***\nMaxAL\n %d\n DailyVar  SinusVar  lLai  lBCCycles lInterc\n f f f f f\n hCritS\n 0\n hdr\n"%len(T)+rows+"end\n***\n")
if __name__=="__main__":
    a=sys.argv
    gen(a[1],int(a[2]),a[3],a[4],a[5],len(a)>6)
