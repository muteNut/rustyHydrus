import os,sys
def gen(d, variant="A"):
    os.makedirs(d,exist_ok=True)
    ns = 1 if variant in("A","C") else 2
    equil = "t" if variant=="A" else "f"
    nl = variant=="B"
    frac = 1.0 if variant=="A" else 0.5
    imm = 0.0
    if variant=="C": frac=1.0; imm=0.03
    if variant=="C": equil="f"
    beta = 0.8 if variant=="A" else 1.0
    sel = f"""Pcp_File_Version=4
*** BLOCK A: BASIC INFORMATION *****************************************
Heading
Test2 solute+heat variant {variant}
LUnit  TUnit  MUnit
cm
days
mmol
lWat   lChem lTemp  lSink lRoot lShort lWDep lScreen lVariabBC lEquil lInverse
 t     t     t      t     f     f      f     f       t         {equil}      f
lSnow  lHP1   lMeteo  lVapor lActiveU lFluxes lIrrig  lDummy  lDummy  lDummy
 f       f       f       f       f       f       f       f       f       f
NMat    NLay  CosAlpha
  2       2       1
*** BLOCK B: WATER FLOW INFORMATION ************************************
MaxIt   TolTh   TolH
  20    0.001     1
TopInf WLayer KodTop InitCond
 t     f      -1       f
BotInf qGWLF FreeD SeepF KodBot DrainF  hSeep
 f     f     t     f      -1     f       0
hTab1   hTabN
    0.001   200
Model   Hysteresis
  0          0
   thr     ths    Alfa      n         Ks       l
   0.078   0.43   0.036   1.56     24.96     0.5
   0.02    0.35   0.145   2.68    712.8      0.5
*** BLOCK C: TIME INFORMATION ******************************************
        dt       dtMin       dtMax     DMul    DMul2  ItMin ItMax  MPL
       0.01      0.01       0.01        1.3     0.7     3     7     3
      tInit        tMax
        0          30
  lPrintD  nPrStep tPrintInterval lEnter
     f       1         1            f
 TPrint(1),TPrint(2),...,TPrint(MPL)
       10   20   30
*** BLOCK E: HEAT TRANSPORT INFORMATION ********************************
Qn Qo Disper B1 B2 B3 Cn Co Cw
   0.57   0.001   5.   1.4e17  7.5e17  3.0e18  1.43e14  1.87e14  3.12e14
   0.60   0.001   5.   1.4e17  7.5e17  3.0e18  1.43e14  1.87e14  3.12e14
Ampl tPeriod iCampbell SnowMF lSnowInput
  0   1   0  0.45  f
kTopT  TTop  kBotT  TBot
  1      20   0     20
*** BLOCK F: SOLUTE TRANSPORT INFORMATION ******************************
Epsi  lUpW  lArtD lTDep  cTolA    cTolR   MaxItC    PeCr    No.Solutes  lTort   iBacter  lFiltr
 0.5    f     t    f     0.001    0.01     20       2.       {ns}           t        0       f
iNonEqual lWatDep lDualNEq lInitM lInitEq lTort
   0        f       f       f      {'t' if variant!='A' else 'f'}       f
Bulk.d.  DisperL.  Frac  Mobile WC
  1.3   2.   {frac}   {imm}
  1.5   1.   {frac}   {imm}
"""
    for j in range(ns):
        mu=0.02 if j==0 else 0.0
        gam=0.0 if j==0 else 0.05
        sel += f"""DifW     DifG  (solute {j+1})
  0.5   0.0
Mat  Ks  Nu  Beta  Henry  SinkL1  SinkS1  SinkG1  SinkL1p SinkS1p SinkG1p SinkL0 SinkS0 SinkG0 Alfa
  {0.3 if j==0 else 0.1}   0.  {beta if j==0 else 1.0}   0.  {mu}  {mu}   0.  {gam}  {gam} 0.  0.  0.  0.  0.15
  {0.2 if j==0 else 0.05}   0.  {beta if j==0 else 1.0}   0.  {mu}  {mu}   0.  {gam}  {gam} 0.  0.  0.  0.  0.15
"""
    kt=" ".join(["-1"]+["1.0"]*ns) 
    sel += "kTopCh cTop kBotCh cBot\n " + " ".join(["-1"]+["1.0"]*ns+["0"]+["0.0"]*ns) + "\ntPulse\n  10\n"
    sel += """*** BLOCK G: ROOT WATER UPTAKE INFORMATION *****************************
     Model   OmegaC
"""
    sel += " 0  " + " ".join(["0.0"]*ns) + "  1\n"
    sel += """P0       P2H      P2L      P3       r2H      r2L
-10     -200     -800     -8000     0.5      0.1
POptm(1),POptm(2),...,POptm(NMat)
-25 -25
lSolRed
 f
*** END OF INPUT FILE 'SELECTOR.IN' ************************************
"""
    open(d+"/Selector.in","w").write(sel)
    lines=["Pcp_File_Version=4","    1","    x       h   Mat  Lay   Beta   Axz   Bxz   Dxz Temp Conc"," 101 0 %d 0.0"%ns]
    for i in range(101):
        x=-i*1.0
        mat=1 if x>-50 else 2
        beta=1.0 if x>-30 else 0.0
        h=-200.0 if mat==1 else -80.0
        c=" ".join(["%.3f"%(5.0 if i<20 and k==0 else 0.0) for k in range(ns)])
        s=""
        if equil=="f": s=" "+" ".join(["%.3f"%0.0]*ns)
        lines.append("%d %.5f %.3f %d %d %.4f 1 1 1 15. %s%s"%(i+1,x,h,mat,mat,beta,c,s))
    lines += ["    3","    21   51   101"]
    open(d+"/Profile.dat","w").write("\n".join(lines)+"\n")
    rec=[]
    T=[(2,3.0,0.3,0.4),(4,0.0,0.5,0.5),(8,0.0,0.5,0.5),(12,8.0,0.3,0.4),(20,0.0,0.6,0.6),(30,0.0,0.6,0.6)]
    hdr="""Pcp_File_Version=4
*** BLOCK I: ATMOSPHERIC INFORMATION  **********************************
MaxAL
  6
 DailyVar  SinusVar  lLai  lBCCycles lInterc
   f       f       f      f        f
 hCritS
  0
 tAtm Prec rSoil rRoot hCritA rB hB hT tTop tBot Ampl cTop cBot
"""
    body=""
    for k,(t,p,e,r) in enumerate(T):
        cs=" ".join(["%.2f 0.0"%(2.0+k) ]*ns)
        body+=f" {t} {p} {e} {r} 10000 0 0 0 {18+k} 20 {3.0} {cs}\n"
    open(d+"/ATMOSPH.IN","w").write(hdr+body+"end\n*** END\n")
if __name__=="__main__":
    gen(sys.argv[1],sys.argv[2])
