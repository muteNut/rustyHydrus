import os,sys
def gen(d, nnode=101, depth=100.0):
    os.makedirs(d,exist_ok=True)
    sel = """Pcp_File_Version=4
*** BLOCK A: BASIC INFORMATION *****************************************
Heading
Test1 infiltration/evaporation with roots
LUnit  TUnit  MUnit
cm
days
mmol
lWat   lChem lTemp  lSink lRoot lShort lWDep lScreen lVariabBC lEquil lInverse
 t     f     f      t     f     f      f     f       t         t      f
lSnow  lHP1   lMeteo  lVapor lActiveU lFluxes lIrrig  lDummy  lDummy  lDummy
 f       f       f       f       f       f       f       f       f       f
NMat    NLay  CosAlpha
  2       2       1
*** BLOCK B: WATER FLOW INFORMATION ************************************
MaxIt   TolTh   TolH       (maximum number of iterations and tolerance for water content and pressure head)
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
       0.001      1e-05       5        1.3     0.7     3     7     3
      tInit        tMax
        0          30
  lPrintD  nPrStep tPrintInterval lEnter
     f       1         1            f
 TPrint(1),TPrint(2),...,TPrint(MPL)
       10   20   30
*** BLOCK G: ROOT WATER UPTAKE INFORMATION *****************************
     Model   (0 - Feddes, 1 - S shape)   OmegaC
       0                                   1
P0       P2H      P2L      P3       r2H      r2L
-10     -200     -800     -8000     0.5      0.1
POptm(1),POptm(2),...,POptm(NMat)
-25 -25
*** END OF INPUT FILE 'SELECTOR.IN' ************************************
"""
    open(d+"/Selector.in","w").write(sel)
    lines=["Pcp_File_Version=4","    1","    x       h   Mat  Lay   Beta   Axz   Bxz   Dxz"," %d 0 0 0.0"%nnode]
    dx=depth/(nnode-1)
    for i in range(nnode):
        x=-i*dx
        mat=1 if x>-50 else 2
        beta=1.0 if x>-30 else 0.0
        h=-200.0 if mat==1 else -80.0
        lines.append("%d %.5f %.3f %d %d %.4f 1 1 1"%(i+1,x,h,mat,mat,beta))
    lines.append("    3")
    lines.append("    21   51   101")
    open(d+"/Profile.dat","w").write("\n".join(lines)+"\n")
    atm = """Pcp_File_Version=4
*** BLOCK I: ATMOSPHERIC INFORMATION  **********************************
MaxAL                    (MaxAL = number of atmospheric data-records)
  6
 DailyVar  SinusVar  lLai  lBCCycles lInterc
   f       f       f      f        f
 hCritS                 (max. allowed pressure head at the soil surface)
  0
    tAtm        Prec       rSoil       rRoot      hCritA          rB          hB          hT
     2           3.0        0.3         0.4        10000          0           0           0
     4           0.0        0.5         0.5        10000          0           0           0
     8           0.0        0.5         0.5        10000          0           0           0
     12          8.0        0.3         0.4        10000          0           0           0
     20          0.0        0.6         0.6        10000          0           0           0
     30          0.0        0.6         0.6        10000          0           0           0
end
*** END OF INPUT FILE 'ATMOSPH.IN' **************************************
"""
    open(d+"/ATMOSPH.IN","w").write(atm)
if __name__=="__main__":
    gen(sys.argv[1])
