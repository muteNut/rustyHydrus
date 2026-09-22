#!/bin/bash
# Builds the ORIGINAL HYDRUS-1D Fortran engine with gfortran (Linux) so that results can be compared.
# Usage: ./build_reference.sh /path/to/fortran_sources   (the *_FOR.txt files or *.FOR files)
# The MS-Fortran specific parts (MSFLIB signal handlers, NARGS/GETARG, getdat/gettim, Windows paths)
# are patched automatically. The physics code is untouched.
set -e
SRC=${1:?directory with HYDRUS.FOR, WATFLOW.FOR, ...}
OUT=${2:-./fortran_ref}
mkdir -p "$OUT"; cd "$OUT"
for f in HYDRUS HYSTER INPUT MATERIAL SINK OUTPUT SOLUTE TIME WATFLOW TEMPER; do
  for cand in "$SRC/${f}_FOR.txt" "$SRC/$f.FOR" "$SRC/$f.for" "$SRC/$f.FOR.txt"; do [ -f "$cand" ] && cp "$cand" "$f.f" && break; done
done
python3 - <<'PY'
import re
def fix_paths(t): return re.sub(r"'\\([A-Za-z_0-9.]+)'", lambda m: "'/"+m.group(1)+"'", t).replace("      use MSFLIB\n","")
s=open('HYDRUS.f',encoding='latin-1').read()
s=re.sub(r"      interface\n.*?      end interface\n","",s,flags=re.S)
s=re.sub(r"      i4ret = SIGNALQQ.*\n","",s)
s=s.replace("iCount = NARGS()","iCount = iargc()+1").replace("call GETARG(i2, cDataPath, status)","call getarg(1, cDataPath)")
i=s.find("*-----------------------------------------------------------------------\n*     Signal handler routines"); j=s.find("      subroutine CloseFiles")
s=s[:i]+s[j:]
open('HYDRUS.f','w',encoding='latin-1').write(fix_paths(s))
for f in ['INPUT','SOLUTE','TIME','OUTPUT','WATFLOW','HYSTER','SINK','MATERIAL','TEMPER']:
    t=open(f+'.f',encoding='latin-1').read().replace("if(iMetHour) then","if(iMetHour.eq.1) then")
    open(f+'.f','w',encoding='latin-1').write(fix_paths(t))
PY
cat > stubs.f <<'F'
      subroutine getdat(a,b,c)
      return
      end
      subroutine gettim(a,b,c,d)
      return
      end
F
gfortran -std=legacy -w -O2 -ffixed-line-length-none -fno-automatic -o hydrus_ref HYDRUS.f HYSTER.f INPUT.f MATERIAL.f SINK.f OUTPUT.f SOLUTE.f TIME.f WATFLOW.f TEMPER.f stubs.f
echo "built $(pwd)/hydrus_ref  -> run:  cd <project dir>; $(pwd)/hydrus_ref . < /dev/null"
echo "NOTE: on Linux the input files must be named exactly Selector.in, Profile.dat, ATMOSPH.IN"
