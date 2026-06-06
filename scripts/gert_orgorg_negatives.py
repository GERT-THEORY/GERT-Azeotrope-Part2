#!/usr/bin/env python3
"""
GERT organic-organic — negative (maximum-boiling) azeotrope criterion.

Parameter-free, mechanistically grounded:
  negative azeotrope  <=>  dC = (a1-a2)(b2-b1) > 0          (donor-acceptor
                            complementarity: cross-association > self-association)
                      AND  |ln r| < (E_HB/RT) * dC          (association overcomes
                            volatility; ln r = (88/R)(Tb2-Tb1)/T via Trouton)

Constants are independent physics, NOT fitted to azeotrope outcomes:
  - Trouton: dSvap = 88 J/mol/K
  - E_HB ~ 25 kJ/mol (H-bond enthalpy scale; result robust over 20-30 kJ/mol)
  - a (alpha), b (beta): tabulated Kamlet-Taft descriptors

Benchmark: ie9b03694 SI (Gmehling Azeotropic Data). On the org-org subset with
full descriptors (n=2141, 68 negatives):
  GERT:  negative recall 75%, precision 46%
  SM12:  negative recall  1%, precision  6%
"""
import pandas as pd, numpy as np
from _kt import AB      # Kamlet-Taft alpha,beta table
from _tb import Tb      # normal boiling points (K)

F="/mnt/user-data/uploads/ie9b03694_si_002.xlsx"
R=8.314; DS_TROUTON=88.0; E_HB=25_000.0     # J/mol

def truth(c):
    c=str(c).strip()
    if '/' in c: return 'amb'
    if c in ('Z','z','Z2','ZC'): return 'none'
    if c.startswith(('AO','AE','AX')) or c=='A': return 'neg' if 'N' in c else 'pos'
    return 'amb'

def load():
    s=pd.read_excel(F,sheet_name="Exp-vs-SMx")
    s['truth']=s.Experimental.map(truth)
    s['SM12']=s['SM12 Prediction'].map(lambda c:{'Z':'none','X':'pos','N':'neg'}.get(str(c).strip(),'amb'))
    nrm=lambda x:str(x).strip().lower()
    s=s[~((s['Component 1'].map(nrm)=='water')|(s['Component 2'].map(nrm)=='water'))]
    s=s[s.truth.isin(['none','pos','neg'])].copy()
    s['ab1']=s['Component 1'].map(lambda x:AB.get(nrm(x)))
    s['ab2']=s['Component 2'].map(lambda x:AB.get(nrm(x)))
    s['Tb1']=s['Component 1'].map(lambda x:Tb.get(nrm(x)))
    s['Tb2']=s['Component 2'].map(lambda x:Tb.get(nrm(x)))
    return s.dropna(subset=['ab1','ab2','Tb1','Tb2']).copy()

def gert_negative(d):
    a1,b1=d.ab1.str[0],d.ab1.str[1]; a2,b2=d.ab2.str[0],d.ab2.str[1]
    dC=(a1-a2)*(b2-b1)                                   # complementarity excess
    Tm=(d.Tb1+d.Tb2)/2
    lnr=(DS_TROUTON/R)*(d.Tb2-d.Tb1).abs()/Tm           # |ln(P1/P2)| via Trouton
    s=E_HB/(R*Tm.mean())                                 # physical scale ~8.6
    return (dC>0) & (lnr < s*dC)

if __name__=="__main__":
    d=load(); d['gneg']=gert_negative(d)
    neg=d[d.truth=='neg']; nn=d[d.truth!='neg']
    rec=d.gneg[neg.index].mean()*100
    prec=d.gneg[neg.index].sum()/max(d.gneg.sum(),1)*100
    print(f"n={len(d)}  negatives={len(neg)}")
    print(f"GERT negative: recall {rec:.0f}%  precision {prec:.0f}%  (false+={d.gneg[nn.index].sum()})")
    print(f"SM12 negative: recall {(d.SM12=='neg')[neg.index].mean()*100:.0f}%  "
          f"precision {(d[d.SM12=='neg'].truth=='neg').mean()*100:.0f}%")
