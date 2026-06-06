#!/usr/bin/env python3
"""Generate the full organic-organic SI table (all 2,141 systems) with component
names, phi_type, Delta C, |ln r|, GERT and SM12 predictions, observed behaviour.
Reproduces GERT raw 73.1% / SM12 73.9%. Output: GERT_Paper2_SI_orgorg_2141.csv
Requires: ie9b03694_si_002.xlsx (Roese et al. 2019) and _kt.py/_tb.py/_hansen.py."""
import pandas as pd
from _kt import AB; from _tb import Tb; from _hansen import HAN
DATA="ie9b03694_si_002.xlsx"; R=8.314; TROUTON=88.0; E_HB=25_000.0; SPINODAL=2.0; A_NA,B_NA=0.10,0.25
APROTIC={'nitromethane','nitroethane','1-nitropropane','2-nitropropane','nitrobenzene',
 '1-methyl-2-nitrobenzene','nitrile acetic acid','nitrile propanoic acid','nitrile butanoic acid','nitrile benzoic acid'}
nrm=lambda x:str(x).strip().lower()
def truth(c):
    c=str(c).strip()
    if '/' in c: return 'amb'
    if c in ('Z','z','Z2','ZC'): return 'none'
    if c.startswith(('AO','AE','AX')) or c=='A': return 'neg' if 'N' in c else 'pos'
    return 'amb'
def classify(a1,b1,a2,b2,dd1,dp1,V1,dd2,dp2,V2,T1,T2,n1,n2):
    if n1 in APROTIC: a1=0.0
    if n2 in APROTIC: a2=0.0
    Tm=(T1+T2)/2.0; Vm=(V1+V2)/2.0
    lnr=(TROUTON/R)*abs(T2-T1)/Tm
    X=a1*b2+a2*b1; S=a1*b1+a2*b2; dC=(a1-a2)*(b2-b1)
    phi=X/(X+S) if (X+S)>0 else float('nan')
    gdisp=Vm*((dd1-dd2)**2+(dp1-dp2)**2)/(R*Tm)
    NA=(a1<A_NA and a2<A_NA and b1<B_NA and b2<B_NA)
    dev=min(gdisp,SPINODAL) if NA else min((E_HB/(R*Tm))*abs(dC),SPINODAL)
    pred='none' if lnr>=dev else ('pos' if (NA or dC<0) else 'neg')
    return pred,X,S,dC,phi,lnr,dev,NA
def main():
    sv=pd.read_excel(DATA,sheet_name="Exp-vs-SMx")
    sv['obs']=sv.Experimental.map(truth)
    sv['SM12']=sv['SM12 Prediction'].map(lambda c:{'Z':'none','X':'pos','N':'neg'}.get(str(c).strip(),'amb'))
    sv=sv[~((sv['Component 1'].map(nrm)=='water')|(sv['Component 2'].map(nrm)=='water'))]
    sv=sv[sv.obs.isin(['none','pos','neg'])].copy()
    rows=[]
    for _,r in sv.iterrows():
        n1,n2=nrm(r['Component 1']),nrm(r['Component 2'])
        ab1,ab2=AB.get(n1),AB.get(n2); h1,h2=HAN.get(n1),HAN.get(n2); t1,t2=Tb.get(n1),Tb.get(n2)
        if None in (ab1,ab2,h1,h2,t1,t2): continue
        pred,X,S,dC,phi,lnr,dev,NA=classify(ab1[0],ab1[1],ab2[0],ab2[1],h1[0],h1[1],h1[2],h2[0],h2[1],h2[2],t1,t2,n1,n2)
        rows.append({'Component 1':r['Component 1'],'Component 2':r['Component 2'],'observed':r.obs,
            'SM12':r.SM12,'GERT':pred,'channel':('disp' if NA else 'assoc'),
            'phi_type':round(phi,3) if phi==phi else '','deltaC':round(dC,3),
            'abs_lnr':round(lnr,3),'W_dev':round(dev,3),
            'GERT_correct':int(pred==r.obs),'SM12_correct':int(r.SM12==r.obs)})
    SI=pd.DataFrame(rows); SI.to_csv('GERT_Paper2_SI_orgorg_2141.csv',index=False)
    print(f"n={len(SI)}  GERT raw={SI.GERT_correct.mean()*100:.1f}%  SM12 raw={SI.SM12_correct.mean()*100:.1f}%")
if __name__=='__main__': main()
