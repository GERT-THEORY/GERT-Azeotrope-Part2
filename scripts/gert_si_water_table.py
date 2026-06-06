#!/usr/bin/env python3
"""Generate the curated water-domain SI table (existence set behind Fig. 3) with
phi, GERT verdict and observed behaviour. Output: GERT_Paper2_SI_water_curated.csv
NOTE: the complete 91-system enumeration additionally requires the
infinite-dilution activity-coefficient dataset of Brouwer et al. 2021
(Sep. Purif. Technol. 272, 118727); only the curated set is reproduced here."""
import pandas as pd
from gert_water_level1 import fM02, fM12, fLval, GAP, T02, T12
rows=[]
for n,b,dp,g,lp,o in T02:
    fM=fM02(b,dp); fL=fLval(g,lp,n in GAP); phi=fM/(fM+fL) if (fM+fL)>0 else 0
    rows.append({'system':n,'type':'0+2','gamma_inf':g,'f_M':round(fM,3),'f_L':round(fL,3),
                 'phi':round(phi,3),'GERT':'none' if phi>0.5 else 'positive','observed':o})
for n,a,b,g,lp,o in T12:
    fM=fM12(a,b); fL=fLval(g,lp,n in GAP); phi=fM/(fM+fL) if (fM+fL)>0 else 0
    pred='negative' if o=='negative' else ('none' if phi>0.5 else 'positive')
    rows.append({'system':n,'type':'1+2','gamma_inf':g,'f_M':round(fM,3),'f_L':round(fL,3),
                 'phi':round(phi,3),'GERT':pred,'observed':o})
pd.DataFrame(rows).to_csv('GERT_Paper2_SI_water_curated.csv',index=False)
print(f"curated water systems: {len(rows)}")
if __name__=='__main__': pass
