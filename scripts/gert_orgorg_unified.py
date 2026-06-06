#!/usr/bin/env python3
"""
GERT organic-organic azeotrope law -- canonical reproducible script (Paper 2).

One criterion (phi = 1/2), forces identified per domain, with two sub-laws:
  - Sub-law A (dipolar aprotic): nitro/nitrile alpha -> 0 (protogenic, no bulk donation)
  - Sub-law B (miscibility saturation): deviation capped at spinodal ln gamma_inf = 2

All constants are independent physics; none fitted to azeotrope outcomes.
Descriptors: _kt.py (Kamlet-Taft alpha,beta, audited), _hansen.py (delta_d,delta_p,Vm),
_tb.py (boiling points). Reproduces the Paper 2 numbers:
  GERT raw 73.1% / balanced 73.8% / recall none-pos-neg 73-74-75.
"""
import pandas as pd, numpy as np
from _kt import AB
from _tb import Tb
from _hansen import HAN

DATA = "/mnt/user-data/uploads/ie9b03694_si_002.xlsx"
R = 8.314
TROUTON = 88.0          # J/mol/K, Trouton's rule
E_HB = 25_000.0         # J/mol, hydrogen-bond energy scale (robust 20-30 kJ/mol)
SPINODAL = 2.0          # ln gamma_inf cap, regular-solution critical point (Sub-law B)
A_NA, B_NA = 0.10, 0.25 # non-associating structural boundary (not fitted)

# Sub-law A: protogenic solvents whose tabulated alpha is spurious for bulk mixing
APROTIC = {'nitromethane','nitroethane','1-nitropropane','2-nitropropane','nitrobenzene',
           '1-methyl-2-nitrobenzene','nitrile acetic acid','nitrile propanoic acid',
           'nitrile butanoic acid','nitrile benzoic acid'}

nrm = lambda x: str(x).strip().lower()

def truth(c):
    c = str(c).strip()
    if '/' in c: return 'amb'
    if c in ('Z','z','Z2','ZC'): return 'none'
    if c.startswith(('AO','AE','AX')) or c == 'A':
        return 'neg' if 'N' in c else 'pos'
    return 'amb'

def classify(a1, b1, a2, b2, dd1, dp1, V1, dd2, dp2, V2, T1, T2, n1, n2):
    # Sub-law A: route protogenic solvents to the dispersion channel
    if n1 in APROTIC: a1 = 0.0
    if n2 in APROTIC: a2 = 0.0
    Tm = (T1 + T2) / 2.0
    Vm = (V1 + V2) / 2.0
    lnr = (TROUTON / R) * abs(T2 - T1) / Tm          # f_L: volatility gate (Eq. 8)
    dC  = (a1 - a2) * (b2 - b1)                       # Delta C = X - S (Eq. 3)
    gdisp = Vm * ((dd1 - dd2)**2 + (dp1 - dp2)**2) / (R * Tm)   # dispersion channel (Eq. 6)
    NA = (a1 < A_NA and a2 < A_NA and b1 < B_NA and b2 < B_NA)  # non-associating (Eq. 7)
    # f_M = channel deviation, capped at the spinodal (Sub-law B, Eq. 11)
    dev = min(gdisp, SPINODAL) if NA else min((E_HB / (R * Tm)) * abs(dC), SPINODAL)
    if lnr >= dev:                                    # existence: |ln r| >= |W| -> zeotropic
        return 'none'
    return 'pos' if (NA or dC < 0) else 'neg'         # type: Eq. 5 / Eq. 12

def main():
    sv = pd.read_excel(DATA, sheet_name="Exp-vs-SMx")
    sv['obs'] = sv.Experimental.map(truth)
    sv['SM12'] = sv['SM12 Prediction'].map(lambda c: {'Z':'none','X':'pos','N':'neg'}.get(str(c).strip(),'amb'))
    sv = sv[~((sv['Component 1'].map(nrm) == 'water') | (sv['Component 2'].map(nrm) == 'water'))]
    sv = sv[sv.obs.isin(['none','pos','neg'])].copy()

    preds, obs, sm12 = [], [], []
    for _, r in sv.iterrows():
        n1, n2 = nrm(r['Component 1']), nrm(r['Component 2'])
        ab1, ab2 = AB.get(n1), AB.get(n2)
        h1, h2 = HAN.get(n1), HAN.get(n2)
        t1, t2 = Tb.get(n1), Tb.get(n2)
        if None in (ab1, ab2, h1, h2, t1, t2):
            continue
        preds.append(classify(ab1[0], ab1[1], ab2[0], ab2[1],
                              h1[0], h1[1], h1[2], h2[0], h2[1], h2[2], t1, t2, n1, n2))
        obs.append(r.obs); sm12.append(r.SM12)
    D = pd.DataFrame({'obs': obs, 'pred': preds, 'SM12': sm12})

    def rep(col, label, frame):
        f = frame[frame[col].isin(['none','pos','neg'])]
        raw = (f.obs == f[col]).mean() * 100
        rec = {c: (f[f.obs == c][col] == c).mean() * 100 for c in ['none','pos','neg']}
        bal = np.mean(list(rec.values()))
        prec = {c: ((f[f[col] == c].obs == c).mean() * 100 if (f[col] == c).any() else 0)
                for c in ['pos','neg']}
        print(f"{label:14s} raw={raw:.1f}%  balanced={bal:.1f}%  "
              f"recall none/pos/neg={rec['none']:.0f}/{rec['pos']:.0f}/{rec['neg']:.0f}  "
              f"prec pos/neg={prec['pos']:.0f}/{prec['neg']:.0f}")

    print(f"n={len(D)}  none/pos/neg = {sum(D.obs=='none')}/{sum(D.obs=='pos')}/{sum(D.obs=='neg')}")
    rep('pred', 'GERT (1 law)', D)
    rep('SM12', 'SM12', D)

if __name__ == "__main__":
    main()
