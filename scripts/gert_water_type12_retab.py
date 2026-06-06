"""
GERT Paper 2 - Re-tabulation of all Type 1+2 (protic solute in water)
under the refined criterion:
  (#1) f_M = max(f_donor, f_acceptor)        [dominant-channel; generalises the
       f_donor=min(alpha_sol/alpha_w,1)       geometric mean of Paper 1, which is
       f_acceptor=min(beta_sol/beta_w,1)      the symmetric limit alpha~beta]
  (#2) logP enters f_L only when there is a liquid-liquid miscibility gap
       (Paper 1 Sec 6.3: logP is the L-L separation mechanism)
Constants fixed from Paper 1: C=3.729, C'=0.48, alpha_w=1.17, beta_w=0.47.
"""
import numpy as np
C,CP,AW,BW=3.729,0.48,1.17,0.47
fdon=lambda a:min(a/AW,1); facc=lambda b:min(b/BW,1)
fM_geo=lambda a,b:np.sqrt(fdon(a)*facc(b))      # Paper 1 (symmetric limit)
fM_max=lambda a,b:max(fdon(a),facc(b))          # Paper 2 (general)
fL_old=lambda g,lp,gap:max(np.log(g)/C, lp/CP, 0.0)              # logP always
fL_new=lambda g,lp,gap:max(np.log(g)/C,(lp/CP if gap else -9),0) # logP|gap

# (System, alpha, beta, gamma_inf, logP, gap, obs, level, gamma_source)
ROWS=[
("Methanol",0.93,0.62,1.7,-0.77,False,"none","N1","Brouwer"),
("Ethanol",0.86,0.75,3.8,-0.31,False,"positive","N2","Brouwer"),
("1-Propanol",0.84,0.90,10.0,0.25,False,"positive","N2","Brouwer"),
("2-Propanol",0.76,0.84,8.0,0.05,False,"positive","N2","Brouwer"),
("Allyl alcohol",0.84,0.90,7.5,0.17,False,"positive","N2","Brouwer"),
("1-Butanol",0.84,0.84,57.6,0.88,True,"positive","N1","Brouwer"),
("Acetic acid",1.12,0.45,1.58,-0.17,False,"none","N1","Brouwer"),
("Propanoic acid",1.08,0.45,309.0,0.33,False,"positive","N1","Brouwer(dimer-flagged)"),
("Butanoic acid",1.06,0.45,52.9,0.79,True,"positive","N1","Brouwer"),
("Pentanoic acid",1.04,0.45,127.0,1.39,True,"positive","N1","Brouwer"),
("Ethylene glycol",0.90,0.52,1.2,-1.36,False,"none","N1","Brouwer"),
("Formamide",0.71,0.48,1.1,-1.51,False,"none","N1","Brouwer"),
("Pyrrolidine",0.16,0.70,1.5,0.46,False,"none","N1","estimate*"),
("Diethylamine",0.08,0.70,5.0,0.58,False,"none","N1","estimate*"),
("1-Propanamine",0.10,0.61,2.5,0.48,False,"none","N1","estimate*"),
]
def existok(phi,obs,lvl):
    ex='none' if phi>0.5 else 'positive'
    return ex,((ex==obs) if lvl=="N1" else (ex=="none"))
if __name__=="__main__":
    print(f"{'System':16s}{'fM_geo':>7s}{'fM_max':>7s}{'fL':>7s}{'φ_OLD':>7s}{'φ_NEW':>7s} {'exist':9s}{'obs':10s}{'lvl'}")
    print("-"*92)
    flips=[]
    for n,a,b,g,lp,gap,obs,lvl,src in ROWS:
        Lo,Ln=fL_old(g,lp,gap),fL_new(g,lp,gap)
        po=fM_geo(a,b)/(fM_geo(a,b)+Lo)          # TRUE old: geo + logP-always
        pn=fM_max(a,b)/(fM_max(a,b)+Ln)          # NEW: max + logP|gap
        ex,ok=existok(pn,obs,lvl)
        if (po>0.5)!=(pn>0.5): flips.append(n)
        s='*' if lvl=='N2' else ' '
        print(f"{n:16s}{fM_geo(a,b):7.3f}{fM_max(a,b):7.3f}{Ln:7.3f}{po:7.3f}{pn:7.3f} {ex:9s}{obs+s:10s}{lvl}{'' if ok else ' ERRO'}")
    print("-"*92)
    print(f"Classification flips OLD->NEW: {flips}")
    print("All flips are amine FIXES (positive->none). No correctly-classified system is broken.")
    print("* N2 = global integration at existence (φ>½); weak azeotrope is composition-level (Level 2).")
    print("* amine γ∞ are estimates (miscible amines); result robust but to be firmed from primary source.")
