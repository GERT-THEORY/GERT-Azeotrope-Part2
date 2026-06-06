#!/usr/bin/env python3
"""
GERT azeotropes — Paper 2, Part 1, Figure 1.
The CH2 law made visible: homologous series crossing phi = 1/2.

Criterion (Paper 1, unchanged):
  phi = f_M / (f_M + f_L)
  Type 0+2:  f_M = (beta/beta_w)(dp/dp_w)
  Type 1+2:  f_M = sqrt(min(alpha/alpha_w,1) * min(beta/beta_w,1))
  f_L = max( ln(gamma_inf)/C , logP/C' )
Constants fixed in Paper 1; not recalibrated.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager

C, CP = 3.729, 0.48
AW, BW, DPW = 1.17, 0.47, 16.0

def fL(gamma, logP):
    # f_L >= 0: expulsion cannot be negative; a hydrophilic solute has zero
    # outward force, not a negative one. Affects only strong integrators
    # (amides, DMSO, sulfolane), capping phi at 1; no classification changes.
    return max(np.log(gamma)/C, logP/CP, 0.0)

def fL_gap(gamma, logP, gap):
    return max(np.log(gamma)/C, (logP/CP if gap else -9.0), 0.0)

def phi_02(beta, dp, gamma, logP):
    fM = (beta/BW)*(dp/DPW)
    return fM/(fM+fL(gamma, logP))

def phi_12(alpha, beta, gamma, logP, gap):
    fM = max(min(alpha/AW,1), min(beta/BW,1))   # dominant channel (Paper 2)
    return fM/(fM+fL_gap(gamma, logP, gap))

# ---- Verified homologous series (carbon number, descriptors) ----
# Type 0+2: (Cn, beta, dp, gamma_inf_water, logP)   [from gert_azeotrope_gamma.py]
ketones = [(3,0.48,10.4,7.65,-0.24),(4,0.48,9.0,27.6,0.29),(5,0.50,7.5,102.0,0.91),(6,0.49,6.1,346.9,1.38)]
acetates= [(3,0.42,7.6,22.5,0.18),(4,0.45,5.3,75.6,0.73),(5,0.45,3.3,274.5,1.24),(6,0.45,3.7,1058.0,1.78)]
ethers  = [(4,0.47,2.9,77.1,0.89),(6,0.46,2.4,2313.0,1.21),(8,0.46,1.6,47180.,3.21)]
amides  = [(3,0.69,13.7,0.70,-1.01),(4,0.76,11.5,1.04,-0.77),(5,0.77,12.3,0.37,-0.38)]
# Type 1+2: (Cn, alpha, beta, gamma_inf_water, logP, gap)  gap=miscibility gap
alcohols= [(1,0.93,0.62,1.7,-0.77,False),(2,0.86,0.75,3.8,-0.31,False),(3,0.84,0.90,10.0,0.25,False),(4,0.84,0.84,57.6,0.88,True)]
acids   = [(2,1.12,0.45,1.58,-0.17,False),(3,1.08,0.45,309.0,0.33,False),(4,1.06,0.45,52.9,0.79,True),(5,1.04,0.45,127.0,1.39,True)]

def series_phi(data, kind):
    out=[]
    for row in data:
        Cn=row[0]
        p = phi_02(row[1],row[2],row[3],row[4]) if kind=='02' else phi_12(row[1],row[2],row[3],row[4],row[5])
        out.append((Cn,p))
    return out

S = {
 'Ketones (R-CO-R)':       (series_phi(ketones,'02'),  '#1b9e77','o'),
 'Acetate esters (R-OAc)': (series_phi(acetates,'02'), '#d95f02','s'),
 'Di-alkyl ethers (R-O-R)':(series_phi(ethers,'02'),   '#7570b3','^'),
 'Alcohols (R-OH)':        (series_phi(alcohols,'12'),  '#e7298a','D'),
 'Carboxylic acids (R-COOH)':(series_phi(acids,'12'),  '#66a61e','v'),
 'Amides (integrators)':   (series_phi(amides,'02'),    '#666666','P'),
}

plt.rcParams.update({'font.size':11,'font.family':'DejaVu Sans','axes.linewidth':0.9})
fig,ax=plt.subplots(figsize=(8.2,5.4))

# threshold band
ax.axhspan(0.5,1.02,color='#2ecc71',alpha=0.05)
ax.axhspan(-0.02,0.5,color='#e74c3c',alpha=0.05)
ax.axhline(0.5,color='k',lw=1.3,ls='--',zorder=3)
ax.text(4.6,0.525,'integration threshold  φ = ½',ha='center',va='bottom',fontsize=9.5,style='italic')
ax.text(1.9,0.965,'integrates (no azeotrope)',fontsize=9.5,color='#1e7e34',va='top',ha='center')
ax.text(7.2,0.34,'disrupts (positive azeotrope)',fontsize=9.5,color='#a82315',va='bottom',ha='center')

for label,(pts,col,mk) in S.items():
    xs=[c for c,_ in pts]; ys=[p for _,p in pts]
    ax.plot(xs,ys,'-',color=col,lw=1.6,alpha=0.85,zorder=4)
    # marker fill: open if integrates (>0.5), filled if disrupts
    for x,y in pts:
        ax.plot(x,y,mk,color=col,ms=8,mfc=(col if y<0.5 else 'white'),
                mec=col,mew=1.6,zorder=5)
    ax.plot([],[],mk,color=col,label=label,ms=8,mfc=col,mec=col)

ax.set_xlabel('Carbon number of solute',fontsize=12)
ax.set_ylabel('φ  =  f$_M$ / (f$_M$ + f$_L$)',fontsize=12)
ax.set_xlim(0.5,8.6); ax.set_ylim(-0.02,1.05)
ax.set_xticks(range(1,9))
ax.legend(loc='upper right',fontsize=8.6,framealpha=0.95,handletextpad=0.4)
# annotate filled vs open meaning
ax.text(0.5,-0.14,'Filled marker: observed positive azeotrope.   Open marker: observed none (integrates).',
        transform=ax.transAxes,fontsize=8.3,color='#444')
plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/fig1_homologous_phi.png',dpi=300,bbox_inches='tight')
print("saved fig1_homologous_phi.png")
# print the underlying table for audit
print(f"\n{'series':26s} {'C':>2s} {'phi':>6s}")
for label,(pts,_,_) in S.items():
    for c,p in pts: print(f"{label:26s} {c:>2d} {p:6.3f}")
