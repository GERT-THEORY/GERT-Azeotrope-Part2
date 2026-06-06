#!/usr/bin/env python3
"""
GERT azeotropes - Paper 2, Part 1, Figure 3 + Level-1 recomputation (REFINED).
Classification plane (f_L, f_M): separatrix phi=1/2 is the diagonal f_M=f_L.
Above -> integrates (none); below -> disrupts (positive). Points coloured by
OBSERVED behaviour; the diagonal shows the PREDICTION.

Refined criterion (Paper 2):
  f_M^(1+2) = max(f_donor, f_acceptor)            [dominant channel; generalises
              f_donor=min(a/AW,1), f_acc=min(b/BW,1) the Paper 1 geometric mean]
  f_L = max( ln(gamma)/C , logP/C' |_gap , 0 )    [logP only if miscibility gap]
The three aqueous amines now integrate (above the line); only the Type-3
negative-deviation case (formic acid) lies off the positive/none axis.
"""
import numpy as np
import matplotlib.pyplot as plt
C, CP, AW, BW, DPW = 3.729, 0.48, 1.17, 0.47, 16.0
def fM02(b, dp):    return (b/BW)*(dp/DPW)
def fM12(a, b):     return max(min(a/AW,1), min(b/BW,1))          # dominant channel
def fLval(g, logP, gap):
    return max(np.log(g)/C, (logP/CP if gap else -9.0), 0.0)     # logP only if gap

# partially-miscible / heteroazeotrope systems (logP term applies)
GAP={"2-Pentanone","3-Pentanone","Cyclopentanone","Cyclohexanone","2-Hexanone",
"Ethyl formate","Ethyl acetate","Propyl acetate","Butyl acetate","Tetrahydropyran",
"Diethyl ether","MTBE","DiisoPr ether","Di-n-propyl ether","Di-n-butyl ether","Anisole",
"Benzonitrile","Butyronitrile","Propanal","Butanal","DCM","Chloroform","1,2-DCE",
"Triethylamine","Thiophene","MEK","1-Butanol","Butanoic acid","Pentanoic acid",
"Furfural","Nitroethane","Methyl formate"}

T02=[("Acetone",0.48,10.4,7.65,-0.24,"none"),("MEK",0.48,9.0,27.6,0.29,"positive"),
("2-Pentanone",0.50,7.5,102,0.91,"positive"),("3-Pentanone",0.45,7.6,113,0.82,"positive"),
("Cyclopentanone",0.52,11.9,16.1,0.40,"positive"),("Cyclohexanone",0.53,8.4,5.14,0.80,"positive"),
("2-Hexanone",0.49,6.1,346.9,1.38,"positive"),
("Methyl formate",0.37,8.4,15.75,0.03,"positive"),("Ethyl formate",0.36,7.2,46.65,0.23,"positive"),
("Methyl acetate",0.42,7.6,22.5,0.18,"positive"),("Ethyl acetate",0.45,5.3,75.6,0.73,"positive"),
("Propyl acetate",0.45,3.3,274.5,1.24,"positive"),("Butyl acetate",0.45,3.7,1058,1.78,"positive"),
("THF",0.55,5.7,17.0,0.46,"positive"),("Tetrahydropyran",0.54,4.5,70.5,0.82,"positive"),
("1,3-Dioxolane",0.45,6.6,9.71,-0.37,"positive"),("1,4-Dioxane",0.37,1.8,5.44,-0.27,"positive"),
("Diethyl ether",0.47,2.9,77.1,0.89,"positive"),("MTBE",0.55,4.3,113,0.94,"positive"),
("DiisoPr ether",0.49,3.4,628,1.52,"positive"),("Di-n-propyl ether",0.46,2.4,2313,1.21,"positive"),
("Di-n-butyl ether",0.46,1.6,47180,3.21,"positive"),("Anisole",0.22,4.1,4000,2.11,"positive"),
("NMP",0.77,12.3,0.37,-0.38,"none"),("DMF",0.69,13.7,0.70,-1.01,"none"),("DMAc",0.76,11.5,1.04,-0.77,"none"),
("Sulfolane",0.39,17.4,2.0,-0.77,"none"),
("Propionitrile",0.37,16.1,36.4,0.16,"positive"),("Benzonitrile",0.37,8.4,1741,1.56,"positive"),
("Butyronitrile",0.45,12.4,115.16,0.53,"positive"),
("Pyridine",0.64,8.8,31.8,0.65,"positive"),("3-MePyridine",0.68,7.8,49.1,1.20,"positive"),
("4-MePyridine",0.67,7.8,42.3,1.22,"positive"),
("Acetaldehyde",0.45,11.3,4.15,-0.34,"none"),("Furfural",0.34,14.9,97.4,0.41,"positive"),
("Propanal",0.41,11.1,24.7,0.59,"positive"),("Butanal",0.41,10.1,37.95,0.88,"positive"),
("DCM",0.00,6.3,250,1.25,"positive"),("Chloroform",0.00,3.1,835,1.97,"positive"),("1,2-DCE",0.00,10.4,641,1.48,"positive"),
("DMSO",0.76,16.4,1.3,-1.35,"none"),("Triethylamine",0.71,0.4,67.5,1.45,"positive"),
("Nitromethane",0.06,18.8,31.6,-0.33,"positive"),("Nitroethane",0.25,15.5,75.0,0.18,"positive"),
("Thiophene",0.15,2.4,1272,1.81,"positive")]
T12=[("Methanol",0.93,0.62,1.7,-0.77,"none"),("Ethanol",0.86,0.75,3.8,-0.31,"positive"),
("1-Propanol",0.84,0.90,10,0.25,"positive"),("2-Propanol",0.76,0.84,8,0.05,"positive"),
("1-Butanol",0.84,0.84,57.6,0.88,"positive"),("Allyl alcohol",0.84,0.90,7.5,0.17,"positive"),
("Acetic acid",1.12,0.45,1.58,-0.17,"none"),("Propanoic acid",1.08,0.45,309,0.33,"positive"),
("Butanoic acid",1.06,0.45,52.9,0.79,"positive"),("Pentanoic acid",1.04,0.45,127,1.39,"positive"),
("Formic acid",1.23,0.38,0.68,-0.54,"negative"),
("Pyrrolidine",0.16,0.70,1.5,0.46,"none"),("Diethylamine",0.08,0.70,5.0,0.58,"none"),
("1-Propanamine",0.10,0.61,2.5,0.48,"none")]

pts=[]
for n,b,dp,g,lp,o in T02: pts.append((n,fM02(b,dp),fLval(g,lp,n in GAP),o))
for n,a,b,g,lp,o in T12: pts.append((n,fM12(a,b),fLval(g,lp,n in GAP),o))

col={'none':'#2ca25f','positive':'#de2d26','negative':'#3182bd','n2':'#e6a817'}
N2={"Ethanol","1-Propanol","2-Propanol","Allyl alcohol"}  # integrate at existence; azeotrope is composition (Fig 2)
plt.rcParams.update({'font.size':11,'font.family':'DejaVu Sans','axes.linewidth':0.9})
fig,ax=plt.subplots(figsize=(7.8,6.4))
lim=2.6
ax.fill_between([0,lim],[0,lim],lim,color='#2ecc71',alpha=0.06)
ax.fill_between([0,lim],0,[0,lim],color='#e74c3c',alpha=0.06)
ax.plot([0,lim],[0,lim],'k--',lw=1.3,zorder=3)
ax.text(2.1,2.25,'phi = 1/2  (f$_M$ = f$_L$)',rotation=45,fontsize=9.5,style='italic',ha='center')
ax.text(0.12,2.35,'integrates\n(none)',color='#1e7e34',fontsize=10,va='top')
ax.text(1.9,0.12,'disrupts\n(positive)',color='#a82315',fontsize=10,va='bottom',ha='center')

lab_map={'none':'none','positive':'positive','negative':'Type-3 (negative)','n2':'composition azeotrope (N2)'}
seen=set()
for n,fM,fL,o in pts:
    cat = 'n2' if n in N2 else o
    lab=lab_map[cat] if cat not in seen else None; seen.add(cat)
    amine = n in ("Pyrrolidine","Diethylamine","1-Propanamine")
    mk = 'D' if amine else ('s' if cat=='n2' else 'o')
    ax.plot(fL,fM,mk,ms=9 if (amine or cat=='n2') else 8,color=col[cat],mfc=col[cat],
            mec='black' if amine else 'white',mew=1.4 if amine else 0.5,zorder=6 if amine else 5,label=lab)
for n,fM,fL,o in pts:
    off={"Pyrrolidine":(6,9),"1-Propanamine":(2,-15),"Diethylamine":(8,9),
         "Acetone":(6,-12),"Sulfolane":(-50,4),"Di-n-butyl ether":(6,4),
         "Chloroform":(6,4),"NMP":(6,4)}
    if n in off:
        ax.annotate(n,(fL,fM),textcoords='offset points',xytext=off[n],fontsize=8,color='#333')

ax.set_xlim(0,lim); ax.set_ylim(0,lim)
ax.set_xlabel('f$_L$  =  max[ ln(gamma$^\\infty$)/C , logP/C$\'$ |$_{gap}$ , 0 ]',fontsize=12)
ax.set_ylabel('f$_M$  (cohesive integration, dominant channel)',fontsize=12)
ax.legend(loc='upper right',fontsize=9.0,framealpha=0.95,title='observed')
plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/fig3_classification_plane.png',dpi=300,bbox_inches='tight')
print("saved fig3_classification_plane.png")
ok=tot=0; wrong=[]
for n,fM,fL,o in pts:
    if o=='negative': continue
    pred='none' if fM>fL else 'positive'
    truth = 'none' if (o=='none' or n in N2) else 'positive'   # N2 integrate at existence
    tot+=1
    if pred==truth: ok+=1
    else: wrong.append(n)
print(f"Existence-level separation (excl. Type-3 negative; N2 scored as integrating): {ok}/{tot}.")
print("Misclassified at existence:", wrong if wrong else "NONE")
