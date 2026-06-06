#!/usr/bin/env python3
"""
GERT azeotropes — Paper 2, Part 1, Figure 2.
Level-2 composition profiles phi(x): the interior minimum that marks a
composition azeotrope, generated WITHOUT new parameters from the Margules
identity W = -R*T*ln(gamma_inf) (same gamma_inf as Level 1).

  f_M(x) = x*dH_w + (1-x)*dH_s + W*x(1-x)
  f_L(x) = T*[x*dS_w + (1-x)*dS_s],  dS = dH/Tb   (Trouton)
  phi(x) = f_M/(f_M+f_L);  azeotrope <-> interior extremum (dphi/dx=0)
"""
import numpy as np
import matplotlib.pyplot as plt

R=8.314; C=3.729; T=358.0
# (dHvap kJ/mol, Tb K, gamma_inf in water)
W_={'water':(40.7,373.15,None)}
sol={'methanol':(35.2,337.70,1.7),'1-propanol':(41.4,370.35,10.0),
     '2-propanol':(39.9,355.41,8.0),'allyl alcohol':(40.0,370.2,7.5)}

def profile(s,n=800):
    hA,TbA,_=W_['water']; hB,TbB,g=sol[s]
    dSA,dSB=hA/TbA,hB/TbB
    Wm=-R*T*np.log(g)/1000.0          # kJ/mol, Margules identity
    x=np.linspace(0.001,0.999,n)       # x = mole fraction of water
    fM=x*hA+(1-x)*hB+Wm*x*(1-x); fL=T*(x*dSA+(1-x)*dSB)
    phi=fM/(fM+fL)
    dp=np.diff(phi)/np.diff(x)
    mins=[(x[i+1],phi[i+1]) for i in range(len(dp)-1) if dp[i]<0 and dp[i+1]>0]
    return x,phi,Wm,mins

cols={'methanol':'#666666','2-propanol':'#d95f02','allyl alcohol':'#7570b3','1-propanol':'#e7298a'}
order=['methanol','2-propanol','allyl alcohol','1-propanol']

plt.rcParams.update({'font.size':11,'font.family':'DejaVu Sans','axes.linewidth':0.9})
fig,ax=plt.subplots(figsize=(8.2,5.4))
for s in order:
    x,phi,Wm,mins=profile(s)
    lab=f"{s}  (γ$^\\infty$={sol[s][2]:g}, W={Wm:.1f})"
    style='--' if not mins else '-'
    ax.plot(x,phi,style,color=cols[s],lw=1.9,label=lab,zorder=4)
    for xm,pm in mins:
        ax.plot(xm,pm,'o',color=cols[s],ms=9,mfc='white',mec=cols[s],mew=1.8,zorder=6)
        ax.annotate(f'min  x={xm:.2f}',(xm,pm),textcoords='offset points',
                    xytext=(6,-12),fontsize=8.5,color=cols[s])

ax.set_xlabel('x  (mole fraction of water)',fontsize=12)
ax.set_ylabel('φ(x)  =  f$_M$(x) / [f$_M$(x) + f$_L$(x)]',fontsize=12)
ax.legend(loc='upper center',fontsize=9,framealpha=0.95)
ax.text(0.015,0.02,'methanol (dashed): monotonic → no interior minimum → no azeotrope (correct)',
        transform=ax.transAxes,fontsize=8.6,color='#444',style='italic')
plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/fig2_phi_profiles.png',dpi=300,bbox_inches='tight')
print("saved fig2_phi_profiles.png")
for s in order:
    x,phi,Wm,mins=profile(s)
    loc=', '.join(f'x={m[0]:.2f}' for m in mins) or 'monotonic'
    print(f"  {s:14s} W={Wm:6.2f}  {loc}  range φ=[{phi.min():.3f},{phi.max():.3f}]")
