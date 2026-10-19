import numpy as np, scipy.io as sio, matplotlib
import os
os.makedirs('figs', exist_ok=True)
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':9,'font.family':'serif','axes.grid':True,'grid.alpha':0.3})
D="C:/Users/Meera Devi Raval/Downloads/SFAS/code/"; cfgs=['icc','jstated','jfig','improved']
L={'icc':'ICC as stated','jstated':'Journal as stated','jfig':'Journal figure-implied','improved':'Improved'}
d={c:sio.loadmat(D+f'ss_{c}.mat',squeeze_me=True) for c in cfgs}
cell=lambda x,k: np.asarray(x[k]).astype(float).ravel()
# Stage-1 spectra
fig,ax=plt.subplots(1,2,figsize=(7,2.4),sharey=True)
for a,c in zip(ax,['icc','jfig']):
    s=d[c]; P=s['P1']; a.plot(s['grid1'],10*np.log10(P/P.max()),lw=1)
    for t in s['th']: a.axvline(t,color='r',ls='--',lw=.7)
    a.set_title(rf"$\alpha_c={s['ac']:g}$ ($d_c={s['ac']/2:g}\lambda$)"); a.set_xlabel(r'$\theta$ (deg)'); a.set_xlim(-60,50)
ax[0].set_ylabel('Normalized MUSIC (dB)'); fig.tight_layout(); fig.savefig('figs/stage1.pdf')
for c in cfgs:
    s=d[c]; rg=s['rg']
    # Fig.2-style
    fig,ax=plt.subplots(1,4,figsize=(7.2,2.0))
    for k in range(4):
        P=cell(s['Pr'],k); Pt=cell(s['Pr_true'],k)
        ax[k].semilogx(rg,10*np.log10(P/P.max()),lw=1,label=r'at $\hat\theta^{(1)}$')
        ax[k].semilogx(rg,10*np.log10(Pt/Pt.max()),lw=.8,ls=':',color='gray',label=r'at true $\theta$')
        ax[k].axvline(s['r'][k],color='r',ls='--',lw=.7)
        ax[k].set_title(rf"$r_{k+1}={s['r'][k]:g}\lambda$, $\hat\theta^{{(1)}}={s['theta1'][k]:.1f}^\circ$",fontsize=7.5)
        ax[k].set_xlabel(r'$r$ ($\lambda$)',fontsize=8); ax[k].tick_params(labelsize=7)
    ax[0].set_ylabel('dB'); ax[0].legend(fontsize=6,loc='lower right'); fig.tight_layout(); fig.savefig(f'figs/fig2_{c}.pdf')
    # Fig.3-style
    fig,ax=plt.subplots(1,4,figsize=(7.2,2.0))
    for k in range(4):
        P2=np.asarray(s['P2'][k],float); tl=cell(s['th_loc'],k); rl=cell(s['r_loc'],k)
        ax[k].pcolormesh(rl,tl,P2/P2.max(),shading='auto',cmap='viridis',rasterized=True)
        ax[k].plot(s['r'][k],s['th'][k],'r+',ms=7,mew=1)
        ax[k].set_title(rf"S{k+1}: $\hat\theta={s['th_hat'][k]:.2f}^\circ$, $\hat r={s['r_hat'][k]:.0f}\lambda$",fontsize=7)
        ax[k].set_xlabel(r'$r$ ($\lambda$)',fontsize=8); ax[k].tick_params(labelsize=6.5); ax[k].grid(False)
    ax[0].set_ylabel(r'$\theta$ (deg)'); fig.tight_layout(); fig.savefig(f'figs/fig3_{c}.pdf')
print('ok')
