import numpy as np, scipy.io as sio
D='/home/claude/sfas/'; names=[('icc','ICC as stated'),('jstated','Journal as stated'),('jfig','Journal fig.-implied'),('improved','Improved')]
rows=[]
for c,n in names:
    m=sio.loadmat(D+f'mc_{c}.mat'); E1,E2,Er=m['E1'],m['E2'],m['Er']; ok=~np.isnan(E1[:,0])
    E1,E2,Er=E1[ok],E2[ok],Er[ok]; T=int(m['t'].item())
    f=lambda v,fmt: ' & '.join(fmt%x for x in v)
    rows.append(f"{n} ({ok.sum()}/{T}) & {f(np.sqrt((E1**2).mean(0)),'%.2f')} & {f(np.sqrt((E2**2).mean(0)),'%.3f')} & "
                f"{f(100*np.median(abs(Er),0),'%.1f')} & {f(100*(abs(Er)<0.05).mean(0),'%.0f')} \\\\")
tex=r"""\begin{table}[h]\centering\scriptsize
\caption{Monte-Carlo summary (SNR 10\,dB, $N=500$; trials with all four Stage-1 peaks resolved / total). Columns
per source: $30\lambda$, $300\lambda$, $1000\lambda$, $5000\lambda$.}\label{tab:mc}
\setlength{\tabcolsep}{2.6pt}
\begin{tabular}{@{}l|cccc|cccc|cccc|cccc@{}}\toprule
 & \multicolumn{4}{c|}{Stage-1 DOA RMSE ($^\circ$)} & \multicolumn{4}{c|}{Final DOA RMSE ($^\circ$)} &
 \multicolumn{4}{c|}{Median $|\Delta r|/r$ (\%)} & \multicolumn{4}{c}{$|\Delta r|/r<5\%$ (\% trials)}\\\midrule
""" + "\n".join(rows) + r"""
\bottomrule\end{tabular}\end{table}"""
open('mc_table.tex','w').write(tex); print(tex)
