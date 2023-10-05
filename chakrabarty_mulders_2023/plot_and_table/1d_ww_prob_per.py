import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from init import *

wflim = 0.3
mlim = 10
data = get_model(models='migrationA', lossmechs={'A': 'photev', 'B': 'all'}, startypes='GM')
# ibest = {'G': 195, 'M': 144}
ibest = {'G': 195, 'M': 249}
# ibest = '1sig'

fig, ax = plt.subplots(1, 2, figsize=(10, 6), sharey='all')
pers = np.linspace(0, 2, 200) #if ibest is None else np.linspace(0, 2, 300)
fplims = 0.8, 3.5

for key in data:
    # df = data[key]['df']
    dfs = data[key]['df']
    Pbwps = []
    Pgwps = []
    if type(ibest) == str:
        for j, df in enumerate(dfs):
            # print(j)
            masks, Ns = bulk_comp_stat_model(df, ('br', 'bw', 'gr', 'gw'), wflim=wflim, mlim=mlim)
            Ps = np.array([N/sum(Ns) for N in Ns])
            # print(Ps)
            kerns = [None for _ in masks]
            for m in range(len(masks)):
                kerns[m] = st.gaussian_kde([np.log10(df[masks[m]]['p']), df[masks[m]]['fp']])
            Pbwp = np.ones(len(pers))
            Pbwb = Pbwp.copy()
            Pgwp = Pbwp.copy()
            for i, p in enumerate(pers):
                fvals = np.array([quad(lambda fp: kern((p, fp)), *fplims)[0] for kern in kerns])
                Pbwp[i] = fvals[1] * Ps[1] / sum(fvals * Ps)
                Pbwb[i] = fvals[1] * Ps[1] / sum(fvals[:2] * Ps[:2])
                Pgwp[i] = fvals[3] * Ps[3] / sum(fvals * Ps)
                # if i==54:
                #     print(fvals*Ps / sum(fvals*Ps))
            if Pbwp[0] < 0.05:
                print(j)
                Pbwps.append(Pbwp)
                Pgwps.append(Pgwp)
        if ibest == 'minmax':
            jmaxs = np.argsort([max(Pbwp) for Pbwp in Pbwps])[[0, 1, -2, -1]]
            print(jmaxs)
            for j, jmax in enumerate(jmaxs):
                if j == 0:
                    pl, = ax[0].plot(10**pers, Pbwps[jmax]*100, lw=2, label=key[-1] + '-type host')
                else:
                    ax[0].plot(10 ** pers, Pbwps[jmax] * 100, lw=2, label=key[-1] + '-type host', c=pl.get_c())
                ax[1].plot(10**pers, Pgwps[jmax]*100, lw=2, c=pl.get_c())
        elif ibest == '1sig':
            Pbwpm = np.mean(Pbwps, axis=0)
            Pbwpe = np.std(Pbwps, axis=0)
            Pgwpm = np.mean(Pgwps, axis=0)
            Pgwpe = np.std(Pgwps, axis=0)
            pl, = ax[0].plot(10 ** pers, Pbwpm * 100, lw=2, label=key[-1] + '-type host')
            ax[0].plot(10 ** pers, (Pbwpm + Pbwpe) * 100, lw=2, c=pl.get_c())
            ax[0].plot(10 ** pers, (Pbwpm - Pbwpe) * 100, lw=2, c=pl.get_c())
            ax[1].plot(10 ** pers, Pgwpm * 100, lw=2, c=pl.get_c())
            ax[1].plot(10 ** pers, (Pgwpm + Pgwpe) * 100, lw=2, c=pl.get_c())
            ax[1].plot(10 ** pers, (Pgwpm - Pgwpe) * 100, lw=2, c=pl.get_c())
    else:
        df = dfs[ibest[key[-1]]]
        masks, Ns = bulk_comp_stat_model(df, ('br', 'bw', 'gr', 'gw'), wflim=wflim, mlim=mlim)
        Ps = np.array([N/sum(Ns) for N in Ns])
        # print(Ps)
        kerns = [None for _ in masks]
        for m in range(len(masks)):
            kerns[m] = st.gaussian_kde([np.log10(df[masks[m]]['p']), df[masks[m]]['fp']])
        Pbwp = np.ones(len(pers))
        Pbwb = Pbwp.copy()
        Pgwp = Pbwp.copy()
        for i, p in enumerate(pers):
            fvals = np.array([quad(lambda fp: kern((p, fp)), *fplims)[0] for kern in kerns])
            Pbwp[i] = fvals[1] * Ps[1] / sum(fvals * Ps)
            Pbwb[i] = fvals[1] * Ps[1] / sum(fvals[:2] * Ps[:2])
            Pgwp[i] = fvals[3] * Ps[3] / sum(fvals * Ps)
        p, = ax[0].plot(10 ** pers, Pbwp * 100, lw=2, label=key[-1] + '-type host')
        # ax[0].plot(10 ** pers, Pbwb * 100, lw=2, ls='--', c=p.get_c())
        ax[1].plot(10 ** pers, Pgwp * 100, lw=2)

log_axis(ax[0])
log_axis(ax[1])
ax[0].set_ylabel('Occurrence among all planets (%)')
ax[0].yaxis.set_label_coords(-0.13, 0.5)
ax[0].set_xlabel('Period (days)')
ax[1].set_xlabel('Period (days)')
ax[0].set_title('Bare water planets', size=20, y=0.92)
ax[1].set_title('Gas-rich water planets', size=20, y=0.92)
ax[0].annotate('', xy=(60, 32), xytext=(4, 32), arrowprops=dict(arrowstyle='<->', linewidth=2))
ax[0].text(14, 34, 'Find more BWPs (RV) or\nVerify absence of\n atm of BWPs (JWST)', size=20, ha='center', va='bottom')
ax[1].annotate('', xy=(90, 40), xytext=(12, 40), arrowprops=dict(arrowstyle='->', linewidth=2))
ax[1].text(11, 40, 'Detect\n atmospheric\n traces of water\nworlds (JWST)', size=20, ha='right', va='center')
legend_colortxt_nosym(ax[0], frameon=False, fontsize=20, bbox_to_anchor=(0.5, 0.89), labelspacing=0.3)
prettify(fig, minorlocators={'x': None, 'y': 5}, margin=dict(left=0.1, right=0.99, bottom=0.1, top=0.99), wspace=0.05, hspace=0.03)
fig.savefig('../results/wwprob1d.jpg')
plt.show()

