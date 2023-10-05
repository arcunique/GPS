import numpy as np
from matplotlib import font_manager
from init import *

lossmech = 'photev'
snowline = {'G': 2.2, 'M': 0.8}
ibest = {'G': 104, 'M': 234}

data = get_model(models='migrationA', lossmechs='photev')

with open('../data/obs.pk', 'rb') as file:
    obsdata = pkl.load(file)

fig, ax = plt.subplots(1, 2, figsize=(20, 8), sharey=True)

for i, key in enumerate(data):
    df = data[key]['df'][ibest[key[-1]]]
    dfobs = obsdata[key[-1]]['df']
    df = df[df['fp'] <= 3]
    # dfobs = dfobs[dfobs['p'] >= df['p'].min()]
    # fpgrid = np.linspace(df['fp'].min(), df['fp'].max(), 13)
    # pgrid = 10
    # fpgrid = np.linspace(0.9, df['fp'].max(), 13)
    fpgrid = np.arange(0.8244, 2.5, 0.164)
    print(fpgrid)
    pgrid = np.logspace(np.log10(0.3), np.log10(df['p'].max()), 13)
    # imin = np.argmin(np.abs(fpgrid - f_mass_radius(0.4)))
    # fpgrid[imin] = f_mass_radius(0.4)
    # imin = np.argmin(np.abs(fpgrid - f_mass_radius(0.7)))
    # fpgrid[imin] = f_mass_radius(0.7)
    print(f_mass_radius(0.3), f_mass_radius(0.7))
    color = 'azure'
    ax[i].axhline(f_mass_radius(0.3), ls='--', lw=3, color=color)
    ax[i].axhline(f_mass_radius(0.7), ls='--', lw=3, color=color)
    wwdet, cont = show_ww_prob2d(df=df, pgrid=pgrid, y='fp', ygrid=fpgrid, ax=ax[i], cax=False, wflim=0.3, cmap='jet_r', alpha=0.6)
    brmask, bwmask, amask = bulk_comp_stat_obs(dfobs, kind=('r', 'w', 'a'), mlim=10)[0]
    ebr = ax[i].errorbar(dfobs[brmask]['p'].values, dfobs[brmask]['fp'].values, dfobs[brmask]['fpe'].values, ls='', marker='o', ms=13, ecolor='k', elinewidth=2, mew=2, markeredgecolor='k', markerfacecolor='r')
    ebw = ax[i].errorbar(dfobs[bwmask]['p'].values, dfobs[bwmask]['fp'].values, dfobs[bwmask]['fpe'].values, ls='', marker='o', ms=13, ecolor='k', elinewidth=2, mew=2, markeredgecolor='k', markerfacecolor='b')
    ebg = ax[i].errorbar(dfobs[amask]['p'].values, dfobs[amask]['fp'].values, dfobs[amask]['fpe'].values, ls='', marker='D', ms=11, color='k', elinewidth=2, mew=2, markeredgecolor='k', markerfacecolor='grey')
    xlim = ax[i].get_xlim()
    offsetfac = 1.2
    ax[i].set_xlim(xlim[0]/offsetfac, xlim[1]*offsetfac)

# font_name = font_manager.FontProperties('Rage Italic').get_name()
ax[0].text(0.4, 1.356, r'70% H$_{\rm 2}0$', size=18, color=color)
ax[0].text(1.09, 1.18, r'30%', size=18, color=color)
ax[0].set_ylim(0.65, 2.5)
cbar_ax = fig.add_axes([0.90, 0.15, 0.03, 0.7])
fig.colorbar(cont, cax=cbar_ax)
cbar_ax.set_ylabel("Probability of water worlds", size=24)
cbar_ax.yaxis.set_label_coords(2.2, 0.5)
cbar_ax.tick_params(labelsize=18)
cbar_ax.set_yticks(np.arange(0, 1.1, 0.1))
kwargs = dict(marker='o', ms=10, ecolor='k', elinewidth=2, mew=2, markeredgecolor='k')
# legend = ax[-1,1].legend([plt.scatter([], [], [], ) for size in sizes], mranges, frameon=True, title=f'Mass ({munit})', facecolor='whitesmoke', framealpha=1, fontsize=12)
leg = ax[1].legend([ebr, ebw, ebg], ['Bare rocky', 'Bare water', 'Gas-rich'], fontsize=24, frameon=True, loc='upper left', bbox_to_anchor=(-0.02, 1), handletextpad=0.1)
plt.text(0.98, 0.02, 'G-type host', size=24, ha='right', va='bottom', transform=ax[0].transAxes)
plt.text(0.98, 0.02, 'M-type host', size=24, ha='right', va='bottom', transform=ax[1].transAxes)
prettify(fig, labelsize=24, ticklabelsize=20, margin=dict(left=0.065, right=0.88, bottom=0.09, top=0.98), wspace=0.02)
ax[0].set_ylabel('$\mathscr{F}$', size=38)
# ax[0].set_facecolor('honeydew')
# ax[1].set_facecolor('gray')
fig.savefig(f'../results/wwprob2d.jpg')
plt.show()