import matplotlib.pyplot as plt

from init import *

mlim = 20
wflim = 0.3
ibest = {'migrationA': {'G': 104, 'M': 234}, 'migrationB': None}
startype = 'GM'

data = get_model(models=('migrationA',), lossmechs={'A': 'photev', 'B': 'all'}, startypes=startype)
# print(data.keys())

with open('../data/obs.pk', 'rb') as file:
    obsdata = pkl.load(file)

ibest = {
    'migrationA_photev_G': 80, 'migrationA_photev_M': 14, 'migrationB_photev_G': 18,
    'migrationB_photev_M': 18, 'migrationB_impact_G': 12, 'migrationB_impact_M': 13
}
print(len(data))

fig, ax = plt.subplots(2, 1, figsize=(9, 13), sharey='all', sharex='all')
ax = ax.ravel()
# ax = [ax]
for i, key in enumerate(data):
    df = data[key]['df'][ibest[key]]
    dfobs = obsdata[key[-1]]['df']
    brmask, grmask = bulk_comp_stat_model(df, ('br', 'gr'), mlim=mlim, wflim=wflim)[0]
    bwmask, gwmask = bulk_comp_stat_model(df, ('bw', 'gw'), mlim=mlim, wflim=wflim)[0]
    brmask &= df['fp'] < f_mass_radius(0.2)
    grmask &= (df['fp'] > f_mass_radius(0.8)) #| (df['m'] > mlim)
    gwmask &= (df['fp'] > f_mass_radius(0.8)) #| (df['m'] > mlim)
    bwmask &= df['fp'] < f_mass_radius(0.7)
    for m, (mask, c) in enumerate(zip((grmask, gwmask, brmask, bwmask,), ('rosybrown', 'deepskyblue', 'firebrick', 'royalblue'))):
        thresh = 0.001 if m > 1 else 0.2
        sns.kdeplot(data=df[mask], x='p', y='fp', fill=True, thresh=thresh, color=c, ax=ax[i], log_scale=(True, False))
        # Scatterplot(df[mask], x='p', y='fp', xerr=False, yerr=False, marker='x', s=40, c=c, ax=ax[i])
    Scatterplot(dfobs, x='p', y='fp', seg=('r', 'w', 'a'), c_r='r', c_w='b', c_a='grey', s=90, ax=ax[i], mlim=mlim, ec='k', linewidth=1.6)
    # fig.suptitle(str(ibest))
    # log_axis(ax[i], y=False)
    ax[i].set_xlim(ax[i].get_xlim()[0], (100 if key[-1] == 'G' else 90))
    ax[i].text(ax[i].get_xlim()[1]-20, 0.83, key[-1]+'-type host', size=21, ha='right', va='bottom')
ax[-2].set_xlabel('Period (days)')
ax[-1].set_xlabel('Period (days)')
textfont = 20
legend = ax[0].legend([plt.plot([], [], lw=0, marker='o', c='w')[0] for _ in range(4)], ['Gas-rich rocky'+' '*5, 'Gas-rich water'+' '*5, 'Bare rocky'+' '*5, 'Bare water'+' '*5],
                      title='KDE from model:', labelspacing=0.3, fontsize=22, loc='upper left', bbox_to_anchor=(-0.08, 1.01), frameon=False)
legend.get_title().set_fontsize(22)
for text, color in zip(legend.get_texts(), ['rosybrown', 'deepskyblue', 'firebrick', 'royalblue']):
    text.set_color(color)
legend = ax[1].legend([ax[-1].errorbar([0], [0], None, xerr=[0.1], marker='o', color='r'), ax[-1].errorbar([0], [0], None, xerr=[0.1], marker='o', color='b'), ax[-1].errorbar([0], [0], None, xerr=[0.1], marker='o', color='grey')], ['Rocky      ', 'Water-rich      ', 'Gas-rich      '],
                      title="From observation:", fontsize=22, frameon=False, loc='upper left', bbox_to_anchor=(-0.01, 1.01),
                      handlelength=1, labelspacing=0.3, )
legend.get_title().set_fontsize(22)
for text, color in zip(legend.get_texts(), ('r', 'b', 'grey')):
    text.set_color(color)

prettify(fig, labelsize=22, minorlocators={'x': None, 'y': None}, margin=dict(left=0.13, right=0.99, bottom=0.058, top=0.99), wspace=0.05, hspace=0.03)
for i in (0, 1):
    ax[i].set_ylabel(r'$\mathscr{F}$', size=32)
    ax[i].set_ylim(0.8, 3)
    ax[i].set_yscale('log', base=2)
    scalar_axis(ax[i], x=False, y=True)
    ax[i].yaxis.set_minor_locator(LogLocator(base=2, subs=np.arange(1.0, 10.0) * 0.1, numticks=10))

fig.savefig('../results/pfp_kde.jpg')

plt.show()
