from init import *
from matplotlib.ticker import LogLocator
mlim = 10
wflim = 0.3
ibest = {'migrationA': {'G': 104, 'M': 234}, 'insitu': None}
# ibest = None

data = get_model(models=('migrationA', 'insitu'), lossmechs='photev')
# data = get_model(models='migrationB', lossmechs='impact')
with open('../data/obs.pk', 'rb') as file:
    obsdata = pkl.load(file)

fig, ax = plt.subplots(2, 2, figsize=(10, 10), sharey='all', sharex='all')
ax = ax.ravel()
for i, key in enumerate(data):
    dfs = data[key]['df']
    mbwps = []
    if 'mig' in key:
        for df in dfs:
            bwmask = bulk_comp_stat_model(df, 'bw', mlim=mlim, wflim=wflim)[0]
            mbwps.append(df[bwmask]['m'].value_counts(bins=(2,20)).values[0])
        print('max mass bwp:', max(mbwps))
    if ibest[key.split('_')[0]]:
        df = dfs[ibest[key.split('_')[0]][key[-1]]]
    elif 'mig' in key:
        df = dfs[np.argmax(mbwps)]
    elif 'ins' in key:
        df = dfs[data[key]['ibest']]
    dfobs = obsdata[key[-1]]['df']
    brmask, grmask = bulk_comp_stat_model(df, ('br', 'gr'), mlim=mlim, wflim=wflim)[0]
    bwmask, gwmask = bulk_comp_stat_model(df, ('bw', 'gw'), mlim=mlim, wflim=wflim)[0]
    brmask &= df['fp'] < f_mass_radius(0.3)
    grmask &= (df['fp'] > f_mass_radius(0.7)) | (df['m'] > mlim)
    gwmask &= (df['fp'] > f_mass_radius(0.7)) | (df['m'] > mlim)
    # print(df[bwmask]['m'].max())
    # dfbw = df[bwmask]
    # print(dfbw[dfbw['m'] > 5][['m', 'a']])
    # print(key[-1])
    # print(df[(df['m']>6.8) & (df['m']<6.81) & (df['a']>0.069) & (df['a']<0.07)][['m', 'rc', 'T', 'a', 'i', 'rp', 'x']])
    model_mr_curves(0.5, 'b', alpha=0.3, lw=3, mmax=10, ax=ax[i])
    model_mr_curves(0, 'r', alpha=0.3, lw=3, mmax=13, ax=ax[i])
    for mask, c in zip((grmask, gwmask, brmask, bwmask,), ('rosybrown', 'deepskyblue', 'firebrick', 'cornflowerblue')):
        # sns.kdeplot(data=df[mask], x='m', y='rp', fill=True, thresh=0.5, color=c, ax=ax[i])
        alpha = 0.5 if c == 'firebrick' else 0.6
        Scatterplot(df[mask], x='m', y='rp', xerr=False, yerr=False, s=100, marker='o', c=c, ax=ax[i], alpha=alpha)
    Scatterplot(dfobs, x='m', y='rp', seg=('r', 'w', 'a'), c_r='r', c_w='b', c_a='grey', marker='D', s=70, ax=ax[i], mlim=10, ec='k', linewidth=2.5)
    log_axis(ax[i], y=False)
ax[0].set_yscale('log', base=2)
scalar_axis(ax[0], x=False, y=True)
ax[0].yaxis.set_minor_locator(LogLocator(base=2, subs=np.arange(1.0, 10.0) * 0.1, numticks=10))
ax[0].set_ylabel(f'Radius ({runit})')
ax[2].set_ylabel(f'Radius ({runit})')
# ax[4].set_ylabel(f'radius ({runit})')
ax[-2].set_xlabel(f'Mass ({munit})')
ax[-1].set_xlabel(f'Mass ({munit})')
# ax[0].set_yticks([1, 2, 3, 4])
ax[0].set_ylim(0.7, ax[0].get_ylim()[1])
ax[0].set_xlim(0.4, ax[0].get_xlim()[1])
legend_colortxt_nosym(ax[0], ['Gas-rich rocky', 'Gas-rich water', 'Bare rocky', 'Bare water'], fontcolors=['rosybrown', 'deepskyblue', 'firebrick', 'royalblue'], labelspacing=0.3, fontsize=17, loc='lower right', bbox_to_anchor=(1.03, -0.02), frameon=False)
legend = ax[1].legend([ax[-1].errorbar([0], [0], None, xerr=[0.1], marker='o', color='r'), ax[-1].errorbar([0], [0], None, xerr=[0.1], marker='o', color='b'), ax[-1].errorbar([0], [0], None, xerr=[0.1], marker='o', color='grey')], ['Rocky', 'Water-rich', 'Gas-rich'],
                      title="From observation", fontsize=17, frameon=False, loc='lower right', bbox_to_anchor=(1.03, -0.02),
                      handlelength=1, labelspacing=0.3, )
legend.get_title().set_fontsize(17)
for text, color in zip(legend.get_texts(), ('r', 'b', 'grey')):
    text.set_color(color)
# ax[1].get_legend().set_title("title")
textfont = 20
plt.text(-0.34, 0.5, ' '*8 + 'Migration-A' + ' '*8, size=textfont, transform=ax[0].transAxes, rotation=90, va='center', ha='center', bbox=dict(facecolor='none', edgecolor='k',  boxstyle='round'))
plt.text(-0.34, 0.5, ' '*10 + 'In-situ' + ' '*10, size=textfont, transform=ax[2].transAxes, rotation=90, va='center', ha='center', bbox=dict(facecolor='none', edgecolor='k',  boxstyle='round'))
# plt.text(-0.38, 0.5, ' '*8 + 'Migration-B +' + ' '*8 + '\nphoto-evaporation', size=textfont, transform=ax[2].transAxes, rotation=90, va='center', ha='center', bbox=dict(facecolor='none', edgecolor='k',  boxstyle='round'))
# plt.text(-0.38, 0.5, ' '*8 + 'Migration-B +' + ' '*8 + '\nimpact', size=textfont, transform=ax[4].transAxes, rotation=90, va='center', ha='center', bbox=dict(facecolor='none', edgecolor='k',  boxstyle='round'))
plt.text(0.5, 1.1, ' '*9 + 'G-type host' + ' '*9, size=textfont, transform=ax[0].transAxes, va='center', ha='center', bbox=dict(facecolor='none', edgecolor='k',  boxstyle='round'))
plt.text(0.5, 1.1, ' '*9 + 'M-type host' + ' '*9, size=textfont, transform=ax[1].transAxes, va='center', ha='center', bbox=dict(facecolor='none', edgecolor='k',  boxstyle='round'))
prettify(fig, minorlocators={'x': None, 'y': None}, margin=dict(left=0.17, right=0.98, bottom=0.07, top=0.925), wspace=0.05, hspace=0.1)
fig.savefig('../results/mr.jpg')
plt.show()

