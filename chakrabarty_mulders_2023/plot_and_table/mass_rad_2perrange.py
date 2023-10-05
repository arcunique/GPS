from matplotlib.colors import to_rgb
from init import *


df = get_model('migrationA', 'photev')['migrationA_photev_M']['df'][234]
with open('../data/obs.pk', 'rb') as file:
    dfobs = pkl.load(file)['M']['df']

colors = {'bw': 'royalblue', 'gw': 'deepskyblue', 'br': 'r', 'gr': 'rosybrown'}
alpha = {'bw': 1, 'gw': 0.8, 'br': 0.6, 'gr': 1}
size = {'bw': 250, 'gw': 150, 'br': 120, 'gr': 120}

fig, ax = plt.subplots(1, 2, sharey='row', figsize=(16, 8))
for comp in colors:
    mask = bulk_comp_stat_model(df, (comp,), mlim=10)[0][0]
    pmask = mask & (df['p'] < 10)
    # size = 150 if 'w' in comp else 100
    # alpha = 0.8 if 'w' in comp else 0.6
    # lw = 2 if 'w' in comp else None
    lw = 0.5
    zorder = 2 if 'w' in comp else 1
    Scatterplot(df[pmask], x='m', y='rp', xerr=False, yerr=False, s=size[comp], c=colors[comp], alpha=alpha[comp], marker='o', ec='w', linewidth=lw, ax=ax[0], zorder=zorder)
    pmask = mask & (df['p'] >= 10)
    Scatterplot(df[pmask], x='m', y='rp', xerr=False, yerr=False, s=size[comp], c=colors[comp], alpha=alpha[comp], marker='o', ec='w', linewidth=lw, ax=ax[1], zorder=zorder)
Scatterplot(dfobs[dfobs['p'] < 10], 'm', 'rp', seg='w', c_w='b', ec='k', marker='D', linewidth=2.5, s=150, ax=ax[0], mlim=10, zorder=3)
Scatterplot(dfobs[dfobs['p'] >= 10], 'm', 'rp', seg='w', c_w='b', ec='k', marker='D', linewidth=2.5, s=150, ax=ax[1], mlim=10, zorder=3)
log_axis(ax[0])
log_axis(ax[1])
handles = [plt.plot([], [], marker='o', ls='', c=colors[comp], ms=10, alpha=alpha[comp])[0] for comp in colors] + [plt.errorbar([], [], [], marker='D', ms=10, ls='', markeredgecolor='k', color='b')]
colors['o'] = 'b'
alpha['o'] = 1
legend = ax[0].legend(handles, ['Bare water', 'Gas-rich water', 'Bare rocky', 'Gas-rich rocky', 'Bare water from\nobservation'],
             frameon=False, fontsize=25, loc='upper left', bbox_to_anchor=(-0.05, 0.98), handletextpad=0)
for text, comp in zip(legend.get_texts(), colors):
    text.set_color(to_rgb(colors[comp]) + (alpha[comp],))
plt.text(0.96, 0.05, 'P < 10 days', size=26, transform=ax[0].transAxes, ha='right', va='bottom')
plt.text(0.96, 0.05, 'P > 10 days', size=26, transform=ax[1].transAxes, ha='right', va='bottom')
ax[0].yaxis.set_label_coords(-0.09, 0.5)
ax[0].set_ylabel(f'Radius ({runit})')
ax[0].set_xlabel(f'Mass ({munit})')
ax[1].set_xlabel(f'Mass ({munit})')
prettify(fig, labelsize=24, ticklabelsize=22, margin=dict(left=0.07, right=0.99, bottom=0.1, top=0.99), wspace=0.02)
fig.savefig('../results/mrww.jpg')
plt.show()
