import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.transforms import IdentityTransform
from matplotlib.gridspec import GridSpec
from matplotlib.colors import to_hex
from init import *
import seaborn as sns
from matplotlib.lines import Line2D

lossmech = 'photev'
wflim = 0.4
mlim = 20

data = get_model(lossmechs=lossmech)

print(data.keys())
hax, hgap = 3, 0.3
hrats = [hax] * len(data)
for i in range(len(data)//2-1):
    hrats.insert(3*i+2, hgap)
print(hrats)
fig, ax = plt.subplots(len(data)+len(data)//2-1, 3, sharex='col', figsize=(15, 15), gridspec_kw={'height_ratios': hrats})
for i in range(len(data)//2-1):
    for j in range(3):
        ax[3*i+2, j].axis('off')
ax = ax[[i for i in range(len(ax)) if i not in 3*np.arange(len(data)//2)+2]]

P = np.logspace(0, 2, 50)  # Period
Rval = 10**(-0.11*np.log10(P)+0.4)  # Radius valley as a func of P

palette = {'bw': 'b', 'br': 'r', 'gw': 'deepskyblue', 'gr': 'rosybrown'}
nmbins = 6
sizelim = (20, 200)
mbins = np.linspace(0, 25, nmbins)
sizes = np.linspace(*sizelim, nmbins-1)
# print(sizes)
mranges = [f'{int(mbins[i])}-{int(mbins[i+1])}' for i in range(nmbins-1)]

for i, key in enumerate(data):
    print(key)
    df = data[key]['df'][data[key]['ibest']]
    if 'impact' in key:
        print(key, sum(df['imx'] >= 0.1)/df.shape[0])
    masks = {}
    for comp in palette:
        df.loc[bulk_comp_stat_model(df, comp, mlim=mlim, wflim=wflim, atmkey='xo')[0], 'ocomp'] = comp
        masks[comp] = bulk_comp_stat_model(df, comp, mlim=mlim, wflim=wflim, atmkey='x')[0]
        df.loc[masks[comp], 'comp'] = comp
    for j in range(nmbins-1):
        df.loc[(df['m'] >= mbins[j]) & (df['m'] < mbins[j+1]), 'mrange'] = j
    sns.scatterplot(data=df, x='p', y='rpo', palette=palette, hue='ocomp', legend=False, ax=ax[i, 0], size='mrange', sizes=sizelim, ec='k')
    sns.scatterplot(data=df, x='p', y='rp', palette=palette, hue='comp', legend=False, ax=ax[i, 1], size='mrange', sizes=sizelim,  ec='k')
    if key[-1] == 'G':
        ax[i, 2].plot(*fulton_rphist(), c='k')
        ax[i, 0].plot(P, Rval, '--k')
        ax[i, 1].plot(P, Rval, '--k')
        plt.text(0.05, 0.05, 'G-type host', size=20, transform=ax[i,0].transAxes)
    else:
        plt.text(0.05, 0.05, 'M-type host', size=20, transform=ax[i,0].transAxes)
    weights = [np.ones(sum(masks[mkey])) / df.shape[0] for mkey in masks]
    ax[i, 2].hist([df[masks[mkey]]['rp'] for mkey in masks], bins=fulton_redges(0.5, 5), weights=weights, stacked=True,
                  histtype='bar', color=[palette[mkey] for mkey in masks], alpha=0.7)
    ax[i, 0].set_xscale('log')
    ax[i, 1].set_xscale('log')
    ax[i, 0].set_ylim(0.5, 5)
    ax[i, 1].set_ylim(0.5, 5)
    scalar_axis(ax[i, 0])
    scalar_axis(ax[i, 1])
    ax[i, 0].set_ylabel(f'Radius ({runit})')
    ax[i, 1].set_ylabel("")
    ax[i, 1].yaxis.set_ticklabels([])
    if i != len(data)-1:
        ax[i, 1].xaxis.set_ticklabels([])
    ax[i, 2].yaxis.tick_right()
    ax[i, 2].yaxis.set_label_position("right")
    ax[i, 2].set_ylabel('Occurrence')
    ax[i, 2].set_xlabel("")
    # ax[i, 2].text(0.92, 0.8, chr(ord('A')+i), size=18, transform=ax[i, 2].transAxes)
ylim, xlim = max([ax[i,2].get_ylim()[1] for i in range(ax.shape[0])]), 5.4
# ylim, xlim = 0.15, 6
for i in range(ax.shape[0]):
    ax[i, 2].set_ylim(ax[i, 2].get_ylim()[0], ylim)
    # ax[i, 2].set_xlim(ax[i, 2].get_xlim()[0], xlim)
    ax[i, 0].set_ylim(ax[i, 0].get_ylim()[0], 4.3)
    ax[i, 2].set_xlim(0.6, 5.4)
legend = ax[-1,1].legend([plt.scatter([], [], color="white", marker='o', s=size, facecolor="whitesmoke", edgecolor='k') for size in sizes], mranges, frameon=True, title=f'Mass ({munit})', facecolor='whitesmoke', framealpha=1, fontsize=12)
legend.get_title().set_fontsize('12')
hpos, vstart, vgap = 0.97, 0.95, 0.15
size = 20
plt.text(hpos, vstart, 'Bare rocky', color=palette['br'], size=size, ha='right', va='top', transform=ax[0, 2].transAxes)
plt.text(hpos, vstart-vgap, 'Bare water', color=palette['bw'], size=size, ha='right', va='top', transform=ax[0, 2].transAxes)
plt.text(hpos, vstart-2*vgap, 'Observed', color='k', size=size, ha='right', va='top', transform=ax[0, 2].transAxes)
plt.text(hpos, vstart, 'Gas-rich rocky', color=palette['gr'], size=size, ha='right', va='top', transform=ax[2, 2].transAxes)
plt.text(hpos, vstart-vgap, 'Gas-rich water', color=palette['gw'], size=size, ha='right', va='top', transform=ax[2, 2].transAxes)
# ax[-1, 0].set_xticks([1, 10, 100])
ax[-1, 1].set_xticks([1, 10, 100])
ax[-1, 0].set_xlabel('Period (days)')
ax[-1, 1].set_xlabel('Period (days)')
ax[-1, 2].set_xlabel(f'Radius ({runit})')
plt.text(-0.27, 0, ' '*11 + 'Migration-A' + ' '*11, size=20, transform=ax[0, 0].transAxes, rotation=90, va='center', bbox=dict(facecolor='none', edgecolor='k',  boxstyle='round'))
plt.text(-0.27, 0, ' '*11 + 'Migration-B' + ' '*11, size=20, transform=ax[2, 0].transAxes, rotation=90, va='center', bbox=dict(facecolor='none', edgecolor='k',  boxstyle='round'))
plt.text(-0.27, 0, ' '*13 + 'In-situ' + ' '*15, size=20, transform=ax[4, 0].transAxes, rotation=90, va='center', bbox=dict(facecolor='none', edgecolor='k', boxstyle='round'))
plt.text(0.5, 1.15, ' '+f'Before atmospheric loss'+' ', size=20, transform=ax[0, 0].transAxes, rotation=0, ha='center', bbox=dict(facecolor='none', edgecolor='k', boxstyle='round'))
if lossmech == 'photev':
    lossmechtxt = ' '*12+f'After photo-evaporation loss and cooling'+' '*14
else:
    lossmechtxt = ' '*22+f'After impact loss and cooling'+' '*22
plt.text(1.03, 1.15, lossmechtxt, size=20, transform=ax[0, 1].transAxes, rotation=0, ha='center', bbox=dict(boxstyle='round', facecolor='none', edgecolor='k'))
# t.set_bbox()
prettify(fig, labelsize=20, ticklabelsize=16, margin=dict(left=0.09, right=0.94, bottom=0.07, top=0.94), wspace=0.05, hspace=0.07)
fig.savefig(f'../results/per_rad_stat_{lossmech}.jpg')
plt.show()