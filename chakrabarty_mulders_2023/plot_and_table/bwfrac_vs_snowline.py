import matplotlib.pyplot as plt
import numpy as np

from init import *

wflim = 0.3, 0.7
mlim = 10

data = get_model(models=('migrationA', 'migrationB'), lossmechs={'A': 'photev', 'B': 'all'}, snowline=None)
# print(data.keys())

with open('../data/obs.pk', 'rb') as file:
    obsdata = pkl.load(file)

fig, ax = plt.subplots(3, 2, figsize=(10, 14), sharey='all', sharex='all')
ax = ax.ravel()

for i, key in enumerate(data):
    print(key)
    snowlines = list(data[key])
    bwfm = np.array([])
    bwfe = np.array([])
    for snowline in snowlines:
        # bwfs.append(data[key][snowline][f'bwp_bp_prob_{wflim[0]}-{wflim[1]}'])
        bwfs = [get_wwprob_using_kde(df, key[-1], wflim=wflim[0], mlim=mlim) for df in data[key][snowline]['df']]
        bwfm = np.append(bwfm, np.mean(bwfs))
        bwfe = np.append(bwfe, np.std(bwfs))
        if snowline == snowline_fiducial['m'][key[-1]]:
            print(key, snowline, bwfm[-1], bwfe[-1])
    # bwfm = np.mean(bwfs, axis=1)
    # bwfe = np.std(bwfs, axis=1)
    fac = 3 if key == 'migrationA_photev_G' else 2
    ax[i].fill_between(snowlines, bwfm - fac * bwfe, bwfm + fac * bwfe, fc='lightgrey')
    ax[i].fill_between(snowlines, bwfm - bwfe, bwfm + bwfe, fc='gray')
    ax[i].plot(snowlines, bwfm, c='k', lw=2)
    Nwo, Nrwo = bulk_comp_stat_obs(obsdata[key[-1]]['df'], err=True, kind=('w', 'rw'), mlim=mlim, cutoff=wflim, nsim=1000)[1]
    bwfobsm, *bwfobse = np.around(calc_mean_eb(Nwo / Nrwo * 100), 1)
    bwfobse = max(bwfobse)
    print('obs', key[-1], bwfobsm, bwfobse)
    # print(key[-1], bwfobsm, bwfobse)
    slobs = 3 if i%2 else 0.7
    eb = ax[i].errorbar([slobs], [bwfobsm], [bwfobse], marker='o', capsize=10, color='b', elinewidth=3, markeredgewidth=3)
    if i == 1:
        ax[i].text(3.3, bwfobsm+bwfobse+4, 'From\nobservation', ha='right', color='b', size=20)
        # ax[i].legend(eb, ['Upper limit from observation'], fontsize=18, loc='upper right', size=18)
        # legend_colortxt_nosym(ax[i], ['Upper limit\nfrom observation\nwith uncertainty'], fontcolors='b', fontsize=16, loc='upper right')
    sl0 = np.interp(bwfobsm + bwfobse, (bwfm - fac * bwfe)[::-1], snowlines[::-1])
    sl1 = np.interp(bwfobsm + bwfobse, (bwfm + fac * bwfe)[::-1], snowlines[::-1])
    sl2 = np.interp(bwfobsm - bwfobse, (bwfm - fac * bwfe)[::-1], snowlines[::-1])
    sl3 = np.interp(bwfobsm - bwfobse, (bwfm + fac * bwfe)[::-1], snowlines[::-1])
    if sl0 >= min(snowlines) and (sl3-sl0) > 0.5:
        ax[i].axvline(sl0, c='deepskyblue')
    if sl3 <= max(snowlines) and (sl3-sl0) > 0.5:
        ax[i].axvline(sl3, c='deepskyblue')
    ax[i].plot([sl0-0.1, sl1+0.1], [bwfobsm + bwfobse, bwfobsm + bwfobse], ls='--', c='b')
    ax[i].plot([sl2-0.1, sl3+0.1], [bwfobsm - bwfobse, bwfobsm - bwfobse], ls='--', c='b')
    # bwfobsm = obsdata[key[-1]][f'bwp_bp_{wflim[0]}-{wflim[1]}'][0]
    # bwfobse = max(obsdata[key[-1]][f'bwp_bp_{wflim[0]}-{wflim[1]}'][1:])
ax[0].set_ylabel('BWP / BP (%)')
ax[2].set_ylabel('BWP / BP (%)')
ax[4].set_ylabel('BWP / BP (%)')
ax[-2].set_xlabel('Snowline location (AU)')
ax[-1].set_xlabel('Snowline location (AU)')
textfont = 20
plt.text(-0.38, 0.5, ' '*8 + 'Migration-A +' + ' '*8 + '\nphoto-evaporation', size=textfont, transform=ax[0].transAxes, rotation=90, va='center', ha='center', bbox=dict(facecolor='none', edgecolor='k',  boxstyle='round'))
plt.text(-0.38, 0.5, ' '*8 + 'Migration-B +' + ' '*8 + '\nphoto-evaporation', size=textfont, transform=ax[2].transAxes, rotation=90, va='center', ha='center', bbox=dict(facecolor='none', edgecolor='k',  boxstyle='round'))
plt.text(-0.38, 0.5, ' '*8 + 'Migration-B +' + ' '*8 + '\nimpact', size=textfont, transform=ax[4].transAxes, rotation=90, va='center', ha='center', bbox=dict(facecolor='none', edgecolor='k',  boxstyle='round'))
plt.text(0.5, 1.1, ' '*9 + 'G-type host' + ' '*9, size=textfont, transform=ax[0].transAxes, va='center', ha='center', bbox=dict(facecolor='none', edgecolor='k',  boxstyle='round'))
plt.text(0.5, 1.1, ' '*9 + 'M-type host' + ' '*9, size=textfont, transform=ax[1].transAxes, va='center', ha='center', bbox=dict(facecolor='none', edgecolor='k',  boxstyle='round'))
prettify(fig, minorlocators={'x': 0.1, 'y': 2}, margin=dict(left=0.205, right=0.98, bottom=0.05, top=0.94), wspace=0.05, hspace=0.1)
fig.savefig('../results/bwf_snowline.jpg')

plt.show()
