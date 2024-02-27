from init import *
migmodel = 'A'
all_snowlines = [0.5, 0.6, 1, 1.1, 1.2, 1.5, 1.8, 2, 2.2, 2.3, 2.4, 2.5, 3, 3.3]

X = lambda: 10**np.random.uniform(-2.5, np.log10(0.1), 1)[0]
Xparametric = False
wflim = 0.3
mlim = 10

repeat = 5
corr_cutoff = None
corr_stop = 0.5
params = {'photev': np.mgrid[0.01:0.05:5j, 0.1:0.5:5j], 'impact': (5, 0.1)}
lossmechs = 'photev', 'impact'
startypes = 'GM'
calc = False
savetype = 'update'
saveto = f'../final_res_for_paper/data/migration{migmodel}1.pk'

snowline = 0.5
show = 'photev_M'
show_params = 0.03, 0.2
dfind = 'best'
plot = True, True, False


def get_model(startype, lossmech, params):
    print(lossmech, startype)
    params = [np.repeat(param, repeat) for param in params]
    print([param for param in zip(*params)])
    if lossmech == 'photev':
        models = ('Gen-M-s22', 'Gen-M-s50') if migmodel == 'A' else ('Gen-M-s10', 'Gen-M-s50')
        kwargs = [{'photevap_eta': lambda vesc: param[0] * (vesc / 25e4) ** param[1]} for param in zip(*params)]
        comment = "eta = params[0] * (Vesc/25e4) ** params[1]"
    elif lossmech == 'impact':
        models = ('Gen-M-s22', 'Gen-M-s50') if migmodel == 'A' else ('Gen-M-s10', 'Gen-M-s50')
        kwargs = [dict(tcol=param[0], imp_mx_cutoff=param[1]) for param in zip(*params)]
        # kwargs = dict(tcol=5, imp_mx_cutoff=0.1, photevap_eta=lambda vesc: params[0] * (vesc / 25e4) ** params[1])
        comment = "collision time: params[0], impact mass-ratio: params[1]"
    elif lossmech == 'photev+impact':
        models = 'Gen-M-s10', 'Gen-M-s50'
        kwargs = [{'photevap_eta': lambda vesc: param[0] * (vesc / 25e4) ** param[1], 'tcol': 5, 'imp_mx_cutoff': 0.1} for param in zip(*params)]
        comment = "eta = params[0] * (Vesc/25e4) ** params[1]"
    star = defstars[startype]
    res = {}
    for snowline in all_snowlines:
        while True:
            star['snowline'] = snowline
            dfs_params = []
            corrs_params = []
            dfs = []
            corrs = []
            final_params = []
            success = False
            for i in range(len(params[0])):
                kwargs = {'photevap_eta': lambda vesc: params[0][i] * (vesc / 25e4) ** params[1][i]}
                df = evolve_genesis_data(models, star, X, Xparametric, rpmax=5, start_age=100, end_age=5000, **kwargs).df
                corr = corr_with_fulton_rphist(df['rp'].values)
                if corr_cutoff:
                    corrs_params.append(corr)
                    dfs_params.append(df)
                    if (i+1) % repeat == 0:
                        # print('\tparams:', i+1, params[0][i], params[1][i], max(corrs_params))
                        if max(corrs_params) >= corr_cutoff:
                            corrs += corrs_params
                            dfs += dfs_params
                            final_params.append([param[i] for param in params])
                            success = True
                        elif max(corrs_params) < corr_stop:
                            print('\tpoor corr:', i+1, max(corrs_params))
                            break
                        corrs_params = []
                        dfs_params = []
                else:
                    dfs.append(df)
                    corrs.append(corr)
            if dfs or not corr_cutoff:
                break
        if not corr_cutoff:
            final_params = params
        res[snowline] = {'df': dfs, 'ibest': np.argmax(corrs), 'corr': corrs, 'params': final_params}
        bwp_bp_probs = [get_wwprob_using_kde(df, startype, iflim, mlim) for df in dfs]
        maxmbw = [df[bulk_comp_stat_model(df, 'bw', mlim=mlim, iflim=iflim)[0]]['m'].max() for df in dfs]
        if isinstance(final_params[0], (list, tuple)):
            print('max, min:', max(final_params), max(final_params, key=lambda x: x[1]), min(final_params), min(final_params, key=lambda x: x[1]))
        print('snowline:', snowline, max(corrs) if corrs else None, max(maxmbw) if maxmbw else None, np.mean(bwp_bp_probs), np.std(bwp_bp_probs), '%')
    res['params'] = params
    res['comment'] = comment
    return res

def save(res, savetype):
    if savetype in ('u', 'update') and os.path.isfile(saveto):
        with open(saveto, 'rb') as fread:
            resold = pkl.load(fread)
        for lossmech in lossmechs:
            for startype in startypes:
                if lossmech + '_' + startype in resold:
                    resold[lossmech + '_' + startype] = {**resold[lossmech + '_' + startype], **res[lossmech + '_' + startype]}
                else:
                    resold[lossmech + '_' + startype] = res[lossmech + '_' + startype]
        res = resold.copy()
    with open(saveto, 'wb') as fwrite:
        pkl.dump(res, fwrite)
    # print()

if calc or not os.path.isfile(saveto):
    res = {}
    for lossmech in lossmechs:
        for startype in startypes:
            res[lossmech + '_' + startype] = get_model(startype, lossmech, params[lossmech] if np.ndim(params[lossmech]) < 2 else [val.ravel() for val in params[lossmech]])
    if savetype and savetype != 'ask':
        save(res, savetype)
else:
    with open(saveto, 'rb') as fread:
        res = pkl.load(fread)

if show:
    # print(res.)
    if calc and len(all_snowlines) == 1:
        snowline = all_snowlines[0]
        show = f'{lossmechs[0]}_{startypes[0]}'
    # print(f'Showing for {show} and snowline={snowline}')
    dfs = res[show][snowline]['df']
    try:
        params = list(zip(*res[show]['params']))
    except:
        params = res[show]['params']
    bwp_bps = [get_wwprob_using_kde(df, show[-1], iflim=wflim, mlim=mlim) for df in dfs]
    if dfind == 'mid':
        dfind = len(dfs)//2
    elif dfind == 'best':
        dfind = res[show][snowline]['ibest']
    elif dfind == 'param':
        dfind = params.index(show_params)
    df = dfs[dfind]
    print('\n')
    print('number:', len(dfs), 'index:', dfind)
    print('bestcorr params:', params[dfind] if isinstance(params[0], (list, tuple)) else params)
    print('corr with fulton hist:', corr_with_fulton_rphist(df['rp'].values))
    # bwp_bps = res[show][snowline]['bwp_bp']
    # print(df.groupby('sy').size())

    with open('../final_res_for_paper/data/obs.pk', 'rb') as fread:
        resobs = pkl.load(fread)
    dfobs = resobs[show[-1]]['df']
    Nwo, Nrwo = bulk_comp_stat_obs(dfobs, err=True, kind=('w', 'rw'), mlim=mlim, cutoff=(wflim, 0.7), nsim=1000)[1]
    bwfobsm, *bwfobse = np.around(calc_mean_eb(Nwo / Nrwo * 100), 1)
    bwfobse = max(bwfobse)
    rfult = fulton_redges(1, 5)

    print('model ww%:', 'mean -', np.mean(bwp_bps), 'std -', np.std(bwp_bps))
    print('Obs ww%:', bwfobsm, bwfobse)

    if plot[1]:
        fig, ax = plt.subplots(1, 2)
        mbwps = []
        for df in dfs:
            bwmask = bulk_comp_stat_model(df, 'bw', mlim=mlim, iflim=wflim)[0]
            try:
                num = df[bwmask]['m'].value_counts(bins=(2, 20)).values[0]
            except:
                num = 0
            mbwps.append(num)
        print('max no of bwp mass>2:', max(mbwps))
        df1 = dfs[np.argmax(mbwps)]
        bwmask, brmask, gmask = bulk_comp_stat_model(df1, ('bw', 'br', 'g'), iflim=wflim, mlim=mlim)[0]
        bmask = bulk_comp_stat_model(df1, 'b', iflim=wflim, mlim=mlim)[0]
        grmask, gwmask = bulk_comp_stat_model(df1, ('gr', 'gw'), iflim=wflim, mlim=mlim)[0]

        print('max mass:', df1[bwmask]['m'].max())

        # sns.kdeplot(data=df[~bmask], x='p', y='fp', fill=True, thresh=0.1, levels=5, color='grey')
        # sns.kdeplot(data=df[bmask], x='p', y='fp', fill=True, thresh=0.1, levels=5, color='magenta')
        # sns.kdeplot(data=df[grmask], x='p', y='fp', fill=True, thresh=0.3, levels=5, color='rosybrown')
        # sns.kdeplot(data=df[gwmask], x='p', y='fp', fill=True, thresh=0.3, levels=5, color='deepskyblue')
        # sns.kdeplot(data=df[brmask], x='p', y='fp', fill=True, thresh=0.3, levels=5, color='r')
        # sns.kdeplot(data=df[bwmask], x='p', y='fp', fill=True, thresh=0.3, levels=5, color='b')

        # sns.kdeplot(data=df[bwmask], x='p', y='fp', fill=True, thresh=0.32, levels=5, color='b')
        # sns.kdeplot(data=df[brmask], x='p', y='fp', fill=True, thresh=0.32, levels=5, color='r')

        Scatterplot(df1[bwmask], x='m', y='rp', xerr=False, yerr=False, c='b', s=20, marker='x', alpha=0.6, ax=ax[0])
        Scatterplot(df1[brmask], x='m', y='rp', xerr=False, yerr=False, c='r', s=20, marker='x', alpha=0.6, ax=ax[0])
        Scatterplot(df1[gmask], x='m', y='rp', xerr=False, yerr=False, c='grey', s=20, marker='x', alpha=0.6, ax=ax[0])
        Scatterplot(dfobs, x='m', y='rp', seg=('r', 'w', 'a'), c_r='r', c_w='b', c_a='grey', s=40,  cutoff=(0.4, 0.7), mlim=mlim, ax=ax[0], ec='k')

        Scatterplot(df1[bwmask], x='p', y='fp', xerr=False, yerr=False, c='b', s=20, marker='x', alpha=0.6, ax=ax[1])
        Scatterplot(df1[brmask], x='p', y='fp', xerr=False, yerr=False, c='r', s=20, marker='x', alpha=0.6, ax=ax[1])
        Scatterplot(df1[gmask], x='p', y='fp', xerr=False, yerr=False, c='grey', s=20, marker='x', alpha=0.6, ax=ax[1])
        Scatterplot(dfobs, x='p', y='fp', seg=('r', 'w', 'a'), c_r='r', c_w='b', c_a='grey', s=40, ax=ax[1],  cutoff=(0.4, 0.7), ec='k', mlim=mlim)


        ax[0].set_xscale('log')
        ax[1].set_xscale('log')
        scalar_axis(ax[0])
        scalar_axis(ax[1])

    if plot[0]:
        df = df1.copy()
        brmask, bwmask, gmask = bulk_comp_stat_model(df, ('br', 'bw', 'g'), iflim=wflim, mlim=mlim)[0]
        # print('br:', sum(brmask), 'bw:', sum(bwmask), 'g:', sum(gmask))
        fig, ax = plt.subplots(figsize=(10, 7))
        # df.hist('rp', bins=rfult, weights=np.ones(df.shape[0]) / df.shape[0], histtype='step', ax=ax)
        weights = [np.ones(sum(bwmask)) / df.shape[0], np.ones(sum(brmask)) / df.shape[0], np.ones(sum(gmask)) / df.shape[0]]
        ax.hist([df[bwmask]['rp'].values, df[brmask]['rp'].values, df[gmask]['rp'].values], bins=rfult, weights=weights, stacked=True, color=['b', 'r', 'grey'])
        ax.plot(*fulton_rphist())
        ax.set_xlabel(f'Radius ({runit})', size=22)
        ax.set_ylabel('Occurrence', size=22)
        ax.tick_params(labelsize=18)
        # fig.savefig('../additional/figures/photev_impact_loweta.jpg')

    if plot[2]:
        fig, ax = plt.subplots(figsize=(10, 7))
        prs = np.array([])
        for sys, dfs in df.groupby('sy'):
            if dfs.size > 1:
                pers = dfs['p'].sort_values().values
                prs = np.append(prs, pers[1:]/pers[:-1])
        ax.hist(prs, bins=500, histtype='step')
        # ax.set_xticks([1, 7/6, 6/5, 5/4, 4/3, 7/5, 3/2, 5/3, 2, 3])
        for x in [1, 7/6, 6/5, 5/4, 4/3, 7/5, 3/2, 5/3, 2, 3]:
            ax.axvline(x, ls='--', c='k')
        ax.set_xlim(1,4)
        # ax.set_xscale('log')

    if plot[0] or plot[1] or plot[2]:
        plt.show()

if savetype == 'ask':
    savetype = input('Savetype? u for update, s for save: ')
    if savetype:
        save(res, savetype)







