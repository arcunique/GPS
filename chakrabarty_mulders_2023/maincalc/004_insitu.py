from init import *

with open(path2data('genesis_all_wfs.pk'), 'rb') as file:
    gendata = pkl.load(file)

all_snowlines = [2.5]
X = lambda: 10**np.random.uniform(-2.5, np.log10(0.1), 1)[0]
Xparametric = True
wflim = 0.4
mlim = 10

repeat = 5
params = {'photev': np.mgrid[0.01:0.05:5j, 0.1:0.5:5j], 'impact': (5, 0.1)}
lossmechs = 'photev', 'impact'
startypes = 'GM'
calc = True
savetype = ''
saveto = '../final_res_for_paper/data/insitu.pk'

snowline = 2.5
show = 'impact_G'
dfind = 'best'
plot = True, True


def calc_final_radii(params, lossmech, Gen):
    if lossmech == 'photev':
        Gen.calc_total_radii_after_losses(True, photevap_loss={'method': "owen+17", 'eta': lambda vesc: params[0]*(vesc/25e4)**params[1], 'rp_calc_meth': 2, 'tmax': 3000, 'tnum': 200, 'Tkh_Myr': 100, 'Tkh_update': True, 'age': 2, 'age_update': True, 'icefrac_cutoff': 0.1})
    elif lossmech == 'impact':
        Gen.calc_total_radii_after_losses(True, impact_loss={'method': "izidoro+22", 'imp_mx_cutoff': params[1], 'tcol': params[0]})
    return Gen.agg_to_df(datalims={'rp': (1, 5)})


def get_model(startype, lossmech, params):
    print(lossmech, startype)
    models = lambda x: '-M' not in x   # and x not in ('Gen-O-s22', 'Gen-P')
    if lossmech == 'photev':
        kwargs = {'photevap_eta': lambda vesc: params[0] * (vesc / 25e4) ** params[1]}
        comment = "eta = params[0] * (Vesc/25e4) ** params[1]"
    elif lossmech == 'impact':
        kwargs = dict(tcol=params[0], imp_mx_cutoff=params[1])
        comment = "collision time: params[0], impact mass-ratio: params[1]"
    elif lossmech == 'photev+impact':
        models = 'Gen-M-s22', 'Gen-M-s50'
        kwargs = {'photevap_eta': lambda vesc: params[0] * (vesc / 25e4) ** params[1], 'tcol': 0, 'imp_mx_cutoff': 0.1}
        comment = "eta = params[0] * (Vesc/25e4) ** params[1]"
    star = defstars[startype]
    # gen = evolve_gendata(models, star=star, X=X, Xsingle=Xsingle, break_reson_chain=True, data=gendata)
    res = {}
    for snowline in all_snowlines:
        star['snowline'] = snowline
        dfs = []
        bwp_bps = []
        gwp_gps = []
        corrs = []
        for i in range(repeat):
            df = evolve_genesis_data(models, star, X, Xparametric, rpmax=5, start_age=100, end_age=5000, **kwargs).df
            df = df[df['rpo'] > 2]
            Nbw, Nb = bulk_comp_stat_model(df, ('bw', 'b'), wflim=wflim, mlim=mlim)[1]
            # print(Nbw, Nb)
            Ngw, Ng = bulk_comp_stat_model(df, ('gw', 'g'), wflim=wflim, mlim=mlim)[1]
            dfs.append(df)
            bwp_bps.append(Nbw / Nb * 100)
            gwp_gps.append(Ngw/Ng*100)
            corrs.append(corr_with_fulton_rphist(df['rp'].values))
        res[snowline] = {'df': dfs, 'bwp_bp': bwp_bps, 'gwp_gp': gwp_gps, 'ibest': np.argmax(corrs)}
        print(snowline, np.mean(bwp_bps), np.std(bwp_bps), '%')
    res['params'] = params
    res['comment'] = comment
    return res

def save(res, savetype):
    if savetype in ('u', 'update') and os.path.isfile(saveto):
        with open(saveto, 'rb') as fread:
            resold = pkl.load(fread)
        for lossmech in lossmechs:
            for startype in startypes:
                res[lossmech + '_' + startype] = {**resold[lossmech + '_' + startype], **res[lossmech + '_' + startype]}
    with open(saveto, 'wb') as fwrite:
        pkl.dump(res, fwrite)
    # print()

if calc or not os.path.isfile(saveto):
    res = {}
    for lossmech in lossmechs:
        for startype in startypes:
            res[lossmech + '_' + startype] = get_model(startype, lossmech, params[lossmech])
    if savetype and savetype != 'ask':
        save(res, savetype)
else:
    with open(saveto, 'rb') as fread:
        res = pkl.load(fread)

if show:
    dfs = res[show][snowline]['df']
    if dfind == 'mid':
        dfind = len(dfs)//2
    elif dfind == 'best':
        dfind = res[show][snowline]['ibest']
    df = dfs[dfind]
    bwp_bps = res[show][snowline]['bwp_bp']

    with open('../final_res_for_paper/data/obs.pk', 'rb') as fread:
        resobs = pkl.load(fread)
    dfobs = resobs[show[-1]]['df']
    bwp_bp_obs = resobs[show[-1]]['bwp_bp']
    rfult = fulton_redges(1, 5)

    print('model ww%:', 'min -', min(bwp_bps), 'max -', max(bwp_bps), 'mean -', np.mean(bwp_bps), 'std -', np.std(bwp_bps))
    print('Obs ww%:', bwp_bp_obs)

    # print(df[(df['m']<2) & (df['rp'] > 1.5)][['m', 'rp', 'rpo', 'rc', 'xo', 'x', 'T']])
    # print(df[(df['m'] > 1.92) & (df['m'] < 1.93)][['m', 'rp', 'rpo', 'rc', 'xo', 'x', 'T', 'a']].values)

    if plot[1]:
        fig, ax = plt.subplots(1, 2)
        bwmask, brmask, gmask = bulk_comp_stat_model(df, ('bw', 'br', 'g'), wflim=wflim, mlim=mlim)[0]
        bmask = bulk_comp_stat_model(df, 'b', wflim=wflim, mlim=mlim)[0]

        print(df[bwmask]['m'].max(), df[bwmask]['m'].min(), df['m'].max())

        sns.kdeplot(data=df[bmask], x='p', y='fp', fill=True, thresh=0.32, levels=5, color='grey')
        sns.kdeplot(data=df[bwmask], x='p', y='fp', fill=True, thresh=0.32, levels=5, color='b')
        sns.kdeplot(data=df[brmask], x='p', y='fp', fill=True, thresh=0.32, levels=5, color='r')

        Scatterplot(df[bwmask], x='m', y='rp', xerr=False, yerr=False, c='b', s=20, marker='x', alpha=0.6, ax=ax[0])
        Scatterplot(df[brmask], x='m', y='rp', xerr=False, yerr=False, c='r', s=20, marker='x', alpha=0.6, ax=ax[0])
        Scatterplot(df[gmask], x='m', y='rp', xerr=False, yerr=False, c='grey', s=20, marker='x', alpha=0.6, ax=ax[0])
        Scatterplot(dfobs, x='m', y='rp', seg=('r', 'w', 'a'), c_r='r', c_w='b', c_a='grey', s=40, cutoff=(0.4, 0.7),
                    mlim=mlim, ax=ax[0], ec='k')

        Scatterplot(df[bwmask], x='p', y='fp', xerr=False, yerr=False, c='b', s=20, marker='x', alpha=0.6, ax=ax[1])
        Scatterplot(df[brmask], x='p', y='fp', xerr=False, yerr=False, c='r', s=20, marker='x', alpha=0.6, ax=ax[1])
        Scatterplot(df[gmask], x='p', y='fp', xerr=False, yerr=False, c='grey', s=20, marker='x', alpha=0.6, ax=ax[1])
        Scatterplot(dfobs, x='p', y='fp', seg=('r', 'w', 'a'), c_r='r', c_w='b', c_a='grey', s=40, ax=ax[1],
                    cutoff=(0.4, 0.7), ec='k', mlim=mlim)

        ax[0].set_xscale('log')
        ax[1].set_xscale('log')
        scalar_axis(ax[0])
        scalar_axis(ax[1])

    if plot[0]:
        brmask, bwmask, gmask = bulk_comp_stat_model(df, ('br', 'bw', 'g'), wflim=wflim, mlim=mlim)[0]
        print('br:', sum(brmask), 'bw:', sum(bwmask), 'g:', sum(gmask))
        fig, ax = plt.subplots(figsize=(10, 7))
        # df.hist('rp', bins=rfult, weights=np.ones(df.shape[0]) / df.shape[0], histtype='step', ax=ax)
        weights = [np.ones(sum(bwmask)) / df.shape[0], np.ones(sum(brmask)) / df.shape[0],
                   np.ones(sum(gmask)) / df.shape[0]]
        # ax.hist([df[bwmask]['rp'].values, df[brmask]['rp'].values, df[gmask]['rp'].values], bins=rfult, weights=weights,
        #         stacked=True, color=['b', 'r', 'grey'])
        # ax.hist(df['rp'].values, bins=rfult, weights=np.ones(df.shape[0])/df.shape[0], color='b', histtype='step')
        rps = []
        weights = []
        mods = []
        for mod, dfg in df.groupby('mn'):
            rps.append(dfg['rp'].values)
            weights.append(np.ones(dfg.shape[0])/df.shape[0])
            mods.append(mod)
        ax.hist(rps, bins=rfult, weights=weights, stacked=True)
        ax.legend(mods)
        ax.plot(*fulton_rphist(), c='k')
        ax.set_xlabel(f'Radius ({runit})', size=22)
        ax.set_ylabel('Occurrence', size=22)
        ax.tick_params(labelsize=18)

    if plot[0] or plot[1]:
        plt.show()

if savetype == 'ask':
    savetype = input('Savetype? u for update, s for save: ')
    if savetype:
        save(res, savetype)







