import matplotlib.pyplot as plt

from init import *

def teq2per(Teq, Teff, rs, ms):
    a = Teff**2/Teq**2*(rs*RS)
    return sma2per(a/AU, ms)

teff = {'M': (2500, 4000), 'K': (4000, 5200), 'G': (5200, 6500)} # These are the ranges of Teff considered for each type of host star
mlim = 10
calc_wflim = 0.3, 0.7
def_wflim = 0.4, 0.7
update = False  # False means calculate from scratch
file = '../final_res_for_paper/data/obs.pk'
updateto = '../final_res_for_paper/data/obs1.pk'

def get_obs_data(startype, wflim):
    res = {}
    if startype == 'M':
        lp = LPStats(which='imp').filter(rp=(1, 4), teff=teff[startype], p=(0, 100), m=(0, 30)).add_more()
        dflporg = lp.df
        dflp = lp.filter(errorpc={'rp': 8, 'm': 25}).df
        dflp['pe'] = get_value_with_error(teq2per, [dflp['T'].values, dflp['teff'].values, dflp['rs'].values, dflp['ms'].values],
                                           [dflp['Te'].values, dflp['teffe'].values, dflp['rse'].values, dflp['mse'].values], get='err')
        dft = TEPCATStats(which='spec').add_more().filter(rp=(1, 4), teff=teff[startype], p=(0, 100), m=(0, 30), errorpc={'rp': 8, 'm': 25}).df
        dft = dft[dft['name'] != 'L168-9 b']
        dfnew = dft[~dft['name'].isin(dflporg['name'].values) & ~dft['name'].isin(dflporg['namealt'].values)]
        dfobs = pd.concat([dflp, dfnew])
        Nwo, Nrwo = bulk_comp_stat_obs(dflp, err=True, kind=('w', 'rw'), mlim=mlim, cutoff=wflim, nsim=1000)[1]
        bwp_bp = np.around(calc_mean_eb(Nwo / Nrwo * 100), 1)
        res = {'df': dfobs, 'dflp': dflp, 'dfnew': dfnew, 'bwp_bp_lp': bwp_bp}
        # print(dfnew[['name', 'ref']])
    elif startype == 'G':
        dfobs = TEPCATStats(which='spec').add_more().filter(rp=(1, 4), teff=teff[startype], p=(0, 100), errorpc={'rp': 8, 'm': 25}).df
        res = {'df': dfobs}
    Nwo, Nrwo = bulk_comp_stat_obs(dfobs, err=True, kind=('w', 'rw'), mlim=mlim, cutoff=wflim, nsim=1000)[1]
    bwp_bp = np.around(calc_mean_eb(Nwo / Nrwo * 100), 1)
    res['bwp_bp'] = bwp_bp
    # print(startype, bwp_bp)
    return res

def update_obs_data(file, calc_wflim, def_wflim=(0.4, 0.7), saveto=''): #
    with open(file, 'rb') as fr:
        res = pkl.load(fr)
    for data in res.values():
        df = data['df']
        if def_wflim:
            for key in ('bwp_bp', 'bwp_bp_lp'):
                if key in data:
                    data[key + f'_{def_wflim[0]}-{def_wflim[1]}'] = data[key]
        if calc_wflim and np.isscalar(calc_wflim[0]):
            calc_wflim = [calc_wflim]
        for wflim in calc_wflim:
            Nwo, Nrwo = bulk_comp_stat_obs(df, err=True, kind=('w', 'rw'), mlim=mlim, cutoff=wflim, nsim=1000)[1]
            bwp_bp = np.around(calc_mean_eb(Nwo / Nrwo * 100), 1)
            data[f'bwp_bp_{wflim[0]}-{wflim[1]}'] = bwp_bp
            if 'dflp' in data:
                df = data['dflp']
                Nwo, Nrwo = bulk_comp_stat_obs(df, err=True, kind=('w', 'rw'), mlim=mlim, cutoff=wflim, nsim=1000)[1]
                bwp_bp = np.around(calc_mean_eb(Nwo / Nrwo * 100), 1)
                data[f'bwp_bp_lp_{wflim[0]}-{wflim[1]}'] = bwp_bp
        print(data['bwp_bp_0.3-0.7'])
    # print(res)
    if saveto:
        if saveto == '.':
            saveto = file
        with open(saveto, 'wb') as fw:
            pkl.dump(res, fw)



if not update:
    res = {}
    for startype in ('G', 'M'):
        res[startype] = get_obs_data(startype, calc_wflim)
        if startype == 'G':
            print(startype, res['G']['bwp_bp'])
        if startype == 'M':
            print(startype, res['M']['bwp_bp'], res['M']['bwp_bp_lp'])


    fig, ax = plt.subplots()
    dflp = res['M']['dflp']
    dfnew = res['M']['dfnew']
    Scatterplot(dflp, x='m', y='rp', seg=('r', 'w', 'a'), mlim=10, c_r='r', c_w='b', c_a='grey', ax=ax)
    Scatterplot(dfnew, x='m', y='rp', c='magenta', ax=ax)
    plt.show()

    save = input('Save? y/n [n]: ')
    if save == 'y':
        with open(file, 'wb') as fwrite:
            pkl.dump(res, fwrite)
else:
    update_obs_data(file, calc_wflim, def_wflim, updateto)


# print(res['G'].keys())
