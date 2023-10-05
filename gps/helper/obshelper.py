import numpy as np
from gps.source import f_mass_radius
import os.path

path2data = lambda *args: os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', *args)


def bulk_comp_stat_obs(df, err=False, kind='', cutoff=(0.4, 0.7), mlim=None, dpperc=(), nsim=100):
    if not err:
        if type(kind) == str:
            rev = False
            mask = df['fp'] > 0
            if kind.startswith('!'):
                kind, rev = kind[1:], True
            if kind == 'r':
                mask = df['fp'] < f_mass_radius(cutoff[0])
            if kind == 'w':
                mask = (df['fp'] >= f_mass_radius(cutoff[0])) & (df['fp'] <= f_mass_radius(cutoff[1]))
                # print('w', sum(mask))
            if kind == 'rw':
                mask = df['fp'] <= f_mass_radius(cutoff[1])
            if kind == 'a':
                mask = df['fp'] > f_mass_radius(cutoff[1])
            if mlim is not None:
                if kind in ('r', 'w', 'rw'):
                    mask &= df['m'] <= mlim
                if kind == 'a':
                    mask |= df['m'] > mlim
            if rev:
                mask = ~mask
            if dpperc:
                dplim = np.percentile(df['dp'][mask], dpperc)
                mask &= (df['dp'] >= dplim[0]) & (df['dp'] <= dplim[1])
            return mask, sum(mask)
        masks = []
        Nkinds = []
        for k in kind:
            res = bulk_comp_stat_obs(df, kind=k, cutoff=cutoff, dpperc=dpperc, mlim=mlim)
            masks.append(res[0])
            Nkinds.append(res[1])
        return masks, Nkinds
    masks = np.ones((nsim, df.shape[0]), dtype=bool) if type(kind) == str else [np.ones((nsim, df.shape[0]), dtype=bool) for _ in kind]
    Nkinds = np.ones(nsim) if type(kind) == str else [np.ones(nsim) for _ in kind]
    fpsim = np.random.uniform(df['fp'].values-df['fpe'].values, df['fp'].values+df['fpe'].values, (nsim, df.shape[0]))
    for i in range(nsim):
        dfi = df.copy()
        dfi['fp'] = np.random.uniform(df['fp'].values-df['fpe'].values, df['fp'].values+df['fpe'].values)
        mask, Nkind = bulk_comp_stat_obs(dfi, kind=kind, cutoff=cutoff, dpperc=dpperc, mlim=mlim)
        if type(kind) == str:
            masks[i] = mask
            Nkinds[i] = Nkind
        else:
            for k in range(len(kind)):
                masks[k][i] = mask[k]
                Nkinds[k][i] = Nkind[k]
    return masks, Nkinds, fpsim


def calc_mean_eb(data, axis=None, percentile=(16, 50, 84), error=True):
    if np.ndim(data) == 1 or axis is None:
        # print(data)
        data = np.ravel(data)
        axis = None
    res = np.percentile(data, percentile, axis=axis)
    if axis in (None, 0):
        lv, mv, uv = res
    elif axis == 1:
        lv, mv, uv = res.T
    if not error:
        return lv, mv, uv
    return mv, mv-lv, uv-mv


def fulton_redges(rpmin, rpmax):
    redges = np.logspace(np.log10(0.5), np.log10(20), 40)
    return redges[(redges < rpmax) & (redges > rpmin)]


def fulton_rphist():
    return np.loadtxt(path2data('fulton18.csv'), unpack=True)


def corr_with_fulton_rphist(rp):
    rf, of = fulton_rphist()
    rbins = fulton_redges(1, 5)
    rc = (rbins[:-1] + rbins[1:]) / 2
    okep = np.interp(rc, rf, of)
    return np.corrcoef(np.histogram(rp, bins=rbins, weights=np.ones(len(rp)) / len(rp))[0], okep)[0,1]