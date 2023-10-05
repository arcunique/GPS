import numpy as np
from .constants import *
from .planet_radius_calc_methods import total_rad_generic

Myr2sec = 86400 * 365.25 * 1e6
Lxuv_decay = lambda t, ms: 1.2e30*ms*(t/100)**-1.5 * Myr2sec
Lxuv = lambda t, ms:  1.2e30*Myr2sec*ms if t<100 else Lxuv_decay(t, ms)
Lsun = 1.2e30*10**3.5


def calc_eta(func, m, rp):
    vesc = (2 * G * m * ME / rp / RE) ** 0.5
    return func(vesc)


def xdot_photevap(t, X, Mp, Rp, a, ms, eta=0.1):
    if np.all(X==0):
        return np.zeros(max([np.shape(v) for v in locals().values()]))
    return -eta*np.pi*Rp**3*Lxuv(t, ms)/(4*np.pi*a**2*G*Mp**2)


def atm_structure_photevap(mc, rc, rp, Teq, X, a, ms=1, eta=0.1, tmax=1000, tnum=200, t=None, rp_calc_meth=0, snapshot=False, **kwargs):
    if np.all(X == 0):
        return rp, X
    mc = np.atleast_1d(mc)
    rc = np.atleast_1d(rc)
    if rp is not None:
        rp = np.atleast_1d(rp)
    else:
        rp = total_rad_generic(rp_calc_meth, mc, X, Teq, rc, kwargs.pop('wf', None), **kwargs)
    X = np.atleast_1d(X)
    Teq = np.atleast_1d(Teq)
    a = np.atleast_1d(a)
    Xf = X.copy()
    rpf = rp.copy()
    ms = ms * np.ones(len(mc))
    mask = Xf > 0.0001
    # rpf[~mask] = rc[~mask]
    if t is None:
        t = np.linspace(0, tmax, tnum)
    dt = np.diff(t)
    # print(len(Xf[Xf < 0.001]) / len(Xf))
    if snapshot:
        rpss = []
        Xss = []
    for i in range(0, len(t) - 1):
        Xf[mask] = Xf[mask] + xdot_photevap(t[i], Xf[mask], mc[mask]*ME, rpf[mask]*RE, a[mask]*AU, ms[mask], eta) * dt[i]
        mask = Xf > 0.0001
        if sum(mask) == 0:
            break
        rpkwargs = kwargs.copy()
        if kwargs.get('age_update', True):
            rpkwargs['age'] = kwargs.get('age', 100) + t[i+1]
        if snapshot:
            rpss.append(rpf.copy())
            Xss.append(Xf.copy())
            tss = t + kwargs.get('age', 0)
    if not snapshot:
        rpf[~mask] = rc[~mask]
        Xf[~mask] = 0.0001
        # print(len(Xf[Xf<0.001])/len(Xf))
        return rpf, Xf
    for i in range(len(rpss)):
        mask = Xss[i] > 0.0001
        rpss[i][~mask] = rc[~mask]
        Xss[i][~mask] = 0.0001
    return tss, rpss, Xss

