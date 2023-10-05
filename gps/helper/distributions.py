import numpy as np
from scipy.integrate import quad
from scipy.interpolate import BPoly


def custom_random_sampler(x0, x1, custDist, size, nControl=10**6):
    samples=[]
    nLoop=0
    while len(samples)<size and nLoop<nControl:
        x=np.random.uniform(low=x0,high=x1)
        prop=custDist(x)
        assert prop>=0 and prop<=1
        if np.random.uniform(low=0,high=1) <= prop:
            samples += [x]
        nLoop+=1
    return np.array(samples)


def broken_power_per_dist(p0, k1, k2, prange=()):
    if not prange:
        return lambda p: 1 / ((p/p0)**-k1 + (p/p0)**-k2)
    return lambda p: broken_power_per_dist(p0, k1, k2)(p) / quad(broken_power_per_dist(p0, k1, k2), *prange)[0]

def broken_power_per_sampler(p0, k1, k2, N, prange=(1, 100)):
    # return custom_random_sampler(*prange, custDist=broken_power_per_dist(p0, k1, k2, prange), size=N)
    p = np.logspace(*np.log10(prange), 100000)
    pdf = broken_power_per_dist(p0, k1, k2)(p)
    return np.random.choice(p, N, p=pdf/sum(pdf))

def bergstein_mass_random_sampler(coeffs, logmrange, N):
    bp = BPoly(np.atleast_2d(coeffs).T, logmrange)
    # func = lambda x: bp(x)/bp.integrate(*logmrange)
    x = np.linspace(*logmrange, 10000)
    # y = 10**x
    # print(y[y<0])
    pdf = bp(x) / bp.integrate(*logmrange)
    y = 10 ** np.random.choice(x, N, p=pdf / sum(pdf))
    # print(y[y<0])
    return y
    # return custom_random_sampler(*logmrange, func, N)


def bergstein_X_random_sampler(coeffs, logxrange, N):
    bp = BPoly(np.atleast_2d(coeffs).T, [0,1])
    cdfs = np.linspace(0, 1, 1000)
    pdfdlogx = np.diff(cdfs)
    logx = (logxrange[1] - logxrange[0]) / (bp(1) - bp(0)) * (bp(cdfs) - bp(0)) + logxrange[0]
    dlogx = np.diff(logx)
    logxc = (logx[1:] + logx[:-1]) / 2
    pdf = pdfdlogx / dlogx
    return 10 ** np.random.choice(logxc, N, p=pdf / sum(pdf))