import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib.colors import to_rgb
from .obshelper import bulk_comp_stat_obs, calc_mean_eb, f_mass_radius
import numpy as np
from matplotlib.ticker import MultipleLocator


def Plot(df, *args, **kwargs):
    xerr, yerr = None, None
    xu, yu = 0, 0
    if len(args) > 25:
        yerr = args[25]
        yu = 25
    elif 'yerr' in kwargs:
        yerr = kwargs['yerr']
    if len(args) > 26:
        xerr = args[26]
        xu = 26
    elif 'xerr' in kwargs:
        xerr = kwargs['xerr']
    if xu or yu:
        args = list(args)
    for var in ('x', 'y'):
        err = locals()[var + 'err']
        # print(var, err)
        if isinstance(err, (tuple, list)) and type(err[0]) == str and type(err[1]) == str:
            # print(var, df)
            err = df.__getitem__(list(err)[:2]).values.T
            # print(err)
            u = locals()[var + 'u']
            if u:
                args[u] = err
            else:
                kwargs[var + 'err'] = err
    xscalar = kwargs.pop('xscalar', True)
    yscalar = kwargs.pop('yscalar', True)
    ax = kwargs.pop('ax', None)
    if ax is None:
        ax = plt.gca()
    # print(kwargs)
    df.plot(*args, ax=ax, **kwargs)
    if xscalar and 'log' in ax.get_xscale():
        ax.xaxis.set_major_formatter(ScalarFormatter())
    if yscalar and 'log' in ax.get_yscale():
        ax.yaxis.set_major_formatter(ScalarFormatter())
    return ax


def Scatterplot(df, x, y, xerr=True, yerr=True, seg=(), cutoff=(0.4, 0.7), mlim=None, dpperc=(), **kwargs):
    if xerr:
        xerr = (x + 'e+', x + 'e-') if x + 'e+' in df and x + 'e-' in df else x + 'e'
    else:
        xerr = None
    if yerr:
        yerr = (y + 'e+', y + 'e-') if y + 'e+' in df and y + 'e-' in df else y + 'e'
    else:
        yerr = None
    if not seg:
        kwargs.setdefault('c', 'k')
        # print(kwargs)
        for key, val in kwargs.items():
            if type(val) == str and val[0] == '@':
                kwargs[key] = df[val[1:]].values * 3
        return Plot(df, x=x, y=y, xerr=xerr, yerr=yerr, kind='scatter', **kwargs)
    # dfw, dfo = seg_water_others(df, cutoff)
    if type(seg) == str:
        seg = (seg,)
    masks = []
    for s in range(len(seg)):
        if seg[s] == 'w':
            masks.append(bulk_comp_stat_obs(df, kind='w', cutoff=cutoff, mlim=mlim, dpperc=dpperc)[0])
        else:
            masks.append(bulk_comp_stat_obs(df, kind=seg[s], cutoff=cutoff, mlim=mlim, dpperc=dpperc)[0])
        _kwargs = kwargs.copy()
        for key in kwargs:
            if key.endswith('_' + seg[s]):
                _kwargs[key.rstrip('_' + seg[s])] = _kwargs.pop(key)
            elif any([key.endswith('_' + sg) for sg in seg]):
                _kwargs.pop(key)
            if type(kwargs[key]) == str and kwargs[key][0] == '@':
                _kwargs[key] = df[masks[s]][kwargs[key][1:]].values
        ax = Plot(df=df[masks[s]], x=x, y=y, xerr=xerr, yerr=yerr, kind='scatter', **_kwargs)
    return masks, ax


def Hist(df, x, bins=10, xscale='lin', normalize=True, weights=None, ax=None, **kwargs):
    if ax is None:
        ax = plt.gca()
    if weights is not None:
        normalize = False
    if not isinstance(df, (list, tuple)):
        if xscale == 'log' and type(bins) == int:
            bins = np.logspace(np.log10(df[x].min()), np.log10(df[x].max()), bins)
        if normalize:
            weights = np.ones(df.shape[0]) / df.shape[0]
        if np.isscalar(weights):
            weights *= np.ones(df.shape[0])
        ax.hist(df[x].values, bins=bins, weights=weights, **kwargs)
    else:
        if normalize:
            weights = [np.ones(dfi.shape[0]) / dfi.shape[0] for dfi in df]
        if np.isscalar(weights):
            weights = [np.ones(dfi.shape[0]) * weights for dfi in df]
        elif isinstance(weights, (list, tuple)):
            weights = list(weights)
            for w in range(len(weights)):
                if np.isscalar(weights[w]):
                    weights[w] = np.ones(df[w].shape[0]) * weights[w]
        X = [dfi[x].values for dfi in df]
        ax.hist(X, bins=bins, weights=weights, **kwargs)


def Histplot(df, x, bins=10, range=None, density=None, weights=None, normalize=False, xerr=None, shade={},
             histtype='step', ax=None, eb_kwargs={}, **kwargs):
    if ax is None:
        ax = plt.gca()
    if histtype == 'step':
        kwargs['where'] = 'mid'
    if type(x) == str:
        x = df[x]
    if normalize and weights is None:
        weights = np.ones(len(x)) / len(x)
    if xerr is None:
        count, binedge = np.histogram(x, bins=bins, range=range, density=density, weights=weights)
    else:
        if type(xerr) == str:
            xerr = df[xerr]
        elif isinstance(xerr, (list, tuple)) and len(xerr) == 2:
            xerr = list(xerr)
            if type(xerr[0]) == str:
                xerr[0] = df[xerr[0]]
            if type(xerr[1]) == str:
                xerr[1] = df[xerr[1]]
            xerr = (xerr[0] + xerr[1]) / 2
        nsim = kwargs.pop('nsim', 100)
        counts, binedges = [], []
        for i, X in enumerate(np.random.normal(x, xerr, (nsim, len(x)))):
            cn, be = np.histogram(X, bins=bins, range=range, density=density, weights=weights)
            counts.append(cn)
            binedges.append(be)
        count, cep, cen = calc_mean_eb(counts, axis=0, error=True)
        count = np.mean(counts, axis=0)
        counte = (cep + cen) / 2
        if all([np.array_equal(binedge, binedges[0]) for binedge in binedges[1:]]):
            binedge = binedges[0]
        else:
            binedge = np.median(binedges, axis=0)
    bincen = (binedge[1:] + binedge[:-1]) / 2
    binw = np.diff(binedge)
    if xerr is not None:
        fmt = eb_kwargs.pop('fmt', '.')
        ebc = eb_kwargs.pop('c', None)
    if histtype == 'bar':
        p, = ax.bar(x=bincen, height=count, width=binw, **kwargs)
    elif histtype == 'step':
        p, = ax.step(x=bincen, y=count, **kwargs)
    # if xerr is not None:
    #     color = ebc if ebc else p.get_color()
    #     ax.errorbar(bincen, count, counte, fmt=fmt, color=color, **eb_kwargs)
    for binrange, binkwargs in shade.items():
        mask = (bincen >= binrange[0]) & (bincen <= binrange[1])
        if histtype == 'bar':
            p, = ax.bar(x=bincen[mask], height=count[mask], width=binw[mask], **{**kwargs, **binkwargs})
        elif histtype == 'step':
            p, = ax.step(x=bincen[mask], y=count[mask], **{**kwargs, **binkwargs})
        if xerr is not None:
            color = ebc if ebc else p.get_color()
            ax.errorbar(bincen, count, counte, fmt=fmt, color=color, **eb_kwargs)
    return ax


def Histstep(x, bins=None, density=None, weights=None, norm=False, ax=None, **kwargs):
    if bins is None:
        bins = 10
    dn, xb = np.histogram(x, bins=bins, weights=weights, density=density)
    xc = (xb[:-1] + xb[1:])/2
    if norm:
        if type(norm) == bool:
            dn = dn / sum(dn)
        else:
            dn = dn / norm
    X = xc
    Y = dn
    if callable(kwargs.get('convertx', None)):
        X = kwargs['convertx'](X)
    if callable(kwargs.get('converty', None)):
        Y = kwargs['converty'](Y)
    if ax is None:
        plt.step(X, Y, where='mid', **kwargs)
    else:
        ax.step(X, Y, where='mid', **kwargs)


def legend_colortxt_nosym(ax=None,*args,**kwargs):
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    if ax is None: ax = plt.gca()
    ha = kwargs.pop('ha', None)
    ebsym = kwargs.pop('ebsym', False)
    fontsizes = kwargs.pop('fontsizes', None)
    fontcolors = kwargs.pop('fontcolors', None)
    alpha = kwargs.pop('alpha', None)
    kwargs.setdefault('handlelength', 0)
    kwargs.setdefault('handletextpad', 0)
    leg = ax.legend(*args,**kwargs,)
    renderer = ax.figure.canvas.get_renderer()
    shifts = [t.get_window_extent(renderer).width for t in leg.get_texts()]
    shift = max(shifts) if shifts else 0
    if fontsizes is None:
        fontsizes = (None,)*len(leg.legendHandles)
    if type(fontcolors)==str:
       fontcolors = (fontcolors,)*len(leg.legendHandles)
    if isinstance(alpha, (int, float, type(None))):
        alpha = (alpha, ) * len(leg.legendHandles)
    for h, handle in enumerate(leg.legendHandles): # leg.get_texts(), fontsizes):
        txt = leg.get_texts()[h]
        if handle is not None:
            if type(handle) != LineCollection or not ebsym:
                handle.set_visible(False)
        if fontcolors is None or fontcolors[h] == '':
            try:
                txt.set_color(handle.get_c())
            except AttributeError:
                try:
                    txt.set_color(handle.get_color()[0])
                except AttributeError:
                    txt.set_color(handle.get_edgecolor())
        elif fontcolors[h] is not None:
            if not alpha[h]:
                txt.set_color(fontcolors[h])
            else:
                txt.set_color(to_rgb(fontcolors[h]) + (alpha[h],))
        if ha is not None:
            txt.set_ha(ha)
            if ha == 'right':
                print(shift)
                txt.set_position((shift, 0))
                txt.set_ha('right')
        if fontsizes[h] is not None: txt.set_fontsize(fontsizes[h])
    return leg


def scalar_axis(ax=None, x=True, y=True):
    if ax is None:
        ax = plt.gca()
    if x:
        ax.xaxis.set_major_formatter(ScalarFormatter())
    if y:
        ax.yaxis.set_major_formatter(ScalarFormatter())


def log_axis(ax=None, x=True, y=False):
    if ax is None:
        ax = plt.gca()
    if x:
        ax.set_xscale('log')
    if y:
        ax.set_yscale('log')
    scalar_axis(ax, x, y)


def minorlocator(ax, xlocator=5, ylocator=None):
    if ax is None:
        ax = plt.gca()
    if xlocator:
        ml = MultipleLocator(xlocator)
        ax.xaxis.set_minor_locator(ml)
    if ylocator:
        ml = MultipleLocator(ylocator)
        ax.yaxis.set_minor_locator(ml)


def model_mr_curves(*args, y='rp', **kwargs):
    mmin = kwargs.pop('mmin', 0.5)
    mmax = kwargs.pop('mmax', 100)
    ax = kwargs.pop('ax', None)
    m = np.logspace(*np.log10([mmin, mmax]), 100)
    newargs = []
    for arg in args:
        if isinstance(arg, (int, float)):
            if y == 'rp':
                yval = f_mass_radius(arg) * m**(1/3.7)
            elif y == 'fp':
                yval = f_mass_radius(arg) * np.ones(100)
            elif y == 'dp':
                yval = m / (f_mass_radius(arg) * m**(1/3.7))**3
            newargs += [m, yval]
        else:
            newargs.append(arg)
    if ax is None:
        ax = plt.gca()
    return ax.plot(*newargs, **kwargs)


def prettify(fig=None, labelsize=20, ticklabelsize=18, ticksides='lrbt', tickdir='in', majortickparams={'length': 7, 'width': 2}, minortickparams={'length': 3, 'width': 1.5}, minorlocators={}, margin={}, wspace=None, hspace=None):
    if fig is None:
        fig = plt.gcf()
    for i, ax in enumerate(np.ravel(fig.axes)):
        if labelsize not in (0, None, [], ()):
            ax.set_xlabel(ax.get_xlabel(), size=labelsize if np.ndim(labelsize) == 0 else np.ravel(labelsize)[i])
            ax.set_ylabel(ax.get_ylabel(), size=labelsize if np.ndim(labelsize) == 0 else np.ravel(labelsize)[i])
        majortickkwargs = {}
        if ticklabelsize not in (0, None, [], ()):
            majortickkwargs["labelsize"] = ticklabelsize if np.ndim(ticklabelsize) == 0 else np.ravel(ticklabelsize)[i]
        ts = ticksides if type(ticksides) == str else np.ravel(ticksides)[i]
        tickdirkwargs = {}
        if 'l' in ts:
            tickdirkwargs["left"] = True
        if 'r' in ts:
            tickdirkwargs["right"] = True
        if 't' in ts:
            tickdirkwargs["top"] = True
        if 'b' in ts:
            tickdirkwargs["bottom"] = True
        if tickdir:
            tickdirkwargs["direction"] = tickdir
        majortickkwargs.update(tickdirkwargs)
        majortickkwargs.update(majortickparams)
        if majortickkwargs:
            ax.tick_params(which='major', **majortickkwargs)
        minortickkwargs = {**tickdirkwargs, **minortickparams}
        if minortickkwargs:
            ax.tick_params(which='minor', labelsize=0, **minortickkwargs)
        xlocator = minorlocators.get('x')
        ylocator = minorlocators.get('y')
        minorlocator(ax, xlocator=xlocator, ylocator=ylocator)
    subplots_adjust_kwargs = margin.copy()
    if wspace not in (None, [], ()):
        subplots_adjust_kwargs["wspace"] = wspace
    if hspace not in (None, [], ()):
        subplots_adjust_kwargs["hspace"] = hspace
    if subplots_adjust_kwargs:
        fig.subplots_adjust(**subplots_adjust_kwargs)
