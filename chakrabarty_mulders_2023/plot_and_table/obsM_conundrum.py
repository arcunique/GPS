import matplotlib.pyplot as plt

from init import *

with open('../data/obs.pk', 'rb') as file:
    res = pkl.load(file)['M']
newcolor = 'forestgreen'
fig, ax = plt.subplots(figsize=(6, 6))
dflp = res['dflp']
dfnew = res['dfnew']
Scatterplot(dflp, x='m', y='rp', seg=('r', 'w', 'a'), mlim=10, c_r='r', c_w='b', c_a='dimgray', ax=ax, ec='k', s=80,)
Scatterplot(dfnew, x='m', y='rp', c=newcolor, marker='x', s=60, ax=ax)
# Scatterplot(dfnew, x='m', y='rp', xerr=False, yerr=False, seg=('w',), mlim=10, c_w='b', marker='x', s=50, ax=ax)
model_mr_curves(0, 'r', 0.5, 'b', ax=ax, alpha=0.6, mmax=29)
model_mr_curves(0.3, 'b--', mmax=22, ax=ax, alpha=0.5)
model_mr_curves(0.7, 'b--', mmax=18, ax=ax, alpha=0.5)
x, yshift = 36, -0.14
ax.text(x, f_mass_radius(0.5)*x**(1/3.7)+yshift, r'50% H$_{\rm 2}$O line', color='b', size=15, ha='right')
x, yshift = 17, 0.01
ax.text(x, f_mass_radius(0.7)*x**(1/3.7)+yshift, r'70%', color='b', size=15, ha='right', va='bottom')
x, yshift = 22, 0.0
ax.text(x, f_mass_radius(0.3)*x**(1/3.7)+yshift, r'30%', color='b', size=15, ha='left', va='top')
x, yshift = 14, -0.03
ax.text(x, f_mass_radius(0)*x**(1/3.7)+yshift, r'Earth-like', color='r', size=15, ha='left', va='top')
ax.set_xlabel(f'Mass ({munit})')
ax.set_ylabel(f'Radius ({runit})')
log_axis(ax)
legend_colortxt_nosym(ax, ['Rocky from LP22', 'Water-rich from LP22', 'Gas-rich from LP22', 'From updated TEPCat'],
                      fontcolors=['r', 'b', 'dimgray', newcolor], fontsize=17, frameon=False, loc='upper left',
                      bbox_to_anchor=(0, 1.02))
prettify(fig, margin=dict(left=0.14, right=0.99, bottom=0.12, top=0.99), wspace=0.03, hspace=0.2)
fig.savefig(f'../results/obsM.jpg')
plt.show()