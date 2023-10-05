import numpy as np
from matplotlib.colors import to_rgb
from init import *

data = get_model(models=('insitu', 'migrationA', 'migrationB'), lossmechs='photev', startypes='G')
print({key: len(val['df']) for key, val in data.items()})

fig, ax = plt.subplots(3, figsize=(6, 10), sharey='all', sharex='all')

titles = ['In-situ', 'Migration-A', 'Migration-B']
colors = [to_rgb('b')+(0.7,), to_rgb('r')+(0.7,)]

for i, key in enumerate(data):
    ind = 50 if 'A' in key else 5
    df = data[key]['df'][ind]
    for mn, dfi in df.groupby('mn'):
        print(mn, dfi['m'].min(), dfi['m'].max(), dfi['p'].min(), dfi['p'].max())
    print()
    print(key, df['m'].min(), df['m'].max(), df['p'].min(), df['p'].max())
    print()
    mask = df['wf'] >= 0.3
    weights = [np.ones(sum(mask))/df.shape[0], np.ones(sum(~mask))/df.shape[0]]
    Hist([df[mask], df[~mask]], 'm', bins=np.logspace(np.log10(df['m'].min()), np.log10(df['m'].max()), 20), weights=1/df.shape[0], stacked=True, color=colors,  ax=ax[i])
    # Histstep(df['p'], bins=np.logspace(0, 2, 20), ax=ax[i])
    log_axis(ax[i])
    ax[i].set_title(titles[i], size=22, y=0.85)
    ax[i].set_ylabel('Occurrence')
    # ax[i].set_ylim(np.array(ax[i].get_ylim()) * np.array([1, 1.2]))
ax[-1].set_xlabel(f'Mass ({munit})')
plt.text(0.96, 0.75, 'Rocky core', transform=ax[-1].transAxes, size=20, c=colors[1], ha='right')
plt.text(0.96, 0.65, 'Water-rich core', transform=ax[-1].transAxes, size=20, c=colors[0], ha='right')
prettify(fig, labelsize=22, ticklabelsize=20, margin=dict(left=0.18, right=0.99, bottom=0.075, top=0.99), wspace=0.02, hspace=0.04, )
fig.savefig('../results/mcif.jpg')
plt.show()

