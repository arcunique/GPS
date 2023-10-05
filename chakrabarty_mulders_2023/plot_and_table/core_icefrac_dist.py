import matplotlib.pyplot as plt

from init import *

dfs = [
    GenPop(datapath('genesis_all_updated.hdf5'), models=lambda x: '-M-' not in x, star='G').get_state('if_2.50', 'p', time=-1).rename(columns={'if_2.50': 'wf'}),
    GenPop(datapath('genesis_all_updated.hdf5'), models=lambda x: '-M-' in x, star='G').get_state('if_2.20', 'p', time=-1).rename(columns={'if_2.20': 'wf'}),
    GenPop(datapath('genesis_all_updated.hdf5'), models=lambda x: '-M-' in x, star='G').get_state('if_0.80', 'p', time=-1).rename(columns={'if_0.80': 'wf'}),
]

models = ['In-situ', 'Migration (G-type host)', 'Migration (M-type host)']

# fig, ax = plt.subplots(1, 4, figsize=(10, 3.5), gridspec_kw={'width_ratios': [3, 3, 3, 2.3]})
# ax[-1].axis('off')
fig, ax = plt.subplots(1, 3, figsize=(10, 4.3))
colors = ['darkred', 'purple', 'darkblue']
# colors = ['#996633', 'purple', '#6666FF']
labels = ['Rocky\nWMF$<$0.1', 'Intermediate\n0.1$\leq$WMF$\leq$0.3', 'Water-rich\nWMF$>$0.3']
# labels = ['< 0.1', '0.1-0.3', '> 0.3']
for i, df in enumerate(dfs):
    df = df[(df['p']<100)].reset_index(drop=True)
    print(df['wf'].min(), df['wf'].max())
    df.loc[df['wf']<0.1, 'wmfcat'] = 0
    df.loc[(df['wf'] >= 0.1) & (df['wf'] < 0.3), 'wmfcat'] = 1
    df.loc[df['wf'] >= 0.3, 'wmfcat'] = 2
    # print(df['wmfcat'].value_counts(normalize=True))
    cdf = df['wmfcat'].value_counts(normalize=True).sort_index()
    print(cdf)
    keys, pcts = cdf.index, cdf.values * 100
    patches, texts, autotexts = ax[i].pie(pcts, colors=[colors[int(key)] for key in keys],
                                          explode=[0.06]*len(pcts),
                                          autopct='%.01f%%', shadow=True, startangle=[0, 155, -120][i], textprops={'fontsize':18})
    ax[i].set_title(models[i], size=20, x=[0.45, 0.43, 0.48][i])
    for at in autotexts:
        at.set_fontsize(15)
        # at.set_color('gold')
        at.set_color('w')
    ax[i].set_facecolor(bgcolor)
legend = ax[1].legend(labels, ncol=3, loc='lower center', bbox_to_anchor=(0.5, -0.22), fontsize=18)
for text, color in zip(legend.get_texts(), colors):
    text.set_color(color)


prettify(fig, margin=dict(left=0, right=1, bottom=0.07, top=0.99), wspace=0.02, hspace=0.01)

fig.savefig('../results/wmfpie.jpg')
plt.show()