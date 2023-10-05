from init import *

a100 = per2sma(100, 1), per2sma(100, 0.35)
print(a100)

fig, ax = plt.subplots(1, 2, figsize=(12, 8), sharey=True)
for i, (model, modtype) in enumerate(zip(('Gen-M-s22', 'Gen-O-s22'), ('migration', 'in-situ'))):
    gp = GenPop(datapath('genesis_all_updated.hdf5'), models=model, star='G')  # 'Gen-M-s50'
    # print(gp.models.keys())
    # print(gp.models['Gen-M-s22']['run_01']['snapshots'].keys())
    df = gp.get_state('a', 'p', f'if_{2.5:0.02f}', 't', bb_only=False, time='all')
    df = df[(df['a'] > 0)]
    # fsids = df['wf'].values
    # gen.agg_to_df('id', 'a', f'if_{2.5:0.02f}', 't', source='ev')
    # gen.statdf = gen.statdf[gen.statdf['id'].isin(fsids)]
    # gen.show_stat_from_df('a-t', ax=ax[i], c=f'@if_{2.5:0.02f}', colorbar=False, alpha=0.6, cmap='coolwarm_r')
    # print(df.columns)
    sns.scatterplot(data=df, x='a', y='t', ax=ax[i], c=df['if_2.50'], alpha=0.6, cmap='coolwarm_r')
    ax[i].set_xlabel('a (AU)')
    # legend_colortxt_nosym(ax[i], [model+' ('+modtype+')'], fontcolors=['k'], fontsize=17, frameon=False, loc='upper center')
    ax[i].text(0.5, 1.01, model+' ('+modtype+')', color='k', size=17, ha='center', transform=ax[i].transAxes)
    ax[i].set_ylim(ax[i].get_ylim()[0], 10.7)
    ax[i].set_xscale('log')
    ax[i].plot([a100[0], a100[0]], [-0.3, 10.6], ls='--', c='k', lw=3)
    # ax[i].plot([a100[1], a100[1]], [0, 10.6], ls='--', c='k', lw=3)
    # ax[i].text(a100[0]+0.05, 10.3, 'P=100 days (G)', ha='left', size=15)
    ax[i].text(a100[0] + 0.05, 10.2, 'P=100 days', ha='left', size=15)
    # ax[i].text(a100[0]-0.15, 10.3, 'P=100 days (M)', ha='right', size=15)
    # ax[i].axvline(a100[1], ymax=11)
    # ax[i].set_xticks([1, 5, 10])
    scalar_axis(ax[i])
ax[1].set_xlim(ax[0].get_xlim())
cbar_ax = fig.add_axes([0.91, 0.1, 0.02, 0.8])
fig.colorbar(ax[1].collections[-1], cax=cbar_ax)
# cbar_ax.yaxis.tick_left()
# cbar_ax.yaxis.set_label_position("left")
cbar_ax.set_ylabel('Water/ice mass-fraction', size=18)
cbar_ax.tick_params(labelsize=18)
ax[0].set_ylabel('Time (Myr)')
prettify(fig, margin=dict(left=0.07, right=0.89, bottom=0.08, top=0.96), wspace=0)
fig.savefig('../results/icefrac_evol.jpg')
plt.show()