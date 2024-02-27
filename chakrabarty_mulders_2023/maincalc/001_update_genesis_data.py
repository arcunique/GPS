from init import *
try:
    import tqdm
    tqdmimp = True
except:
    print("Install tqdm to see progressbar")
    tqdmimp = False


''' 
Uncomment this section to update the masses of impactors and the water mass-fraction evolution for different snowline.
Once saved this section better be commented to run the next section    
'''

# all_snowlines = [0.5, 0.6, 1, 1.1, 1.2, 1.5, 1.8, 2, 2.2, 2.3, 2.4, 2.5, 3, 3.3]
# genfile = loadhdf2dict(path2data('genesis_all.hdf5'))
# gp = GenPop(genfile)
# print(gp.models['Gen-HM']['run_01']['collisions'].keys())
# gpsys = tqdm.tqdm(gp.systems.values()) if not tqdmimp else tqdm.tqdm(gp.systems.values())
# for gs in gpsys:
#     gs.update_impact_phases()
#     for snowline in [0.5]:
#         gs.calc_wmf(snowline=snowline)
# print(gp.models['Gen-HM']['run_01']['collisions'].keys())
# print(gp.models['Gen-HM']['run_01']['snapshots'].keys())
# gp.save_genesis(path2data('genesis_all_updated.hdf5'))

''' 
Uncomment this section to check broken resonance chain condition. Here we are ensuring 95% broken resonance chain and 5 % stable.    
'''

genfile = path2data('genesis_all_updated.hdf5')
for star in ('G', 'M'):
    gp = GenPop(genfile, star=star)
    gp.filter_runs_by_res_chain_break(0.95)
    df = gp.get_state()
    df = df[(df['a'] > 0) & (df['p'] < 100)]
    df.to_csv(path2data(f'genesis_all_{star}.csv'))



