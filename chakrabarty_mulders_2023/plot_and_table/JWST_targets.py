import warnings
from itertools import product

import numpy as np
import pandas as pd

from init import *
from astropy.modeling.models import BlackBody
import astropy.units as U
try:
    from pathos.multiprocessing import Pool
    from itertools import product
    parallelimport = True
except:
    parallelimport = False
import tqdm
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

snowline = {'G': 2.2, 'M': 0.8}
iflim = 0.3, 0.7
mlim = 10

calc = False
tables = True

data = get_model(models=('migrationA', 'migrationB'), lossmechs={'A': 'photev', 'B': 'all'})
modelkeys = []
for key in data:
    if key not in modelkeys:
        modelkeys.append('_'.join(key.split('_')[:2]))
modeldata = {key: {'G': data[key+'_G'], 'M': data[key+'_M']} for key in modelkeys}

with open('data/obs.pk', 'rb') as file:
    obsdata = pkl.load(file)
    for key in obsdata:
        obsdata[key]['df']['Host type'] = key


def update_pe(df):
    planets = '(' + ', '.join(["'" + pl + "'" for pl in df['name'].values]) + ')'
    df1 = NasaExoplanetArchive.query_criteria(table="ps", select=('pl_name', 'pl_orbpererr1', 'pl_orbpererr2'), where=f"pl_name IN {planets}").to_pandas()
    df1 = df1.groupby('pl_name').mean().reset_index()
    df1 = df1.rename(columns={'pl_name': 'name'})
    df1['pe'] = df1.apply(lambda row: np.max([row['pl_orbpererr1'], row['pl_orbpererr2']]), axis=1)
    dfu = pd.DataFrame({'name': ['K2-314 d', 'GJ 3053 b', 'GJ 3053 c', 'L98-59 c', 'L98-59 d', 'CD-60 8051 b'],
                        'pe': [0.005, 0.00041, 0.00003, 2.5e-6, 6e-6, 0.00018]})
    df = pd.merge(df[[col for col in df.columns if col != 'pe']], df1[['name', 'pe']], on='name', how='left')
    # print(df[df['name'].isnull()])
    mask = dfu['name'].isin(df['name'].values)
    df.loc[df['name'].isin(dfu[mask]['name'].values), 'pe'] = dfu[mask]['pe'].values
    # print('upd', df[df['pe'].isnull()]['name'])
    return df


def func(arg):
    kdeobj, p, fp = arg
    return kdeobj.comp_pdf_given_detected('w', p, fp)


def calc_ww_prob_of_indiv_planet(kdeobjs, rowobs, nsim=(30, 30), obserr=True, npool=1):
    if type(kdeobjs) == period_F_KDE:
        kdeobjs = [kdeobjs]
    if obserr:
        prange = np.random.normal(rowobs['p'], rowobs['pe'], nsim[0])
        fprange = np.random.normal(rowobs['fp'], rowobs['fpe'], nsim[1])
    else:
        prange = rowobs[['p']]
        fprange = rowobs[['fp']]
    args = list(product(kdeobjs, prange[prange > 0], fprange[fprange > 0]))
    wwprob = []
    if npool < 2:
        kdeobj: period_F_KDE
        for kdeobj, p, fp in tqdm.tqdm(args, desc=rowobs['name']):
            wwprob.append(kdeobj.comp_pdf_given_detected('w', p, fp))
    else:
        if not parallelimport:
            raise ImportError('Either install Pathos for multiprocessing or use nppol=1.')
        pool = Pool(npool)
        with tqdm.tqdm(total=len(args), leave=False, position=0, desc=rowobs['name']) as pbar:
            for prob in pool.imap_unordered(func, args):
                wwprob.append(prob)
                pbar.update()
        pool.close()
    v = np.round(np.percentile(wwprob, [16, 50, 84]), 1) * 100
    return v[1], v[1]-v[0], v[2]-v[1], np.round(np.mean(wwprob)*100, 1), np.round(np.std(wwprob)*100, 1)
    # return np.round(np.mean(wwprob)*100, 1), np.round(np.std(wwprob)*100, 1)


def transpec_metric(row):
    scalefac = 0.19 if row['rp'] < 1.5 else 1.26 if 1.5 <= row['rp'] < 2.75 else 1.28 if 2.75 < row['rp'] < 4 else 1.15
    return scalefac * row['rp']**3 * row['T'] / row['m'] / row['rs']**2 * 10**(-row['jmag']/5)


def str2tuple(x):
    try:
        return eval(x)
    except Exception as e:
        # print(e, type(e))
        if "name 'nan' is not defined" in str(e):
            return (np.nan, np.nan)
        return

def rename_table_columns(df):
    return df.rename(columns={'name': 'planet', 'p': 'P', 'pe': 'Pe', 'm': 'mp', 'me+': 'mpe+', 'me-': 'mpe-', 'fp': 'Fp', 'fpe': 'Fpe',
                              'T': 'teq', 'Te': 'teqe', 'wwprob_migrationA_photev': 'migrationA_photev_wwprob',
                              'wwprob_migrationB_photev': 'migrationB_photev_wwprob',
                              'wwprob_migrationB_impact': 'migrationB_impact_wwprob'})


def calculate(obsdata, npool, comp='', oberr=True, models='', startypes='', planets='', saveto=''):
    dfobs = {}
    if not startypes:
        startypes = list(obsdata)
    for startype in startypes:
        dfobs[startype] = obsdata[startype]['df']
        brmask, bwmask, gmask = bulk_comp_stat_obs(dfobs[startype], kind=('r', 'w', 'a'), mlim=10, cutoff=iflim)[0]
        dfobs[startype].loc[brmask, 'Likely nature'] = 'bare rocky'
        dfobs[startype].loc[bwmask, 'Likely nature'] = 'bare water'
        dfobs[startype].loc[gmask, 'Likely nature'] = 'gas-rich'
        if comp == 'br':
            dfobs[startype] = dfobs[startype][brmask].reset_index(drop=True)
        elif comp == 'bw':
            dfobs[startype] = dfobs[startype][bwmask].reset_index(drop=True)
        elif comp == 'g':
            dfobs[startype] = dfobs[startype][gmask].reset_index(drop=True)
        elif comp == 'b':
            dfobs[startype] = dfobs[startype][brmask|bwmask].reset_index(drop=True)
        if planets:
            if type(planets) == str:
                planets = [planets]
            dfobs[startype] = dfobs[startype][dfobs[startype]['name'].isin(planets)].reset_index(drop=True)
        if not models:
            models = list(modeldata)
        elif type(models) == str:
            models = [models]
        for key in models:
            print(startype, key)
            if dfobs[startype].empty:
                continue
            kdeobjs = [period_F_KDE(df, iflim[0], mlim) for df in modeldata[key][startype]['df']]
            dfobs[startype]['wwprob_' + key] = dfobs[startype].apply(
                lambda row: calc_ww_prob_of_indiv_planet(kdeobjs, row, obserr=oberr, npool=npool), axis=1)
    dfobs = pd.concat([dfobs[startype] for startype in ('G', 'M')], ignore_index=True)
    if saveto:
        dfobs.to_csv(f'interim_results/{saveto}', index=False)
    return dfobs


def save_full_table(df, saveto=''):
    # print(df.columns)
    df: pd.DataFrame = rename_table_columns(df)
    # print(df[['migrationA_photev_wwprob', 'migrationB_photev_wwprob', 'migrationB_photev_wwprob']])
    for col in ('migrationA_photev_wwprob', 'migrationB_photev_wwprob', 'migrationB_impact_wwprob'):
        df[col] = df[col].apply(str2tuple).values
        # print(df[col].apply(type))
        df[col+'e'] = df[col].apply(lambda x: x[1] if hasattr(x, '__len__') else None).values
        df[col] = df[col].apply(lambda x: x[0] if hasattr(x, '__len__') else None).values
    df = df[['planet', 'mp', 'mpe+', 'mpe-', 'rp', 'rpe+', 'rpe-', 'P', 'Pe', 'teq', 'teqe', 'dp', 'dpe',
             'Fp', 'Fpe', 'ms', 'mse+', 'mse-', 'rs', 'rse+', 'rse-', 'teff', 'teffe', 'logg', 'met',
             'Host type', 'Likely nature', 'migrationA_photev_wwprob', 'migrationA_photev_wwprobe',
             'migrationB_photev_wwprob', 'migrationB_photev_wwprobe',
             'migrationB_impact_wwprob', 'migrationB_impact_wwprobe',
             'TSM', 'K_RV', 'ref']]
    header = '''# Title: Where are the Water Worlds? Identifying the Exo-water-worlds Using Models of Planet Formation and Atmospheric Evolution
# Authors: Aritra Chakrabarty & Gijs D. Mulders
# Created On: 29-08-2023
# Table: Properties of the sample list of observed planets with precise radius (8%) and mass (25%) measurements along with probabilities of being water worlds
# Columns:
#    planet: Planet name
#    mp: Planet mass in Earth-unit
#    mpe+: Upper unc. in planet mass
#    mpe-: Lower unc. in planet mass
#    rp: Planet radius in Earth-unit
#    rpe+: Upper unc. in planet radius
#    rpe-: Lower unc. in planet radius
#    P: Period (days)
#    Pe: Unc. im period
#    teq: Equilibrium temperature (K) of planet
#    teqe: Unc. in equilibrium temperature
#    dp: Planet density in Earth-unit
#    dpe: Unc. in planet density
#    Fp: F-value of planet
#    Fpe: Unc. in F-value
#    ms: Stellar mass in solar unit
#    mse+: Upper unc. in stellar mass
#    mse-: Lower unc. in stellar mass
#    rs: stellar radius in solar unit
#    rse+: Upper unc. in stellar radius
#    rse-: Lower unc. in stellar radius
#    teff: Effectove temperature (K) of host star
#    teffe: Unc. in effective temperature
#    logg: logg of host star
#    met: metallicity of hist star
#    Host type: Type of host star used in models (G or M)
#    Likely nature: Likely nature of planet (bare rocky, bare water, or gas-rich)
#    migrationA_photev_wwprob: Probability of being a water world (%) according to migration-A + photo-evaporation model
#    migrationA_photev_wwprobe: Unc. in probability
#    migrationB_photev_wwprob: Probability of being a water world (%) according to migration-B + photo-evaporation model
#    migrationB_photev_wwprobe: Unc. in probability
#    migrationB_imact_wwprob: Probability of being a water world (%) according to migration-B + impact model
#    migrationB_imact_wwprobe: Unc. in probability
#    TSM: Transmission spectroscopic metric for JWST/NIRISS (Kempton et al. 2018)
#    K_RV: RV semi-amplitude (m/s) 
#    ref: Reference
'''
    # print(df[['planet', 'ref']])
    with open(f'../results/{saveto}', 'w') as fw:  # allplanets_wwprob_tsm_K.csv
        fw.write(header)
        df.to_csv(fw, encoding='utf-8', mode='a', index=False)
    # df.to_csv('results/allplanets_wwprob_tsm_K.csv', index=False)



def latex_table_for_paper(df):
    dfobs['wwprob_migrationA_photev'] = dfobs['wwprob_migrationA_photev'].apply(str2tuple)
    dfobs['wwprob_migrationB_photev'] = dfobs['wwprob_migrationB_photev'].apply(str2tuple)
    # print(dfobs['TSM']>10)
    mask = df.apply(
        lambda row: (((row['Likely nature'] == 'bare water') and (row['wwprob_migrationA_photev'][0] >= 40)) or
                    ((row['Likely nature'] == 'gas-rich') and (row['wwprob_migrationA_photev'][0] >= 60) and (row['wwprob_migrationA_photev'][1] < 10))), axis=1)  #
    # mask = df['wwprob_migrationA_photev'].apply(lambda x: sum(x) > 75)
    df = df.rename(columns={'name': 'Planet', 'wwprob_migrationA_photev': 'ww prob (\%)', 'K_RV': 'K (m/s)', 'Likely nature': 'Nature'})
    df.loc[df['Nature'].apply(lambda x: 'bare' in x), 'Nature'] = 'bare'
    df['Period (days)'] = df.apply(lambda row: (f"{row['p']:0.02f}", (f"{row['pe']:0.04f}" if row['pe']>=1e-4 else f"{row['pe']}")), axis=1)
    df = df[mask][['Planet', 'Host type', 'Nature', 'Period (days)', 'ww prob (\%)', 'TSM', 'K (m/s)']]
    df = df.sort_values('Nature')
    df.fillna('--', inplace=True)
    # print(df['Period (days)'].apply(type))
    caption = r"List of planets with mean  of being a water world (ww) $> 40$\%, J-mag of the host stars $<$ 10, and TSM of the planets $>$ 10"
    label = r"tab:jwst"
    comments = r"Probabilities of being a water world (ww prob) are calcualted with the help of our migration-A + photo-evapoartion and migration-B + photo-evaporation models. " \
               r"The transmission spectroscopic metric (TSM) values are calculated by following \cite{kempton18}."
    print(create_aastex_table(df, caption=caption, label=label, comments=comments))


if __name__ == '__main__':
    if calc:
        dfobs = calculate(obsdata, 1, oberr=True, startypes='', saveto='allplanets_with_wwprob.csv')
        # dfobs = calculate(obsdata, 1, oberr=True, planets='Kepler-138 d', saveto='kepler138d_with_wwprob.csv')  # allplanets_with_wwprob.csv
    else:
        dfobs = pd.read_csv('../interim_results/allplanets_with_wwprob.csv')

    if tables:
        df_jmag_k = pd.read_csv('../interim_results/allplanets_jmag_K.csv')
        df_jmag_k['K_RV'] = df_jmag_k['K_RV'].apply(str2tuple)
        dfobs = pd.merge(dfobs, df_jmag_k, on='name', how='left')
        mask = dfobs['T'].isnull()
        dfobs.loc[mask, 'T'] = dfobs[mask]['teff'] * (dfobs[mask]['rs'] * RS / dfobs[mask]['a'] / AU) ** 0.5
        dfobs['TSM'] = dfobs.apply(lambda row: np.round(transpec_metric(row), 1), axis=1).values
        save_full_table(dfobs, 'allplanets_wwprob_tsm_K.csv')
        print(dfobs)
        latex_table_for_paper(dfobs)
        # print(pd.read_csv('../results/allplanets_wwprob_tsm_K.csv', comment='#'))
    else:
        dfobs = dfobs.rename(columns={'wwprob_migrationB_impact': 'wwI', 'wwprob_migrationB_photev': 'wwB',
                                      'wwprob_migrationA_photev': 'wwA'})
        print(dfobs[[col for col in dfobs if col in ('name', 'host', 'p', 'fp', 'fpe') or 'ww' in col]])





