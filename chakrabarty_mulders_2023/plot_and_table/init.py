import os.path
import sys
sys.path.append('../..')

from gps import *
from utils import *

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import pickle as pkl
from matplotlib.ticker import LogLocator, ScalarFormatter

matplotlib.use('Qt5Agg')

bgcolor = 'lightyellow' #'floralwhite'
plt.rcParams['axes.facecolor'] = bgcolor
# plt.rcParams['figure.facecolor'] = bgcolor
# matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{calrsfs}']

snowline_fiducial = {'m': {'G': 2.2, 'M': 0.8}, 'i': 2.5}
def get_model(models='all', lossmechs='all', startypes='all', snowline='fid', order='', path=r'../data'):
    mkeys = {'migrationA': 'a', 'migrationB': 'b', 'insitu': 'c'}
    lkeys = {'photev': 'a', 'impact': 'b'}
    skeys = {'G': 'a', 'M': 'b'}
    if models == 'all':
        models = list(mkeys)
    elif type(models) == str:
        models = [models]
    if lossmechs == 'all':
        lossmechs = {model[-1]: list(lkeys) for model in models}
    elif type(lossmechs) == str:
        lossmechs = {model[-1]: [lossmechs] for model in models}
    elif type(lossmechs) == dict:
        for key in lossmechs:
            if lossmechs[key] == 'all':
                lossmechs[key] = list(lkeys)
            elif type(lossmechs[key]) == str:
                lossmechs[key] = [lossmechs[key]]
    if startypes == 'all':
        startypes = list(skeys)
    if type(snowline) == str and snowline == 'fid':
        snowline = snowline_fiducial
    data = {}
    for model in models:
        filename = model+'.pk' if 'migration' in model else model+'.pk'
        print(os.path.abspath(path))
        with open(os.path.join(path, filename), 'rb') as fr:
            dat = pkl.load(fr)
        # print(model, lossmechs[model[-1]])
        for lossmech in lossmechs[model[-1]]:
            for startype in startypes:
                key = lossmech+'_'+startype
                # if model == 'insitu':
                #     print(key, dat.keys())
                if key in dat:
                    dkey = model+'_'+key
                    if not snowline:
                        data[dkey] = {sl: dat[key][sl] for sl in dat[key] if isinstance(sl, (int, float))}
                    if np.isscalar(snowline):
                        data[dkey] = dat[key][snowline]
                    if isinstance(snowline, (list, tuple)):
                        data[dkey] = {sl: dat[key][sl] for sl in snowline}
                    if type(snowline) == dict:
                        try:
                            data[dkey] = dat[key][snowline[startype]]
                        except KeyError:
                            try:
                                data[dkey] = dat[key][snowline[model[0]][startype]]
                            except TypeError:
                                data[dkey] = dat[key][snowline[model[0]]]
    if order:
        orders = lambda x: np.array([mkeys[x[0]], lkeys[x[1]], skeys[x[2]]])[[order.index('m'), order.index('l'), order.index('s')]]
        data = dict(sorted(data.items(), key=lambda x: ''.join(orders(x[0].split('_')))))
    return data


