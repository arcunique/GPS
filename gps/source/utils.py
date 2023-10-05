import warnings

import h5py
import os
import numpy as np
import json
from .constants import *

datapath = lambda *args: os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', *args)
with open(datapath('defstars.json')) as f:
    defstars = json.load(f)

def sma2per(sma, ms=1):
    if not isinstance(sma, (int, float, np.ndarray)):
        sma = np.array(sma)
    if not isinstance(ms, (int, float, np.ndarray)):
        ms = np.array(ms)
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        # try:
        #     2 * np.pi * ((sma * AU) ** 3 / G / (ms * MS)) ** 0.5 / 86400
        # except:
        #     print(sma, ms, (sma**3 * AU**3))
        #     raise
    return 2 * np.pi * ((sma**3 * AU**3) / G / (ms * MS)) ** 0.5 / 86400


def per2sma(p, ms=1):
    if not isinstance(p, (int, float, np.ndarray)):
        p = np.array(p)
    if not isinstance(ms, (int, float, np.ndarray)):
        ms = np.array(ms)
    return ((p * 86400 / 2 / np.pi) ** 2 * G * ms * MS) ** (1 / 3) / AU

def savedict2hdf(data: dict, file):
    if type(file) == str:
        file = h5py.File(file, 'w')
    for key, values in data.items():
        if not isinstance(values, dict):
            try:
                file.create_dataset(key, data=values)
            except TypeError:
                print(key, values, data)
                raise
        else:
            # testprint('dict', key, values.keys())
            group = file.create_group(key)
            savedict2hdf(values, group)
    if type(file) == h5py.File:
        file.close()

def loadhdf2dict(file):
    if type(file) == str:
        file = h5py.File(file)
    data = {}
    for key, values in file.items():
        if isinstance(values, h5py.Dataset):
            data[key] = values[()]
        elif isinstance(values, h5py.Group):
            data[key] = loadhdf2dict(values)
    if isinstance(file, h5py.File):
        file.close()
    return data

def extreme(q):
    return min(q), max(q)