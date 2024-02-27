import numpy as np
import random
import pandas as pd
from pathos.multiprocessing import Pool
from .utils import *
from .constants import *

with open(datapath('defstars.json')) as f:
    defstars = json.load(f)


class GenSys:
    '''
    A class to define planetary systems
    '''

    def __init__(self, system: (dict, str), star: (dict, str)):
        '''
        :param system: dictionary or file name from which the dict will be read. It needs to be in the format
                      {"snapshots": {"mass: [], "sma": [], "t": [], ...}, "collisions":{"collision_time": [],
                      "particle_i": [], "particle_j": []}, 'planets': {"plaent 1": {"bb": [], "id": 1}, ...}}
                      type: dict or str

        :param star: dictionary in the same format as each star type defined in "defstars" or a type (F, G, or M) given in defstars
                     type: dict or str
        '''
        self.sys = None
        self.imp_sum = []
        if system:
            self. load(system)
        if star:
            self.add_star(star)
            
    def load(self, system):
        if isinstance(system, dict):
            self.sys = system
        elif isinstance(system, str):
            self.sys = loadhdf2dict(system)

    def save(self, saveto):
        savedict2hdf(self.sys, saveto)
        
    def add_star(self, star):
        if type(star) == str:
            self.star = defstars[star]
        elif isinstance(star, dict):
            self.star = star
        self._if = f'if_{self.star["snowline"]:0.2f}'

    def remove_vars(self, *args):
        for arg in args:
            self.sys['snapshots'].pop(arg, None)

    def rename_vars(self, rename: dict):
        for old, new in rename.items():
            try:
                self.sys['snapshots'][new] = self.sys['snapshots'].pop(old)
            except Exception as e:
                print(e)

    def copy_vars_to(self, obj, *args):
        if not args:
            args = self.sys['snapshots'].keys()
        for arg in args:
            try:
                obj.sys['snapshots'][arg] = self.sys['snapshots'][arg].copy()
            except Exception as e:
                print(e)

    def update_impact_phases(self):
        ''' Updates masses of impacting bodies right before collision in the 'collisions' section '''
        # print('here')
        ss, co = self.sys['snapshots'], self.sys['collisions']
        t0 = ss['t'] == min(ss['t'])
        mbb = {bb: mbb for bb, mbb in zip(ss['id'][t0], ss['mass'][t0])}
        co['mass_i'] = np.zeros(len(co['collision_time']))
        co['mass_j'] = co['mass_i'].copy()
        ss['imp_mr_max'] = np.zeros(len(ss['t']))
        for ct in np.sort(np.unique(co['collision_time'])):
            mask = co['collision_time'] == ct
            co['mass_i'][mask] = np.array([mbb[i] for i in co['particle_i'][mask]])
            co['mass_j'][mask] = np.array([mbb[j] for j in co['particle_j'][mask]])
            mbb.update(zip(co['particle_i'][mask], co['mass_i'][mask] + co['mass_j'][mask]))
            # mask = (ss['t'] > ct) & np.isin(ss['id'], co['particle_i'])
            # uids = np.unique(ss['id'][mask])
            # for uid in uids:
            #     mask1 = (co['particle_i'] == uid) & (co['mass_j'] != 0)
            #     ss['imp_mr_max'][mask & (ss['id'] == uid)] = np.max(co['mass_j'][mask1]/co['mass_i'][mask1])

    def calc_wmf(self, check=True, wmf0_in=0, wmf0_out=0.5, snowline=None, wmf0_func=None, key=None):
        f'''
        Calculates evolution (snapshots) of water mass-fraction that changes due to collisions and saves in the 
        "snapshots" section with a new key which is either named as the argument "key" if provided or as f"{snowline:0.2f}"
        :param check: if True checks whether the snowline key already exists in the "snapshots" section and skips if already exists, if False then calculates always 
                      type: bool
        :param wmf0_in: Inner wmf value in case a step function is assumed, "wmf0_func" should be None
                        type: NoneType or int or float
        :param wmf0_out: Outer wmf value in case a step function is assumed, "wmf0_func" should be None
                         type: NoneType or int or float
        :param snowline: Snowline distance from host star in AU in case a step function is assumed, "wmf0_func" should be None
                         type: NoneType or int or float
        :param wmf0_func: wmf as a function of the distance from host star in AU
                          type: NoneType or function
        :param key: if provided the value is stored under this key in the "snapshots" section else the key used is f"{snowline:0.2f}"
        '''
        # print('here')
        if snowline is None:
            snowline = self.star['snowline']
        # print(check)
        ss, pb, co = self.sys['snapshots'], self.sys['planets'], self.sys['collisions']
        if check and ss.get(f'if_{snowline:0.2f}') is not None:
            return
        tss = ss['t']
        t0 = tss == min(tss)
        cts = np.insert(np.insert(np.unique(co['collision_time']), 0, min(tss)), len(co['collision_time']) + 1, max(tss) + 1)
        bbs, sma0 = ss['id'][t0], ss['sma'][t0]
        if wmf0_func and callable(wmf0_func):
            ifc = dict(zip(bbs, [wmf0_func(a) for a in sma0]))
        else:
            ifc = dict(zip(bbs, np.where(sma0 > snowline, wmf0_out, wmf0_in)))
        mc = {bbs[i]: ss['mass'][t0][i] for i in range(len(bbs))}
        ifss = np.zeros(len(tss))
        for i in range(len(cts) - 1):
            ti = (tss >= cts[i]) & (tss < cts[i + 1])
            for bb in bbs:
                ifss[ti & (ss['id'] == bb)] = ifc[bb]
            if cts[i + 1] in co['collision_time']:
                for bbi, bbj in zip(co['particle_i'][co['collision_time'] == cts[i + 1]],
                                    co['particle_j'][co['collision_time'] == cts[i + 1]]):
                    ifc[bbi] = (mc[bbi] * ifc[bbi] + mc[bbj] * ifc[bbj]) / (mc[bbi] + mc[bbj])
                    mc[bbi] += mc[bbj]
        if key is None:
            key = f'if_{snowline:0.2f}'
        ss[key] = ifss

    @property
    def available_symbols(self):
        return 'm', 'T', 'a', 'p', 'wf'

    def decode_symbols(self, *syms):
        arrs = []
        if type(syms) == str:
            syms = [syms]
        for sym in syms:
            # print(sym)
            if sym in self.sys['snapshots']:
                arrs.append(self.sys['snapshots'][sym])
            if sym == 'm':
                # print('here')
                arrs.append(self.sys['snapshots']['mass'])
            if sym == 'T':
                arrs.append(self.star['Teq@1AU'] / self.sys['snapshots']['sma'] ** 0.5)
            if sym == 'a':
                arrs.append(self.sys['snapshots']['sma'])
            if sym == 'p':
                arrs.append(sma2per(self.sys['snapshots']['sma'], self.star['mass']))
            if sym == 'wf':
                arrs.append(self.sys['snapshots'].get(self._if, np.full(len(self.sys['snapshots']['t']), np.nan)))
            # if sym == 'ims':
            #     arrs.append(self.imp_sum)
        return arrs

    def get_state(self, *args, bb_only=True, time=-1):
        '''
        Returns a dataFrame with states of "snapshots" variables as each column, e,g. states of "t", "mass", "sma", etc. at a point of time
        :param args: The keys from the "snapshots" section or their shorthands: "m" for "mass", "T" for Teq, "a" for "sma", "p" for period, "wf" for default wmf
        :param bb_only: if True, considers only those building block ids that make up the final planets and ignores other ids
        :param time: at what time the states are calculated (default=-1), 0 and -1 imply the initial and final states respectively
        :return: Pandas DataFrame
        '''
        if not args:
            args = self.available_symbols
        if time != 'all':
            if time == -1 or time >= max(self.sys['snapshots']['t']):
                tc = self.sys['snapshots']['t'] == max(self.sys['snapshots']['t'])
            else:
                t = np.unique(self.sys['snapshots']['t'])
                tc = self.sys['snapshots']['t'] == t[t >= time][0]
        else:
            tc = Ellipsis
        if bb_only:
            tc &= np.isin(self.sys['snapshots']['id'], np.array([val['id'] for val in self.sys['planets'].values()]))
        if 'ims' in args:
            self.get_impact_phase_summary(bb_only=False)
        for arr in self.decode_symbols(*args):
            try:
                arr[tc]
            except:
                print(len(arr), len(tc))
                raise
        return pd.DataFrame(dict(zip(args, [arr[tc] for arr in self.decode_symbols(*args)])))

    def get_impact_phase_summary(self, bb_only=True):
        ss, pb, co = self.sys['snapshots'], self.sys['planets'], self.sys['collisions']
        tf = ss['t'] == max(ss['t'])
        if bb_only:
            tf &= np.isin(self.sys['snapshots']['id'], np.array([val['id'] for val in self.sys['planets'].values()]))
        return np.array([{key: val[co['particle_i'] == bb] for key, val in co.items()} for bb in ss['id'][tf]])

    def get_resonance_status(self, time=-1, bb_only=False):
        ''' Calculates if the planets are in resonance chain or broken'''
        ss, pb = self.sys['snapshots'], self.sys['planets']
        prres = np.array([1, 7/6, 6/5, 5/4, 4/3, 7/5, 3/2, 5/3, 2, 3])
        if time == -1 or time >= max(ss['t']):
            tc = ss['t'] == max(ss['t'])
        else:
            t = np.unique(ss['t'])
            tc = ss['t'] == t[t >= time][0]
        if bb_only:
            tc &= np.array([ss['id'] in val['id'] for val in pb.values()])
        per = sma2per(ss['sma'][tc], self.star['mass'])
        per = np.sort(per)
        if len(per) > 1:
            perrat = per[1:] / per[:-1]
            if np.all([np.min(np.abs(pr - prres)) < 0.05 * prres[np.argmin(np.abs(pr - prres))] for pr in perrat]):
                return True
            return False


class GenPop:
    '''
    A class to define a population of planetary systems
    '''
    def __init__(self, genfile, models='all', star=''):
        '''

        :param genfile: dictionary or file name from which the dict will be read.
                        type: dict or str
        :param models: model name or a list of model names or "all" models (default 'all') to be chosen from the dict.
                       It can also be a function that returns True/False for each model in the dict. In that case, the
                       models for which the function returns True will be selected.
                       type: str or list or function
        :param star: dictionary in the same format as each star type defined in "defstars" or a type (F, G, or M) given in defstars
                     type: dict or str
        '''
        self.models = {}
        if genfile:
            self.load_genesis(genfile, models, star)

    def load_genesis(self, data, models='all', star=''):
        if type(data) == str:
            self.models = loadhdf2dict(data)
        else:
            self.models = data
        if not (type(models) == str and models == 'all'):
            if type(models) == str:
                models = [models]
            elif callable(models):
                models = [mod for mod in self.models if models(mod)]
            self.models = {key: val for key, val in self.models.items() if key in models}
        self.get_systems(star)

    def save_genesis(self, file):
        savedict2hdf(self.models, file)

    def get_systems(self, star):
        self.systems = {(mod, run): GenSys(self.models[mod][run], star) for mod in self.models for run in self.models[mod]}

    def filter_runs_by_res_chain_break(self, f=0.95, time=-1, bb_only=False):
        '''
        Selects systems in resonance and in broken resonance in the given ratio
        :param f: Fraction of all systems that are supposed to have broken resonance chains
                  type: float
        :param time: Time at which the resonance status is checked. 0 and -1 imply initial and final times resp.
                     type: int or float
        :param bb_only: If the resonance will be checked for building blocks making up the final planets (True or False)
        :return: status in dict format
        '''
        reson = []
        Nbreak = 0
        Nsingle = 0
        for sys, gs in self.systems.items():
            status = gs.get_resonance_status(time=time, bb_only=bb_only)
            if status is None:
                Nsingle += 1
            elif status:
                reson.append(sys)
            else:
                Nbreak += 1
        Nreson_org = len(reson)
        Nreson_req = int(Nbreak * (1 - f) / f)
        if Nreson_req < Nreson_org:
            reson_keep = random.sample(reson, Nreson_req)
            reson_rem = [val for val in reson if val not in reson_keep]
        else:
            reson_rem = []
            Nreson_req = Nreson_org
        for mod, run in reson_rem:
            self.models[mod].pop(run)
            if len(self.models[mod]) == 0:
                self.models.pop(mod)
        return {'broken': Nbreak, 'in-reson_old': Nreson_org, 'in-reson-new': Nreson_req, 'Nsingle': Nsingle}

    def apply(self, func, npool=1):
        if npool < 2:
            for gs in self.systems.values():
                func(gs)
        else:
            pool = Pool(npool)
            # print(self.systems.values())
            for _ in pool.map(func, list(self.systems.values())):
                pass
            pool.close()

    @property
    def available_symbols(self):
        return 'm', 'T', 'a', 'p', 'wf', 'ims', 'sys'

    def get_state(self, *args, bb_only=True, time=-1):
        ''' Same as get_state function of GenSys '''
        dfs = []
        if not args:
            args = self.available_symbols
        for sys, gs in self.systems.items():
            df = gs.get_state(*args, bb_only=bb_only, time=time)
            # print(df)
            if 'sys' in args:
                df['sys'] = [sys for _ in range(df.shape[0])]
            if 'ims' in args:
                df['ims'] = gs.get_impact_phase_summary(bb_only=bb_only)
            dfs.append(df)
        return pd.concat(dfs, axis=0).reset_index(drop=True)

