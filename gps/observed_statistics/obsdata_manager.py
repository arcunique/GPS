from numpy import inf
from .utils import get_value_with_error
import pandas as pd

class ObsStats:

    columns = ['name', 'rp', 'm', 'p', 'T', 'teff', 'rs', 'ms', 'logg']

    def __init__(self, df=None):
        if df is not None:
            self.df = pd.DataFrame(df)
        else:
            self.df = pd.DataFrame()
        self._df = pd.DataFrame()

    def filter(self, inplace=False, errorpc={}, **kwargs):
        df = self.df.copy()
        for key, val in kwargs.items():
            if key in df:
                if isinstance(val, (int, float, str)):
                    df = df[self.df[key] == val]
                else:
                    val = list(val)
                    if val[0] is None:
                        val[0] = -inf
                    if val[1] is None:
                        val[1] = inf
                    df = df[(df[key] >= val[0]) & (df[key] <= val[1])]
        if inplace:
            self.df = df.copy()
            if errorpc:
                self.filter_by_errorpc(**errorpc, inplace=True)
        else:
            OS = ObsStats(df=df.reset_index(drop=True))
            if errorpc:
                OS = OS.filter_by_errorpc(**errorpc)
            return OS

    def filter_by_errorpc(self, inplace=False, **kwargs):
        df = self.df.copy()
        for key, val in kwargs.items():
            if key in df:
                epc = (df[key+'e+'] + df[key+'e-']) / 2 / df[key] * 100
                df = df[epc <= val]
        if inplace:
            self.df = df.copy()
        else:
            OS = ObsStats(df=df.reset_index(drop=True))
            return OS

    def add_more(self, inplace=False):
        df = self.df.copy()
        df['rpe'] = (df['rpe+'] + df['rpe-']) / 2
        df['me'] = (df['me+'] + df['me-']) / 2
        df['rpe%'] = df['rpe'] / df['rp'] * 100
        df['me%'] = df['me'] / df['m'] * 100
        if 'dp' not in df:
            mask = df['rp'] > 0
        else:
            mask = df['dp'].isnull()
        if not mask.empty:
            func = lambda m, rp: m / rp**3
            df.loc[mask, 'dp'] = df[mask].apply(lambda x: get_value_with_error(func, x[['m', 'rp']], get='val'), axis=1)
            dpe = df[mask].apply(lambda x: get_value_with_error(func, x[['m', 'rp']].values, x[['me', 'rpe']].values, get='err'), axis=1)
            if 'dpe' in df:
                df.loc[mask, 'dpe'] = dpe
            else:
                df.loc[mask, 'dpe+'] = dpe
                df.loc[mask, 'dpe-'] = dpe
        if 'dpe' in df:
            df['dpe+'] = df['dpe']
            df['dpe-'] = df['dpe']
        else:
            df['dpe'] = (df['dpe+'] + df['dpe-']) / 2
        if 'fp' not in df:
            func = lambda m, rp: rp / m ** (1/ 3.7)
            df['fp'] = df.apply(lambda x: get_value_with_error(func, [x['m'], x['rp']], get='val'), axis=1)
            df['fpe'] = df.apply(lambda x: get_value_with_error(func, x[['m', 'rp']].values, x[['me', 'rpe']].values, get='err'), axis=1)
            df['fpe+'] = df['fpe']
            df['fpe-'] = df['fpe']
        for col in ('p', 'T', 'rs', 'ms', 'teff'):
            if col+'e+' in df and col+'e' not in df:
                df[col+'e'] = (df[col+'e+'].values + df[col+'e-'].values) / 2
        if inplace:
            self.df = df.copy()
        else:
            OS = ObsStats(df=df.reset_index(drop=True))
            return OS

    @property
    def dffull(self):
        if self._df.empty:
            self._df = self.add_more().df
        return self._df



