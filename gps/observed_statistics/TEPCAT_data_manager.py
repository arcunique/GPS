import os
import pandas as pd
from .obsdata_manager import ObsStats
from .utils import get_value_with_error
from source.gps.utils import per2sma, sma2per
from source.gps.constants import RJ, RE, MJ, ME

TEPCATfilepath = r"D:\Modeling\codes\GPS\source\data\TEPCAT"

def remove_pre0_from_name(name):
    nameparts = name.split()
    for i in range(len(nameparts)):
        if nameparts[i].isdecimal():
            nameparts[i] = str(int(nameparts[i]))
        elif '-' in nameparts[i]:
            nameparts[i] = '-'.join([str(int(part)) if part.isdecimal() else part for part in nameparts[i].split('-')])
    return ' '.join(nameparts)

class TEPCATStats(ObsStats):

    def __init__(self, df=None, which='spec', loadfrom='', suffix='', check_consistency=True, rearrange=True, rpmax=6, keepall=False, pd_kwargs={}, consistency_kwargs={}):
        if df is None:
            if loadfrom:
                df = pd.read_csv(loadfrom, **pd_kwargs)
            elif which:
                if which in ('spec', 'specific'):
                    df = pd.read_csv(os.path.join(TEPCATfilepath, f"allplanets-csv{suffix}.csv"))
                    df1 = pd.read_csv(os.path.join(TEPCATfilepath, f"allinfo-csv{suffix}.csv"))
                    df = df.merge(df1[['System', 'Period(day)', 'Perioderr']], how='left', on='System')
                elif which == 'all':
                    df = pd.read_csv(os.path.join(TEPCATfilepath, f"allinfo-csv{suffix}.csv"))
                    df1 = pd.read_csv(os.path.join(TEPCATfilepath, f"allplanets-csv{suffix}.csv"))
                    df = df.merge(df1[['System', 'Recent_reference']], how='left', on='System')
        super().__init__(df=df)
        self.rename_columns()
        if not self.df.empty:
            if rpmax is not None:
                self.df = self.df[self.df["R_b"] <= rpmax]
            if check_consistency:
                for key in ('null', 'errors'):
                    consistency_kwargs.setdefault(key, True)
                self.check_consistency(**consistency_kwargs)
            if rearrange:
                self.rearrange(keepall)
            self.df.reset_index(drop=True, inplace=True)

    def check_consistency(self, null=True, errors=True):
        mcols = ["R_b", "M_b", "Period(day)"]
        mcolsep = ['R_berru', 'M_berru', 'Period(day)err']
        mcolsen = ['R_berrd', 'M_berrd']
        for col in self.df:
            if self.df[col].dtype.kind in 'iuf':
                 self.df.loc[self.df[col] == -1, col] = None
            if self.df[col].dtype == object:
                 self.df.loc[self.df[col] == -1, col] = ''
        if null:
            rows = self.df['a(AU)'].isnull()
            if sum(rows) > 0:
                self.df.loc[rows, 'a(AU)'], self.df.loc[rows, 'a(AU)erru'] = get_value_with_error(per2sma,
                                                (self.df['Period(day)'][rows].values, self.df['M_A'][rows].values),
                                                (self.df['Period(day)err'][rows].values,
                                                 (self.df['M_Aerru'][rows].values + self.df['M_Aerrd'][rows].values) / 2))
                self.df.loc[rows, 'a(AU)errd'] = self.df.loc[rows, 'a(AU)erru'].values
            rows = self.df['Period(day)'].isnull()
            if sum(rows) > 0:
                self.df.loc[rows, 'Period(day)'], self.df.loc[rows, 'Period(day)err'] = get_value_with_error(sma2per,
                                                 (self.df['a(AU)'][rows].values, self.df['M_A'][rows].values),
                                                 ((self.df['a(AU)erru'][rows].values + self.df['a(AU)errd'][rows].values) / 2,
                                                  (self.df['M_Aerru'][rows].values + self.df['M_Aerrd'][rows].values) / 2))
            self.df = self.df[(self.df[mcols] != 0).all(axis=1) & ~self.df[mcols].isnull().any(axis=1)]
        if errors:
            self.df = self.df[~self.df[mcolsep].isnull().any(axis=1) | ~self.df[mcolsen].isnull().any(axis=1)]

    def rename_columns(self):
        columns = self.df.columns.values
        renmap = {}
        i = 0
        while i < len(columns):
            if 'erru' in columns[i] or 'err' in columns[i]:
                renmap[columns[i]] = columns[i - 1] + 'err'
                if i<len(columns)-1 and ('errd' in columns[i+1] or 'err' in columns[i+1]):
                    renmap[columns[i]] = columns[i - 1] + 'erru'
                    renmap[columns[i + 1]] = columns[i - 1] + 'errd'
                    i += 1
            i += 1
        self.df.rename(columns=renmap, inplace=True)

    def rearrange(self, keepall=False):
        renmap = {'Period(day)': 'p', 'a(AU)': 'a', 'M_b': 'm', 'R_b': 'rp', 'g_b': 'g', 'rho_b': 'dp', 'Teq': 'T',
                  'Teff': 'teff', 'M_A': 'ms', 'R_A': 'rs', 'loggA': 'logg', '[Fe/H]': 'met', 'RA(deg)': 'ra',
                  'Dec(deg)': 'dec', 'Recent_reference': 'ref'}
        self.df['name'] = self.df['System'].str.replace('_', ' ')
        self.df['name'] = self.df['name'].apply(lambda x: x + ' b' if not x[-1].isalpha() or not x[-1].islower()
                            else x[:-1] + ' ' + x[-1] if x[-2] != ' ' else x)
        self.df['host'] = self.df['name'].str[:-2]
        for col in list(renmap):
            if col+'erru' in self.df:
                renmap[col+'erru'] = renmap[col] + 'e+'
            if col+'errd' in self.df:
                renmap[col+'errd'] = renmap[col] + 'e-'
            if col+'err' in self.df:
                renmap[col+'err'] = renmap[col] + 'e'
        self.df = self.df.rename(columns=renmap)
        for key in ('m', 'me+', 'me-'):
            self.df[key] *= MJ / ME
        for key in ('rp', 'rpe+', 'rpe-'):
            self.df[key] *= RJ/ RE
        for key in ('dp', 'dpe+', 'dpe-'):
            self.df[key] /= 5.513 # ME / (4/3 * pi * RE**3)
        if not keepall:
            self.df = self.df[['name', 'host'] + [col for col in renmap.values() if col in self.df]]
        self.df['name'] = self.df['name'].apply(remove_pre0_from_name)




if __name__ == '__main__':
    pd.set_option('display.max_columns', None)
    TS = TEPCATStats(which='spec')
    print(TS.df[['name', 'm']][TS.df['m'] == 0])


