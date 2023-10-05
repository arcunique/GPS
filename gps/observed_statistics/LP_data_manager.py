import os
import pandas as pd
from .obsdata_manager import ObsStats


LPfilepath = r"D:\Modeling\codes\Genesis\data\luque&palle22\science.abl7164_data_s1"

class LPStats(ObsStats):

    def __init__(self, df=None, which='imp', loadfrom='', rearrange=True, keepall=False, pd_kwargs={}):
        if df is None:
            if loadfrom:
                df = pd.read_csv(loadfrom, **pd_kwargs)
            elif which:
                if which in ('imp', 'improved'):
                    df = pd.read_csv(os.path.join(LPfilepath, "STPM_improved_210721.csv"))
                elif which in ('org', 'original'):
                    df = pd.read_csv(os.path.join(LPfilepath, "STPM_original_210721.csv"))
        super().__init__(df=df)
        if not self.df.empty:
            if rearrange:
                self.rearrange(keepall)
            self.df.reset_index(drop=True, inplace=True)

    def rearrange(self, keepall=False):
        renmap = {'Star': 'host', 'R_Rterra': 'rp', 'M_Mterra': 'm', "Porb_d": 'p', "rho_gcm-3": 'dp',
                  'Teq_K': 'T', "Teff_K": 'teff', "R_Rsol": 'rs', "M_Msol": 'ms', "SpT": 'sp', "ParameterRef": 'ref'}
        self.df['name'] = self.df['Star'] + ' ' + self.df['Planet']
        self.df['namealt'] = self.df['AltName'] + ' ' + self.df['Planet']
        for col in list(renmap):
            if 'eu' + col in self.df:
                renmap['eu' + col] = renmap[col] + 'e+'
            if 'ed' + col in self.df:
                renmap['ed' + col] = renmap[col] + 'e-'
            if 'e' + col in self.df:
                renmap['e' + col] = renmap[col] + 'e'
        self.df = self.df.rename(columns=renmap)
        for key in ('dp', 'dpe+', 'dpe-'):
            self.df[key] /= 5.513 # ME / (4/3 * pi * RE**3)
        if not keepall:
            self.df = self.df[['name', 'namealt'] + list(renmap.values())]
