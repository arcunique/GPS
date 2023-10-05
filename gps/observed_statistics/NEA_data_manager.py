import os.path

try:
    from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
    astroquery_exists = True
except ImportError:
    astroquery_exists = False
import pandas as pd
from .obsdata_manager import ObsStats
from datetime import datetime, timezone, timedelta
from html2text import html2text

cache_path = r"D:\Modeling\codes\Genesis\data\NEA"

def html2title(html_str, sep='[]'):
    s = html2text(html_str)
    return ' '.join( s[s.find(sep[0])+1: s.find(sep[1])].split()).replace('https://ui.adsabs.harvard.edu/abs/', '')

class NEAStats(ObsStats):

    def __init__(self, df=None, cached=True, query=True, loadfrom='', cache_res=True, true_mass=True, rpmax=6, check_consistency=True, rearrange=True, keepall=False, pd_kwargs={}, consistency_kwargs={}, reflink=False):
        if df is None:
            if loadfrom and os.path.exists(loadfrom):
                pd_kwargs.setdefault('comment', '#')
                df = pd.read_csv(loadfrom, **pd_kwargs)
            elif cached and self.cached_file:
                df = pd.read_csv(self.cached_file)
            elif query:
                mcols = ["pl_rade", "pl_bmasse", "pl_orbper", "pl_orbsmax", "pl_eqt", "st_teff", "st_rad", "st_mass", "st_logg"]
                columns = ["pl_name", "hostname", "default_flag", "pl_bmassprov", "pl_controv_flag", "pl_refname"]
                for col in mcols:
                    columns += [col, col+'err1', col+'err2', col+'lim']
                qkwargs = {}
                if true_mass:
                    qkwargs["where"] = "pl_bmassprov like 'Mass'"
                df = NasaExoplanetArchive.query_criteria(table="ps", select=columns, **qkwargs).to_pandas()
                if cache_res:
                    dtnow = datetime.now(timezone(timedelta(hours=-7))).strftime("%Y.%m.%d_%H.%M.%S")
                    df.to_csv(os.path.join(cache_path, f'PS_{dtnow}.csv'), index=False)
        super().__init__(df=df)
        if not self.df.empty:
            if rpmax is not None:
                self.df = self.df[self.df["pl_rade"] <= rpmax]
            if check_consistency:
                for key in ('errors', 'lim', 'default'):
                    consistency_kwargs.setdefault(key, True)
                consistency_kwargs.setdefault('overlap', 1.5)
                self.check_consistency(**consistency_kwargs)
            if rearrange:
                self.rearrange(keepall, reflink)
            self.df.reset_index(drop=True, inplace=True)

    @property
    def cached_file(self):
        files = [os.path.join(cache_path, file) for file in os.listdir(cache_path)]
        if files:
            return max(files, key=os.path.getctime)

    def check_consistency(self, errors=True, lim=True, overlap=1.5, default=True, true_mass=True):
        mcols = {"rp": "pl_rade", "m": "pl_bmasse", "p": "pl_orbper"}
        self.df = self.df[self.df["pl_controv_flag"] == 0]
        if true_mass:
            self.df = self.df[self.df['pl_bmassprov'] == 'Mass']
        if errors:
            if errors in (True, 'all'):
                mcolsep = [col + 'err1' for col in mcols.values()]
                mcolsen = [col + 'err2' for col in mcols.values()]
            else:
                if type(errors) == str:
                    errors = [errors]
                mcolsep = [mcols[key] + 'err1' for key in errors]
                mcolsen = [mcols[key] + 'err2' for key in errors]
            self.df = self.df[~self.df[mcolsep].isnull().any(axis=1) | ~self.df[mcolsen].isnull().any(axis=1)]
        if lim:
            if lim in (True, 'all'):
                mcolslim = [col + 'lim' for col in mcols.values()]
            else:
                if type(lim) == str:
                    lim = [lim]
                mcolslim = [mcols[key] + 'lim' for key in lim]
            self.df = self.df[(self.df[mcolslim] == 0).all(axis=1)]
        if overlap is not None:
            mcols = ["pl_rade", "pl_bmasse"]
            mcolsep = [col + 'err1' for col in mcols]
            mcolsen = [col + 'err2' for col in mcols]
            allowed_planets = []
            # print(self.df)
            for pl, df in self.df.groupby('pl_name'):
                if df["default_flag"].any():
                    defval = df[df["default_flag"] == 1].iloc[0][mcols].values
                    dfndef = df[df["default_flag"] == 0]
                    if not dfndef.empty:
                        allowp = (dfndef[mcols].values + overlap * dfndef[mcolsep].values) >= defval
                        allown = (dfndef[mcols].values + overlap * dfndef[mcolsen].values) <= defval
                        allow = allowp.all() & allown.all()
                    else:
                        allow = True
                else:
                    allow = False
                if allow:
                    allowed_planets.append(pl)
            self.df = self.df[self.df['pl_name'].isin(allowed_planets)]
        if default:
            self.df = self.df[self.df["default_flag"] == 1]

    def rearrange(self, keepall=False, reflink=False):
        renmap = {'pl_name': 'name', 'hostname': 'host', 'pl_rade': 'rp', 'pl_bmasse': 'm', "pl_orbper": 'p', "pl_orbsmax": 'a',
                  'pl_eqt': 'T', "st_teff": 'teff', "st_rad": 'rs', "st_mass": 'ms', "st_logg": 'logg', "pl_refname": 'ref'}
        if reflink:
            self.df['pl_refname'] = self.df['pl_refname'].apply(html2title, args=('()',))
        else:
            self.df['pl_refname'] = self.df['pl_refname'].apply(html2title, args=('[]',))
        for col in list(renmap):
            if col + 'err1' in self.df:
                renmap[col + 'err1'] = renmap[col] + 'e+'
            if col + 'err2' in self.df:
                self.df[col + 'err2'] = - self.df[col + 'err2']
                renmap[col + 'err2'] = renmap[col] + 'e-'
        self.df = self.df.rename(columns=renmap)
        if not keepall:
            self.df = self.df[list(renmap.values())]



if __name__ == '__main__':
    pd.set_option('display.max_columns', None)
    df = NEAStats().filter(rp=(1,5), teff=(4000, 5000)).filter_by_errorpc(rp=8, m=20).dffull
    print(df.head())
    # mcols = ["pl_rade", "pl_bmasse", "pl_orbper", "pl_eqt", "st_teff", "st_rad", "st_mass", "st_logg"]
    # columns = ["pl_name", "hostname", "default_flag", "pl_bmassprov", "pl_controv_flag"]
    # for col in mcols:
    #     columns += [col, col + 'err1', col + 'err2', col + 'lim']
    # df = NasaExoplanetArchive.query_criteria(table="ps", where="pl_bmassprov like 'Mass' AND hostname like 'AU Mic'").to_pandas()
    # print(df)

