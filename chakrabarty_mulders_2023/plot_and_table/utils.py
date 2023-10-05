import numpy as np
import pandas as pd

def create_aastex_table(df: pd.DataFrame, loc='!h', ha='c', vbarpos=None, hbarpos=None, caption='', label='', tablenum=None, tablewidth='0pt', comments=None, cellsuff={}):
    s = r'''
\begin{{deluxetable}}{{{coltype}}}[{loc}]
{c_caption}\tablecaption{{{caption}}}
{c_label}\label{{{label}}}
{c_tablenum}\tablenum{{{tablenum}}}
{c_tablewidth}\tablewidth{{{tablewidth}}}
\tablehead{{
{head}
}}
\startdata
{data}
\enddata
{c_comments}\tablecomments{{{comments}}}
\end{{deluxetable}}
    '''
    params = {key: val for key, val in locals().items() if key not in ('s', 'df', 'ha', 'vbarpos', 'hbarpos', 'cellsuff')}
    nspace = s.split('\n')[1].find(r'\begin')
    for key in list(params):
        if params[key] is None:
            params['c_'+key] = '%'
            params[key] = ''
        else:
            params['c_'+key] = ''
    if len(ha) != df.shape[1]:
        ha = ha[0] * df.shape[1]
    if not vbarpos:
        params['coltype'] = ha
    elif vbarpos == 'o':
        params['coltype'] = '|' + ha + '|'
    elif vbarpos == 'lo':
        params['coltype'] = '|' + ha
    elif vbarpos == 'ro':
        params['coltype'] = ha + '|'
    elif vbarpos == 'i':
        params['coltype'] = '|'.join(ha)
    elif vbarpos == 'a':
        params['coltype'] = '|' + '|'.join(ha) + '|'
    elif hasattr(vbarpos, '__len__'):
        if type(vbarpos) == str and vbarpos.isdecimal():
            vbarpos = [int(ch) for ch in vbarpos]
        params['coltype'] = ''.join(['|'*vbarpos[i]+ha[i] for i in range(df.shape[1])]) + '|'*vbarpos[-1]
    nhead = max([col.count('\n')+1 for col in df.columns])
    head = []
    for i in range(nhead):
        head.append(' & '.join([r"\colhead{{{}}}".format(col.split('\n')[i] if len(col.split('\n')) > i else '') for col in df.columns]))
    params['head'] = (r' \\' + '\n' + ' '*nspace).join(head)
    data = []
    for r, row in df.iterrows():
        values = list(row.values)
        for v in range(len(values)):
            if values[v] and not np.isscalar(values[v]): # and all([('int' in str(type(val)) or 'float' in str(type(val))) for val in values[v]])
                if len(values[v]) == 2:
                    values[v] = f'${values[v][0]} \pm {values[v][1]}$'
                elif len(values[v]) == 3:
                    values[v] = f'${values[v][0]}_{{{values[v][1]}}}^{{+{values[v][2]}}}$'
            else:
                values[v] = str(values[v])
            if (r, df.columns[v]) in cellsuff:
                values[v] += cellsuff[(r, df.columns[v])]
            if r in cellsuff:
                values[v] += cellsuff[r]
            if df.columns[v] in cellsuff:
                values[v] += cellsuff[df.columns[v]]
        # print([type(val) for val in values])
        print(values)
        data.append(' & '.join(values))
    params['data'] = (r' \\' + '\n' + ' '*nspace).join(data)
    return s.format(**params)


if __name__ == '__main__':
    df = pd.DataFrame({'a\nbc': [(1, 0.5),2,3], 'xyz': [(10, 1, 2), 20, 30]})
    print(create_aastex_table(df, cellsuff={'a\nbc': r' \%'}, vbarpos=[0,1,2]))