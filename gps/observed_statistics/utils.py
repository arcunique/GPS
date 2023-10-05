import numpy as np


def get_value_with_error(func, x, xe=None, get='both'):
    val = func(*x)
    if xe is None or get == 'val':
        return val
    # print(val.shape)
    err = val * ((np.array(xe)/np.array(x))**2).sum(0)**0.5
    if get == 'err':
        return err
    return val, err
