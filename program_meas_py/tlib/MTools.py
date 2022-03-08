# -*- coding: utf-8 -*-

import os, sys
import numpy as np
import types


def get_limits(indg, diff=1):

    len_indg = len(indg)

    if len_indg == 0:
        keg = 0
        beg = np.zeros((0, 2), dtype=np.int64)
    elif len_indg == 1:
        keg = 1
        beg = indg*np.ones((1, 2), dtype=np.int64)
    else:
        kk = indg[1:]-indg[:-1]
        kk = np.where(kk>diff)[0]
        if len(kk)==0:
            keg = 1
            beg = np.array([indg[0], indg[-1]]).reshape((1,2))
        else:
            yy1 = list(indg[kk+1])
            yy1.insert(0, indg[0])
            yy2 = list(indg[kk])
            yy2.append(indg[-1])
            keg = len(yy1)
            beg = np.zeros((keg, 2), dtype=np.int64)
            beg[:,0] = yy1
            beg[:,1] = yy2

    return keg, beg


def full_record(data1, vempty=88888):

#    print(type(vempty))

    dl = len(data1)
    data2 = data1[:]

    if var_len(vempty) == 1:
        ind = np.where(data1 >= vempty)[0]
        ke, be = get_limits(ind)
    elif type(vempty) == list:
        be = np.array(vempty)
    else:
        be = vempty
    sz_be = be.shape
    if sz_be[1] != 2:
        return
    ke = sz_be[0]

    if ke==1 and be[0,0]==0 and be[0,1]==len(data1):
        data2 = np.zeros_like(data1)
        return

    for i in range(ke):
        if i==0:
            if be[i, 0] == 0:
                tmp = data1[be[i, 1] + 1]
                data2[be[i, 0]:be[i, 1]+1] = tmp
                continue
        if i==ke-1:
            if be[i, 1] == dl-1:
                tmp = data1[be[i, 0] - 1]
                data2[be[i, 0]:be[i, 1]+1] = tmp
                continue

        i1 = be[i, 0] - 1
        i2 = be[i, 1] + 1
        data2[i1:i2+1] = np.linspace(data1[i1], data1[i2], i2 - i1+1)

    return data2
    
    

def var_class(x):

    return x.__class__.__name__


def var_len(x):

    name = var_class(x)
    if (name == 'int') or (name == 'float') or (name == 'float64') or (name == 'complex'):
        ln = 1
    elif name == 'list':
        ln = len(x)
    elif name == 'ndarray':
        ln = x.size
    else:
        ln = 0
    
    return ln