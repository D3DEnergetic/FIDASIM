#!/bin/sh
"exec" "$FIDASIM_DIR/deps/python" "$0" "$@"
# -*- coding: utf-8 -*-

"""
Adapted from pyFIDASIM/vmec_read.py
    @author: micha
"""

import numpy as np

from scipy.io import netcdf

def read_vmec(file_name):
    """
    """
    try:
        f = netcdf.netcdf_file(file_name, mode='r', version=4)
    except IOError:
        print(f'File does not exist: {file_name}')
        raise

    ns = f.variables['ns'].data

    xm = f.variables['xm'].data
    xn = f.variables['xn'].data
    md = len(xm)

    xm_nyq = f.variables['xm_nyq'].data
    xn_nyq = f.variables['xn_nyq'].data
    md_nyq = len(xm_nyq)

    ds = 1. / (ns - 1)

    fourierAmps = {
            'R'        : f.variables['rmnc'].data,
            'Z'        : f.variables['zmns'].data,
            'Lambda'   : f.variables['lmns'].data,
            'Bu'       : f.variables['bsubumnc'].data,
            'Bv'       : f.variables['bsubvmnc'].data,
            'Bs'       : f.variables['bsubsmns'].data,
            'Bmod'     : f.variables['bmnc'].data,
            'Jacobian' : f.variables['gmnc'].data,
            'dR_ds'    : np.gradient(f.variables['rmnc'].data, ds, axis=0),
            'dR_du'    : -f.variables['rmnc'].data * xm,
            'dR_dv'    : f.variables['rmnc'].data * xn,
            'dZ_ds'    : np.gradient(f.variables['zmns'].data, ds, axis=0),
            'dZ_du'    : f.variables['zmns'].data * xm,
            'dZ_dv'    : -f.variables['zmns'].data * xn
            }
    
    if md == md_nyq:
        nyq_limit = False

        cos_keys = ['R', 'Jacobian', 'Bu', 'Bv', 'Bmod']
        sin_keys = ['Z', 'Lambda', 'Bs']

        cos_nyq_keys = []
        sin_nyq_keys = []
    else:
        nyq_limit = True

        cos_keys = ['R', 'Jacobian']
        sin_keys = ['Z', 'Lambda']

        cos_nyq_keys = ['Bu', 'Bv', 'Bmod']
        sin_nyq_keys = ['Bs']

    cos_keys.extend(['dR_ds', 'dZ_du', 'dZ_dv'])
    sin_keys.extend(['dR_du', 'dR_dv', 'dZ_ds'])

    f.close()

    return {'fourierAmps':fourierAmps,
            'ns':ns, 'xm':xm, 'xn':xn, 'md':md,
            'xm_nyq':xm_nyq, 'xn_nyq':xn_nyq, 'md_nyq':md_nyq,
            'nyq_limit':nyq_limit, 'cos_keys':cos_keys, 'sin_keys':sin_keys,
            'cos_nyq_keys':cos_nyq_keys, 'sin_nyq_keys':sin_nyq_keys,
            'keys':['R', 'Z', 'Bs', 'Bv', 'Bu', 'dR_ds', 'dR_dv', 'dR_du', 'dZ_ds', 'dZ_dv', 'dZ_du']}
