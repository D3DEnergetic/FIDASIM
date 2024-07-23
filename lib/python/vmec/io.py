#!/bin/sh
"exec" "$FIDASIM_DIR/deps/python" "$0" "$@"
# -*- coding: utf-8 -*-

"""
read_vmec adapted from pyFIDASIM/vmec_read.py (author: micha) and matlabVMEC/VMEC/read_vmec.m (author: S. Lazerson)
"""

import numpy as np

from scipy.io import netcdf

def read_vmec(file_name):
    """
    """
    if file_name[-3:] == '.nc': # netCDF file
        f = netcdf.netcdf_file(file_name, mode='r', version=4)
    else: # text file
        f = read_vmec_txt(file_name)

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


def read_vmec_txt(file_name):
    """
    """
    with open(file_name, 'r') as file:
        flines = file.readlines()

    version = float(flines[0].strip().split()[-1])
    if ',' in flines[1]:
        sep = ','
        csv = True
    else:
        sep = None
        csv = False

    data = [version]
    for line in flines[1:]:
        for val in line.strip().split(sep=sep):
            if not val.isalpha():
                data.append(float(val))
    data = np.array(data) # numpy array will be easier to manipulate

    f = {}
    cursor = 1
    #if version > 6.5 and version <= 6.95:
    f['lfreeb'] = 0
    keys = ['wb', 'wp', 'gamma', 'pfac', 'rmax_surf', 'rmin_surf', 'zmax_surf', 'nfp', 'ns', 'mpol', 'ntor', 'mnmax', 'itfsq', 'niter', 'iasym', 'ireconstruct', 'ierr_vmec', 'imse', 'itse', 'nbsets', 'nobd', 'nextcur', 'nstore_seq']
    for ikey, key in enumerate(keys):
        f[key] = data[ikey+1]
        cursor += 1

    if f['ierr_vmec'] and f['ierr_vmec'] != 4:
        raise Exception(f'ierr_vmec: {f["ierr_vmec"]}')

    if f['nbsets'] > 0:
        f['nbfld'] = data[cursor:cursor+f['nbsets']]
        cursor += f['nbsets']

    if f['iasym'] > 0:
        d1_x, d2_x = 16, 14
    else:
        d1_x, d2_x = 13, 11
    data1 = data[cursor:cursor+d1_x*f['mnmax']].reshape((d1_x, f['mnmax']))
    cursor += d1_x*f['mnmax']
    data2 = data[cursor:cursor+d2_x*f['mnmax']*(f['ns']-1)].reshape((d2_x, f['mnmax']*(f['ns']-1)))
    cursor += d2_x*f['mnmax']*(f['ns']-1)

    f['xm'] = data1[0, :].astype(int)
    f['xn'] = data1[1, :].astype(int)

    keys = ['rmnc', 'zmns', 'lmns', 'bmnc', 'gmnc', 'bsubumnc', 'bsubvmnc', 'bsubsmns', 'bsupumnc', 'bsupvmnc', 'currvmnc' ]
    if f['iasym'] > 0:
        keys.extend(['rmns', 'zmnc', 'lmnc'])
    for ikey, key in enumerate(keys):
        f[key] = np.vstack([data1[ikey+2, :], data2[ikey, :].reshape((f['mnmax'], f['ns'] - 1)).T]).T

    data3 = data[cursor:cursor+13*(f['ns']-1)].reshape((13, f['ns']-1))
    cursor += 13*(f['ns']-1)
    keys = ['iotas', 'mass', 'pres', 'beta_vol', 'phip', 'buco', 'bvco', 'phi', 'vp', 'overr', 'jcuru', 'jcurv', 'specw']
    for ikey, key in enumerate(keys):
        f[key] = data3[ikey, :]

    curs0 = cursor
    keys = ['aspect', 'betatot', 'betapol', 'betator', 'betaxis', 'b0', 'isigna', 'input_extension', 'IonLarmor', 'VolAvgB', 'RBtor0', 'RBtor', 'Itor', 'Aminor', 'Rmajor', 'Volume']
    for ikey, key in enumerate(keys):
        f[key] = data[curs0+ikey]
        cursor += 1
    
    data4 = data[cursor:cursor+6*(f['ns']-2)].reshape((6, f['ns']-2))
    cursor += 6*(f['ns']-2)
    keys = ['Dmerc', 'Dshear', 'Dwell', 'Dcurr', 'Dgeod', 'equif']
    for ikey, key in enumerate(keys):
        f[key] = data4[ikey, :]

    f['curlabel'] = np.zeros(f['nextcur'])
    if f['nextcur'] > 0:
        f['lfreeb'] = 1
        f['extcur'] = data[cursor:cursor+f['nextcur']]
        cursor += f['nextcur']
        """
        Original matlabVMEC code:
        fscanf(fid, '\n');
        rem=f.nextcur;
        j-0;
        while rem > 0
            line=fgetl(fid);
            fscanf(fid, '\n')l
            test=line(1);
            index=findstr(line,test);
            for i=1:size(index,2)/2;
                f.curlabel{i+j}=strtrim(line(index(2*i-1)+1:index(2*i)-1));
            end
            j=j+size(index,2)/2;
            rem=rem-size(index,2)/2;
        end
        """
    data5 = data[cursor:cursor+4*f['nstore_seq']].reshape((4, f['nstore_seq']))
    cursor += 4*f['nstore_seq']
    keys = ['sqt', 'wdot', 'jdotb', 'bdotgradv']
    for ikey, key in enumerate(keys):
        f[key] = data5[ikey, :]

    if f['ireconstruct'] > 0:
        if f['imse'] >= 2 or f['itse'] > 0:
            curs0 = cursor
            keys = ['twsgt', 'msewgt', 'isnodes']
            for ikey, key in enumerate(keys):
                f[key] = data[curs0+ikey]
                cursor += 1

            data5_1 = data[cursor:cursor+3*f['isnodes']].reshape((3, f['isnodes']))
            cursor += 3*f['isnodes']
            keys = ['sknots', 'ystark', 'y2stark']
            for ikey, key in enumerate(keys):
                f[key] = data5_1[ikey, :]

            f['ipnodes'] = data[cursor]
            cursor += 1

            data5_2 = data[cursor:cursor+3*f['ipnodes']].reshape((3, f['ipnodes']))
            cursor += 3*f['ipnodes']
            keys = ['pknots', 'ythom', 'y2thom']
            for key, ikey in enumerate(keys):
                f[key] = data5_2[ikey, :]

            data5_3 = data[cursor:cursor+7*(2*f['ns']-1)].reshape((7, 2*f['ns']-1))
            cursor += 7*(2*f['ns']-1)
            keys = ['anglemse', 'rmid', 'qmid', 'shear', 'presmid', 'alfa', 'curmid']
            for key, ikey in enumerate(keys):
                f[key] = data5_3[ikey, :]

            data5_4 = data[cursor:cursor+3*f['imse']].reshape((3, f['imse']))
            cursor += 3*f['imse']
            keys = ['rstark', 'datastark', 'qmeas']
            for ikey, key in enumerate(keys):
                f[key] = data5_4[ikey, :]

            data5_5 = data[cursor:cursor+2*f['itse']].reshape((2, f['itse']))
            cursor += 2*f['itse']
            keys = ['trhom', 'datathom']
            for ikey, key in enumerate(keys):
                f[key] = data5_5[ikey, :]

        if f['nobd'] > 0:
            data5_6 = data[cursor:cursor+3*f['nobd']].reshape((3, f['nobd']))
            cursor += 3*f['nobd']
            keys = ['dsiext', 'plflux', 'fsiobt']
            for ikey, key in enumerate(keys):
                f[key] = data5_6[ikey, :]

            f['flmwgt'] = data[cursor]
            cursor += 1

        nbfldn = f['nbfld'].sum()
        if nbfldn > 0:
            f['bcoil'] = []
            f['plbfld'] = []
            f['bbc'] = []
            for n in range(f['nbsets']):
                data5_7 = data[cursor:cursor+3*f['nbfld'][n]].reshape((3, f['nbfld'][n]))
                cursor += 3*f['nbfld'][n]
                keys = ['bcoil', 'plbfld', 'bbc']
                for ikey, key in enumerate(keys):
                    f[key].append(data5_7[ikey, :])

            f['bcwgt'] = data[cursor]
            cursor += 1

        curs0 = cursor
        keys = ['phidiam', 'delphid', 'nsets', 'nparts', 'nlim']
        for ikey, key in enumerate(keys):
            f[key] = data[curs0+ikey]
            cursor += 1

        f['nsetsn'] = data[cursor:cursor+f['nsets']]
        cursor += f['nsets']

        f['pfcspec'] = np.zeros((f['nparts'], max(f['nsetsn']), f['nsets']))
        for k in range(f['nsets']):
            for j in range(f['nsetsn'][k]):
                for i in range(f['nparts']):
                        f['pfcspec'][i, j, k] = data[cursor]
                        cursor += 1

        f['limitr'] = data[cursor:cursor+f['nlim']]
        cursor += f['nlim']
        keys = ['rlim', 'zlim']
        for key in keys:
            f[key] = np.zeros((max(f['limitr']), f['nlim']))
        for j in range(f['nlim']):
            for i in range(f['limitr'][j]):
                f['rlim'][i, j] = data[cursor]
                cursor += 1
                f['zlim'][i, j] = data[cursor]
                cursor += 1

        curs0 = cursor
        keys = ['nrgrid', 'nzgrid', 'tokid', 'rx1', 'rx2', 'zy1', 'zy2', 'conif', 'imatch_phiedge']
        for ikey, key in enumerate(keys):
            f[key] = data[curs0+ikey]
            cursor += 1

    return f
