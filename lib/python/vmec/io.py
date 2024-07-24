#!/bin/sh
"exec" "$FIDASIM_DIR/deps/python" "$0" "$@"
# -*- coding: utf-8 -*-

"""
read_vmec adapted from pyFIDASIM/vmec_read.py (author: micha) and matlabVMEC/VMEC/read_vmec.m (author: S. Lazerson)
"""

import numpy as np

from scipy.io import netcdf_file

def read_vmec(file_name):
    """
    """
    if file_name[-3:] == '.nc': # netCDF file
        f = read_vmec_nc(file_name)
        nc_flag = True
    else: # text file
        f = read_vmec_txt(file_name)
        nc_flag = False

    ns = f['ns']
    xm = f['xm']
    xn = f['xn']
    md = len(xm)

    xm_nyq = f['xm_nyq']
    xn_nyq = f['xn_nyq']
    md_nyq = len(xm_nyq)

    ds = 1. / (ns - 1)

    fourier_amps = {
            'R'        : f['rmnc'],
            'Z'        : f['zmns'],
            'Lambda'   : f['lmns'],
            'Bu'       : f['bsubumnc'],
            'Bv'       : f['bsubvmnc'],
            'Bs'       : f['bsubsmns'],
            'Bmod'     : f['bmnc'],
            'Jacobian' : f['gmnc'],
            'dR_ds'    : np.gradient(f['rmnc'], ds, axis=0),
            'dR_du'    : -f['rmnc'] * xm,
            'dR_dv'    : f['rmnc'] * xn,
            'dZ_ds'    : np.gradient(f['zmns'], ds, axis=0),
            'dZ_du'    : f['zmns'] * xm,
            'dZ_dv'    : -f['zmns'] * xn
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

    """
    msize = f['xm'].max()
    nsize = max(f['xn']/f['nfp']) - min(f['xn']/f['nfp']) + 1
    f['rbc'] = np.zeros((msize+1, nsize, f['ns']))
    f['zbs'] = np.zeros((msize+1, nsize, f['ns']))
    f['lbs'] = np.zeros((msize+1, nsize, f['ns']))

    f['rsc'] = np.zeros((msize+1, nsize, f['ns']))
    f['rus'] = np.zeros((msize+1, nsize, f['ns']))
    f['rvs'] = np.zeros((msize+1, nsize, f['ns']))
    f['zss'] = np.zeros((msize+1, nsize, f['ns']))
    f['zuc'] = np.zeros((msize+1, nsize, f['ns']))
    f['zvc'] = np.zeros((msize+1, nsize, f['ns']))
    """

    return {'fourierAmps':fourier_amps,
            'ns':ns, 'xm':xm, 'xn':xn, 'md':md,
            'xm_nyq':xm_nyq, 'xn_nyq':xn_nyq, 'md_nyq':md_nyq,
            'nyq_limit':nyq_limit, 'cos_keys':cos_keys, 'sin_keys':sin_keys,
            'cos_nyq_keys':cos_nyq_keys, 'sin_nyq_keys':sin_nyq_keys,
            'keys':['R', 'Z', 'Bs', 'Bv', 'Bu', 'dR_ds', 'dR_dv', 'dR_du', 'dZ_ds', 'dZ_dv', 'dZ_du']}


def read_vmec_nc(file_name):
    """
    """
    wout_netcdf = netcdf_file(file_name, mode='r', version=4)
    f = {}
    for key, var in wout_netcdf.variables.items():
        print(key)
        if var.data.ndim == 0:
            val = np.array([var.data])[0]
            if float(val).is_integer():
                type_id = int
            else:
                type_id = float
            f[key] = type_id(val)
        elif isinstance(var.data[0], np.bytes_):
            f[key] = ''.join([val.decode('utf-8') for val in var.data]).strip()
        elif var.data.ndim > 1 and isinstance(var.data[0][0], np.bytes_):
            f[key] = [''.join([val2.decode('utf-8') for val2 in val1]).strip() for val1 in var.data]
        else:
            arr_val = var.data
            if all([float(val).is_integer() for val in arr_val.flatten()]):
                type_id = int
            else:
                type_id = float
            f[key] = arr_val.astype(type_id)
    wout_netcdf.flush()
    wout_netcdf.fp.close()
    wout_netcdf.close()
    return f


def read_vmec_txt(file_name):
    """
    """
    data = read_vmec_txt_data(file_name)
    
    f = {}
    cursor = 1
    #if version > 6.5 and version <= 6.95:
    f['lfreeb'] = 0
    keys = {'wb':float, 'wp':float, 'gamma':float, 'pfac':float, 'rmax_surf':float, 'rmin_surf':float, 'zmax_surf':float, 'nfp':int, 'ns':int, 'mpol':int, 'ntor':int, 'mnmax':int, 'itfsq':int, 'niter':int, 'iasym':int, 'ireconstruct':int, 'ierr_vmec':int, 'imse':int, 'itse':int, 'nbsets':int, 'nobd':int, 'nextcur':int, 'nstore_seq':int}
    for ikey, (key, type_id) in enumerate(keys.items()):
        f[key] = type_id(data[ikey+1])
        cursor += 1

    if f['ierr_vmec'] and f['ierr_vmec'] != 4:
        raise Exception(f'ierr_vmec: {f["ierr_vmec"]}')

    if f['nbsets'] > 0:
        f['nbfld'] = data[cursor:cursor+f['nbsets']].astype(int)
        cursor += f['nbsets']

    if f['iasym'] > 0:
        d1_x, d2_x = 16, 14
    else:
        d1_x, d2_x = 13, 11
    data1 = data[cursor:cursor+d1_x*f['mnmax']].reshape((f['mnmax'], d1_x)).T
    cursor += d1_x*f['mnmax']
    data2 = data[cursor:cursor+d2_x*f['mnmax']*(f['ns']-1)].reshape((f['mnmax']*(f['ns']-1), d2_x)).T
    cursor += d2_x*f['mnmax']*(f['ns']-1)

    f['xm'] = data1[0, :].astype(int)
    f['xn'] = data1[1, :].astype(int)

    keys = ['rmnc', 'zmns', 'lmns', 'bmnc', 'gmnc', 'bsubumnc', 'bsubvmnc', 'bsubsmns', 'bsupumnc', 'bsupvmnc', 'currvmnc']
    if f['iasym'] > 0:
        keys.extend(['rmns', 'zmnc', 'lmnc'])
    for ikey, key in enumerate(keys):
        f[key] = np.vstack([data1[ikey+2, :], data2[ikey, :].reshape((f['ns']-1, f['mnmax']))])

    data3 = data[cursor:cursor+13*(f['ns']-1)].reshape((f['ns']-1, 13)).T
    cursor += 13*(f['ns']-1)
    keys = ['iotas', 'mass', 'pres', 'beta_vol', 'phip', 'buco', 'bvco', 'phi', 'vp', 'overr', 'jcuru', 'jcurv', 'specw']
    for ikey, key in enumerate(keys):
        f[key] = data3[ikey, :]

    curs0 = cursor
    keys = ['aspect', 'betatot', 'betapol', 'betator', 'betaxis', 'b0', 'isigna', 'IonLarmor', 'VolAvgB', 'RBtor0', 'RBtor', 'Itor', 'Aminor', 'Rmajor', 'Volume']
    for ikey, key in enumerate(keys):
        f[key] = data[curs0+ikey]
        cursor += 1
    
    data4 = data[cursor:cursor+6*(f['ns']-2)].reshape((f['ns']-2, 6)).T
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
    data5 = data[cursor:cursor+2*f['nstore_seq']].reshape((f['nstore_seq'], 2)).T
    cursor += 2*f['nstore_seq']
    keys = ['sqt', 'wdot']
    for ikey, key in enumerate(keys):
        f[key] = data5[ikey, :]

    data6 = data[cursor:cursor+2*f['ns']].reshape((f['ns'], 2)).T
    cursor += 2*f['ns']
    keys = ['jdotb', 'bdotgradv']
    for ikey, key in enumerate(keys):
        f[key] = data6[ikey, :]

    if f['ireconstruct'] > 0:
        if f['imse'] >= 2 or f['itse'] > 0:
            curs0 = cursor
            keys = {'twsgt':float, 'msewgt':float, 'isnodes':int}
            for ikey, (key, type_id) in enumerate(keys.items()):
                f[key] = type_id(data[curs0+ikey])
                cursor += 1

            data7_1 = data[cursor:cursor+3*f['isnodes']].reshape((f['isnodes'], 3)).T
            cursor += 3*f['isnodes']
            keys = ['sknots', 'ystark', 'y2stark']
            for ikey, key in enumerate(keys):
                f[key] = data7_1[ikey, :]

            f['ipnodes'] = float(data[cursor])
            cursor += 1

            data7_2 = data[cursor:cursor+3*f['ipnodes']].reshape((f['ipnodes'], 3)).T
            cursor += 3*f['ipnodes']
            keys = ['pknots', 'ythom', 'y2thom']
            for key, ikey in enumerate(keys):
                f[key] = data7_2[ikey, :]

            data7_3 = data[cursor:cursor+7*(2*f['ns']-1)].reshape((2*f['ns']-1, 7)).T
            cursor += 7*(2*f['ns']-1)
            keys = ['anglemse', 'rmid', 'qmid', 'shear', 'presmid', 'alfa', 'curmid']
            for key, ikey in enumerate(keys):
                f[key] = data7_3[ikey, :]

            data7_4 = data[cursor:cursor+3*f['imse']].reshape((f['imse'], 3)).T
            cursor += 3*f['imse']
            keys = ['rstark', 'datastark', 'qmeas']
            for ikey, key in enumerate(keys):
                f[key] = data7_4[ikey, :]

            data7_5 = data[cursor:cursor+2*f['itse']].reshape((f['itse'], 2)).T
            cursor += 2*f['itse']
            keys = ['trhom', 'datathom']
            for ikey, key in enumerate(keys):
                f[key] = data7_5[ikey, :]

        if f['nobd'] > 0:
            data7_6 = data[cursor:cursor+3*f['nobd']].reshape((f['nobd'], 3)).T
            cursor += 3*f['nobd']
            keys = ['dsiext', 'plflux', 'fsiobt']
            for ikey, key in enumerate(keys):
                f[key] = data7_6[ikey, :]

            f['flmwgt'] = data[cursor]
            cursor += 1

        nbfldn = f['nbfld'].sum()
        if nbfldn > 0:
            f['bcoil'] = []
            f['plbfld'] = []
            f['bbc'] = []
            for n in range(f['nbsets']):
                data7_7 = data[cursor:cursor+3*f['nbfld'][n]].reshape((f['nbfld'][n], 3)).T
                cursor += 3*f['nbfld'][n]
                keys = ['bcoil', 'plbfld', 'bbc']
                for ikey, key in enumerate(keys):
                    f[key].append(data7_7[ikey, :])

            f['bcwgt'] = data[cursor]
            cursor += 1

        curs0 = cursor
        keys = {'phidiam':float, 'delphid':float, 'nsets':int, 'nparts':int, 'nlim':int}
        for ikey, (key, type_id) in enumerate(keys.items()):
            f[key] = type_id(data[curs0+ikey])
            cursor += 1

        f['nsetsn'] = data[cursor:cursor+f['nsets']].astype(int)
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
        keys = {'nrgrid':int, 'nzgrid':int, 'tokid':float, 'rx1':float, 'rx2':float, 'zy1':float, 'zy2':float, 'conif':float, 'imatch_phiedge':float}
        for ikey, (key, type_id) in enumerate(keys.items()):
            f[key] = type_id(data[curs0+ikey])
            cursor += 1

    if 'lrfplogical' not in f:
        f['lrfplogical'] = 0

    if 'xm_nyq' not in f:
        f['xm_nyq'] = f['xm']

    if 'xn_nyq' not in f:
        f['xn_nyq'] = f['xn']

    if 'mnmax_nyq' not in f:
        f['mnmax_nyq'] = f['mnmax']

    """
    f = half2fullmesh(f)
    """
    
    return f


def read_vmec_txt_data(file_name):
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

    return data

