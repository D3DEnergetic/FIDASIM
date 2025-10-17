#!/bin/sh
"exec" "$FIDASIM_DIR/deps/python" "$0" "$@"
# -*- coding: utf-8 -*-

"""
read_vmec adapted from:
    pyFIDASIM/vmec_read.py (author: micha)
    matlabVMEC/VMEC/read_vmec.m (author: S. Lazerson)
    STELLOPT/LIBSTELL/vmec.py and STELLOPT/LIBSTELL/libstell.py
    STELLOPT/LIBSTELL/Sources/Modules/read_wout_mod.f90
"""

import os
import numpy as np

from scipy.io import netcdf_file

def read_vmec(file_name):
    """
    #+#read_vmec
    #+ Read in vmec data from `file_name`
    #+***
    #+##Input Arguments
    #+    **file_name**: Name of the VMEC .wout file
    #+
    #+##Output Arguments
    #+    **wout**: VMEC dictionary
    #+
    #+##Example Usage
    #+```python
    #+>>> wout = read_vmec(file_name)
    #+```
    """
    if file_name[-3:] == '.nc': # netCDF file
        f = read_vmec_nc(file_name)
        nc_flag = True
    else: # text file
        f = read_vmec_txt(file_name)
        nc_flag = False

    if 'lrfplogical' not in f:
        f['lrfplogical'] = 0

    if 'xm_nyq' not in f:
        f['xm_nyq'] = f['xm']

    if 'xn_nyq' not in f:
        f['xn_nyq'] = f['xn']

    if 'mnmax_nyq' not in f:
        f['mnmax_nyq'] = f['mnmax']

    for key in ['buco', 'bvco', 'vp', 'specw']:
        f[key] = h2f(f[key])
    
    try:
        f['overr'] = h2f(f['overr'])
        lasym = f['lasym']
    except:
        f['over_r'] = h2f(f['over_r'])
        lasym = f['lasym__logical__']
    
    for mn in range(f['mnmax']):
        f['lmns'][:, mn] = h2f(f['lmns'][:, mn])
    for mn in range(f['mnmax_nyq']):
        for key in ['bmnc', 'gmnc', 'bsupumnc', 'bsupvmnc', 'bsubsmns', 'bsubumnc', 'bsubvmnc']:
            f[key][:, mn] = h2f(f[key][:, mn])

    if lasym:
        for mn in range(f['mnmax']):
            f['lmnc'][:, mn] = h2f(f['lmnc'][:, mn])
        for mn in range(f['mnmax_nyq']):
            for key in ['bmns', 'gmns', 'bsupumns', 'bsupvmns', 'bsubsmnc', 'bsubumns', 'bsubvmns']:
                f[key][:, mn] = h2f(f[key][:, mn])

    #f['eplasma'] = 1.5 * 4 * (np.pi ** 2) * sum(f['vp'] * f['pres'])

    f['mn00'] = None
    for mn in range(f['mnmax']):
        if f['xm'][mn] == 0 and f['xn'][mn] == 0:
            f['mn00'] = mn


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

    return {'f':f, 'fourier_amps':fourier_amps,
            'ns':ns, 's_dom':np.linspace(0, 1, ns),
            'xm':xm, 'xn':xn, 'md':md,
            'xm_nyq':xm_nyq, 'xn_nyq':xn_nyq, 'md_nyq':md_nyq,
            'nyq_limit':nyq_limit, 'cos_keys':cos_keys, 'sin_keys':sin_keys,
            'cos_nyq_keys':cos_nyq_keys, 'sin_nyq_keys':sin_nyq_keys,
            'data_source':os.path.abspath(file_name),
            'keys':['R', 'Z', 'Bs', 'Bv', 'Bu', 'dR_ds', 'dR_dv', 'dR_du', 'dZ_ds', 'dZ_dv', 'dZ_du']}

def read_vmec_nc(file_name):
    """
    #+#read_vmec_nc
    #+ Read in vmec data from `file_name` NETCDF file
    #+***
    #+##Input Arguments
    #+    **file_name**: Name of the VMEC .wout NETCDF file
    #+
    #+##Output Arguments
    #+    **wout**: VMEC dictionary
    #+
    #+##Example Usage
    #+```python
    #+>>> wout = read_vmec_nc(file_name)
    #+```
    """
    wout_netcdf = netcdf_file(file_name, mode='r', version=4)
    f = {}
    for key, var in wout_netcdf.variables.items():
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
    #+#read_vmec_txt
    #+ Read in vmec data from `file_name` text file
    #+***
    #+##Input Arguments
    #+    **file_name**: Name of the VMEC .wout text file
    #+
    #+##Output Arguments
    #+    **wout**: VMEC dictionary
    #+
    #+##Example Usage
    #+```python
    #+>>> wout = read_vmec_txt(file_name)
    #+```
    """
    print(f'Reading VMEC data from text file {os.path.abspath(file_name)}')
    data = read_vmec_txt_data(file_name)
    vers = data[0]
    f = {'version':vers, 'lfreeb':0}
    del data[0] # I don't think this is the best way to do this - WHayashi
    
    mu0 = 4e-7 * np.pi

    keys = {'wb':float,
            'wp':float,
            'gamma':float,
            'pfac':float,
            'rmax_surf':float,
            'rmin_surf':float,
            'zmax_surf':float,
            'nfp':int,
            'ns':int,
            'mpol':int,
            'ntor':int,
            'mnmax':int,
            'mnmax_nyq':int,
            'itfsq':int,
            'niter':int,
            'iasym':int,
            'ireconstruct':int,
            'ierr_vmec':int,
            'imse':int,
            'itse':int,
            'nbsets':int,
            'nobd':int,
            'nextcur':int,
            'nstore_seq':int
            }
    
    if vers <= 8.00:
        del keys['mnmax_nyq']
    if vers <= 6.50:
        del keys['zmax_surf']
    if vers <= 6.20:
        del keys['nstore_seq']
        f['nstore_seq'] = 100
    if vers <= 5.10:
        for key in ['ierr_vmec', 'rmax_surf', 'rmin_surf']:
            del keys[key]

    for key, type_id in keys.items():
        f[key] = type_id(data[0])
        del data[0]

    f['lasym'] = bool(f['iasym'] > 0)

    if f['ierr_vmec'] and f['ierr_vmec'] != 4:
        raise Exception(f'ierr_vmec: {f["ierr_vmec"]}')

    if f['nbsets'] > 0:
        f['nbfld'] = np.array(data[:f['nbsets']], dtype=int)
        del data[:f['nbsets']]

    f['mgrid_file'] = str(data[0])
    del data[0]

    if vers <= 8.00:
        if f['lasym']:
            d1_x, d2_x = 16, 14
        else:
            d1_x, d2_x = 13, 11
        data1 = np.array(data[:d1_x*f['mnmax']]).reshape((f['mnmax'], d1_x)).T
        del data[:data1.size]
        data2 = np.array(data[:+d2_x*f['mnmax']*(f['ns']-1)]).reshape((f['mnmax']*(f['ns']-1), d2_x)).T
        del data[:data2.size]

        f['xm'] = np.array(data1[0, :], dtype=int)
        f['xn'] = np.array(data1[1, :], dtype=int)

        keys = ['rmnc', 'zmns', 'lmns', 'bmnc', 'gmnc', 'bsubumnc', 'bsubvmnc', 'bsubsmns', 'bsupumnc', 'bsupvmnc', 'currvmnc']
        if f['lasym'] and vers > 6.50:
            keys.extend(['rmns', 'zmnc', 'lmnc'])
        for ikey, key in enumerate(keys):
            f[key] = np.vstack([data1[ikey+2, :], np.array(data2[ikey, :]).reshape((f['ns']-1, f['mnmax']))])
    else:
        xyshape1, xyshape2 = (f['ns'], f['mnmax']), (f['ns'], f['mnmax_nyq'])
        f['xm'], f['xn'] = np.zeros(f['mnmax']), np.zeros(f['mnmax'])
        f['rmnc'], f['zmns'], f['lmns'] = np.zeros(xyshape1), np.zeros(xyshape1), np.zeros(xyshape1)
        f['xm_nyq'], f['xn_nyq'] = np.zeros(f['mnmax_nyq']), np.zeros(f['mnmax_nyq'])
        f['bmnc'], f['gmnc'] = np.zeros(xyshape2), np.zeros(xyshape2)
        f['bsubumnc'], f['bsubvmnc'], f['bsubsmns'], f['bsupumnc'], f['bsupvmnc'] = np.zeros(xyshape2), np.zeros(xyshape2), np.zeros(xyshape2), np.zeros(xyshape2), np.zeros(xyshape2)
        if f['lasym']:
            f['rmns'], f['zmnc'], f['lmnc'] = np.zeros(xyshape1), np.zeros(xyshape1), np.zeros(xyshape1)
            f['bmns'], f['gmns'] = np.zeros(xyshape2), np.zeros(xyshape2)
            f['bsubumns'], f['bsubvmns'], f['bsubsmnc'], f['bsupumns'], f['bsupvmns'] = np.zeros(xyshape2), np.zeros(xyshape2), np.zeros(xyshape2), np.zeros(xyshape2), np.zeros(xyshape2)
        for i in range(f['ns']):
            for j in range(f['mnmax']):
                if i == 0:
                    for key in ['xm', 'xn']:
                        f[key][j] = data[0]
                        del data[0]
                for key in ['rmnc', 'zmns', 'lmns']:
                    f[key][i, j] = data[0]
                    del data[0]
                if f['lasym']:
                    for key in ['rmns', 'zmnc', 'lmnc']:
                        f[key][i, j] = data[0]
                        del data[0]
            for j in range(f['mnmax_nyq']):
                if i == 0:
                    for key in ['xm_nyq', 'xn_nyq']:
                        f[key][j] = data[0]
                        del data[0]
                for key in ['bmnc', 'gmnc', 'bsubumnc', 'bsubvmnc', 'bsubsmns', 'bsupumnc', 'bsupvmnc']:
                    f[key][j] = data[0]
                    del data[0]
                if f['lasym']:
                    for key in ['bmns', 'gmns', 'bsubumns', 'bsubvmns', 'bsubsmnc', 'bsupumns', 'bsupvmns']:
                        f[key][i, j] = data[0]
                        del data[0]
        f['mnyq'] = max(f['xm_nyq'])
        f['nnyq'] = max(f['xn_nyq'])

        f['currvmnc'] = np.zeros((f['ns'], f['mnmax_nyq']))
        f['currumnc'] = np.zeros((f['ns'], f['mnmax_nyq']))
        ohs = f['ns'] - 1
        hs = 1. / ohs
        ns = f['ns']
        shalf, sfull = np.zeros(ns), np.zeros(ns)
        for i in range(f['ns']):
            shalf[i] = np.sqrt(hs * (i - 1.5))
            sfull[i] = np.sqrt(hs * (i - 1))
        js1 = np.arange(2, ns)
        js = np.arange(1, ns-1)
        for mn in range(f['mnmax_nyq']):
            if f['xm_nyq'] % 2 == 1:
                t1 = 0.5 * (shalf[js1] * f['bsubsmns'][js1, mn] + shalf[js] * f['bsubsmns'][hs, mn]) / sfull[js]
                bu0 = f['bsubumnc'][js, mn] / shalf[js]
                bu1 = f['bsubumnc'][js1, mn] / shalf[js1]
                t2  = ohs * (bu1-bu0) * sfull[js] + 0.25 * (bu0+bu1) / sfull[js]
                bv0 = f['bsubvmnc'][js, mn] / shalf[js]
                bv1 = f['bsubvmnc'][js1, mn] / shalf[js1]
                t3  = ohs * (bv1-bv0) * sfull[js] + 0.25 * (bv0+bv1) / sfull[js]
            else:
                t1  = 0.5 * (f['bsubsmns'][js1, mn] + f['bsubsmns'][js, mn])
                t2  = ohs * (f['bsubumnc'][js1, mn] + f['bsubumnc'][js, mn])
                t3  = ohs * (f['bsubvmnc'][js1, mn] + f['bsubvmnc'][js, mn])
            f['currumnc'][js, mn] = -1 * float[f['xn_nyq'][mn]] * t1 - t3
            f['currvmnc'][js, mn] = -1 * float[f['xn_nyq'][mn]] * t1 + t2

        f['currumnc'][0, :] = 0.
        f['currvmnc'][0, :] = 0.
        for i in range(f['mnmax_nyq']):
            if f['xm_nyq'][i] == 0:
                f['currumnc'][0, i] = 2 * f['currumnc'][1, i] - f['currumnc'][2, i]
                f['currvmnc'][0, i] = 2 * (f['ns'] - 1) * f['bsubumnc'][1, i]
        f['currumnc'][-1, :] = 2 * f['currumnc'][-2, :] - f['currumnc'][-3, :]
        f['currvmnc'][-1, :] = 2 * f['currvmnc'][-2, :] - f['currvmnc'][-3, :]
        f['currumnc'] = f['currumnc'] / mu0
        f['currvmnc'] = f['currvmnc'] / mu0

        if f['lasym']:
            f['currvmns'] = np.zeros((f['ns'], f['mnmax_nyq']))
            f['currumns'] = np.zeros((f['ns'], f['mnmax_nyq']))
            for mn in range(f['mnmax_nyq']):
                if f['xm_nyq'] % 2 == 1:
                    t1  = 0.5 * (shalf[js1] * f['bsubsmnc'][js1, mn] + shalf[js] * f['bsubsmnc'][js, mn]) / sfull[js]
                    bu0 = f['bsubumns'][js, mn] / shalf[js]
                    bu1 = f['bsubumns'][js1, mn] / shalf[js1]
                    t2  = ohs * (bu1-bu0) * sfull[js] + 0.25 * (bu0+bu1) / sfull[js]
                    bv0 = f['bsubvmns'][js, mn] / shalf[js]
                    bv1 = f['bsubvmns'][js1, mn] / shalf[js1]
                    t3  = ohs * (bv1-bv0) * sfull[js] + 0.25 * (bv0+bv1) / sfull[js]
                else:
                    t1  = 0.5 * (f['bsubsmnc'][js1, mn] + f['bsubsmnc'][js, mn])
                    t2  = ohs * (f['bsubumns'][js1, mn] + f['bsubumns'][js, mn])
                    t3  = ohs * (f['bsubvmns'][js1, mn] + f['bsubvmns'][js, mn])
                f['currumns'][js, mn] = -1 * float(f['xn_nyq'][mn]) * t1 - t3
                f['currvmns'][js, mn] = -1 * float(f['xn_nyq'][mn]) * t1 + t2

            f['currumns'][0, :] = 0
            f['currvmns'][0, :] = 0
            for i in range(f['mnmax_nyq']):
                if f['xm_nyq'][i] == 0:
                    f['currumns'][0, i] = 2 * f['currumns'][1, i] - f['currumns'][3, i]
                    f['currvmns'][0, i] = 2 * (f['ns']-1) * f['bsubumns'][1, i]
            f['currumns'][-1, :] = 2 * f['currumns'][-2, :] - f['currumns'][-3, :]
            f['currvmns'][-1, :] = 2 * f['currvmns'][-2, :] - f['currvmns'][-3, :]
            f['currumns'] = f['currumns'] / mu0
            f['currvmns'] = f['currvmns'] / mu0

    if vers > 6.95:
        data2_1 = np.array(data[:f['ns']*6]).reshape((f['ns'], 6)).T
        del data[:data2_1.size]
        keys = ['iotaf', 'presf', 'phipf', 'phi', 'jcuru', 'jcurv']
        for ikey, key in enumerate(keys):
            f[key] = data2_1[ikey, :]

    data3_shape = [0, 0]
    if vers <= 6.05 or (vers > 6.20 and vers <= 6.50):
        data3_shape[0] = int(f['ns']/2)
    else:
        data3_shape[0] = f['ns'] - 1
    if vers <= 6.05:
        data3_shape[1] = 12
    elif vers > 6.05 and vers <= 6.95:
        data3_shape[1] = 13
    elif vers > 6.95:
        data3_shape[1] = 10
    data3 = np.array(data[:np.prod(data3_shape)]).reshape(data3_shape).T
    del data[:data3.size]
    keys = ['iotas', 'mass', 'pres', 'beta_vol', 'phip', 'buco', 'bvco', 'phi', 'vp', 'overr', 'jcuru', 'jcurv', 'specw']
    if vers <= 6.05:
        keys.remove('beta_vol')
    elif vers > 6.95:
        for key in ['phi', 'jcuru', 'jcurv']:
            keys.remove(key)
    for ikey, key in enumerate(keys):
        f[key] = data3[ikey, :]

    keys = ['aspect', 'betatot', 'betapol', 'betator', 'betaxis', 'b0', 'isigna', 'input_extension', 'IonLarmor', 'VolAvgB', 'RBtor0', 'RBtor', 'Itor', 'Aminor', 'Rmajor', 'Volume']
    if vers <= 6.05:
        for key in ['isigna', 'input_extension', 'IonLarmor', 'VolAvgB', 'RBtor0', 'RBtor', 'Itor', 'Aminor', 'Rmajor', 'Volume']:
            keys.remove(key)
    for key in keys:
        f[key] = data[0]
        del data[0]
    
    data4 = np.array(data[:6*(f['ns']-2)]).reshape((f['ns']-2, 6)).T
    del data[:data4.size]
    keys = ['Dmerc', 'Dshear', 'Dwell', 'Dcurr', 'Dgeod', 'equif']
    for ikey, key in enumerate(keys):
        f[key] = data4[ikey, :]

    f['curlabel'] = np.zeros(f['nextcur'])
    if f['nextcur'] > 0:
        f['lfreeb'] = 1
        f['extcur'] = data[:f['nextcur']]
        del data[:f['nextcur']]
        """
        Original matlabVMEC code:
        fscanf(fid, '\n')
        rem=f.nextcur
        j-0
        while rem > 0
            line=fgetl(fid)
            fscanf(fid, '\n')l
            test=line(1)
            index=findstr(line,test)
            for i=1:size(index,2)/2
                f.curlabel{i+j}=strtrim(line(index(2*i-1)+1:index(2*i)-1))
            end
            j=j+size(index,2)/2
            rem=rem-size(index,2)/2
        end
        """
    data5 = np.array(data[:2*f['nstore_seq']]).reshape((f['nstore_seq'], 2)).T
    del data[:data5.size]
    keys = ['fsqt', 'wdot']
    for ikey, key in enumerate(keys):
        f[key] = data5[ikey, :]

    if vers > 6.05:
        data6 = np.array(data[:2*f['ns']]).reshape((f['ns'], 2)).T
        del data[:data6.size]
        keys = ['jdotb', 'bdotgradv']
        for ikey, key in enumerate(keys):
            f[key] = data6[ikey, :]

    if f['ireconstruct'] > 0:
        if f['imse'] >= 2 or f['itse'] > 0:
            keys = {'tswgt':float, 'msewgt':float, 'isnodes':int}
            for key, type_id in keys.items():
                f[key] = type_id(data[0])
                del data[0]

            data7_1 = np.array(data[:3*f['isnodes']]).reshape((f['isnodes'], 3)).T
            del data[:3*f['isnodes']]
            keys = ['sknots', 'ystark', 'y2stark']
            for ikey, key in enumerate(keys):
                f[key] = data7_1[ikey, :]

            f['ipnodes'] = float(data[0])
            del data[0]

            data7_2 = np.array(data[:3*f['ipnodes']]).reshape((f['ipnodes'], 3)).T
            del data[:3*f['ipnodes']]
            keys = ['pknots', 'ythom', 'y2thom']
            for ikey, key in enumerate(keys):
                f[key] = data7_2[ikey, :]

            data7_3 = np.array(data[:7*(2*f['ns']-1)]).reshape((2*f['ns']-1, 7)).T
            del data[:7*(2*f['ns']-1)]
            keys = ['anglemse', 'rmid', 'qmid', 'shear', 'presmid', 'alfa', 'curmid']
            for ikey, key in enumerate(keys):
                f[key] = data7_3[ikey, :]

            data7_4 = np.array(data[:3*f['imse']]).reshape((f['imse'], 3)).T
            del data[:3*f['imse']]
            keys = ['rstark', 'datastark', 'qmeas']
            for ikey, key in enumerate(keys):
                f[key] = data7_4[ikey, :]

            data7_5 = np.array(data[:2*f['itse']]).reshape((f['itse'], 2)).T
            del data[:2*f['itse']]
            keys = ['rthom', 'datathom']
            for ikey, key in enumerate(keys):
                f[key] = data7_5[ikey, :]

        if f['nobd'] > 0:
            data7_6 = np.array(data[:3*f['nobd']]).reshape((f['nobd'], 3)).T
            del data[:3*f['nobd']]
            keys = ['dsiext', 'plflux', 'fsiobt']
            for ikey, key in enumerate(keys):
                f[key] = data7_6[ikey, :]

            f['flmwgt'] = data[0]
            del data[0]

        nbfldn = f['nbfld'].sum()
        if nbfldn > 0:
            f['bcoil'] = []
            f['plbfld'] = []
            f['bbc'] = []
            for n in range(f['nbsets']):
                data7_7 = np.array(data[:3*f['nbfld'][n]]).reshape((f['nbfld'][n], 3)).T
                del data[:3*f['nbfld'][n]]
                keys = ['bcoil', 'plbfld', 'bbc']
                for ikey, key in enumerate(keys):
                    f[key].append(data7_7[ikey, :])

            f['bcwgt'] = data[0]
            del data[0]

        keys = {'phidiam':float, 'delphid':float, 'nsets':int, 'nparts':int, 'nlim':int}
        for key, type_id in keys.items():
            f[key] = type_id(data[0])
            del data[0]

        f['nsetsn'] = np.array(data[:f['nsets']], dtype=int)
        del data[:f['nsets']]

        f['pfcspec'] = np.zeros((f['nparts'], max(f['nsetsn']), f['nsets']))
        for k in range(f['nsets']):
            for j in range(f['nsetsn'][k]):
                for i in range(f['nparts']):
                        f['pfcspec'][i, j, k] = data[0]
                        del data[0]

        f['limitr'] = data[:f['nlim']]
        del data[:f['nlim']]
        keys = ['rlim', 'zlim']
        for key in keys:
            f[key] = np.zeros((max(f['limitr']), f['nlim']))
        for j in range(f['nlim']):
            for i in range(f['limitr'][j]):
                f['rlim'][i, j] = data[0]
                del data[0]
                f['zlim'][i, j] = data[0]
                del data[0]

        keys = {'nrgrid':int, 'nzgrid':int, 'tokid':float, 'rx1':float, 'rx2':float, 'zy1':float, 'zy2':float, 'condif':float, 'imatch_phiedge':float}
        for key, type_id in keys.items():
            f[key] = type_id(data[0])
            del data[0]

    if len(data) == 0:
        print(f'All data read in')
    else:
        print(f'Partial data read, # of remaining values: {len(data)}')

    if f['version'] <= 6.05:
        for key in ['mass', 'pres', 'jcuru', 'jcurv', 'jdotb']:
            if key in f:
                f[key] = f[key] / mu0
        f['phi'] = -1 * f['phi']
        
    return dict(sorted(f.items(), key=lambda x:x[0]))

def read_vmec_txt_data(file_name, decode_bytes='utf-8'):
    """
    #+#read_vmec_txt_data
    #+ Helper function for `read_vmec_txt`
    #+ Reads data from file_name into list object
    #+***
    #+##Input Arguments
    #+    **file_name**: Name of the VMEC .wout text file
    #+
    #+##Keyword Arguments
    #+     **decode_bytes**: Encoding standard
    #+
    #+##Output Arguments
    #+    **data**: 1-D list object
    #+
    #+##Example Usage
    #+```python
    #+>>> data = read_vmec_txt_data(file_name, decode_bytes='utf-8')
    #+```
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
            elif isinstance(val, (bytes, np.bytes_)):
                data.append(val.decode(decode_bytes))
            elif isinstance(val, str):
                data.append(val)

    return data

def h2f(var_half):
    """
    #+#h2f
    #+ Helper function for `read_vmec`
    #+ Converts half mesh objects to full mesh objects
    #+***
    #+##Input Arguments
    #+    **var_half**: Half mesh data array
    #+
    #+##Output Arguments
    #+    **temp**: Full mesh data array
    #+
    #+##Example Usage
    #+```python
    #+>>> temp = h2f(var_half)
    #+```
    """
    temp = np.zeros(var_half.size)
    temp[0] = 1.5 * var_half[1] - 0.5 * var_half[2]
    temp[1:-1] = 0.5 * (var_half[1:-1] + var_half[2:])
    temp[-1] = 2.0 * var_half[-1] - 1.0 * var_half[-2]
    return temp
