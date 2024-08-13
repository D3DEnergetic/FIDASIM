#!/bin/sh
"exe" "$FIDASIM_DIR/deps/python" "$0" "$@"
# -*- coding: utf-8 -*-

"""
fourier_transform_3D and Brzp_transform adapted from pyFIDASIM/vmec_read.py
    @author: micha
"""

import numpy as np

from scipy.interpolate import griddata

def fourier_transform_3D(wout, ntheta=16, nphi=20, thetamin=None, thetamax=None, phimin=None, phimax=None):
    """
    #+#fourier_transform_3D
    #+ Performs 3D Fourier transform on objects with key in `wout['keys']`
    #+***
    #+##Input Arguments
    #+    **wout**: VMEC dictionary
    #+
    #+##Keyword Arguments
    #+    **ntheta**: Number of theta points
    #+
    #+    **nphi**: Number of phi points
    #+
    #+    **thetamin**: Minimum theta value
    #+
    #+    **thetamax**: Maximum theta value
    #+
    #+    **phimin**: Minimum phi value
    #+
    #+    **phimax**: Maximum phi value
    #+
    #+##Output Arguments
    #+    **new_wout**: VMEC dictionary
    #+
    #+##Example Usage
    #+```python
    #+>>> wout = read_vmec(file_name)
    #+>>> new_wout = fourier_transform_3D(wout, ntheta=16, nphi=20, thetamin=0, thetamax=2*pi, phimin=0, phimax=2*pi)
    #+```
    """
    new_wout = wout.copy()
    if (thetamin is None and thetamax is None) or np.isclose([thetamin, thetamax], [0, 2*np.pi], atol=0.5*np.pi/180).all():
        thetamin = 0
        thetamax = 2*np.pi * (1 - 1/ntheta) # Removes double counting 0 (2*pi)
    if (phimin is None and phimax is None) or np.isclose([phimin, phimax], [0, 2*np.pi], atol=0.5*np.pi/180).all():
        phimin = 0
        phimax = 2*np.pi * (1 - 1/nphi)

    theta = np.linspace(thetamin, thetamax, ntheta)
    phi = np.linspace(phimin, phimax, nphi)

    pol, tor = np.meshgrid(theta, phi)
    pol_xm = np.dot(new_wout['xm'].reshape(new_wout['md'], 1), pol.reshape(1, nphi * ntheta))
    tor_xn = np.dot(new_wout['xn'].reshape(new_wout['md'], 1), tor.reshape(1, nphi * ntheta))

    cos_mu_nv = np.cos(pol_xm - tor_xn)
    sin_mu_nv = np.sin(pol_xm - tor_xn)

    if new_wout['nyq_limit']:
        for key in new_wout['keys']:
            if key in new_wout['cos_nyq_keys'] or key in new_wout['sin_nyq_keys']:
                pol_nyq_xm = np.dot(new_wout['xm_nyq'].reshape(new_wout['md_nyq'], 1), pol.reshape(1, nphi * ntheta))
                tor_nyq_xn = np.dot(new_wout['xn_nyq'].reshape(new_wout['md_nyq'], 1), pol.reshape(1, nphi * ntheta))

                cos_nyq_mu_nv = np.cos(pol_nyq_xm - tor_nyq_xn)
                sin_nyq_mu_nv = np.sin(pol_nyq_xm - tor_nyq_xn)

    inverse_fourier_amps = {}
    for key in new_wout['keys']:
        if key in new_wout['cos_keys']:
            inverse_fourier_amps[key] = np.dot(new_wout['fourier_amps'][key], cos_mu_nv).reshape(new_wout['ns'], nphi, ntheta)
        elif key in new_wout['sin_keys']:
            inverse_fourier_amps[key] = np.dot(new_wout['fourier_amps'][key], sin_mu_nv).reshape(new_wout['ns'], nphi, ntheta)
        elif new_wout['nyq_limit'] and key in new_wout['cos_nyq_keys']:
            inverse_fourier_amps[key] = np.dot(new_wout['fourier_amps'][key], cos_nyq_mu_nv).reshape(new_wout['ns'], nphi, ntheta)
        elif new_wout['nyq_limit'] and key in new_wout['sin_nyq_keys']:
            inverse_fourier_amps[key] = np.dot(new_wout['fourier_amps'][key], sin_nyq_mu_nv).reshape(new_wout['ns'], nphi, ntheta)
        else:
            raise NameError(f'key = "{key}" is not available')

    new_wout.update({
        'inverse_fourier_amps':inverse_fourier_amps,
        'ntheta':ntheta, 'thetamin':thetamin, 'thetamax':thetamax,
        'nphi':nphi, 'phimin':phimin, 'phimax':phimax,
        'theta':theta, 'phi':phi
        })

    return new_wout

def Brzp_transform(wout, cc_in_out=None, nrgrid=41, nzgrid=25):
    """
    #+#Brzp_transform
    #+ Converts B(s,v,u) to B(R,Z,tor)
    #+***
    #+##Input Arguments
    #+    **wout**: VMEC dictionary
    #+
    #+##Keyword Arguments
    #+    **cc_in_out**: Pair of COCOS indices for field coordinate conversion
    #+
    #+    **nrgrid**: Number of R points
    #+
    #+    **nzgrid**: Number of Z points
    #+
    #+##Output Arguments
    #+    **new_wout**: VMEC dictionary
    #+
    #+##Example Usage
    #+```python
    #+>>> wout = read_vmec(file_name)
    #+>>> new_wout = fourier_transform_3D(wout)
    #+>>> new_new_wout = Brzp_transform(new_wout, cc_in_out=(3,5), nrgrid=61, nzgrid=25)
    #+```
    """
    new_wout = wout.copy()
    inverse_fourier_amps = new_wout['inverse_fourier_amps']

    Bs = inverse_fourier_amps['Bs']
    Bv = inverse_fourier_amps['Bv']
    Bu = inverse_fourier_amps['Bu']

    R = inverse_fourier_amps['R']
    dR_ds = inverse_fourier_amps['dR_ds']
    dR_dv = inverse_fourier_amps['dR_dv']
    dR_du = inverse_fourier_amps['dR_du']

    dZ_ds = inverse_fourier_amps['dZ_ds']
    dZ_dv = inverse_fourier_amps['dZ_dv']
    dZ_du = inverse_fourier_amps['dZ_du']

    B_norm = np.zeros_like(R)
    denom = dR_ds * dZ_du - dR_du * dZ_ds
    nz = np.where(denom != 0)
    B_norm[nz] = 1. / denom[nz]

    # Calculates Br, Bz, Btor in (s,v,u) domain
    Br_svu = (dZ_du * Bs - dZ_ds * Bu) * B_norm
    Bz_svu = (dR_ds * Bu - dR_du * Bs) * B_norm
    Bt_svu = ( ( ( Bs * (dR_du * dZ_dv - dR_dv * dZ_du) + Bu * (dR_dv * dZ_ds - dR_ds * dZ_dv) ) * B_norm ) + Bv ) / R

    new_wout.update({
        'Br_svu':Br_svu,
        'Bz_svu':Bz_svu,
        'Bt_svu':Bt_svu
        })

    Z = inverse_fourier_amps['Z']
    ns = new_wout['ns']

    S = np.zeros([ns, new_wout['ntheta']])
    for i in range(new_wout['ntheta']):
        S[:, i] = new_wout['s_dom']

    idx = np.arange(0, ns-1)

    for key in ['Br', 'Bz', 'Bt']:
        new_wout[key] = np.zeros([nrgrid, nzgrid, new_wout['nphi']])
    rho_grid = np.zeros([nrgrid, nzgrid, new_wout['nphi']])

    rmin, rmax = R.min(), R.max()
    zmin, zmax = Z.min(), Z.max()
    rgrid, zgrid = np.mgrid[rmin:rmax:complex(nrgrid), zmin:zmax:complex(nzgrid)]
    rarr, zarr = rgrid[:, 0], zgrid[0, :]
    for iphi in range(new_wout['nphi']):
        rzdata = np.array([R[:, iphi, :].flatten(), Z[:, iphi, :].flatten()]).T
        rho_grid[:, :, iphi] = np.sqrt(griddata(rzdata, S.flatten(), (rgrid, zgrid), fill_value=0))

        idx = np.arange(1, ns, 1)
        rzdata = np.array([R[idx, iphi, :].flatten(), Z[idx, iphi, :].flatten()]).T
        for bkey in ['Br', 'Bz', 'Bt']:
            new_wout[key][:, :, iphi] = griddata(rzdata, new_wout[bkey][idx, iphi, :].flatten(), (rgrid, zgrid), fill_value=0)

    if cc_in_out:
        if isinstance(cc_in_out, (tuple, list, np.ndarray)) and len(cc_in_out) == 2:
            new_wout = convert_COCOS_VMEC(new_wout, cc_in_out)
        else:
            raise ValueError(f'cc_in_out must be a list-like object of length 2\n{cc_in_out=}')

    new_wout.update({
        'rho_grid':rho_grid,
        'R':rarr, 'Z':zarr,
        'nr':rarr.size, 'nz':zarr.size
        })

    return new_wout

def convert_COCOS_VMEC(wout, cc_in_out):
    """
    #+#convert_COCOS_VMEC
    #+ Converts wout COCOS cc_in_out[0] to COCOS cc_in_out[1]
    #+***
    #+##Input Arguments
    #+    **wout**: VMEC dictionary
    #+
    #+##Keyword Arguments
    #+    **cc_in_out**: Pair of COCOS indices for field coordinate conversion
    #+
    #+##Output Arguments
    #+    **new_wout**: VMEC dictionary
    #+
    #+##Example Usage
    #+```python
    #+>>> wout = read_vmec(file_name)
    #+>>> new_wout = convert_COCOS_VMEC(wout, cc_in_out=(3,5))
    #+```
    """
    new_wout = {key:val for key, val in wout.items() if key not in ['Br', 'Bz', 'Bt']}

    cc_in, cc_out = cc_in_out
    sigma_RphZ_out = -1 if cc_in % 2 == 0 else 1
    sigma_RphZ_in = -1 if cc_in % 2 == 0 else 1
    sigma_RphZ_eff = sigma_RphZ_out * sigma_RphZ_in

    new_wout['Br'] = wout['Br'] * sigma_RphZ_eff
    new_wout['Bz'] = wout['Bz'] * sigma_RphZ_eff
    new_wout['Bt'] = wout['Bt'] * sigma_RphZ_eff

    return new_wout

def fields_from_wout(wout, time=None, er=None, ez=None, et=None):
    """
    #+#fields_from_wout
    #+ Extracts field components from `wout` into FIDASIM-readable fields object
    #+***
    #+##Input Arguments
    #+    **wout**: VMEC dictionary
    #+
    #+##Keyword Arguments
    #+    **time**: Experiment time
    #+
    #+    **er**: E-field radial component
    #+
    #+    **ez**: E-field vertical component
    #+
    #+    **et**: E-field toroidal component
    #+
    #+##Output Arguments
    #+    **fields**: FIDASIM fields dictionary
    #+
    #+##Example Usage
    #+```python
    #+>>> wout = read_vmec(file_name)
    #+>>> new_wout = fourier_transform_3D(wout)
    #+>>> new_new_wout = Brzp_transform(new_wout, cc_in_out=(3,5), nrgrid=61, nzgrid=25)
    #+>>> fields = fields_from_wout(new_new_wout, time=0, er=None, ez=None, et=None)
    #+```
    """
    if time is None:
        time = 0
    if 'data_source' not in wout:
        data_source = ''
    else:
        data_source = wout['data_source']

    br, bz, bt = wout['br'], wout['bz'], wout['bt']
    btot =  np.sqrt(br**2 + bz**2 + bt**2)
    zeros = np.where(btot == 0)
    bmask = np.zeros(btot.shape)
    bmask[zeros] = 1

    if er is None:
        er = np.zeros(btot.shape)
    if ez is None:
        ez = np.zeros(btot.shape)
    if et is None:
        et = np.zeros(btot.shape)

    return {'time':time, 'mask':bmask,
            'data_source':data_source,
            'br':br, 'bz':bz, 'bt':bt,
            'er':er, 'ez':ez, 'et':et}

def grid_from_wout(wout):
    """
    #+#grid_from_wout
    #+ Extracts interpolation grid from `wout` into FIDASIM-readable igrid object
    #+***
    #+##Input Arguments
    #+    **wout**: VMEC dictionary
    #+
    #+##Output Arguments
    #+    **igrid**: FIDASIM igrid dictionary
    #+
    #+##Example Usage
    #+```python
    #+>>> wout = read_vmec(file_name)
    #+>>> new_wout = fourier_transform_3D(wout)
    #+>>> new_new_wout = Brzp_transform(new_wout, cc_in_out=(3,5), nrgrid=61, nzgrid=25)
    #+>>> igrid = grid_from_wout(new_new_wout)
    #+```
    """
    r, z, phi = wout['R'], wout['Z'], wout['phi']
    nr, nz, nphi = wout['nr'], wout['nz'], wout['nphi']
    z2d, r2d = np.mesgrid(z, r)
    return {'r':r, 'z':z, 'phi':phi,
            'nr':nr, 'nz':nz, 'nphi':nphi,
            'r2d':r2d, 'z2d':z2d}
