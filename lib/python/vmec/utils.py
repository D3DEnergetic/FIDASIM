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

def Brzp_transform(wout, cc_in_out=None):
    """
    """
    new_wout = wout.copy()
    if cc_in_out:
        if isinstance(cc_in_out, (tuple, list, np.ndarray)) and len(cc_in_out) == 2:
            new_wout = convert_COCOS_VMEC(new_wout, cc_in_out)
        else:
            raise ValueError(f'cc_in_out must be a list-like object of length 2\n{cc_in_out=}')
        
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

    Br = (dZ_du * Bs - dZ_ds * Bu) * B_norm
    Bz = (dR_ds * Bu - dR_du * Bs) * B_norm
    Bt = ( ( ( Bs * (dR_du * dZ_dv - dR_dv * dZ_du) + Bu * (dR_dv * dZ_ds - dR_ds * dZ_dv) ) * B_norm ) + Bv ) / R

    new_wout.update({
        'Br':Br,
        'Bz':Bz,
        'Bt':Bt
        })

    return new_wout

def map_B_and_rho(wout, nrgrid=41, nzgrid=25):
    """
    """
    new_wout = wout.copy()
    inverse_fourier_amps = new_wout['inverse_fourier_amps']

    R = inverse_fourier_amps['R']
    Z = inverse_fourier_amps['Z']
    ns = new_wout['ns']

    S = np.zeros([ns, new_wout['ntheta']])
    for i in range(new_wout['ntheta']):
        S[:, i] = new_wout['s_dom']

    idx = np.arange(0, ns-1)
    brzphi_grid = np.zeros([nrgrid, nzgrid, 3, new_wout['nphi']])
    rho_grid = np.zeros([nrgrid, nzgrid, new_wout['nphi']])

    rmin = R.min()
    rmax = R.max()
    zmin = Z.min()
    zmax = Z.max()
    rgrid, zgrid = np.mgrid[rmin:rmax:complex(nrgrid), zmin:zmax:complex(nzgrid)]
    rarr, zarr = rgrid[:, 0], zgrid[0, :]
    #rarr, zarr = np.zeros((nrgrid, new_wout['nphi'])), np.zeros((nzgrid, new_wout['nphi']))


    for iphi in range(new_wout['nphi']):
        """
        rmin = R[:, iphi, :].min()
        rmax = R[:, iphi, :].max()
        zmin = Z[:, iphi, :].min()
        zmax = Z[:, iphi, :].max()
        rgrid, zgrid = np.mgrid[rmin:rmax:complex(nrgrid), zmin:zmax:complex(nzgrid)]
        rarr[:, iphi] = rgrid[:, 0]
        zarr[:, iphi] = zgrid[0, :]
        """

        rzdata = np.array([R[:, iphi, :].flatten(), Z[:, iphi, :].flatten()]).T
        rho_grid[:, :, iphi] = np.sqrt(griddata(rzdata, S.flatten(), (rgrid, zgrid), fill_value=0))

        idx = np.arange(1, ns, 1)
        rzdata = np.array([R[idx, iphi, :].flatten(), Z[idx, iphi, :].flatten()]).T
        for ib, bkey in enumerate(['Br', 'Bz', 'Bt']):
            brzphi_grid[:, :, ib, iphi] = griddata(rzdata, new_wout[bkey][idx, iphi, :].flatten(), (rgrid, zgrid), fill_value=0)

    new_wout.update({
        'brzphi_grid':brzphi_grid,
        'rho_grid':rho_grid,
        'rarr':rarr, 'zarr':zarr
        })

    return new_wout

def convert_COCOS_VMEC(wout, cc_in_out):
    """
    """
    inverse_fourier_amps = wout['inverse_fourier_amps']
    new_inverse_fourier_amps = {key:val for key, val in inverse_fourier_amps.items() if key not in ['Bs', 'Bv', 'Bu']}
    new_wout = {key:val for key, val in wout.items() if key != 'inverse_fourier_amps'}

    cc_in, cc_out = cc_in_out
    sigma_RphZ_out = -1 if cc_in % 2 == 0 else 1
    sigma_RphZ_in = -1 if cc_in % 2 == 0 else 1
    sigma_RphZ_eff = sigma_RphZ_out * sigma_RphZ_in

    new_inverse_fourier_amps['Bs'] = inverse_fourier_amps['Bs'] * sigma_RphZ_eff # B_rho, radial-like component
    new_inverse_fourier_amps['Bv'] = inverse_fourier_amps['Bv'] * sigma_RphZ_eff # B_phi, toroidal-like component
    new_inverse_fourier_amps['Bu'] = inverse_fourier_amps['Bu'] * sigma_RphZ_eff # B_theta, poloidal-like component

    new_wout['inverse_fourier_amps'] = new_inverse_fourier_amps
    return new_wout
