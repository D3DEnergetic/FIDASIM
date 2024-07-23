#!/bin/sh
"exe" "$FIDASIM_DIR/deps/python" "$0" "$@"
# -*- coding: utf-8 -*-

import numpy as np

def fourier_transform_3D(vmec, ntheta=5, nphi=10, thetamin=0, thetamax=2*np.pi, phimin=0, phimax=2*np.pi):
    theta = np.linspace(thetamin, thetamax, ntheta)
    phi = np.linspace(phimin, phimax, nphi)

    pol, tor = np.meshgrid(theta, phi)
    pol_xm = np.dot(vmec['xm'].reshape(vmec['md'], 1), pol.reshape(1, nphi * ntheta))
    tor_xn = np.dot(vmec['xn'].reshape(vmec['md'], 1), tor.reshape(1, nphi * ntheta))

    cos_mu_nv = np.cos(pol_xm - tor_xn)
    sin_mu_nv = np.sin(pol_xm - tor_xn)

    if vmec['nyq_limit']:
        for key in vmec['keys']:
            if key in vmec['cos_nyq_keys'] or key in vmec['sin_nyq_keys']:
                pol_nyq_xm = np.dot(vmec['xm_nyq'].reshape(vmec['md_nyq'], 1), pol.reshape(1, nphi * ntheta))
                tor_nyq_xn = np.dot(vmec['xn_nyq'].reshape(vmec['md_nyq'], 1), pol.reshape(1, nphi * ntheta))

                cos_nyq_mu_nv = np.cos(pol_nyq_xm - tor_nyq_xn)
                sin_nyq_my_nv = np.sin(pol_nyq_xm - tor_nyq_xn)

    invFourAmps = {}
    for key in vmec['keys']:
        if key in vmec['cos_keys']:
            invFourAmps[key] = np.dot(vmec['fourierAmps'][key], cos_mu_nv).reshape(vmec['ns'], nphi, ntheta)
        elif key in vmec['sin_keys']:
            invFourAmps[key] = np.dot(vmec['fourierAmps'][key], sin_mu_nv).reshape(vmec['ns'], nphi, ntheta)
        elif vmec['nyq_limit'] and key in vmec['cos_nyq_keys']:
            invFourAmps[key] = np.dot(vmec['fourierAmps'][key], cos_nyq_mu_nv).reshape(vmec['ns'], nphi, ntheta)
        elif vmec['nyq_limit'] and key in vmec['sin_nyq_keys']:
            invFourAmps[key] = np.dot(vmec['fourierAmps'][key], sin_nyq_mu_nv).reshape(vmec['ns'], nphi, ntheta)
        else:
            raise NameError(f'key = "{key}" is not available')

    return invFourAmps

def Brzp_transform(invFourAmps):
    Bs = invFourAmps['Bs']
    Bv = invFourAmps['Bv']
    Bu = invFourAmps['Bu']

    R= invFourAmps['R']
    dR_ds = invFourAmps['dR_ds']
    dR_dv = invFourAmps['dR_dv']
    dR_du = invFourAmps['dR_du']

    dZ_ds = invFourAmps['dZ_ds']
    dZ_dv = invFourAmps['dZ_dv']
    dZ_du = invFourAmps['dZ_du']

    B_norm = 1. / (dR_ds * dZ_du - dR_du * dZ_ds)
    Br = (dZ_du * Bs - dZ_ds * Bu) * B_norm
    Bp = ( ( ( Bs * (dR_du * dZ_dv - dR_dv * dZ_du) + Bu * (dR_dv * dZ_ds - dR_ds * dZ_dv) ) * B_norm ) + Bv ) / R
    Bz = (dR_ds * Bu - dR_du * Bs) * B_norm
    
    return {'Br':Br, 'Bz':Bz, 'Bp':Bp}
