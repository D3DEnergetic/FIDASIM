#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy.interpolate
import matplotlib._cntr as cntr

def fluxmap(g):
    npts = g['nw']
    dpsi = g['ssibry'] - g['ssimag']

    # Create phi array excluding boundary
    psi_eqdsk = np.linspace(0,1,npts-1,endpoint=False)
    q_eqdsk = g['qpsi'][0:(npts-1)]

    # Create a fine psi grid
    nint = 100
    psi = np.linspace(0,psi_eqdsk[-1],100)
    q = scipy.interpolate.interp1d(psi_eqdsk, q_eqdsk,'cubic')(psi)
    flux = 2*scipy.integrate.cumtrapz(q,psi,initial=0.0)*np.abs(dpsi)

    pts = 101
    theta = np.linspace(0,2*np.pi,pts,endpoint=False)
    # Find r,theta of psi before boundary
    psi_i = g['ssimag'] + psi_eqdsk[-1]*dpsi # unnormalized
    R,Z = np.meshgrid(g['r'],g['z'])
    psi_contour = cntr.Cntr(R,Z,g['psirz']).trace(psi_i,psi_i,0)[0]
    x_c = psi_contour[:,0] - g['rmaxis']
    y_c = psi_contour[:,1] - g['zmaxis']
    r_c = np.sqrt(x_c**2 + y_c**2)
    theta_c = np.arctan2(y_c,x_c)
    theta_c = np.where(theta_c < 0, theta_c + 2*np.pi, theta_c)
    sw = np.argsort(theta_c)
    theta_c = theta_c[sw]
    r_c = r_c[sw]
    theta_c, sw = np.unique(theta_c,return_index=True)
    r_c = r_c[sw]
    r_i = scipy.interpolate.interp1d(theta_c,r_c,"cubic",fill_value='extrapolate')(theta)

    # Find r, theta of boundary
    x_b = g['bdry'][:,0] - g['rmaxis']
    y_b = g['bdry'][:,1] - g['zmaxis']
    r_b = np.sqrt(x_b**2 + y_b**2)
    theta_b = np.arctan2(y_b,x_b)
    theta_b = np.where(theta_b < 0, theta_b + 2*np.pi, theta_b) #[0,2pi]
    sw = np.argsort(theta_b)
    theta_b = theta_b[sw]
    r_b = r_b[sw]
    theta_b, sw = np.unique(theta_b,return_index=True)
    r_b = r_b[sw]
    r_b = scipy.interpolate.interp1d(theta_b,r_b,"cubic",fill_value='extrapolate')(theta)

    #Integrate to find flux at boundary
    eps = r_i/g['rmaxis']
    y1 = (eps - np.log(1 + eps*np.cos(theta))/np.cos(theta))/np.cos(theta)
    eps = r_b/g['rmaxis']
    y2 = (eps - np.log(1 + eps*np.cos(theta))/np.cos(theta))/np.cos(theta)
    fpsi = 0.5*(g['fpol'][-2] + g['fpol'][-1])
    r_integral = np.abs(fpsi)*g['rmaxis']*(y2 - y1)
    dflux = 2*np.sum(r_integral)/pts
    flux_b = flux[-1] + dflux

    # Find rho and psi
    bcentr = np.abs(g['bcentr'])
    rho = np.linspace(0.0, np.sqrt(flux_b/bcentr), 101)

    flux_new = bcentr*rho**2
    psi = scipy.interpolate.interp1d(np.append(flux,flux_b),np.append(psi, 1.0),'cubic',fill_value='extrapolate')(flux_new)
    psi = g['ssimag'] + dpsi*psi

    return {'rho':rho,'psi':psi}

def rho_rz(g,r_pts,z_pts,norm=True):

    r = g['r']
    z = g['z']

    r_pts = np.array(r_pts)
    z_pts = np.array(z_pts)
    x_pts = r_pts.flatten()
    y_pts = z_pts.flatten()

    psirz = scipy.interpolate.interp2d(r,z,g['psirz'],'cubic')
    psi_pts = np.array([psirz(x,y) for (x,y) in zip(x_pts,y_pts)])

    a = fluxmap(g)
    rho = a['rho']
    psi = a['psi']

    rho_itp = scipy.interpolate.interp1d(psi,rho,'cubic',fill_value=(np.min(rho),np.nan),bounds_error=False)
    rho_pts = rho_itp(psi_pts)
    w = np.isnan(rho_pts)
    rho_pts[w] = np.max(rho)
    if norm:
        rho_pts = rho_pts/np.max(rho)
        rho_pts[w] = ((-psi_pts[w] + g['ssimag'])/(g['ssimag'] - g['ssibry']))**0.5

    return rho_pts.reshape(r_pts.shape)

