#!/bin/sh
"exec" "$FIDASIM_DIR/deps/python" "$0" "$@"
# -*- coding: utf-8 -*-

import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

from scipy.integrate import cumtrapz
from skimage.measure import find_contours

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
    flux = 2*cumtrapz(q,psi,initial=0.0)*dpsi #DLiu 10/13 to avoid negative flux
    
    pts = 101
    theta = np.linspace(0,2*np.pi,pts,endpoint=False)
    # Find r,theta of psi before boundary
    psi_i = g['ssimag'] + psi_eqdsk[-1]*dpsi # unnormalized
    R,Z = np.meshgrid(g['r'],g['z'])
    psi_v_candidates = find_contours(g['psirz'], psi_i)    #find contours that match the psi_i value
    psi_v_candidates = [
        p for p in psi_v_candidates 
        if (p[0,0] == p[-1,0])  # Only select closed contours
        and (
            # make sure the contour has both positive and negative z values (avoid closed contours in the diverter)
            np.any(np.interp(p[:, 0], range(0, len(g['z'])), g['z']) - g['zmaxis'] > 0) and 
            np.any(np.interp(p[:, 0], range(0, len(g['z'])), g['z']) - g['zmaxis'] < 0)
        )
    ]
    if not psi_v_candidates:
        raise ValueError(f"No suitable closed contour found at psi_i = {psi_i}")
    else:
        psi_v = psi_v_candidates[0]
        
    x_c = np.interp(psi_v[:,1], range(0,len(g['r'])), g['r']) - g['rmaxis']
    y_c = np.interp(psi_v[:,0], range(0,len(g['z'])), g['z']) - g['zmaxis']    
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

def rho_rz_orig(g,r_pts,z_pts,norm=True):

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

def unique_sort_index(sequence):
    unique = sorted(set(sequence))
    nindex=len(unique)
    index=np.zeros(nindex,dtype=np.int32)
    for idx in np.arange(nindex):
        index[idx]=np.where(sequence == unique[idx])[0][0]
    return index

def rho_rz(g,r_pts,z_pts,norm=True,psi_pts=None, do_linear=False):
    #Reference 4dlib/EFITLIB/kupfer/idl_geqdsk/rho_rz.pro
    #DLiu on 10/14/2023    
     
    if len(r_pts) != len(z_pts):
       print('RHO_RZ:  ERROR - # x_pts NE # y_pts ', len(r_pts), len(z_pts))
       
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

    if do_linear:
        rho_itp = scipy.interpolate.interp1d(psi,rho,'linear',fill_value=(np.min(rho),np.nan),bounds_error=False)    
    else:
        rho_itp = scipy.interpolate.interp1d(psi,rho,'cubic',fill_value=(np.min(rho),np.nan),bounds_error=False)

    rho_pts = rho_itp(psi_pts)
    w = np.isnan(rho_pts)
    rho_pts[w] = np.max(rho)
 
    #
    # -- now find all points that are outside the last closed flux surface 
    # -- and set rho_pts = max(rho) for all such points.  
    #
    # -- obtain boundary points from GEQDSK --
    #
    nbdry=np.int32(g['nbdry'])

    x_bdry = g['bdry'][0:nbdry-1,0]
    y_bdry = g['bdry'][0:nbdry-1,1]

    # -- express boundary points in (r,theta) coordinates --
    # Here (r,theta) denotes ordinary polar coordinates,
    # where the origin (r = 0) is at the magnetic axis.
    # Boundary points are (r,theta) = (r_bdry,theta_bdry).
    x_bdry = x_bdry - g['rmaxis']
    y_bdry = y_bdry - g['zmaxis']
    rr     = np.sqrt(x_bdry**2 + y_bdry**2)

    theta = rr * 0.  # Give theta right size = 0.

    i = np.where( rr != 0.0)[0]
    if i[0] != -1: 
        i_count=len(i)
    else:
        i_count=0
	   
    # Print('Rho_rz: i for where (_bdry.Ne. 0.0) = ', inum)
    if (i_count != 0):
        theta[i] = np.arccos(x_bdry[i]/rr[i])
	

    j = np.where(y_bdry < 0.)[0]
    if j[0] != -1: 
        in_y_bdry=len(j)
    else:
        in_y_bdry=0	
	
    if (in_y_bdry > 0):
        theta[j] = 2.*np.pi - theta[j]


    #
    # -- sort and make ends periodic --
    #
    k = unique_sort_index(theta)
    #Print, 'Rho_rz: N_elements(uniq(i)) = ', len(k)


    theta_uniq      = theta[k]
    rr_uniq         = rr[k]
    nrr             = len(rr_uniq)
    r_bdry            = np.zeros(nrr+2)
    theta_bdry        = np.zeros(nrr+2)
    theta_bdry[1:nrr+1] = theta_uniq
    r_bdry[1:nrr+1]     = rr_uniq
    r_bdry[0]         = rr_uniq[nrr-1]
    theta_bdry[0]     = theta_uniq[nrr-1] - 2.*np.pi
    r_bdry[nrr+1]     = rr_uniq[0]
    theta_bdry[nrr+1]   = theta_uniq[0] + 2.*np.pi


    #
    # -- express (x_pts,y_pts) in (r,theta) coordinates --
    #

    xx_pts = x_pts - g['rmaxis']
    yy_pts = y_pts - g['zmaxis']
    rr_pts = np.sqrt(xx_pts**2 + yy_pts**2)

    theta_pts = rr_pts * 0.   # give theta_pts the size of rr_pts and = 0.

    ii = np.where(rr_pts != 0.0)[0]
    if ii[0] != -1: 
        ii_count=len(ii)
    else:
        ii_count=0

    if (ii_count != 0):
        theta_pts[ii] = np.arccos(xx_pts[ii]/rr_pts[ii])
	
    jj = np.where(yy_pts < 0.)[0]
    if jj[0] != -1: 
        jj_count=len(jj)
    else:
        jj_count=0	
	
    if (jj_count > 0):
        theta_pts[jj] = 2.*np.pi - theta_pts[jj]	


    #
    # -- interpolate to get the radial position of the boundary 
    #    evaluated at theta = theta_pts --	
    r_boundary_itp = scipy.interpolate.interp1d(theta_bdry, r_bdry, 'cubic')
    r_boundary = r_boundary_itp(theta_pts)

    index = np.where(rr_pts > r_boundary)[0]
    if index[0] != -1: 
       index_count=len(i)
    else:
       index_count=0    

    #Print, 'Rho_rz: i_count for where rr_pts Gt r_bdry = ', i_count

    rhobnd = np.max(rho)
    
    if (index_count > 0):
        rho_pts[index] = rhobnd
	       
    if norm:
        rho_pts = rho_pts/np.max(rho)
        if index_count > 0 :
            rho_pts[index] = ((-psi_pts[index] + g['ssimag'])/(g['ssimag'] - g['ssibry']))**0.5

    return rho_pts.reshape(r_pts.shape)
