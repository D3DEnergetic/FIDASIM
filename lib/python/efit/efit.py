#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import re
import scipy.interpolate
import matplotlib.pyplot as plt
from matplotlib.pyplot import contour, clf

def file_numbers(fp):
    """Generator to get numbers from a text file"""
    toklist = []
    while True:
        line = fp.readline()
        if not line: break
        # Match numbers in the line using regular expression
        pattern = r'[+-]?\d*[\.]?\d+(?:[Ee][+-]?\d+)?'
        toklist = re.findall(pattern, line)
        for tok in toklist:
            yield tok

def readg(f):
    """ Reads a G-EQDSK file

    Parameters
    ----------

    f = Input file. Can either be a file-like object,
        or a string. If a string, then treated as a file name
        and opened.

    Returns
    -------

    """
    #TODO: Make this into a class and precalculate psirz spline
    if isinstance(f, str):
        # If the input is a string, treat as file name
        with open(f) as fh: # Ensure file is closed
            return readg(fh) # Call again with file object

    # Read the first line, which should contain the mesh sizes
    desc = f.readline()
    if not desc:
        raise IOError("Cannot read from input file")

    s = desc.split() # Split by whitespace
    if len(s) < 3:
        raise IOError("First line must contain at least 3 numbers")

    idum = int(s[-3])
    nw = int(s[-2])
    nh = int(s[-1])

    # Use a generator to read numbers
    token = file_numbers(f)

    rdim   = float(next(token))
    zdim   = float(next(token))
    rcentr = float(next(token))
    rleft  = float(next(token))
    zmid   = float(next(token))

    rmaxis = float(next(token))
    zmaxis = float(next(token))
    simag  = float(next(token))
    sibry  = float(next(token))
    bcentr = float(next(token))

    current= float(next(token))
    simag  = float(next(token))
    xdum   = float(next(token))
    rmaxis = float(next(token))
    xdum   = float(next(token))

    zmaxis = float(next(token))
    xdum   = float(next(token))
    sibry  = float(next(token))
    xdum   = float(next(token))
    xdum   = float(next(token))

    # Read arrays
    def read_array(n, name="Unknown"):
        data = np.zeros([n])
        try:
            for i in np.arange(n):
                data[i] = float(next(token))
        except:
            raise IOError("Failed reading array '"+name+"' of size ", n)
        return data

    # read 2d array
    def read_2d(nw, nh, name="Unknown"):
        data = np.zeros([nw, nh])
        for j in np.arange(nh):
            for i in np.arange(nw):
                data[i,j] = float(next(token))
        return data

    fpol   = read_array(nw, "fpol")
    pres   = read_array(nw, "pres")
    ffprim = read_array(nw, "ffprim")
    pprime = read_array(nw, "pprime")
    psirz  = read_2d(nw, nh, "psirz")
    qpsi   = read_array(nw, "qpsi")

    # Read boundary and limiters, if present
    nbbbs  = int(next(token))
    limitr = int(next(token))

    if nbbbs > 0:
        rbbbs = np.zeros([nbbbs])
        zbbbs = np.zeros([nbbbs])
        for i in range(nbbbs):
            rbbbs[i] = float(next(token))
            zbbbs[i] = float(next(token))
    else:
        rbbbs = [0]
        zbbbs = [0]

    if limitr > 0:
        rlim = np.zeros([limitr])
        zlim = np.zeros([limitr])
        for i in range(limitr):
            rlim[i] = float(next(token))
            zlim[i] = float(next(token))
    else:
        rlim = [0]
        zlim = [0]

    # Construct R-Z mesh
    r = np.linspace(rleft, rleft + rdim, nw)
    z = np.linspace(zmid - 0.5*zdim, zmid + 0.5*zdim, nh)

    # Create dictionary of values to return
    result = {'nw': nw, 'nh':nh,        # Number of horizontal and vertical points
              'r':r, 'z':z,                     # Location of the grid-poinst
              'rdim':rdim, 'zdim':zdim,         # Size of the domain in meters
              'rcentr':rcentr, 'bcentr':bcentr, # Reference vacuum toroidal field (m, T)
              'rleft':rleft,                  # R of left side of domain
              'zmid':zmid,                      # Z at the middle of the domain
              'rmaxis':rmaxis, 'zmaxis':zmaxis,     # Location of magnetic axis
              'ssimag':simag, # Poloidal flux at the axis (Weber / rad)
              'ssibry':sibry, # Poloidal flux at plasma boundary (Weber / rad)
              'current':current,
              'psirz':psirz.T,    # Poloidal flux in Weber/rad on grid points
              'fpol':fpol,  # Poloidal current function on uniform flux grid
              'ffprim':ffprim, # derivative of poloidal flux
              'pres':pres,  # Plasma pressure in nt/m^2 on uniform flux grid
              'pprime':pprime, # derivative of pressure
              'qpsi':qpsi,  # q values on uniform flux grid
              'nbdry':nbbbs, 'bdry':np.vstack((rbbbs,zbbbs)).T, # Plasma boundary
              'limitr':limitr, 'lim':np.vstack((rlim,zlim)).T} # Wall boundary

    return result

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
    psi_contour = contour(g['r'],g['z'],g['psirz'], levels=[psi_i]).collections[0].get_paths()[0].vertices
    clf()
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

