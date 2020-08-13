#!/bin/sh
"exec" "$FIDASIM_DIR/deps/python" "$0" "$@"
# -*- coding: utf-8 -*-

import numpy as np
import re

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
