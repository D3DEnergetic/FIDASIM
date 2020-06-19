#!/usr/bin/env python
# -*- coding: utf-8 -*-

#+#FIDASIM Utilities
#+This file contains useful FIDASIM utilities
#+***
from __future__ import print_function
import os
from os.path import dirname
import subprocess
import platform
import numpy as np
import copy
import h5py
import efit
from scipy.io import netcdf
from scipy.interpolate import interp1d, interp2d, NearestNDInterpolator
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt

def get_fidasim_dir():
    """
    #+#get_fidasim_dir
    #+ Gets FIDASIM install directory
    #+***
    #+##Output Arguments
    #+     **directory**: FIDASIM install directory.
    #+##Example Usage
    #+```python
    #+>>> fida_dir = get_fidasim_dir()
    #+```
    """

    directory = dirname(dirname(dirname(dirname(os.path.abspath(__file__)))))

    return directory

def get_version(fidasim_dir):
    """
    #+#get_version
    #+ Gets FIDASIM version number from git.
    #+ Falls back to reading VERSION file when git is not available
    #+***
    #+##Input Arguments
    #+    **fidasim_dir**: FIDASIM install directory
    #+
    #+##Output Arguments
    #+     **version**: FIDAIM version number.
    #+
    #+##Example Usage
    #+```python
    #+>>> version = get_version(get_fidasim_dir())
    #+```
    """
    version = ''
    alt = False

    if platform.system() == 'Windows':
        alt = True
    else:
        # Location of .git folder
        git_dir = r'{}{}.git'.format(fidasim_dir, os.path.sep)

        # git is installed if git_file is a file
        proc = subprocess.Popen('command -v git', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        git_file = proc.communicate()[0].decode('utf-8')
        git_file = git_file.replace('\n', '')

        # Check that .git folder is present and git is installed
        if os.path.isfile(git_file) and os.path.isdir(git_dir):
            try:
                version = subprocess.check_output(['git', '--git-dir={}'.format(git_dir), 'describe', '--tags', '--always', '--dirty'])
                version = version.replace('\n', '')
            except:
                alt = True
        else:
            alt = True

    # If above didn't work, read version file
    if alt:
        # Git 'version' filepath
        ver_file = '{}{}VERSION'.format(fidasim_dir, os.path.sep)

        if os.path.isfile(ver_file):
            with open(ver_file) as f:
                version = f.read()

    return version

def aabb_intersect(rc, dr, r0, d0):
    """
    #+#aabb_intersect
    #+Calculates intersection length of a ray and an axis aligned bounding box (AABB)
    #+***
    #+##Input Arguments
    #+     **rc**: Center of AABB
    #+
    #+     **dr**: [length, width, height] of AABB
    #+
    #+     **r0**: starting point of ray
    #+
    #+     **d0**: direction of ray
    #+
    #+##Output Arguments
    #+     **intersect**: Intersection length of ray and AABB
    #+
    #+     **ri**: Optional, ray enterence point
    #+
    #+     **rf**: Optional, ray exit point
    #+
    #+##Example Usage
    #+```python
    #+>>> intersect, r_enter, r_exit = aabb_intersect([0,0,0], [1,1,1], [-1,0,0], [1,0,0])
    #+>>> print(intersect)
    #+    1.0
    #+>>> print(r_enter)
    #+    -0.5  0.0  0.0
    #+>>> print(r_exit)
    #+     0.5  0.0  0.0
    #+```
    """
    v0 = d0 / np.sqrt(np.sum(d0 ** 2.))

    # There are 6 sides to a cube/grid
    side_inter = np.zeros(6)

    # Intersection points of ray with planes defined by grid
    ipnts = np.zeros((3, 6))

    # Find whether ray intersects each side
    for i in range(6):
        j = int(np.floor(i / 2))
        ind = np.arange(3, dtype=int)
        ind = ind[ind != j]
        if np.abs(v0[j]) > 0.:   # just v0[j] != 0 right?
            # Intersection point with plane
            ipnts[:, i] = r0 + v0 * (((rc[j] + (np.mod(i, 2) - 0.5) * dr[j]) - r0[j]) / v0[j])

            # Check if point on plane is within grid side
            if (np.abs(ipnts[ind[0], i] - rc[ind[0]]) <= 0.5 * dr[ind[0]]) and \
               (np.abs(ipnts[ind[1], i] - rc[ind[1]]) <= 0.5 * dr[ind[1]]):
                side_inter[i] = 1

    intersect = 0.0
    r_enter = copy.deepcopy(r0)
    r_exit = copy.deepcopy(r0)
    ind = np.arange(side_inter.size)
    ind = ind[side_inter != 0]
    nw = side_inter[ind].size
    if nw >= 2:
        #Find two unique intersection points
        nunique = 0
        for i in range(nw - 1):
            if np.sum(ipnts[:, ind[0]] == ipnts[:, ind[i + 1]]) != 3:
                ind = [ind[0], ind[i + 1]]
                nunique = 2
                break

        if nunique == 2:
            vi = ipnts[:, ind[1]] - ipnts[:, ind[0]]
            vi = vi / np.sqrt(np.sum(vi ** 2.))
            dot_prod = np.sum(v0 * vi)
            if dot_prod > 0.0:
                r_enter = ipnts[:, ind[0]]
                r_exit = ipnts[:, ind[1]]
            else:
                r_enter = ipnts[:, ind[1]]
                r_exit = ipnts[:, ind[0]]

            # Calculate intersection length
            intersect = np.sqrt(np.sum((r_exit - r_enter) ** 2.))

    return intersect, r_enter, r_exit

def tb_zyx(alpha, beta, gamma):
    """
    #+#tb_zyx
    #+Calculates Tait-Bryan z-y'-x" active rotation matrix given rotation angles `alpha`,`beta`,`gamma` in radians
    #+***
    #+##Arguments
    #+     **alpha**: rotation angle about z [radians]
    #+
    #+     **beta**: rotation angle about y' [radians]
    #+
    #+     **gamma**: rotation angle about x" [radians]
    #+
    #+##Return Value
    #+     Rotation Matrix [prefida](|url|/sourcefile/prefida.pro.html)
    #+
    #+##Example Usage
    #+```python
    #+ >>> rot_mat = tb_zyx(np.pi/2, 0.0, np.pi/3)
    #+```
    """
    sa = np.sin(alpha)
    ca = np.cos(alpha)
    sb = np.sin(beta)
    cb = np.cos(beta)
    sg = np.sin(gamma)
    cg = np.cos(gamma)

    r = np.zeros((3, 3))

    r[0, 0] = ca * cb
    r[0, 1] = ca * sb * sg - cg * sa
    r[0, 2] = sa * sg + ca * cg * sb
    r[1, 0] = cb * sa
    r[1, 1] = ca * cg + sa * sb * sg
    r[1, 2] = cg * sa * sb - ca * sg
    r[2, 0] = -sb
    r[2, 1] = cb * sg
    r[2, 2] = cb * cg

    return r

def uvw_to_xyz(alpha, beta, gamma, uvw, origin=np.zeros(3)):
    """
    #+#uvw_to_xyz
    #+ Express non-rotated coordinate `uvw` in rotated `xyz` coordinates
    #+***
    #+##Arguments
    #+     **alpha**: Rotation angle about z [radians]
    #+
    #+     **beta**: Rotation angle about y' [radians]
    #+
    #+     **gamma**: Rotation angle about x" [radians]
    #+
    #+     **uvw**: Point in rotated coordinate system, (3, n)
    #+
    #+##Keyword Arguments
    #+     **origin**: Origin of rotated coordinate system in non-rotated (uvw) coordinates, (3)
    #+
    #+##Output Arguments
    #+     **xyz**: 'uvw' in 'xyz' coordinates
    #+
    #+##Example Usage
    #+```python
    #+>>> xyz = uvw_to_xyz(np.pi/2., 0.0, np.pi/3., uvw, origin=[.1, .2, 0.])
    #+```
    """

    # Make np arrays
    uvw = np.array(uvw, dtype=float)
    origin = np.array(origin, dtype=float)

    # Do checks as this code does not allow multiple points to be entered (yet)
    if uvw.ndim == 2:
        s = uvw.shape
        if s[0] != 3:
            raise ValueError('uvw must be (3, n), but it has shape {}'.format(uvw.shape))
        n = s[1]
    elif uvw.ndim == 1:
        if uvw.size != 3:
            raise ValueError('uvw must have length 3, but it has length {}'.format(uvw.size))
        n = 1
    else:
        raise ValueError('uvw must be (3) or (3, n)')

    if origin.ndim != 1:
        raise ValueError('origin must be 1D, but it has shape {}'.format(origin.shape))

    if origin.size != 3:
        raise ValueError('origin must have length 3, but it has length {}'.format(origin.size))

    # Shift origin
    uvw_shifted = uvw - np.squeeze(np.tile(origin, (n, 1)).T)

    # Get rotation matrix
    r = tb_zyx(alpha, beta, gamma)

    # Apply rotation matrix
    xyz = np.dot(r.T, uvw_shifted)

    return xyz

def xyz_to_uvw(alpha, beta, gamma, xyz, origin = np.zeros(3)):
    """
    #+##`xyz_to_uvw(alpha, beta, gamma, xyz, origin=[0,0,0])`
    #+Express rotated coordinate `xyz` in non-rotated `uvw` coordinates
    #+###Arguments
    #+     **alpha**: Rotation angle about z [radians]
    #+
    #+     **beta**: Rotation angle about y' [radians]
    #+
    #+     **gamma**: Rotation angle about x" [radians]
    #+
    #+     **xyz**: Point in rotated coordinate system
    #+
    #+###Keyword Arguments
    #+     **origin**: Origin of rotated coordinate system in non-rotated (uvw) coordinates.
    #+
    #+###Example Usage
    #+```python
    #+>>> uvw = xyz_to_uvw(np.pi/2,0.0,np.pi/3,xyz)
    #+```
    """
    xyz = np.array(xyz)

    # Do checks as this code does not allow multiple points to be entered (yet)
    if xyz.ndim == 2:
        s = xyz.shape
        if s[0] != 3:
            raise ValueError('xyz must be (3, n), but it has shape {}'.format(uvw.shape))
        n = s[1]
    elif xyz.ndim == 1:
        if xyz.size != 3:
            raise ValueError('xyz must have length 3, but it has length {}'.format(uvw.size))
        n = 1
    else:
        raise ValueError('xyz must be (3) or (3, n)')

    if origin.ndim != 1:
        raise ValueError('origin must be 1D, but it has shape {}'.format(origin.shape))

    if origin.size != 3:
        raise ValueError('origin must have length 3, but it has length {}'.format(origin.size))

    R = tb_zyx(alpha,beta,gamma)

    uvw = np.dot(R, xyz)

    return uvw + np.squeeze(np.tile(origin, (n, 1)).T)

def line_basis(r0, v0):
    """
    #+#line_basis
    #+Calculates basis from a line with +x in the direction of line
    #+***
    #+##Arguments
    #+    **r0**: Starting point of line [cm]
    #+
    #+    **v0**: Direction of line
    #+
    #+##Example Usage
    #+```python
    #+>>> basis = line_basis([0,0,0],[0,-1,0])
    #+>>> x = np.dot(basis,np.array([1,1,0])) ;Transforms a point in line-space ([1,1,0]) to real space
    #+>>> x
    #+    [1, -1, 0]
    #+```
    """
    r0 = np.array(r0)
    v0 = np.array(v0)
    rf = r0 + v0
    dis = np.sqrt(np.sum(v0**2))
    beta = np.arcsin((r0[2] - rf[2])/dis)
    alpha = np.arctan2((rf[1] - r0[1]),(rf[0]-r0[0]))

    R = tb_zyx(alpha,beta,0.0)
    return R

def rz_grid(rmin, rmax, nr, zmin, zmax, nz, phimin=0.0, phimax=0.0, nphi=1):
    """
    #+#rz_grid
    #+Creates interpolation grid
    #+***
    #+##Arguments
    #+    **rmin**: Minimum radius [cm]
    #+
    #+    **rmax**: Maximum radius [cm]
    #+
    #+    **nr**: Number of radii
    #+
    #+    **zmin**: Minimum Z value [cm]
    #+
    #+    **zmax**: Maximum Z value [cm]
    #+
    #+    **nz**: Number of Z values
    #+
    #+    **phimin**: Minimum Phi value [rad]
    #+
    #+    **phimax**: Maximum Phi value [rad]
    #+
    #+    **nphi**: Number of Phi values 
    #+
    #+##Return Value
    #+Interpolation grid dictionary
    #+
    #+##Example Usage
    #+```python
    #+>>> grid = rz_grid(0,200.0,200,-100,100,200,phimin=4*np.pi/3,phimax=5*np.pi/3,nphi=5)
    #+```
    """
    dr = (rmax - rmin) / nr
    dz = (zmax - zmin) / nz
    dphi = (phimax - phimin) / nphi
    r = rmin + dr * np.arange(nr, dtype=np.float64)
    z = zmin + dz * np.arange(nz, dtype=np.float64)
    phi = phimin + dphi * np.arange(nphi, dtype=np.float64)

    r2d = np.tile(r, (nz, 1)).T
    z2d = np.tile(z, (nr, 1))

    grid = {'r2d': r2d,
            'z2d': z2d,
            'r': r,
            'z': z,
            'phi': phi,
            'nr': nr,
            'nz': nz,
            'nphi': nphi}

    return grid

def colored(text, color): #, on_color=None, attrs=None):
    """
    #+#colored
    #+ Return text string formatting for color in terminal
    #+***
    #+##Input Arguments
    #+     **text**: String to be colored
    #+
    #+     **color**: Desired color of string. Red, green, yellow, blue, magenta, cyan, or white.
    #+
    #+##Output Arguments
    #+     **text**: Text formated to have "color" in terminal.
    #+##Example Usage
    #+```python
    #+>>> text = colored("Text to be red", 'red')
    #+>>> print(text)
    #+```
    """
    # Copyright (c) 2008-2011 Volvox Development Team
    #
    # Permission is hereby granted, free of charge, to any person obtaining a copy
    # of this software and associated documentation files (the "Software"), to deal
    # in the Software without restriction, including without limitation the rights
    # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    # copies of the Software, and to permit persons to whom the Software is
    # furnished to do so, subject to the following conditions:
    #
    # The above copyright notice and this permission notice shall be included in
    # all copies or substantial portions of the Software.
    #
    # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    # THE SOFTWARE.
    #
    # Author: Konstantin Lepa <konstantin.lepa@gmail.com>
    COLORS = dict(list(zip(['grey',
                            'red',
                            'green',
                            'yellow',
                            'blue',
                            'magenta',
                            'cyan',
                            'white',],
                            list(range(30, 38)))))

    RESET = '\033[0m'

    if os.getenv('ANSI_COLORS_DISABLED') is None:
        fmt_str = '\033[%dm%s'

        text = fmt_str % (COLORS[color], text)

        text += RESET

    return text

def info(string):
    """
    #+#info
    #+Print a informational message
    #+***
    #+##Arguments
    #+     **str**: message
    #+
    #+##Example Usage
    #+```python
    #+>>> info("This is an informative message")
    #+```
    """
    print(colored('INFO: ' + string, 'cyan'))

def warn(string):
    """
    #+#warn
    #+Print a warning message
    #+***
    #+##Arguments
    #+     **string**: message
    #+
    #+##Example Usage
    #+```python
    #+>>> warn("This may be a problem")
    #+```
    """
    print(colored('WARNING: ' + string, 'magenta'))

def error(string, halt=False):
    """
    #+#error
    #+Print a error message
    #+***
    #+##Arguments
    #+     **string**: message
    #+
    #+##Keyword Arguments
    #+     **halt**: Halt program execution
    #+
    #+##Example Usage
    #+```python
    #+>>> error("Error message")
    #+```
    """
    print(colored('ERROR: {}'.format(string), 'red'))

    if halt:
        raise Exception()

def success(string):
    """
    #+#success
    #+Print a success message
    #+***
    #+##Arguments
    #+     **string**: message
    #+
    #+##Example Usage
    #+```python
    #+>>> success("Yay!!!")
    #+```
    """
    print(colored('SUCCESS: ' + string, 'green'))

def beam_grid(nbi, rstart,
              nx=None, ny=None, nz=None, dv=8.0,
              length=100.0, width=80.0, height=80.0):
    """
    #+#beam_grid
    #+ Calculates settings for a grid that aligns with the neutral beam.
    #+***
    #+##Arguments
    #+    **nbi**: [Neutral beam geometry structure](|url|/page/03_technical/01_prefida_inputs.html#neutral-beam-geometry-structure)
    #+
    #+    **rstart**: Radial start position of beam grid [cm]
    #+
    #+##Keyword Arguments
    #+    **dV**: Cell volume [\(cm^3\)]: Defaults to 8.0
    #+
    #+    **nx**: Number of cells in length: Default determined by `dV`
    #+
    #+    **ny**: Number of cells in width: Default determined by `dV`
    #+
    #+    **nz**: Number of cells in height: Default determined by `dV`
    #+
    #+    **length**: Length of grid along beam sightline. [cm]: Defaults to 100 cm
    #+
    #+    **width**: Width of grid [cm]: Defaults to 100 cm
    #+
    #+    **height**: Height of grid [cm]: Defaults 80 cm
    #+
    #+##Return Value
    #+    Structure containing beam grid settings suitable for the Namelist File
    #+
    #+##Example Usage
    #+```python
    #+>>> grid = beam_grid(nbi,200.0,nx=100,ny=50,nz=50,length=100,width=50,height=50)
    #+```
    """

    if width < nbi['widy']:
        warn("Grid width is smaller then the source width")
        print("width: {}".format(width))
        print("source width: {}".format(nbi['widy']))

    if height < nbi['widz']:
        warn("Grid height is smaller then the source height")
        print("height: {}".format(height))
        print("source height: {}".format(nbi['widz']))

    dv3 = dv ** (1. / 3.)

    if nx is None:
        nx = round(length / dv3)

    if ny is None:
        ny = round(width / dv3)

    if nz is None:
        nz = round(height / dv3)

    xmin = 0.
    xmax = length
    ymin = -width / 2.
    ymax = width / 2.
    zmin = -height / 2.
    zmax = height / 2.

    src = nbi['src']
    axis = nbi['axis'] / np.sqrt(np.sum(nbi['axis'] ** 2))
    pos = src + 100. * axis

    if np.sqrt(src[0] ** 2 + src[1] ** 2) < rstart:
        error("Source radius cannot be less then rstart", halt=True)

    dis = np.sqrt(np.sum((src - pos) ** 2.0))
    beta = np.arcsin((src[2] - pos[2]) / dis)
    alpha = np.arctan2((pos[1] - src[1]), (pos[0] - src[0]))
    gamma = 0.
    a = axis[0] ** 2 + axis[1] ** 2
    b = 2. * (src[0] * axis[0] + src[1] * axis[1])
    c = src[0] ** 2 + src[1] ** 2 - rstart ** 2
    t = (-b - np.sqrt(b ** 2 - 4. * a * c)) / (2. * a)
    origin = src + t * axis

    beam_grid = {'nx': nx,
                 'ny':  ny,
                 'nz':  nz,
                 'xmin': xmin,
                 'xmax': xmax,
                 'ymin': ymin,
                 'ymax': ymax,
                 'zmin': zmin,
                 'zmax': zmax,
                 'alpha': alpha,
                 'beta': beta,
                 'gamma': gamma,
                 'origin': origin}

    return beam_grid

def write_data(h5_obj, dic, desc=dict(), units=dict(), name=''):
    """
    #+#write_data
    #+ Write h5 datasets with attributes 'description' and 'units'
    #+***
    #+##Arguments
    #+     **h5_obj**: An h5 file or group object from h5py
    #+
    #+     **dic**: Dict of data to save as h5 datasets
    #+
    #+##Keyword Arguments
    #+     **name**: Name/description of dic for clarity in raising errors
    #+
    #+     **desc**: Dict with same keys as dic describing each item in dic
    #+
    #+     **units**: Dict with same keys as dic providing units of data in dic, doesn't have to be all keys of dic.
    #+
    #+##Example Usage
    #+```python
    #+>>> write_data(h5_obj, dic, desc, units)
    #+```
    """
    for key in dic:
        if isinstance(dic[key], dict):
            h5_grp = h5_obj.create_group(key)
            write_data(h5_grp, dic[key])
            continue

        # Transpose data to match expected by Fortran and historically provided by IDL
        if isinstance(dic[key], np.ndarray):
            if dic[key].ndim >= 2:
                dic[key] = dic[key].T

        # Make strings of fixed length as required by Fortran.
        # See http://docs.h5py.org/en/latest/strings.html#fixed-length-ascii
        if isinstance(dic[key], str):
            dic[key] = np.string_(dic[key])

        # Create dataset
        ds = h5_obj.create_dataset(key, data = dic[key])

        # Add descrption attribute
        if key in desc:
            ds.attrs['description'] = desc[key]

        # Add units attribute (if present)
        if key in units:
            ds.attrs['units'] = units[key]

def read_geqdsk(filename, grid, poloidal=False):
    """
    #+#read_geqdsk
    #+Reads an EFIT GEQDSK file
    #+***
    #+##Arguments
    #+    **filename**: GEQDSK file
    #+
    #+    **grid**: Interpolation grid
    #+
    #+##Keyword Arguments
    #+    **poloidal**: Return rho_p (sqrt(normalized poloidal flux)) instead of rho (sqrt(normalized toroidal flux))
    #+
    #+##Return Value
    #+Electronmagnetic fields structure, rho, btipsign
    #+
    #+##Example Usage
    #+```python
    #+>>> fields, rho, btipsign = read_geqdsk("./g133223.00200",grid)
    #+```
    """
    dims = grid['r2d'].shape
    r_pts = grid['r2d'].flatten()/100
    z_pts = grid['z2d'].flatten()/100
    g = efit.readg(filename)
    btipsign = np.sign(g["current"]*g["bcentr"])

    fpol = g["fpol"]
    psiaxis = g["ssimag"]
    psiwall = g["ssibry"]
    r = g["r"]
    z = g["z"]

    psi_arr = np.linspace(psiaxis, psiwall, len(fpol))
    fpol_itp = interp1d(psi_arr, fpol, 'cubic', fill_value=fpol[-1],bounds_error=False)
    psirz_itp = interp2d(r, z, g["psirz"], 'cubic')

    if poloidal:
        rhogrid = np.array([psirz_itp(rr,zz) for (rr,zz) in zip(r_pts,z_pts)]).reshape(dims)
        rhogrid = np.sqrt((rhogrid - g["ssimag"])/(g["ssibry"] - g["ssimag"]))
    else:
        rhogrid=efit.rho_rz(g,r_pts,z_pts,norm=True).reshape(dims)

    br = np.array([psirz_itp(rr,zz,dy=1)/rr for (rr,zz) in zip(r_pts,z_pts)]).reshape(dims)
    bz = np.array([-psirz_itp(rr,zz,dx=1)/rr for (rr,zz) in zip(r_pts,z_pts)]).reshape(dims)
    bt = np.array([fpol_itp(psirz_itp(rr,zz))/rr for (rr,zz) in zip(r_pts,z_pts)]).reshape(dims)

    er = br*0
    ez = bz*0
    et = bt*0

    mask = np.ones(dims,dtype=np.int32)

    equil = {"time":0.0,"data_source":os.path.abspath(filename), "mask":mask,
             "br":br,"bt":bt,"bz":bz,"er":er,"et":et,"ez":ez}

    return equil, rhogrid, btipsign

def read_ncdf(filename, vars=None):
    '''
    #+#read_ncdf
    #+Reads a flat NetCDF file
    #+***
    #+##Arguments
    #+    **filename**: NetCDF file
    #+
    #+##Keyword Arguments
    #+    **vars**: List of variables to read
    #+
    #+##Return Value
    #+Structure containing NetCDF variables
    #+
    #+##Example Usage
    #+```python
    #+>>> a = read_ncdf("./123324H01_fi_1.cdf")
    #+```
    '''

    d = dict()
    d['err'] = 1
    if os.path.isfile(filename):
        d['err'] = 0
        f = netcdf.netcdf_file(filename, 'r', mmap=False)
        variables = f.variables
        if vars != None:
            for k in vars:
                # need to check case sensitibity
                if k in variables.keys():
                    v = variables[k]
                    if tuple() == v.shape:
                        d[k] = v.getValue()
                    else:
                        d[k] = v[:]
        else:
            for k,v in variables.items():
                if tuple() == v.shape:
                    d[k] = v.getValue()
                else:
                    d[k] = v[:]
        f.close()
    else:
        error('FILE DOES NOT EXIST: '+filename)

    return d

def extract_transp_plasma(filename, intime, grid, rhogrid,
                          dn0out=None, scrapeoff=None,rho_scrapeoff=0.1):
    '''
    #+#extract_transp_plasma
    #+Extracts `plasma` structure from a TRANSP run
    #+***
    #+##Arguments
    #+    **filename**: TRANSP output file e.g. [TRANSP_RUNID].CDF
    #+
    #+    **intime**: Time of interest [s]
    #+
    #+    **grid**: Interpolation grid
    #+
    #+    **rhogrid**: sqrt(normalized torodial flux) mapped onto the interpolation grid
    #+
    #+##Keyword Arguments
    #+    **dn0out**: Wall Neutral density value `dn0out` variable in transp namelist
    #+
    #+    **scrapeoff**: scrapeoff decay length
    #+
    #+    **rho_scrapeoff**: scrapeoff length, default = 0.1
    #+
    #+##Example Usage
    #+```python
    #+>>> plasma = extract_transp_plasma("./142332H01.CDF", 1.2, grid, rho)
    #+```
    '''

    var_list = ["X","TRFLX","TFLUX","TIME","NE","NH","ND","NT","NIMP","TE","TI","ZEFFI","OMEGA","DN0WD","XZIMP"]

    zz = read_ncdf(filename, vars=var_list)

    t = zz['TIME']
    idx = np.argmin(abs(t-intime))
    time = t[idx].astype('float64')

    print(' * Selecting profiles at :', time, ' s') #pick the closest timeslice to TOI

    impurity_charge = np.max(zz["XZIMP"]).astype("int16")
    transp_ne = zz['NE'][idx,:] #cm^-3
    transp_nimp = zz['NIMP'][idx,:] #cm^-3
    transp_nn = zz['DN0WD'][idx,:] #cm^-3

    if 'NH' in zz:
        transp_nh = zz['NH'][idx,:] #cm^-3
    else:
        transp_nh = 0*transp_ne

    if 'ND' in zz:
        transp_nd = zz['ND'][idx,:] #cm^-3
    else:
        transp_nd = 0*transp_ne

    if 'NT' in zz:
        transp_nt = zz['NT'][idx,:] #cm^-3
    else:
        transp_nt = 0*transp_ne

    transp_te = zz['TE'][idx,:]*1.e-3 # kev
    transp_ti = zz['TI'][idx,:]*1.e-3 # kev
    transp_zeff = zz['ZEFFI'][idx,:]
    rho_cb = np.sqrt(zz['TRFLX'][idx,:]/zz['TFLUX'][idx])
    # center each rho b/c toroidal flux is at cell boundary
    rho = 0.e0*rho_cb
    rho[0] = 0.5*rho_cb[0]
    for i in range(len(rho_cb)-1):
        rho[i+1] = rho_cb[i+1] - 0.5*(rho_cb[i+1] - rho_cb[i])

    if 'OMEGA' not in zz.keys():
        error('OMEGA not found in TRANSP file. Assuming no plasma rotation')
        transp_omega=0.0*transp_te
    else:
        transp_omega = zz['OMEGA'][idx,:] # rad/s

    if dn0out == None:
        dn0out = transp_nn[-1]
    if scrapeoff == None:
        scrapeoff = 0.0

    if scrapeoff > 0.0:
        drho = abs(rho[-1] - rho[-2])
        rho_sc = rho[-1] + drho*(range(np.ceil(rho_scrapeoff/drho)) + 1)
        sc = np.exp(-(rho_sc - rho[-1])/scrapeoff)
        transp_ne = np.append(transp_ne,transp_ne[-1]*sc)
        transp_nimp = np.append(transp_nimp,transp_nimp[-1]*sc)
        transp_nh = np.append(transp_nh,transp_nh[-1]*sc)
        transp_nd = np.append(transp_nd,transp_nd[-1]*sc)
        transp_nt = np.append(transp_nt,transp_nt[-1]*sc)
        transp_te = np.append(transp_te,transp_te[-1]*sc)
        transp_ti = np.append(transp_ti,transp_ti[-1]*sc)
        transp_nn = np.append(transp_nn,0*sc + dn0out)
        transp_zeff = np.append(transp_zeff, (transp_zeff[-1]-1)*sc + 1)
        transp_omega = np.append(transp_omega,transp_omega[-1]*sc)
        rho = np.append(rho, rho_sc)

    profiles = {"rho":rho,
                "dene":np.where(transp_ne > 0, transp_ne, 0.0),
                "denimp":np.where(transp_nimp > 0, transp_nimp, 0.0),
		"denn":np.where(transp_nn > 0, transp_nn, 0.0),
                "te":np.where(transp_te > 0, transp_te, 0.0),
                "ti":np.where(transp_ti > 0, transp_ti, 0.0),
                "zeff":np.where(transp_zeff > 1.0, transp_zeff, 1.0),
                "omega":transp_omega}
    if 'NH' in zz:
        profiles['denh'] = np.where(transp_nh > 0, transp_nh, 0.0)
    if 'ND' in zz:
        profiles['dend'] = np.where(transp_nd > 0, transp_nd, 0.0)
    if 'NT' in zz:
        profiles['dent'] = np.where(transp_nt > 0, transp_nt, 0.0)


    # Interpolate onto r-z grid
    dims = rhogrid.shape
    f_dene = interp1d(rho,transp_ne,fill_value='extrapolate')
    dene = f_dene(rhogrid)
    dene = np.where(dene > 0.0, dene, 0.0).astype('float64')

    f_denimp = interp1d(rho,transp_nimp,fill_value='extrapolate')
    denimp = f_denimp(rhogrid)
    denimp = np.where(denimp > 0.0, denimp, 0.0).astype('float64')

    f_denh = interp1d(rho,transp_nh,fill_value='extrapolate')
    denh = f_denh(rhogrid)
    denh = np.where(denh > 0.0, denh, 0.0).astype('float64')

    f_dend = interp1d(rho,transp_nd,fill_value='extrapolate')
    dend = f_dend(rhogrid)
    dend = np.where(dend > 0.0, dend, 0.0).astype('float64')

    f_dent = interp1d(rho,transp_nt,fill_value='extrapolate')
    dent = f_dent(rhogrid)
    dent = np.where(dent > 0.0, dent, 0.0).astype('float64')

    f_denn = interp1d(rho,np.log(transp_nn),fill_value=np.nan,bounds_error=False)
    log_denn = f_denn(rhogrid)
    denn = np.where(~np.isnan(log_denn), np.exp(log_denn), 0.0).astype('float64')

    f_te = interp1d(rho,transp_te,fill_value='extrapolate')
    te = f_te(rhogrid)
    te = np.where(te > 0, te, 0.0).astype('float64')

    f_ti = interp1d(rho,transp_ti,fill_value='extrapolate')
    ti = f_ti(rhogrid)
    ti = np.where(ti > 0, ti, 0.0).astype('float64')

    f_zeff = interp1d(rho,transp_zeff, fill_value=1.0, bounds_error=False)
    zeff = f_zeff(rhogrid)
    zeff = np.where(zeff > 1, zeff, 1.0).astype('float64')

    f_omega = interp1d(rho,transp_omega,fill_value='extrapolate')
    vt = grid['r2d']*f_omega(rhogrid).astype('float64')
    vr = np.zeros(dims,dtype='float64')
    vz = np.zeros(dims,dtype='float64')

    max_rho = max(abs(rho))

    mask = np.zeros(dims,dtype='int')
    w = np.where(rhogrid <= max_rho) #where we have profiles
    mask[w] = 1

    deni = np.concatenate((denh.reshape(1,dims[0],dims[1]),
                           dend.reshape(1,dims[0],dims[1]),
                           dent.reshape(1,dims[0],dims[1])),axis=0)

    ai = np.array([1.007276466879e0, 2.013553212745e0,3.01550071632e0])
    w_ai = [a in zz for a in ['NH','ND','NT']]

    # SAVE IN PROFILES STRUCTURE
    plasma={"data_source":os.path.abspath(filename),"time":time,"impurity_charge":impurity_charge,
            "nthermal":np.sum(w_ai).astype('int16'), "species_mass":ai[w_ai], "deni":deni[w_ai,:,:],"profiles":profiles,
            "mask":mask,"dene":dene,"denimp":denimp,"denn":denn,"te":te,"ti":ti,
            "vr":vr,"vt":vt,"vz":vz,"zeff":zeff}

    return plasma

def read_nubeam(filename, grid, e_range=(), p_range=(), btipsign=-1, species=1):
    """
    #+#read_nubeam
    #+Reads NUBEAM fast-ion distribution function
    #+***
    #+##Arguments
    #+    **filename**: NUBEAM guiding center fast-ion distribution function file e.g. 159245H01_fi_1.cdf
    #+
    #+    **grid**: Interpolation grid
    #+
    #+##Keyword Arguments
    #+    **btipsign**: Sign of the dot product of the magnetic field and plasma current
    #+
    #+    **e_range**: Energy range to consider
    #+
    #+    **p_range**: Pitch range to consider
    #+
    #+    **species**: Fast-ion species number. Defaults to 1
    #+
    #+##Return Value
    #+Distribution structure
    #+
    #+##Example Usage
    #+```python
    #+>>> dist = read_nubeam("./159245H02_fi_1.cdf",grid,btipsign=-1)
    #+```
    """

    species_var = "SPECIES_{}".format(species)
    sstr = read_ncdf(filename,vars=[species_var])[species_var].tostring().decode('UTF-8')
    print("Species: "+sstr)
    var = read_ncdf(filename, vars=["TIME","R2D","Z2D","E_"+sstr,"A_"+sstr,"F_"+sstr,"RSURF","ZSURF","BMVOL"])

    ngrid = len(var["R2D"])

    try:
        time = var["TIME"][0]
    except:
        time = var["TIME"]

    r2d = var["R2D"]
    z2d = var["Z2D"]
    rsurf = var["RSURF"].T
    zsurf = var["ZSURF"].T
    bmvol = var["BMVOL"]
    pitch = var["A_"+sstr]
    energy = var["E_"+sstr]*1e-3
    fbm = var["F_"+sstr].T*1e3
    fbm = np.where(fbm > 0.0, 0.5*fbm, 0.0) #0.5 to convert to pitch instead of solid angle d_omega/4pi

    if btipsign < 0:
        fbm = fbm[:,::-1,:] #reverse pitch elements

    if not e_range:
        e_range = (np.min(energy), np.max(energy))

    if not p_range:
        p_range = (np.min(pitch), np.max(pitch))

    # Trim distribution according to e/p_range
    we = np.logical_and(energy >= e_range[0], energy <= e_range[1])
    wp = np.logical_and(pitch >= p_range[0], pitch <= p_range[1])
    energy = energy[we]
    nenergy = len(energy)
    pitch = pitch[wp]
    npitch = len(pitch)
    fbm = fbm[we,:,:]
    fbm = fbm[:,wp,:]
    dE = np.abs(energy[1] - energy[0])
    dp = np.abs(pitch[1] - pitch[0])
    emin, emax = np.maximum(np.min(energy) - 0.5*dE, 0.0), np.max(energy) + 0.5*dE
    pmin, pmax = np.maximum(np.min(pitch) - 0.5*dp, -1.0), np.minimum(np.max(pitch)+0.5*dp, 1.0)

    print('Energy min/max: ', emin, emax)
    print('Pitch min/max: ',pmin, pmax)

    nr = grid["nr"]
    nz = grid["nz"]
    r = grid["r"]
    z = grid["z"]
    rgrid = grid["r2d"]
    zgrid = grid["z2d"]
    dr = np.abs(r[1] - r[0])
    dz = np.abs(z[1] - z[0])

    fdens = np.sum(fbm,axis=(0,1))*dE*dp
    ntot = np.sum(fdens*bmvol)
    print('Ntotal in phase space: ',ntot)

    tri = Delaunay(np.vstack((r2d,z2d)).T) # Triangulation for barycentric interpolation
    pts = np.array([xx for xx in zip(r2d,z2d)])
    itp = NearestNDInterpolator(pts,np.arange(ngrid)) #to find indices outside simplices

    points = np.array([xx for xx in zip(rgrid.flatten(),zgrid.flatten())])
    t = tri.find_simplex(points)

    denf = np.zeros((nr,nz))
    fbm_grid = np.zeros((nenergy,npitch,nr,nz))
    for (ind,tt) in enumerate(t):
        i,j = np.unravel_index(ind,(nr,nz))
        if tt == -1:
            ii = int(itp(r[i],z[j]))
            denf[i,j] = fdens[ii]
            fbm_grid[:,:,i,j] = fbm[:,:,ii]
        else:
            b = tri.transform[tt,:2].dot(np.transpose(points[ind] - tri.transform[tt,2]))
            s = tri.simplices[tt,:]
            #perform barycentric linear interpolation
            denf[i,j] = b[0]*fdens[s[0]] + b[1]*fdens[s[1]] + (1 - np.sum(b))*fdens[s[2]]
            fbm_grid[:,:,i,j] = b[0]*fbm[:,:,s[0]] + b[1]*fbm[:,:,s[1]] + (1-np.sum(b))*fbm[:,:,s[2]]

    denf[denf < 0] = 0

    # Correct for points outside of seperatrix
    rmaxis = np.mean(rsurf[:,0])
    zmaxis = np.mean(zsurf[:,0])
    r_sep = rsurf[:,-1]
    z_sep = zsurf[:,-1]

    #plt.triplot(r2d,z2d,tri.simplices.copy())
    #plt.plot(r2d,z2d,'o')
    #plt.plot(r_sep,z_sep)
    #plt.show()
    x_bdry = r_sep - rmaxis
    y_bdry = z_sep - zmaxis
    r_bdry = np.sqrt(x_bdry**2 + y_bdry**2)
    theta_bdry = np.arctan2(y_bdry,x_bdry)
    theta_bdry = np.where(theta_bdry < 0.0, theta_bdry + 2*np.pi, theta_bdry) #[0,2pi]
    w = np.argsort(theta_bdry)
    theta_bdry = theta_bdry[w]
    r_bdry = r_bdry[w]
    theta_bdry, w = np.unique(theta_bdry,return_index=True)
    r_bdry = r_bdry[w]
    itp = interp1d(theta_bdry,r_bdry,'cubic',fill_value='extrapolate')

    x_pts = grid["r2d"] - rmaxis
    y_pts = grid["z2d"] - zmaxis
    r_pts = np.sqrt(x_pts**2 + y_pts**2)
    theta_pts = np.arctan2(y_pts,x_pts)
    theta_pts = np.where(theta_pts < 0.0, theta_pts + 2*np.pi, theta_pts) #[0,2pi]
    r_bdry_itp = itp(theta_pts)

    w = r_pts >= r_bdry_itp + 2
    denf[w] = 0.0
    fbm_grid[:,:,w] = 0.0

    # enforce correct normalization
    ntot_denf = 2*np.pi*dr*dz*np.sum(r*np.sum(denf,axis=1))
    denf = denf*(ntot/ntot_denf)
    ntot_fbm = (2*np.pi*dE*dp*dr*dz)*np.sum(r*np.sum(fbm_grid,axis=(0,1,3)))
    fbm_grid = fbm_grid*(ntot/ntot_denf)


    fbm_dict={"type":1,"time":time,"nenergy":nenergy,"energy":energy,"npitch":npitch,
              "pitch":pitch,"f":fbm_grid,"denf":denf,"data_source":os.path.abspath(filename)}

    return fbm_dict

def nubeam_geometry(nubeam, angle=0.0, verbose=False):
    """
    #+#nubeam_geometry
    #+Calculates the FIDASIM beam geometry from the beam geometry variables in the TRANSP/NUBEAM namelist
    #+***
    #+##Arguments
    #+     **NUBEAM**: Dictionary containing the following
    #+
    #+     **NUBEAM["NAME"]**: Ion source name
    #+
    #+     **NUBEAM["NBSHAP"]**: Ion source shape 1=rectangular, 2=circular
    #+
    #+     **NUBEAM["FOCLZ"]**: Vertical focal length [cm]
    #+
    #+     **NUBEAM["FOCLR"]**: Horizontal focal length [cm]
    #+
    #+     **NUBEAM["DIVZ"]**: Vertical divergence [rad]
    #+
    #+     **NUBEAM["DIVR"]**: Horizontal divergence [rad]
    #+
    #+     **NUBEAM["BMWIDZ"]**: Ion source half height [cm]
    #+
    #+     **NUBEAM["BMWIDR"]**: Ion source half width [cm]
    #+
    #+     **NUBEAM["RTCENA"]**: Radius of tangency point [cm]
    #+
    #+     **NUBEAM["XLBTNA"]**: Distance from center of beam source grid to tangency point [cm]
    #+
    #+     **NUBEAM["XBZETA"]**: Torodial angle [deg] Positive angles defined to be in the counter-clockwise direction
    #+
    #+     **NUBEAM["XYBSCA"]**: Elevation above/below vacuum vessel midplane of center of beam source grid [cm]
    #+
    #+     **NUBEAM["NLJCCW"]**: Orientation of Ip. 1 for True/Counter-clockwise current, 0 or -1 for False/Clock-wise current
    #+
    #+     **NUBEAM["NLCO"]**: 1 for Co-beam, 0 or -1 for Counter-beam
    #+
    #+     **NUBEAM["NBAPSHA"]**: Vector of aperture shapes 1=rectangular, 2=circular
    #+
    #+     **NUBEAM["XLBAPA"]**: Vector of distances from center of beam source grid to the aperture plane [cm]
    #+
    #+     **NUBEAM["XYBAPA"]**: Vector of elevation above/below vacuum vessel midplane of beam centerline at aperture [cm]
    #+
    #+     **NUBEAM["RAPEDGA"]**: Vector of aperture half-widths [cm]
    #+
    #+     **NUBEAM["XZPEDGA"]**: Vector of aperture half-heights [cm]
    #+
    #+     **NUBEAM["XRAPOFFA"]**: Vector of horizontal (y) offsets relative to the +x aligned beam centerline [cm]
    #+
    #+     **NUBEAM["XZAPOFFA"]**: Vector of vertical (z) offsets relative to the +x aligned beam centerline [cm]
    #+
    #+##Keyword Arguments
    #+     **angle**: Angle to add to XBZETA to rotate the beams into correct coordinates [deg]
    #+
    #+     **verbose**: Print out positions
    #+
    #+##Return Value
    #+ Neutral beam structure
    #+
    #+##Example Usage
    #+```python
    #+>>> nbi = nubeam_geometry(nubeam)
    #+```
    """

    if nubeam["NLCO"] == 0:
        nubeam["NLCO"] = -1

    if "NLJCCW" in nubeam:
        if nubeam["NLJCCW"] == 0:
            nubeam["NLJCCW"] = -1
    else:
        warn("Current orientation not specified. Assuming Counter-clockwise.")
        nubeam["NLJCCW"] = 1

    phi_s = (nubeam["XBZETA"] + angle)*np.pi/180.0
    zs = nubeam["XYBSCA"]
    za = nubeam["XYBAPA"][0]
    alpha = np.arcsin((zs-za)/nubeam["XLBAPA"][0])
    pdst = nubeam["XLBTNA"]*np.cos(alpha)
    rs = np.sqrt(nubeam["RTCENA"]**2 + pdst**2)
    dat = nubeam["XLBTNA"] - nubeam["XLBAPA"][0]
    pdat = dat*np.cos(alpha)
    ra = np.sqrt(nubeam["RTCENA"]**2 + pdat**2.0)
    beta_s = np.arccos(nubeam["RTCENA"]/rs)
    beta_a = np.arccos(nubeam["RTCENA"]/ra)
    phi_a = phi_s + nubeam["NLJCCW"]*nubeam["NLCO"]*(beta_s-beta_a)

    src = np.array([rs*np.cos(phi_s), rs*np.sin(phi_s),zs])
    aper_src = np.array([ra*np.cos(phi_a), ra*np.sin(phi_a),za])
    axis = (aper_src - src)
    axis = axis/np.sqrt(np.sum(axis**2))
    pos = src + axis*nubeam["XLBTNA"]

    if verbose:
        print('Source position: ',src)
        print('1st Aperture position: ',aper_src)
        print('Tangency position: ', pos)

    nbi = {"data_source":"TRANSP/NUBEAM namelist","name":nubeam["NAME"],
           "shape":nubeam["NBSHAP"],"src":src,"axis":axis,
           "focy":nubeam["FOCLR"],"focz":nubeam["FOCLZ"],
           "divy":np.repeat(nubeam["DIVR"],3),
           "divz":np.repeat(nubeam["DIVZ"],3),
           "widy":nubeam["BMWIDR"], "widz":nubeam["BMWIDZ"],
           "naperture":len(nubeam["NBAPSHA"]),"ashape":nubeam["NBAPSHA"],
           "awidy":nubeam["RAPEDGA"],"awidz":nubeam["XZPEDGA"],
           "aoffy":nubeam["XRAPOFFA"],"aoffz":nubeam["XZAPOFFA"],
           "adist":nubeam["XLBAPA"] }

    return nbi
