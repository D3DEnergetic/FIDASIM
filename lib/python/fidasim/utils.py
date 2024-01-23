#!/bin/sh
"exec" "$FIDASIM_DIR/deps/python" "$0" "$@"
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
import matplotlib.path as mplPath

is_sorted = lambda a: np.all(a[:-1] <= a[1:])

FIDASIM_default_COCOS = 5 
class COCOS:
    '''
    #+COCOS
    #+COCOS class object
    #+***
    '''
    def __init__(self, index=FIDASIM_default_COCOS):
        self.cocos = index
        self.exp_Bp = 0 if index < 10 else 1
        self.sigma_Bp = 1 if index in [1,2,5,6,11,12,15,16] else -1
        self.sigma_RphZ = 1 if index%2 != 0 else -1
        self.sigma_rhothph = 1 if index in [1,2,7,8,11,12,17,18] else -1
        self.sign_q = self.sigma_rhothph
        self.sign_pprime = -1*self.sigma_Bp 

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

def point_in_polygon(poly,points):
    path = mplPath.Path(poly)
    inside = path.contains_points(points)
    return inside

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
            raise ValueError('xyz must be (3, n), but it has shape {}'.format(xyz.shape))
        n = s[1]
    elif xyz.ndim == 1:
        if xyz.size != 3:
            raise ValueError('xyz must have length 3, but it has length {}'.format(xyz.size))
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
            'nr': nr,
            'nz': nz}

    if nphi > 1:
        grid.setdefault('nphi',nphi)
        grid.setdefault('phi', phi)

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

def detector_aperture_geometry(g,wn):
    #+#detector_aperture_geometry
    #+ Defines charged fusion product diagnostic geometry
    #+ Sources: detector_aperture_geometry.dat and MAST_29908_221_g_ensemble.nml
    #+ NOTE: This file is hardcoded to be useful for developers. General users could modify parameters under
    #+       wn = 0 or 1, and correspondingly set the wn flag when calling detector_aperture_geometry.pro
    #+       Alternatively, the user could manually define the geometry, see run_tests.pro and test_cfpd.pro
    #+***
    #+##Arguments
    #+    **g**: GEQDSK file
    #+
    #+    **wn**: indicates Werner (1) or Netepneko (0) diagnostic geometry
    #+
    #+##Outputs
    #+    **rdist**: Radial coordinates of detector [m]
    #+
    #+    **zdist**: Vertical coordinates of detector [m]
    #+
    #+    **v**: Detector orientations [m/s, rad/s, m/s]
    #+
    #+    **d**: Collimator length [m]
    #+
    #+    **rc**: Outer collimator spacing [m]
    #+
    #+##Example Usage
    #+```python
    #+>>> g = 'g000001.01000'
    #+>>> geometry = detector_aperture_geometry(g,0)
    #+```

    if wn == 1: # Source: MAST_29908_221_g_ensemble.nml
        print('Using MAST geometry 1 (Werner)')
        ## Detector orientation
        theta_port = np.radians(np.repeat(40.0, 4))

        #In MAST_29908_221_g_ensemble.nml, phi_port < 0, but it is made positive here for
        #comparison with detector_aperture_geometry.dat
        phi_port = np.radians([30.0, 35.0, 40.0, 45.0])

        v = np.zeros((3,4))
        v[0,:] = - np.sin(theta_port)*np.cos(phi_port)#    radial direction
        v[1,:] = - np.sin(phi_port)                #    toroidal direction
        v[2,:] = np.cos(theta_port)*np.cos(phi_port)  #    vertical direction

        ## Detector position
        ZDist = np.array([0.030014, 0.038311, 0.013981, 0.024414])

        RDist = np.array([1.668328, 1.661375, 1.648754, 1.64])

        PHDangle = np.array([1.38277674641, 1.39294543453, 1.37658569418, 1.3855051501])

        D = 0.0357      #                        detector-collimator spacing (m)

        RC = np.repeat(2.5e-3, 4) # outer collimator radius (m)

        RCD = np.repeat(2.5e-3, 4) # inner collimator radius (m)

    else: # Source: detector_aperture_geometry.dat
        print('Using MAST geometry 0 (Netepenko)')

        ## Detector orientation
        theta_port = np.radians([46.6875370324, 47.6648339458, 50.031360382, 51.5837100275])

        phi_port = np.radians([38.8895519328, 44.2509710868, 48.3160975078, 53.6006103875])

        v = np.zeros((3,4))
        v[0,:] = - np.sin(theta_port)*np.cos(phi_port)#    radial direction
        v[1,:] = - np.sin(phi_port)                #    toroidal direction
        v[2,:] = np.cos(theta_port)*np.cos(phi_port)  #    vertical direction

        ## Detector position
        ZDist = np.array([0.036701870576, 0.0409530530375, 0.0232888146809, 0.0301116448993])

        RDist = np.array([1.66830563971, 1.67508495819, 1.68813419386, 1.69658132076])

        PHDangle = np.radians([79.8774705956, 79.2421358615, 80.3144072462, 79.7395308235])

        D = 0.04      #                        detector-collimator spacing (m)

        RC = np.array([1.288e-3, 1.294e-3, 1.318e-3, 1.343e-3])# outer collimator radius (m)

        RCD = np.array([1.288e-3, 1.294e-3, 1.318e-3, 1.343e-3])# inner collimator radius (detector radius) (m)

    return {'rdist':RDist,'zdist':ZDist,'v':v,'d':D,'rc':RC}

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

def read_geqdsk(filename, grid, poloidal=False, ccw_phi=True, exp_Bp=0, **convert_COCOS_kw):
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
    #+    **ccw_phi**: Argument for identify_COCOS. Toroidal direction from top view, True if counter-clockwise, False if clockwise
    #+
    #+    **exp_Bp**: Argument for identify_COCOS. 0 if poloidal flux divided by 2 pi, 1 if using effective poloidal flux
    #+
    #+    **convert_COCOS_kw**: Keyword arguments to pass to convert_COCOS via transform_COCOS_from_geqdsk
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
    g = transform_COCOS_from_geqdsk(efit.readg(filename), ccw_phi=ccw_phi, exp_Bp=exp_Bp, **convert_COCOS_kw)
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

def transform_COCOS_from_geqdsk(g, ccw_phi=True, exp_Bp=0, **convert_COCOS_kw):
    '''
    #+#transform_COCOS_from_geqdsk
    #+Identifies the COCOS index of a GEQDSK dictionary object and converts to fidasim COCOS
    #+Reference:
    #+    O. Sauter and S. Yu. Medvedev, Tokamak Coordinate Conventions: COCOS, 
    #+    Computer Physics Communications 184, 293 (2013).
    #+***
    #+##Arguments
    #+    **g**: GEQDSK dictionary object
    #+
    #+##Keyword Arguments
    #+    **ccw_phi**: Toroidal direction from top view, True if counter-clockwise, False if clockwise
    #+
    #+    **exp_Bp**: 0 if poloidal flux divided by 2 pi, 1 if using effective poloidal flux
    #+
    #+    **conver_COCOS_kw**: Keyword arguments for convert_COCOS 
    #+
    #+##Return Value
    #+geqdsk dict with converted COCOS
    #+
    #+##Example Usage
    #+```python
    #+>>> g = efit.readg(filename)
    #+>>> g = read_COCOS_from_geqdsk(g)
    #+```
    '''
    cc_out = COCOS() # cc_out = FIDASIM COCOS
    
    cc_in = identify_COCOS(g, ccw_phi=ccw_phi, exp_Bp=exp_Bp)
    
    if cc_in.cocos != cc_out.cocos:
        return convert_COCOS(g.copy(), cc_in, cc_out, **convert_COCOS_kw)
    else:
        return g.copy()

def identify_COCOS(g, ccw_phi=True, exp_Bp=0):
    '''
    #+#identify_COCOS
    #+Identifies the COCOS index of a GEQDSK dictionary object
    #+Reference:
    #+    O. Sauter and S. Yu. Medvedev, Tokamak Coordinate Conventions: COCOS, 
    #+    Computer Physics Communications 184, 293 (2013).
    #+***
    #+##Arguments
    #+    **g**: GEQDSK dictionary object
    #+
    #+##Keyword Arguments
    #+    **ccw_phi**: Toroidal direction from top view, True if counter-clockwise, False if clockwise
    #+
    #+    **exp_Bp**: 0 if poloidal flux divided by 2 pi, 1 if using effective poloidal flux
    #+
    #+##Return Value
    #+COCOS object
    #+
    #+##Example Usage
    #+```python
    #+>>> g = efit.readg(filename)
    #+>>> g_cocos = identify_COCOS(g)
    #+```
    '''
    # Sauter, eq. 22    
    sigma_Bp_in = -1 * np.sign(g['pprime'][0] * g['current'])
    sigma_RphZ_in = 1 if ccw_phi else -1
    sigma_rhothph_in = np.sign(g['qpsi'][0] * g['current'] * g['bcentr'])
    
    sigmas = [sigma_Bp_in, sigma_RphZ_in, sigma_rhothph_in]
 
    if sigmas == [1, 1, 1]:
        index = 1
    elif sigmas == [1, -1, 1]:
        index = 2
    elif sigmas == [-1, 1, -1]:
        index = 3
    elif sigmas == [-1, -1, -1]:
        index = 4
    elif sigmas == [1, 1, -1]:
        index = 5
    elif sigmas == [1, -1, -1]:
        index = 6
    elif sigmas == [-1, 1, 1]:
        index = 7
    elif sigmas == [-1, -1, 1]:
        index = 8
    else:
        index = FIDASIM_default_COCOS
    
    return COCOS(index+(10*exp_Bp))

def convert_COCOS(g, cc_in, cc_out, sigma_Ip=None, sigma_B0=None, l_d=[1,1], l_B=[1,1], exp_mu0=[0,0]):
    '''
    #+#convert_COCOS
    #+Converts a GEQDSK dictionary according to cc_in --> cc_out
    #+Reference:
    #+    O. Sauter and S. Yu. Medvedev, Tokamak Coordinate Conventions: COCOS, 
    #+    Computer Physics Communications 184, 293 (2013).
    #+***
    #+##Arguments
    #+    **g**: GEQDSK dictionary object
    #+
    #+    **cc_in**: COCOS object, input
    #+
    #+    **cc_out**: COCOS object, output
    #+
    #+##Keyword Arguments
    #+    **sigma_Ip**: Tuple of current sign, (in, out)
    #+
    #+    **sigma_B0**: Tuple of toroidal field sign, (in, out)
    #+
    #+    **l_d**: Tuple of length scale, (in, out)
    #+
    #+    **l_B**: Tuple of field magnitude scale, (in, out)
    #+
    #+    **exp_mu0**: Tuple of exponents for mu0, (in, out)
    #+
    #+##Return Value
    #+GEQDSK dictionary object
    #+
    #+##Example Usage
    #+```python
    #+>>> g = efit.readg(filename)
    #+>>> cc_out = COCOS(3)
    #+>>> cc_in = COCOS(1)
    #+>>> converted_g = convert_COCOS(g, cc_in, cc_out)
    #+```
    '''
    # Sauter, Appendix C
    print(f'CONVERT_COCOS: cocos_in ({cc_in.cocos}) != cocos_out ({cc_out.cocos}), applying COCOS conversion.')
    mu0 = 4*np.pi*1e-7
        
    l_d_eff = l_d[1]/l_d[0]
    l_B_eff = l_B[1]/l_B[0]
    exp_mu0_eff = exp_mu0[1] - exp_mu0[0]

    exp_Bp_eff = cc_out.exp_Bp - cc_in.exp_Bp
    sigma_Bp_eff = cc_out.sigma_Bp * cc_in.sigma_Bp
    sigma_RphZ_eff = cc_out.sigma_RphZ * cc_in.sigma_RphZ
    sigma_rhothph_eff = cc_out.sigma_rhothph * cc_in.sigma_rhothph

    if sigma_Ip is None:
        sigma_Ip_eff = sigma_RphZ_eff
    else:
        sigma_Ip_eff = np.prod(sigma_Ip)

    if sigma_B0 is None:
        sigma_B0_eff = sigma_RphZ_eff
    else:
        sigma_B0_eff = np.prod(sigma_B0)

    for key, val in g.items():
        if key in ['r', 'rdim', 'rleft', 'rbbbs', 'rlim', 'rcentr', 'rmaxis', 'z', 'zdim', 'zmid', 'zbbbz', 'zlim', 'zmaxis', 'nbdry', 'lim']:
            g[key] = np.array(val) * l_d_eff
        elif key == 'bcentr':
            g[key] = np.array(val) * l_B_eff * sigma_B0_eff
        elif key in ['simag', 'ssimag', 'sibry', 'ssibry', 'psi', 'psirz']:
            g[key] = np.array(val) * sigma_Ip_eff * sigma_Bp_eff * ((2*np.pi)**exp_Bp_eff) * (l_d_eff**2) * (l_B_eff)
        elif key == 'current':
            g[key] = np.array(val) * sigma_Ip_eff * l_d_eff * l_B_eff / (mu0**exp_mu0_eff)
        elif key == 'fpol':
            g[key] = np.array(val) * sigma_B0_eff * l_d_eff * l_B_eff
        elif key == 'pres':
            # Sauter, end of Sec. 4
            g[key] = np.array(val) * (l_d_eff**2) / (mu0**exp_mu0_eff)
        elif key == 'ffprim':
            g[key] = np.array(val) * sigma_Ip_eff * sigma_Bp_eff / ((2*np.pi)**exp_Bp_eff) * l_B_eff
        elif key == 'pprime':
            g[key] = np.array(val) * sigma_Ip_eff * sigma_Bp_eff / ((2*np.pi)**exp_Bp_eff) * l_B_eff / ((mu0**exp_mu0_eff) * (l_d_eff**2))
        elif key == 'qpsi':
            g[key] = np.array(val) * sigma_Ip_eff * sigma_B0_eff * sigma_rhothph_eff
    
    return g

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
    plasma={"data_source":os.path.abspath(filename),"time":time,"impurity_charge":int(impurity_charge),
            "nthermal":int(np.sum(w_ai)), "species_mass":ai[w_ai], "deni":deni[w_ai,:,:],"profiles":profiles,
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

def read_imas_nbi(imas_file, beam_index, time_index):
    g = h5py.File(imas_file,'r')

    path = 'nbi/unit/{}/'.format(beam_index)

    ab = g[path+'species/a'][()]
    Zb = g[path+'species/z_n'][()]
    einj = g[path+'energy/data'][()][time_index]*1e-3 #keV
    pinj = g[path+'power_launched/data'][()][time_index]*1e-6 #MW
    cfracs = g[path+'beam_current_fraction/data'][()][:,time_index]

    shape = g[path+'source/geometry_type'][()]
    if shape == 3:
        shape = 1
    elif shape == 2
        shape == 2
    else
        shape = 1

    src_r = g[path+'source/centre/r'][()]
    src_z = g[path+'source/centre/z'][()]
    src_phi = g[path+'source/centre/phi'][()]

    src_u = src_r*np.cos(src_phi)
    src_v = src_r*np.sin(src_phi)
    uvw_src = np.array([src_u,src_v,src_z])
    uvw_axis = np.array([g[path+'source/x3_unit_vector/{}'.format(i)][()] for i in ['x','y','z'])

    # FIDASIM doesnt support beamlets yet so just average
    nbeamlets = len(g[path+'beamlets_group'].keys())
    widz = np.mean([g[path+'beamlets_group/{}/width_vertical'.format(i)][()] for i in range(nbeamlets)])
    widy = np.mean([g[path+'beamlets_group/{}/width_horizontal'.format(i)][()] for i in range(nbeamlets)])

    focz = np.mean([g[path+'beamlets_group/{}/focus/focal_length_vertical'.format(i)][()] for i in range(nbeamlets)])
    focy = np.mean([g[path+'beamlets_group/{}/focus/focal_length_horizontal'.format(i)][()] for i in range(nbeamlets)])

    divy = np.array([np.mean([g[path+'beamlets_group/{}/divergence_component/{}/horizontal'.format(i,j)][()] for i in range(nbeamlets)]) for j in range(3)])
    divz = np.array([np.mean([g[path+'beamlets_group/{}/divergence_component/{}/vertical'.format(i,j)][()] for i in range(nbeamlets)]) for j in range(3)])


    naper = len(g[path+'aperture'].keys())
    aper_shape = [g[path+'aperture/{}/geometry_type'.format(i)][()] for i in range(naper)]
    ashape = np.array([1 if a == 3 else 2 for a in aper_shape]) # ignores outline shape

    aper_rs = [g[path+'aperture/{}/centre/r'.format(i)][()] for i in range(naper)]
    aper_zs = [g[path+'aperture/{}/centre/z'.format(i)][()] for i in range(naper)]
    aper_phis = [g[path+'aperture/{}/centre/phi'.format(i)][()] for i in range(naper)]

    aper_xs = [r*np.cos(p) for (r,p) in zip(aper_rs,aper_phis)]
    aper_ys = [r*np.sin(p) for (r,p) in zip(aper_rs,aper_phis)]

    sx = src_u
    sy = src_v
    sz = src_z
    adist = np.array([np.sqrt((ax - sx)**2 + (ay - sy)**2 + (az - sz)**2) for (ax,ay,az) in zip(aper_xs,aper_ys,aper_zs)])

    #assume zero offset
    aoffz = np.zeros(naper)
    aoffy = np.zeros(naper)

    awidy = np.array([g[path+'aperture/{}/x1_width'.format(i)][()] for i in range(naper)]
    awidz = np.array([g[path+'aperture/{}/x2_width'.format(i)][()] for i in range(naper)]

    nbi = {"name":"beam: {}".format(beam_index),"shape":shape,"data_source":imas_file,
           "src":uvw_src, "axis":uvw_axis, "widy":widy, "widz":widz,
           "divy":divy, "divz":divz, "focy":focy, "focz":focz,
           "naperture":naper, "ashape":ashape, "adist":adist,
           "awidy":awidy, "awidz":awidz, "aoffy":aoffy, "aoffz":aoffz}

    g.close()

    return nbi, {"ab":ab,"Zb":Zb,"einj":einj,"pinj":pinj,"current_fractions":cfracs}

def read_imas_equilbrium(imas_file, time_index; impurity_charge=None):
    g = h5py.File(imas_file)

    path = "core_profiles/profiles_1d/{}/".format(time_index)

    psi = g[path+"grid/psi"][()]
    te = g[path+"electrons/temperature"][()]/1000 #keV
    dene = g[path+"electrons/density"][()]*1e-6 # cm-3
    zeff = g[path+"zeff"][()]
    denn = g[path+"neutral/0/density"][()]*1e-6 # cm-3
    vtor = g[path+"ion/0/velocity/toroidal"][()]*100 #cm/s

    # different number of ions
    ti = [g[path+"ion/{}/temperature".format(i)][()]/1000 for i in g[path+"ion"].keys()]#keV
    deni = [g[path+"ion/{}/density".format(i)][()]*1e-6 for i in g[path+"ion"].keys()]#keV
    species_mass = np.array([g[path+"ion/{}/element/0/a".format(i)][()]*1.007 for i in g[path+"ion"].keys()]#keV
    species_charge = np.array([g[path+"ion/{}/element/0/z_n".format(i)][()] for i in g[path+"ion"].keys()]#keV
    nthermal = len(species_mass)

    # get fields
    path = "equilibrium/time_slice/{}/profiles_2d/0/".format(time_index)

    r = g[path+"grid/dim1"][()]*100
    nr = len(r)
    z = g[path+"grid/dim2"][()]*100
    nz = len(z)
    r2d = np.tile(r, (nz, 1)).T
    z2d = np.tile(z, (nr, 1))
    br = g[path+"b_field_r"][()]
    bz = g[path+"b_field_z"][()]
    bt = g[path+"b_field_tor"][()]
    psi_rz = g[path+"psi"][()]
    bdry_r = g[path+"boundary/outline/r"][()][:-1]*100
    bdry_z = g[path+"boundary/outline/z"][()][:-1]*100
    bdry_points = [(r,z) for (r,z) in zip(bdry_r,bdry_z)]
    mask = np.zeros(psi_rz.shape,dtype='int')
    r2d_flat = r2d.flatten()
    z2d_flat = z2d.flatten()
    points = np.concatenate([r2d_flat.reshape(1,len(r2d_flat)),z2d_flat.reshape(1,len(z2d_flat))],axis=0)
    w = point_in_polygon(bdry_points,points)
    mask[w] = 1

    if !is_sorted(psi):
        psi = psi[::-1]
        te = te[::-1]
        dene = dene[::-1]
        zeff = zeff[::-1]
        denn = denn[::-1]
        ti = [ti_i[::-1] for ti_i in ti]
        deni = [deni_i[::-1] for deni_i in deni]
        vtor = vtor[::-1]

    te_rz = interp1d(psi,te,'cubic',fill_value='extrapolate')(psi_rz)
    dene_rz = interp1d(psi,dene,'cubic',fill_value='extrapolate')(psi_rz)
    zeff_rz = interp1d(psi,zeff,'cubic',fill_value='extrapolate')(psi_rz)
    denn_rz = interp1d(psi,denn,'cubic',fill_value='extrapolate')(psi_rz)
    vtor_rz = interp1d(psi,vtor,'cubic',fill_value='extrapolate')(psi_rz)

    # for ti assume the same temperature (using average over all ions)
    ti_rz = np.mean(np.concatenate([interp1d(psi,ti_i,'cubic',fill_value='extrapolate')(psi_rz).reshape((1,) + psi_rz.shape) for ti_i in ti],axis=0),axis=0)
    deni_rz = np.concatenate([interp1d(psi,deni_i,'cubic',fill_value='extrapolate')(psi_rz).reshape((1,) + psi_rz.shape) for deni_i in deni],axis=0)

    vt_rz = r2d*vtor_rz
    vr_rz = 0*vt_rz
    vz_rz = 0*vt_rz

    if impurity_charge is None:
        impurity_charge = np.max(species_charge) # assume impurity charge is largest 

    denimp_rz = dene_rz * (zeff_rz - 1.0) / (impurity_charge*(impurity_charge-1))

    plasma = {"nr":nr, "nz":nz,"r":r,"z":z,"time":0.0,"r2d":r2d,"z2d":z2d,
              "data_source":imas_file, impurity_charge:impurity_charge,"nthermal":nthermal, "species_mass":species_mass,"mask":mask,
              "te":te_rz,"ti":ti_rz,"dene":dene_rz,"denn":denn_rz,"zeff":zeff_rz,"deni":deni_rz,"denimp":denimp_rz,
              "vr":vr_rz,"vz":vz_rz,"vt":vt_rz,"description":"plasma parameters","coordinate_system":"cylindrical"}
    
    equil = {"nr":nr,"nz":nz,"r":r,"z":z,"time":0.0,"r2d":r2d,"z2d":z2d,
             "data_source":imas_file, "mask":mask, "br":br, "bz":bz, "bt":bt, "er":0*bt, "ez":0*bt, "et":0*bt,
             "description":"fields","coordinate_system":"cylindrical"}

    g.close()

    return plasma, equil

