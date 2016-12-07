#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from lib.python_prefida.tb_zyx import tb_zyx


def uvw_to_xyz(alpha, beta, gamma, uvw, origin=None):
    #+##`uvw_to_xyz(alpha, beta, gamma, uvw, origin=[0,0,0])`
    #+ Express non-rotated coordinate `uvw` in rotated `xyz` coordinates
    #+###Arguments
    #+     **alpha**: Rotation angle about z [radians]
    #+
    #+     **beta**: Rotation angle about y' [radians]
    #+
    #+     **gamma**: Rotation angle about x" [radians]
    #+
    #+     **uvw**: Point in rotated coordinate system, (3)
    #+
    #+###Keyword Arguments
    #+     **origin**: Origin of rotated coordinate system in non-rotated (uvw) coordinates, (3)
    #+
    #+###Example Usage
    #+```python
    #+>>> xyz = uvw_to_xyz(np.pi/2., 0.0, np.pi/3., uvw, origin=[.1, .2, 0.])
    #+```
    if origin is None:
        origin = [0., 0., 0.]

    # Make np arrays
    uvw = np.array(uvw, dtype=float)
    origin = np.array(origin, dtype=float)

    # Do checks as this code does not allow multiple points to be entered (yet)
    if uvw.ndim != 1:
        raise ValueError('uvw must be 1D, but it has shape {}'.format(uvw.shape))
    if origin.ndim != 1:
        raise ValueError('origin must be 1D, but it has shape {}'.format(origin.shape))
    if uvw.size != 3:
        raise ValueError('uvw must have length 3, but it has length {}'.format(uvw.size))
    if origin.size != 3:
        raise ValueError('origin must have length 3, but it has length {}'.format(origin.size))

#    s = uvw.shape  # size(uvw,/dim)

#    if s.size != 2:
#        s = [s, 1]

#    if uvw.ndim != 2:

    # Shift origin
#    uvw_shifted = transpose(uvw - tile_array(origin,1,s[1]))
#    uvw_shifted = np.transpose(uvw - np.tile(origin, (1, s[1])))
    uvw_shifted = uvw - origin

    # Get rotation matrix
    r = tb_zyx(alpha, beta, gamma)

    # Apply rotation matrix
#    xyz = R ## uvw_shifted
    xyz = np.dot(r, uvw_shifted)

    return np.transpose(xyz)
###############################################################################
if __name__ == "__main__":

    a = uvw_to_xyz(1., 0.2, 0.1, [1.2, 2.3, 2.1], origin=[.1, .2, 0.])
    print(a)
