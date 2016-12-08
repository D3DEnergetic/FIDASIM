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
    #+     **uvw**: Point in rotated coordinate system, (3, n)
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
    xyz = np.dot(r, uvw_shifted)

    return xyz
###############################################################################
if __name__ == "__main__":

    uvw = np.tile(np.array([1.2, 2.3, 2.1]), (5, 1)).T      # (3, 5)
    a = uvw_to_xyz(1., 0.2, 0.1, uvw, origin=[.1, .2, 0.])
    print(a.shape)
