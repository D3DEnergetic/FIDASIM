#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


def uvw_to_xyz(alpha, beta, gamma, uvw, origin=None):
    """
    ##`uvw_to_xyz(alpha, beta, gamma, uvw, origin=[0,0,0])`
     Express non-rotated coordinate `uvw` in rotated `xyz` coordinates
    ###Arguments
         **alpha**: Rotation angle about z [radians]

         **beta**: Rotation angle about y' [radians]

         **gamma**: Rotation angle about x" [radians]

         **xyz**: Point in rotated coordinate system

    ###Keyword Arguments
         **origin**: Origin of rotated coordinate system in non-rotated (uvw) coordinates.

    ###Example Usage
    ```idl
    IDL> xyz = uvw_to_xyz(!DPI/2,0.0,!DPI/3,uvw)
    ```
    """
    if origin is None:
        origin = [0., 0., 0.]

    s = uvw.shape  # size(uvw,/dim)

    if s.size != 2:
        s = [s, 1]

#    uvw_shifted = transpose(uvw - tile_array(origin,1,s[1]))
    uvw_shifted = np.transpose(uvw - np.tile(origin, (1, s[1])))

    R = np.transpose(tb_zyx(alpha, beta, gamma))

#    xyz = R ## uvw_shifted
    xyz = np.dot(R, uvw_shifted)

    return np.transpose(xyz)
