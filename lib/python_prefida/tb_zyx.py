#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


def tb_zyx(a, b, g):
    """
    #tb_zyx
    Calculates Tait-Bryan z-y'-x" active rotation matrix given rotation angles `alpha`,`beta`,`gamma` in radians
    ##Arguments
         **alpha**: rotation angle about z [radians]

         **beta**: rotation angle about y' [radians]

         **gamma**: rotation angle about x" [radians]

    ##Return Value
         Rotation Matrix [prefida](|url|/sourcefile/prefida.pro.html)

    ##Example Usage
    ```idl
     IDL> rot_mat = tb_zyx(!DPI/2, 0.0, !DPI/3)
    ```
    """
    sa = np.sin(a)
    ca = np.cos(a)
    sb = np.sin(b)
    cb = np.cos(b)
    sg = np.sin(g)
    cg = np.cos(g)

    R = np.zeros((3, 3))

    R[0, 0] = ca * cb
    R[1, 0] = ca * sb * sg - cg * sa
    R[2, 0] = sa * sg + ca * cg * sb
    R[0, 1] = cb * sa
    R[1, 1] = ca * cg + sa * sb * sg
    R[2, 1] = cg * sa * sb - ca * sg
    R[0, 2] = -sb
    R[1, 2] = cb * sg
    R[2, 2] = cb * cg

    return R
