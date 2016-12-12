#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


def tb_zyx(alpha, beta, gamma):
    #+#tb_zyx
    #+Calculates Tait-Bryan z-y'-x" active rotation matrix given rotation angles `alpha`,`beta`,`gamma` in radians
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
    #+```dist
    #+ >>> rot_mat = tb_zyx(!DPI/2, 0.0, !DPI/3)
    #+```
    sa = np.sin(alpha)
    ca = np.cos(alpha)
    sb = np.sin(beta)
    cb = np.cos(beta)
    sg = np.sin(gamma)
    cg = np.cos(gamma)

    r = np.zeros((3, 3))

    r[0, 0] = ca * cb
    r[1, 0] = ca * sb * sg - cg * sa
    r[2, 0] = sa * sg + ca * cg * sb
    r[0, 1] = cb * sa
    r[1, 1] = ca * cg + sa * sb * sg
    r[2, 1] = cg * sa * sb - ca * sg
    r[0, 2] = -sb
    r[1, 2] = cb * sg
    r[2, 2] = cb * cg

    return r
