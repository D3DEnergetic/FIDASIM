#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import copy


def aabb_intersect(rc, dr, r0, d0):
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
    v0 = d0 / np.sqrt(np.sum(d0 ** 2.))

    # There are 6 sides to a cube/grid
    side_inter = np.zeros(6)

    # Intersection points of ray with planes defined by grid
    ipnts = np.zeros((3, 6))

    # Find whether ray intersects each side
    for i in range(6):
        j = np.floor(i / 2)
        ind = ([0, 1, 2] != j)
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
    w = (side_inter != 0)
    nw = side_inter[w].size
    if nw >= 2:
        #Find two unique intersection points
        nunique = 0
        for i in range(nw - 1):
            if np.sum(ipnts[:, w[0]] == ipnts[:, w[i + 1]]) != 3:
                w = [w[0], w[i + 1]]
                nunique = 2
                break

        if nunique == 2:
            vi = ipnts[:, w[1]] - ipnts[:, w[0]]
            vi = vi / np.sqrt(np.sum(vi ** 2.))
            dot_prod = np.sum(v0 * vi)
            if dot_prod > 0.0:
                r_enter = ipnts[:, w[0]]
                r_exit = ipnts[:, w[1]]
            else:
                r_enter = ipnts[:, w[1]]
                r_exit = ipnts[:, w[0]]

            # Calculate intersection length
            intersect = np.sqrt(np.sum((r_exit - r_enter) ** 2.))

    return intersect, r_enter, r_exit
