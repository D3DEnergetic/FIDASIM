#!/bin/sh
"exec" "$FIDASIM_DIR/deps/python" "$0" "$@"

import argparse
import numpy as np
from scipy.interpolate import interp1d
import fidasim as fs

fida_dir = fs.utils.get_fidasim_dir()
test_dir = fida_dir + '/test'
grid = fs.utils.rz_grid(100.0, 240.0, 70, -100.0, 100.0, 100)
equil, rho, btipsign = fs.utils.read_geqdsk(test_dir+'/g000001.01000', grid)
#equil["ez"] = equil["ez"]/equil["ez"]
print(equil["bz"])
print(equil["br"])
print(equil["bt"])
