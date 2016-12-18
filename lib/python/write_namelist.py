#!/usr/bin/env python
# -*- coding: utf-8 -*-

from lib.python_prefida.info import info
from lib.python_prefida.get_version import get_version
from lib.python_prefida.get_fidasim_dir import get_fidasim_dir
from lib.python_prefida.success import success
import datetime


def write_namelist(filename, inputs):
    #+##`write_namelist, filename, inputs`
    #+Writes namelist file
    #+
    #+###Input Arguments
    #+     **filename**: Name of the namelist file
    #+
    #+     **inputs**: Input dictionary
    #+
    #+###Example Usage
    #+```python
    #+>>> write_namelist(filename, inputs)
    #+```
    info("Writing namelist file...")

    fidasim_version = get_version(get_fidasim_dir())

    with open(filename, "w") as f:
        f.write("!! Created: {}\n".format(datetime.datetime.now()))
        f.write("!! FIDASIM version: {}\n".format(fidasim_version))
        f.write("!! Comment: {}\n".format(inputs['comment']))
        f.write("&fidasim_inputs\n\n")

        f.write("!! Shot Info\n")
        f.write("shot = {:d}    !! Shot Number\n".format(inputs['shot']))
        f.write("time = {:f}    !! Time [s]\n".format(inputs['time']))
        f.write("runid = '{}'   !! runID\n".format(inputs['runid']))
        f.write("result_dir = '{}'    !! Result Directory\n\n".format(inputs['result_dir']))

        f.write("!! Input Files\n")
        f.write("tables_file = '{}'   !! Atomic Tables File\n".format(inputs['tables_file']))
        f.write("equilibrium_file = '" + inputs['equilibrium_file'] + "'    !! File containing plasma parameters and fields\n")
        f.write("geometry_file = '" + inputs['geometry_file'] + "'    !! File containing NBI and diagnostic geometry\n")
        f.write("distribution_file = '" + inputs['distribution_file'] + "'    !! File containing fast-ion distribution\n\n")

        f.write("!! Simulation Switches\n")
        f.write("calc_bes = {:d}    !! Calculate Beam Emission and Halo Spectra\n".format(inputs['calc_bes']))
        f.write("calc_brems = {:d}    !! Calculate Bremsstrahlung\n".format(inputs['calc_brems']))
        f.write("calc_fida = {:d}    !! Calculate FIDA Spectra\n".format(inputs['calc_fida']))
        f.write("calc_npa = {:d}   !! Calculate NPA\n".format(inputs['calc_npa']))
        f.write("calc_birth = {:d}    !! Calculate Birth Profile\n".format(inputs['calc_birth']))
        f.write("calc_fida_wght = {:d}    !! Calculate FIDA weights\n".format(inputs['calc_fida_wght']))
        f.write("calc_npa_wght = {:d}    !! Calculate NPA weights\n".format(inputs['calc_npa_wght']))
        f.write("dump_dcx = {:d}    !! Dump DCX neutrals and spectra\n\n".format(inputs['dump_dcx']))

        f.write("!! Debugging Switches\n")
        f.write("no_flr = {:d}    !! Turn off Finite Larmor Radius effects\n".format(inputs['no_flr']))
        f.write("load_neutrals = {:d}    !! Load neutrals from neutrals file\n".format(inputs['load_neutrals']))
        f.write("neutrals_file = '" + inputs['neutrals_file'] + "'    !! File containing the neutral density\n")
        f.write("verbose = {:d}    !! Verbose\n\n".format(inputs['verbose']))

        f.write("!! Monte Carlo Settings\n")
        f.write("n_fida = {:d}    !! Number of FIDA mc particles\n".format(inputs['n_fida']))
        f.write("n_npa = {:d}    !! Number of NPA mc particles\n".format(inputs['n_npa']))
        f.write("n_nbi = {:d}    !! Number of NBI mc particles\n".format(inputs['n_nbi']))
        f.write("n_halo = {:d}    !! Number of HALO mc particles\n".format(inputs['n_halo']))
        f.write("n_dcx = {:d}     !! Number of DCX mc particles\n".format(inputs['n_dcx']))
        f.write("n_birth = {:d}    !! Number of BIRTH mc particles\n\n".format(inputs['n_birth']))

        f.write("!! Neutral Beam Settings\n")
        f.write("ab = {:f}     !! Beam Species mass [amu]\n".format(inputs['ab']))
        f.write("pinj = {:f}     !! Beam Power [MW]\n".format(inputs['pinj']))
        f.write("einj = {:f}     !! Beam Energy [keV]\n".format(inputs['einj']))
        f.write("current_fractions(1) = {:f} !! Current Fractions (Full component)\n".format(inputs['current_fractions'][0]))
        f.write("current_fractions(2) = {:f} !! Current Fractions (Half component)\n".format(inputs['current_fractions'][1]))
        f.write("current_fractions(3) = {:f} !! Current Fractions (Third component)\n\n".format(inputs['current_fractions'][2]))

        f.write("!! Plasma Settings\n")
        f.write("ai = {:f}     !! Ion Species mass [amu]\n".format(inputs['ai']))
        f.write("impurity_charge = {:d}     !! Impurity Charge\n\n".format(inputs['impurity_charge']))

        f.write("!! Beam Grid Settings\n")
        f.write("nx = {:d}    !! Number of cells in X direction (Into Plasma)\n".format(inputs['nx']))
        f.write("ny = {:d}    !! Number of cells in Y direction\n".format(inputs['ny']))
        f.write("nz = {:d}    !! Number of cells in Z direction\n".format(inputs['nz']))
        f.write("xmin = {:f}     !! Minimum X value [cm]\n".format(inputs['xmin']))
        f.write("xmax = {:f}     !! Maximum X value [cm]\n".format(inputs['xmax']))
        f.write("ymin = {:f}     !! Minimum Y value [cm]\n".format(inputs['ymin']))
        f.write("ymax = {:f}     !! Maximum Y value [cm]\n".format(inputs['ymax']))
        f.write("zmin = {:f}     !! Minimum Z value [cm]\n".format(inputs['zmin']))
        f.write("zmax = {:f}     !! Maximum Z value [cm]\n\n".format(inputs['zmax']))

        f.write("!! Tait-Bryan Angles for z-y`-x`` rotation\n")
        f.write("alpha = {:f}     !! Rotation about z-axis [rad]\n".format(inputs['alpha']))
        f.write("beta  = {:f}     !! Rotation about y`-axis [rad]\n".format(inputs['beta']))
        f.write("gamma = {:f}     !! Rotation about x``-axis [rad]\n\n".format(inputs['gamma']))

        f.write("!! Beam Grid origin in machine coordinates (cartesian)\n")
        f.write("origin(1) = {:f}     !! U value [cm]\n".format(inputs['origin'][0]))
        f.write("origin(2) = {:f}     !! V value [cm]\n".format(inputs['origin'][1]))
        f.write("origin(3) = {:f}     !! W value [cm]\n\n".format(inputs['origin'][2]))

        f.write("!! Wavelength Grid Settings\n")
        f.write("nlambda = {:d}    !! Number of Wavelengths\n".format(inputs['nlambda']))
        f.write("lambdamin = {:f}    !! Minimum Wavelength [nm]\n".format(inputs['lambdamin']))
        f.write("lambdamax = {:f}    !! Maximum Wavelength [nm]\n\n".format(inputs['lambdamax']))

        f.write("!! Weight Function Settings\n")
        f.write("ne_wght = {:d}    !! Number of Energies for Weights\n".format(inputs['ne_wght']))
        f.write("np_wght = {:d}    !! Number of Pitches for Weights\n".format(inputs['np_wght']))
        f.write("nphi_wght = {:d}    !! Number of Gyro-angles for Weights\n".format(inputs['nphi_wght']))
        f.write("emax_wght = {:f}    !! Maximum Energy for Weights [keV]\n".format(inputs['emax_wght']))
        f.write("nlambda_wght = {:d}    !! Number of Wavelengths for Weights \n".format(inputs['nlambda_wght']))
        f.write("lambdamin_wght = {:f}    !! Minimum Wavelength for Weights [nm]\n".format(inputs['lambdamin_wght']))
        f.write("lambdamax_wght = {:f}    !! Maximum Wavelength for Weights [nm]\n\n".format(inputs['lambdamax_wght']))
        f.write("/\n\n")

    success("Namelist file created: {}\n".format(filename))
