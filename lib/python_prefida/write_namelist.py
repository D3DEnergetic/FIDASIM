#!/usr/bin/env python
# -*- coding: utf-8 -*-

from python_prefida import info
from python_prefida import get_version
from python_prefida import get_fidasim_dir
from python_prefida import success
import datetime


def write_namelist(filename, inputs):
    """
    ;+##`write_namelist, filename, inputs`
    ;+Writes namelist file
    ;+
    ;+###Input Arguments
    ;+     **filename**: Name of the namelist file
    ;+
    ;+     **inputs**: Input structure
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> write_namelist, filename, inputs
    ;+```
    """
    info("Writing namelist file...")

    fidasim_version = ''  # get_version(get_fidasim_dir())

    with open(filename, "w") as f:

        f.write("!! Created: {}".format(datetime.datetime.now()))
        f.write("!! FIDASIM version: {}".format(fidasim_version))
        f.write("!! Comment: {}".format(inputs['comment']))
        f.write("&fidasim_inputs")
        f.write("")

        f.write("!! Shot Info")
        f.write("shot = {:d}    !! Shot Number".format(inputs['shot']))
        f.write("time = {:f}    !! Time [s]".format(inputs['time']))
        f.write("runid = {}   !! runID".format(inputs['runid']))
        f.write("result_dir = '{}'    !! Result Directory".format(inputs['result_dir']))
        f.write("")

        f.write("!! Input Files")
        f.write("tables_file = '{}'   !! Atomic Tables File".format(inputs['tables_file']))
        f.write("equilibrium_file = '" + inputs['equilibrium_file'] + "'    !! File containing plasma parameters and fields")
        f.write("geometry_file = '" + inputs['geometry_file'] + "'    !! File containing NBI and diagnostic geometry")
        f.write("distribution_file = '" + inputs['distribution_file'] + "'    !! File containing fast-ion distribution")
        f.write("")

        f.write("!! Simulation Switches")
        f.write("calc_bes = {:d}    !! Calculate Beam Emission and Halo Spectra".format(inputs['calc_bes']))
        f.write("calc_brems = {:d}    !! Calculate Bremsstrahlung".format(inputs['calc_brems']))
        f.write("calc_fida = {:d}    !! Calculate FIDA Spectra".format(inputs['calc_fida']))
        f.write("calc_npa = {:d}   !! Calculate NPA".format(inputs['calc_npa']))
        f.write("calc_birth = {:d}    !! Calculate Birth Profile".format(inputs['calc_birth']))
        f.write("calc_fida_wght = {:d}    !! Calculate FIDA weights".format(inputs['calc_fida_wght']))
        f.write("calc_npa_wght = {:d}    !! Calculate NPA weights".format(inputs['calc_npa_wght']))
        f.write("dump_dcx = {:d}    !! Dump DCX neutrals and spectra".format(inputs['dump_dcx']))
        f.write("")

        f.write("!! Debugging Switches")
        f.write("no_flr = {:d}    !! Turn off Finite Larmor Radius effects".format(inputs['no_flr']))
        f.write("load_neutrals = {:d}    !! Load neutrals from neutrals file".format(inputs['load_neutrals']))
        f.write("neutrals_file = '" + inputs['neutrals_file'] + "'    !! File containing the neutral density")
        f.write("verbose = {:d}    !! Verbose".format(inputs['verbose']))
        f.write("")

        f.write("!! Monte Carlo Settings")
        f.write("n_fida = {:d}    !! Number of FIDA mc particles".format(inputs['n_fida']))
        f.write("n_npa = {:d}    !! Number of NPA mc particles".format(inputs['n_npa']))
        f.write("n_nbi = {:d}    !! Number of NBI mc particles".format(inputs['n_nbi']))
        f.write("n_halo = {:d}    !! Number of HALO mc particles".format(inputs['n_halo']))
        f.write("n_dcx = {:d}     !! Number of DCX mc particles".format(inputs['n_dcx']))
        f.write("n_birth = {:d}    !! Number of BIRTH mc particles".format(inputs['n_birth']))
        f.write("")

        f.write("!! Neutral Beam Settings")
        f.write("ab = {:f}     !! Beam Species mass [amu]".format(inputs['ab']))
        f.write("pinj = {:f}     !! Beam Power [MW]".format(inputs['pinj']))
        f.write("einj = {:f}     !! Beam Energy [keV]".format(inputs['einj']))
        f.write("current_fractions(1) = {:f} !! Current Fractions (Full component)".format(inputs['current_fractions'][0]))
        f.write("current_fractions(2) = {:f} !! Current Fractions (Half component)".format(inputs['current_fractions'][1]))
        f.write("current_fractions(3) = {:f} !! Current Fractions (Third component)".format(inputs['current_fractions'][2]))
        f.write("")

        f.write("!! Plasma Settings")
        f.write("ai = {:f}     !! Ion Species mass [amu]".format(inputs['ai']))
        f.write("impurity_charge = {:d}     !! Impurity Charge".format(inputs['impurity_charge']))
        f.write("")

        f.write("!! Beam Grid Settings")
        f.write("nx = {:d}    !! Number of cells in X direction (Into Plasma)".format(inputs['nx']))
        f.write("ny = {:d}    !! Number of cells in Y direction".format(inputs['ny']))
        f.write("nz = {:d}    !! Number of cells in Z direction".format(inputs['nz']))
        f.write("xmin = {:f}     !! Minimum X value [cm]".format(inputs['xmin']))
        f.write("xmax = {:f}     !! Maximum X value [cm]".format(inputs['xmax']))
        f.write("ymin = {:f}     !! Minimum Y value [cm]".format(inputs['ymin']))
        f.write("ymax = {:f}     !! Maximum Y value [cm]".format(inputs['ymax']))
        f.write("zmin = {:f}     !! Minimum Z value [cm]".format(inputs['zmin']))
        f.write("zmax = {:f}     !! Maximum Z value [cm]".format(inputs['zmax']))
        f.write("!! Tait-Bryan Angles for z-y`-x`` rotation")
        f.write("alpha = {:f}     !! Rotation about z-axis [rad]".format(inputs['alpha']))
        f.write("beta  = {:f}     !! Rotation about y`-axis [rad]".format(inputs['beta']))
        f.write("gamma = {:f}     !! Rotation about x``-axis [rad]".format(inputs['gamma']))
        f.write("!! Beam Grid origin in machine coordinates (cartesian)")
        f.write("origin(1) = {:f}     !! U value [cm]".format(inputs['origin'][0]))
        f.write("origin(2) = {:f}     !! V value [cm]".format(inputs['origin'][1]))
        f.write("origin(3) = {:f}     !! W value [cm]".format(inputs['origin'][2]))
        f.write("")

        f.write("!! Wavelength Grid Settings")
        f.write("nlambda = {:d}    !! Number of Wavelengths".format(inputs['nlambda']))
        f.write("lambdamin = {:f}    !! Minimum Wavelength [nm]".format(inputs['lambdamin']))
        f.write("lambdamax = {:f}    !! Maximum Wavelength [nm]".format(inputs['lambdamax']))
        f.write("")

        f.write("!! Weight Function Settings")
        f.write("ne_wght = {:d}    !! Number of Energies for Weights".format(inputs['ne_wght']))
        f.write("np_wght = {:d}    !! Number of Pitches for Weights".format(inputs['np_wght']))
        f.write("nphi_wght = {:d}    !! Number of Gyro-angles for Weights".format(inputs['nphi_wght']))
        f.write("emax_wght = {:f}    !! Maximum Energy for Weights [keV]".format(inputs['emax_wght']))
        f.write("nlambda_wght = {:d}    !! Number of Wavelengths for Weights ".format(inputs['nlambda_wght']))
        f.write("lambdamin_wght = {:f}    !! Minimum Wavelength for Weights [nm]".format(inputs['lambdamin_wght']))
        f.write("lambdamax_wght = {:f}    !! Maximum Wavelength for Weights [nm]".format(inputs['lambdamax_wght']))
        f.write("")
        f.write("/")
        f.write("")

#    f.close()

    success("Namelist file created: {}".format(filename))
