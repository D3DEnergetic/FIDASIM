#!/bin/sh
"exec" "$FIDASIM_DIR/deps/python" "$0" "$@"

import argparse
import numpy as np
from scipy.interpolate import interp1d
import fidasim as fs

def test_chords():
    ulens = np.zeros(3)
    vlens = np.full(3,-170.0)
    wlens = np.full(3,100.0)
    lens = np.vstack((ulens,vlens,wlens))

    ulos = np.zeros(3)
    vlos = np.array([-200.0,-170.0,-140.0])
    wlos = np.zeros(3)
    los = np.vstack((ulos,vlos,wlos))
    axis = los - lens
    for i in range(3):
        axis[:,i] = axis[:,i]/np.linalg.norm(axis[:,i])

    sigma_pi = np.ones(3)
    spot_size = np.zeros(3)
    radius = np.sqrt(ulos**2 + vlos**2)
    id = np.array([b"f1",b"f2",b"f3"])

    chords = {"nchan":3, "system":"SPECTRAL","id":id, "data_source":"run_nobeamtests:test_chords",
              "lens":lens, "axis":axis, "spot_size":spot_size, "sigma_pi":sigma_pi,
              "radius":radius}

    return chords

def test_npa():
    nchan = 3
    ulens = np.zeros(nchan)
    vlens = np.repeat(-170.0,nchan)
    wlens = np.repeat(100.0,nchan)
    lens = np.vstack((ulens,vlens,wlens))

    ulos = np.zeros(nchan)
    vlos = np.array([-200.0,-170.0,-140.0])
    wlos = np.zeros(nchan)
    radius = np.sqrt(ulos**2 + vlos**2)
    id = np.array([b"c1",b"c2",b"c3"])

    a_cent  = np.zeros((3,nchan))
    a_redge = np.zeros((3,nchan))
    a_tedge = np.zeros((3,nchan))
    d_cent  = np.zeros((3,nchan))
    d_redge = np.zeros((3,nchan))
    d_tedge = np.zeros((3,nchan))

    ac = np.array([0.0, 0.0, 0.0])
    ar = np.array([0.0, 3.0, 0.0])
    at = np.array([0.0, 0.0, 3.0])

    dc = np.array([-50.0, 0.0, 0.0])
    dr = np.array([-50.0, 3.0, 0.0])
    dt = np.array([-50.0, 0.0, 3.0])

    for i in range(nchan):
        r0 = np.array([ulens[i],vlens[i],wlens[i]])
        rf = np.array([ulos[i],vlos[i],wlos[i]])
        R = fs.utils.line_basis(r0,rf-r0)
        a_cent[:,i]  = np.dot(R, ac) + r0
        a_redge[:,i] = np.dot(R, ar) + r0
        a_tedge[:,i] = np.dot(R, at) + r0

        d_cent[:,i]  = np.dot(R, dc) + r0
        d_redge[:,i] = np.dot(R, dr) + r0
        d_tedge[:,i] = np.dot(R, dt) + r0

    npa_chords = {"nchan":nchan, "system":"NPA", "data_source":"run_nobeamtests.py:test_npa",
                  "id":id, "a_shape":np.repeat(2,nchan), "d_shape":np.repeat(2,nchan),
                  "a_cent":a_cent, "a_redge":a_redge, "a_tedge":a_tedge,
                  "d_cent":d_cent, "d_redge":d_redge, "d_tedge":d_tedge,"radius":radius}

    return npa_chords

def test_nc():
    """
    Generates test data for Neutron Collimator geometry.
    NC system is below the plasma, looking upward through plasma and NBI.
    All channels share a common aperture, with different detector positions
    to create a fan. All channels are in a single R,Z plane (phi = 0 degrees).

    Returns:
        nc_chords: A dictionary representing the Neutron Collimator geometry.
    """

    # Chords - fan of sightlines in R,Z plane at phi=0
    nchan = 10

    # All apertures at the same location (shared aperture)
    u_ap = 0.0  # At phi=0 plane
    v_ap = -150.0  # Radial position
    w_ap = -80.0  # Below plasma
    ulens = np.full(nchan, u_ap)
    vlens = np.full(nchan, v_ap)
    wlens = np.full(nchan, w_ap)
    lens = np.vstack((ulens, vlens, wlens))

    # Detectors spread out to create fan of sightlines through shared aperture
    # For a true fan, we need to vary the ANGLE, not maintain parallel sightlines
    # Place all detectors at same distance below aperture, but spread in R-Z
    # All at same phi (u=0)

    det_distance = 100.0  # Distance below aperture (cm)

    # Fan angles: spread vertically from more horizontal to more vertical
    # Angles relative to horizontal in the R-Z plane
    angles = np.linspace(30, 60, nchan)  # degrees from horizontal

    ulos = np.zeros(nchan)
    vlos = np.zeros(nchan)
    wlos = np.zeros(nchan)

    for i in range(nchan):
        angle_rad = np.radians(angles[i])
        # Place detector at fixed distance from aperture, at varying angles
        # Pointing downward and outward (negative v, negative w from aperture)
        vlos[i] = v_ap - det_distance * np.cos(angle_rad)
        wlos[i] = w_ap - det_distance * np.sin(angle_rad)
    radius = np.sqrt(ulos**2 + vlos**2)
    id = np.array([b"nc" + str(i+1).encode('utf-8') for i in range(nchan)])  # Generate IDs

    # Initialize arrays for aperture and detector data
    a_cent = np.zeros((3, nchan))
    a_redge = np.zeros((3, nchan))
    a_tedge = np.zeros((3, nchan))
    d_cent = np.zeros((3, nchan))
    d_redge = np.zeros((3, nchan))
    d_tedge = np.zeros((3, nchan))

    # Aperture and detector half-widths (cm)
    ap_width = 3.0
    det_width = 3.0

    # Loop through each channel and calculate positions
    for i in range(nchan):
        # Aperture position (same for all channels)
        a_cent[:, i] = [ulens[i], vlens[i], wlens[i]]

        # Detector position (different for each channel)
        d_cent[:, i] = [ulos[i], vlos[i], wlos[i]]

        # Calculate line-of-sight direction
        los_dir = d_cent[:, i] - a_cent[:, i]
        los_dir = los_dir / np.linalg.norm(los_dir)

        # Create orthogonal basis for aperture/detector planes
        # Since all channels are in phi=0 plane (x=0), we can use standard basis
        # Radial direction (in y-z plane)
        radial = np.array([0.0, vlens[i], wlens[i]])
        radial = radial / np.linalg.norm(radial)

        # Toroidal direction (perpendicular to r,z plane)
        toroidal = np.array([1.0, 0.0, 0.0])

        # Vertical direction (in plane, perpendicular to radial)
        vertical = np.cross(toroidal, radial)

        # Define aperture edges
        a_redge[:, i] = a_cent[:, i] + ap_width * radial
        a_tedge[:, i] = a_cent[:, i] + ap_width * toroidal

        # Define detector edges (same orientation as aperture)
        d_redge[:, i] = d_cent[:, i] + det_width * radial
        d_tedge[:, i] = d_cent[:, i] + det_width * toroidal

    # Create the Neutron Collimator geometry dictionary
    nc_chords = {
        'nchan': nchan,
        'system': "NC",
        'data_source': "run_nobeamtests.py:test_nc",
        'id': id,
        'a_shape': np.repeat(2,nchan),  # All apertures have shape 2
        'd_shape': np.repeat(2,nchan),  # All detectors have shape 2
        'a_cent': a_cent,
        'a_redge': a_redge,
        'a_tedge': a_tedge,
        'd_cent': d_cent,
        'd_redge': d_redge,
        'd_tedge': d_tedge,
        'radius': radius
    }

    return nc_chords

def test_profiles(filename, grid, rhogrid):
    prof = fs.utils.read_ncdf(filename)

    impurity_charge = 6
    nthermal = 1
    species_mass = np.array([2.01410178e0])
    rho = prof["rho"]
    dene_rho = prof["dene"]*1e-6
    ti_rho = prof["ti"]*1e-3
    te_rho = prof["te"]*1e-3
    zeff_rho = prof["zeff"]*1e0
    omega_rho = prof["omega"]*1e0

    dene = interp1d(rho, dene_rho,fill_value='extrapolate')(rhogrid)
    dene = np.where(dene > 0.0, dene, 0.0).astype('float64')

    ti = interp1d(rho, ti_rho,fill_value='extrapolate')(rhogrid)
    ti = np.where(ti > 0.0, ti, 0.0).astype('float64')

    te = interp1d(rho, te_rho,fill_value='extrapolate')(rhogrid)
    ti = np.where(te > 0.0, te, 0.0).astype('float64')

    zeff = interp1d(rho, zeff_rho,fill_value='extrapolate')(rhogrid)
    zeff = np.where(te > 1.0, te, 1.0).astype('float64')

    denimp = dene*(zeff - 1)/(impurity_charge*(impurity_charge-1))
    deni = (dene - impurity_charge*denimp).reshape(1,grid['nr'],grid['nz'])
    deni = np.where(deni > 0.0, deni, 0.0).astype('float64')

    vt = grid["r2d"]*interp1d(rho, omega_rho,fill_value='extrapolate')(rhogrid)
    vr = 0*vt
    vz = 0*vt
    denn = zeff*0 + 1e8
    max_rho = np.max(np.abs(rho))

    mask = np.where(rhogrid <= max_rho, np.int64(1), np.int64(0))

    profiles = {"time":1.0, "data_source":filename, "mask":mask,
                "deni":deni,"denimp":denimp,"species_mass":species_mass,
                "nthermal":nthermal,"impurity_charge":impurity_charge,
                "te":te, "ti":ti, "vr":vr, "vt":vt, "vz":vz,
                "dene":dene, "zeff":zeff, "denn":denn,"profiles":prof}

    return profiles

def run_test(args):
    fida_dir = fs.utils.get_fidasim_dir()
    test_dir = fida_dir + '/test'

    # Basic inputs - PASSIVE ONLY (no beam)
    # No beam parameters needed for passive-only mode

    basic_inputs = {"device":"test", "shot":1, "time":1.0,
                    # No beam parameters - this is a passive-only run
                    "ab":2.01410178e0,  # Fast-ion mass still needed for distribution
                    "lambdamin":647e0, "lambdamax":667e0, "nlambda":2000,
                    # Particle counts for Monte Carlo (only passive needed)
                    "n_pfida":50000000, "n_pnpa":50000000,
                    "n_fida":0, "n_npa":0, "n_nbi":0,
                    "n_halo":0, "n_dcx":0, "n_birth":0,
                    # Weight function parameters (required but not used for passive-only)
                    "ne_wght":50, "np_wght":50,"nphi_wght":100,"emax_wght":100e0,
                    "nlambda_wght":1000,"lambdamin_wght":647e0,"lambdamax_wght":667e0,
                    # Neutron collimator parameters
                    "ne_nc":40,"np_nc":40,"emax_nc_wght":100e0,
                    "nenergy_nc":100,"emin_nc":1000e0,"emax_nc":3000e0,
                    # ALL ACTIVE MEASUREMENTS DISABLED
                    "calc_npa":0, "calc_brems":0,"calc_fida":0,
                    "calc_bes":0, "calc_dcx":0, "calc_halo":0, "calc_cold":0,
                    "calc_birth":0, "calc_fida_wght":0,"calc_npa_wght":0,
                    "calc_cfpd":0, "calc_res":0,"calc_nc_wght":0,
                    # PASSIVE MEASUREMENTS ENABLED
                    "calc_pfida":1, "calc_pnpa":2, "calc_neutron":1, "calc_neut_spec":1,
                    "result_dir":args.path, "tables_file":fida_dir+'/tables/atomic_tables.h5'}

    # No beam grid parameters needed for passive-only mode
    # The preprocessing module now handles this automatically

    inputs = basic_inputs.copy()
    # No beam grid parameters needed
    inputs["comment"] = "Passive-only test: NO BEAM GEOMETRY, passive FIDA/NPA/neutron only"
    inputs["runid"] = args.runid

    # NO BEAM GEOMETRY - test_beam() not called
    # nbi variable is NOT defined - this is the key test!

    spec = test_chords()   # Needed for passive FIDA spectroscopy
    npa = test_npa()       # Needed for passive NPA detector geometry
    nc = test_nc()         # Needed for neutron collimator

    grid = fs.utils.rz_grid(100.0, 240.0, 70, -100.0, 100.0, 100)

    equil, rho, btipsign = fs.utils.read_geqdsk(test_dir+'/g000001.01000', grid, ccw_phi=True, exp_Bp=0)
    fbm = fs.utils.read_nubeam(test_dir+'/test_fi_1.cdf', grid, btipsign = btipsign)

    plasma = test_profiles(test_dir+'/test_profiles.cdf',grid,rho)
    plasma['deni'] = plasma['deni'] - fbm['denf'].reshape(1,grid['nr'],grid['nz'])
    plasma['deni'] = np.where(plasma['deni'] > 0.0, plasma['deni'], 0.0).astype('float64')

    # Call prefida with nbi=None for passive-only mode - this is the critical test
    # Pass None as 3rd argument to maintain backward compatibility
    # geometry.h5 will be created WITHOUT any /nbi group
    fs.prefida(inputs, grid, None, plasma, equil, fbm, spec=spec, npa=npa, nc=nc)

    return 0

def main():
    parser = argparse.ArgumentParser(description="Creates a FIDASIM passive-only test case (no beam)")

    parser.add_argument('path', help='Result directory')
    parser.add_argument('-r', '--runid', default = 'nobeamtests', help='Test runid')

    args = parser.parse_args()

    run_test(args)

if __name__=='__main__':
    main()
