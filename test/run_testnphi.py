#!/bin/sh
"exec" "$FIDASIM_DIR/deps/python" "$0" "$@"

import argparse
import numpy as np
from scipy.interpolate import interp1d
import fidasim as fs

def test_beam(beta=0.0):
    uvw_src = np.array([0.0, -530.0 - 2.0*np.cos(beta), 2.0*np.sin(beta)])
    uvw_pos = np.array([0.0, -530.0, 0.0])
    uvw_axis = uvw_pos - uvw_src
    uvw_axis = uvw_axis/np.linalg.norm(uvw_axis)

    focy = 999999.9e0
    focz = 1000e0

    divy = np.full(3,8.73e-3)
    divz = np.full(3,2.27e-2)

    widy = 6.0
    widz = 24.0

    naperture = 1
    ashape = np.array([1])
    awidy = np.array([8.85])
    awidz = np.array([24.0])
    aoffy = np.array([0.0])
    aoffz = np.array([0.0])
    adist = np.array([186.1])

    nbi = {"name":"test_beam","shape":1,"data_source":'run_tests:test_beam',
           "src":uvw_src, "axis":uvw_axis, "widy":widy, "widz":widz,
           "divy":divy, "divz":divz, "focy":focy, "focz":focz,
           "naperture":naperture, "ashape":ashape, "adist":adist,
           "awidy":awidy, "awidz":awidz, "aoffy":aoffy, "aoffz":aoffz}

    return nbi

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

    chords = {"nchan":3, "system":"SPECTRAL","id":id, "data_source":"run_tests:test_chords",
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

    npa_chords = {"nchan":nchan, "system":"NPA", "data_source":"run_tests.py:test_npa",
                  "id":id, "a_shape":np.repeat(2,nchan), "d_shape":np.repeat(2,nchan),
                  "a_cent":a_cent, "a_redge":a_redge, "a_tedge":a_tedge,
                  "d_cent":d_cent, "d_redge":d_redge, "d_tedge":d_tedge,"radius":radius}

    return npa_chords

def test_nc():
    """
    Generates test data for Neutron Collimator geometry.

    Returns:
        nc_chords: A dictionary representing the Neutron Collimator geometry.
    """

    # Chords
    nchan = 47
    ulens = np.arange(nchan) * 3 + 101
    vlens = np.zeros(nchan)
    wlens = np.full(nchan, -100.0)
    lens = np.vstack((ulens, vlens, wlens))

    ulos = np.arange(nchan) * 3 + 101
    vlos = np.zeros(nchan)
    wlos = np.zeros(nchan)
    radius = np.sqrt(ulos**2 + vlos**2)
    id = np.array([b"c" + str(i).encode('utf-8') for i in range(nchan)])  # Generate IDs

    # Initialize arrays for aperture and detector data
    a_cent = np.zeros((3, nchan))
    a_redge = np.zeros((3, nchan))
    a_tedge = np.zeros((3, nchan))
    d_cent = np.zeros((3, nchan))
    d_redge = np.zeros((3, nchan))
    d_tedge = np.zeros((3, nchan))

    # Define base points for aperture and detector
    ac = np.array([0.0, 0.0, 0.0])
    ar = np.array([0.0, 3.0, 0.0])
    at = np.array([0.0, 0.0, 3.0])
    dc = np.array([-50.0, 0.0, 0.0])
    dr = np.array([-50.0, 3.0, 0.0])
    dt = np.array([-50.0, 0.0, 3.0])

    # Loop through each channel and calculate positions
    for i in range(nchan):
        r0 = np.array([ulens[i], vlens[i], wlens[i]])
        rf = np.array([ulos[i], vlos[i], wlos[i]])

        R = fs.utils.line_basis(r0,rf-r0)

        # Calculate aperture and detector positions using rotation and translation
        a_cent[:, i] = np.dot(R, ac) + r0
        a_redge[:, i] = np.dot(R, ar) + r0
        a_tedge[:, i] = np.dot(R, at) + r0

        d_cent[:, i] = np.dot(R, dc) + r0
        d_redge[:, i] = np.dot(R, dr) + r0
        d_tedge[:, i] = np.dot(R, dt) + r0

    # Create the Neutron Collimator geometry dictionary
    nc_chords = {
        'nchan': nchan,
        'system': "NC",
        'data_source': "run_tests.py:test_nc",
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

    te = interp1d(rho, te_rho,fill_value='extrapolate')(rhogrid)
    te = np.where(te > 0.0, te, 0.0).astype('float64')

    ti = interp1d(rho, ti_rho,fill_value='extrapolate')(rhogrid)
    ti = np.where(ti > 0.0, ti, 0.0).astype('float64')

    zeff = interp1d(rho, zeff_rho,fill_value=0.0,bounds_error=False)(rhogrid)
    zeff = np.where(zeff >= 1.0, zeff, 1.0).astype('float64')

    omega = interp1d(rho, omega_rho,fill_value='extrapolate')(rhogrid)

    denimp = dene*(zeff-1.0)/impurity_charge**2.0
    denimp = np.where(denimp > 0.0, denimp, 0.0).astype('float64')

    deni = (dene*(1.0-(zeff-1.0)/impurity_charge))
    deni = np.where(deni > 0.0, deni, 0.0).astype('float64').reshape(1,grid['nr'],grid['nz'])

    vt = grid["r2d"]*omega
    vr = 0*vt
    vz = 0*vt

    # Create mask for valid plasma region - needs to be 2D (nr, nz)
    max_rho = np.max(np.abs(rho))
    mask = np.where(rhogrid <= max_rho, np.int64(1), np.int64(0))

    # For denn, use shape (nr, nz) not (3, nr, nz)
    denn = zeff*0 + 1e8

    profiles = {"time":1.0, "data_source":filename, "mask":mask,
                "deni":deni,"denimp":denimp,"species_mass":species_mass,
                "nthermal":nthermal,"impurity_charge":impurity_charge,
                "te":te, "ti":ti, "vr":vr, "vt":vt, "vz":vz,
                "dene":dene, "zeff":zeff, "denn":denn}

    return profiles

def run_test(args):
    import os
    output_dir = args.path.rstrip('/')+'/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    fida_dir = fs.utils.get_fidasim_dir()
    test_dir = fida_dir + '/test'

    einj = 72.5 #keV
    pinj = 1.7  #MW

    cgfitf = np.array([-0.109171,0.0144685,-7.83224e-5])
    cgfith = np.array([0.0841037,0.00255160,-7.42683e-8])
    ffracs = cgfitf[0] + cgfitf[1]*einj + cgfitf[2]*einj**2
    hfracs = cgfith[0] + cgfith[1]*einj + cgfith[2]*einj**2
    tfracs = 1.0-ffracs-hfracs
    current_fractions = np.array([ffracs,hfracs,tfracs])

    basic_inputs = {"device":"test", "shot":1, "time":1.0,
                    "einj":einj, "pinj":pinj, "current_fractions":current_fractions,
                    "ab":2.01410178e0,
                    "lambdamin":647.0, "lambdamax":667.0, "nlambda":2000,
                    "n_fida":5000000, "n_npa":5000000, "n_nbi":50000,
                    "n_pfida":50000000, "n_pnpa":50000000,
                    "n_halo":500000, "n_dcx":500000, "n_birth":10000,
                    "ne_wght":50, "np_wght":50,"nphi_wght":100,"emax_wght":100e0,
                    "nlambda_wght":1000,"lambdamin_wght":647e0,"lambdamax_wght":667e0,
                    "ne_nc":40,"np_nc":40,"emax_nc_wght":100e0,
                    "nenergy_nc":100,"emin_nc":1000e0,"emax_nc":3000e0,
                    "calc_npa":2, "calc_brems":1,"calc_fida":1,"calc_neutron":1,
                    "calc_cfpd":0, "calc_neut_spec":1, "calc_res":0,
                    "calc_bes":1, "calc_dcx":1, "calc_halo":1, "calc_cold":1,
                    "calc_birth":1, "calc_fida_wght":1,"calc_npa_wght":1,"calc_nc_wght":1,
                    "calc_pfida":1, "calc_pnpa":2,
                    "result_dir":output_dir, "tables_file":fida_dir+'/tables/atomic_tables.h5'}

    basic_bgrid = {"nx":50, "ny":60, "nz":70,
                   "xmin":-50.0, "xmax":50.0,
                   "ymin":-230.0,"ymax":-110.0,
                   "zmin":-70.0, "zmax":70.0,
                   "alpha":0.0, "beta":0.0, "gamma":0.0,
                   "origin":np.zeros(3)}

    inputs = basic_inputs.copy()
    inputs.update(basic_bgrid)
    inputs["comment"] = "Testing passive grid nphi=50 without phi array"
    inputs["runid"] = args.runid

    nbi = test_beam()
    spec = test_chords()
    npa = test_npa()
    nc = test_nc()

    # Create standard 2D grid for now
    grid = fs.utils.rz_grid(100.0, 240.0, 70, -100.0, 100.0, 100)

    # MODIFIED: Manually add nphi to grid without phi array for passive grid only
    # This will be written to equilibrium file and detected by Fortran
    grid['nphi'] = 50

    print("\n=== Testing passive grid nphi functionality ===")
    print("Grid keys:", grid.keys())
    print("Grid has nphi:", 'nphi' in grid)
    print("Grid has phi array:", 'phi' in grid)
    if 'nphi' in grid:
        print("Grid nphi value:", grid['nphi'])
    print("This should create an equilibrium file with nphi=50 but no phi array")
    print("Fortran will detect this and use nphi=50 for passive grid only\n")

    equil, rho, btipsign = fs.utils.read_geqdsk(test_dir+'/g000001.01000', grid, ccw_phi=True, exp_Bp=0)
    fbm = fs.utils.read_nubeam(test_dir+'/test_fi_1.cdf', grid, btipsign = btipsign)

    plasma = test_profiles(test_dir+'/test_profiles.cdf',grid,rho)
    plasma['deni'] = plasma['deni'] - fbm['denf'].reshape(1,grid['nr'],grid['nz'])
    plasma['deni'] = np.where(plasma['deni'] > 0.0, plasma['deni'], 0.0).astype('float64')

    # No need to expand arrays - preprocessing now correctly handles nphi without phi as 2D
    fs.prefida(inputs, grid, nbi, plasma, equil, fbm, spec=spec, npa=npa, nc=nc)

    return 0

def main():
    parser = argparse.ArgumentParser(description="Creates a FIDASIM test case with passive grid nphi")

    parser.add_argument('path', help='Result directory')
    parser.add_argument('-r', '--runid', default = 'test_nphi', help='Test runid')

    args = parser.parse_args()

    run_test(args)

if __name__=='__main__':
    main()