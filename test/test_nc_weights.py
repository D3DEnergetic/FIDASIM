#!/usr/bin/env python
"""
Test script for Neutron Collimator Weight Functions
====================================================

This script demonstrates how to:
1. Run FIDASIM with NC weight functions enabled
2. Read the NC weight function output
3. Analyze and visualize the results
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Add FIDASIM Python modules to path
sys.path.insert(0, os.path.join(os.environ['FIDASIM_DIR'], 'lib/python'))
import fidasim as fs

def test_nc_weights():
    """Run a test case with NC weight functions enabled"""
    
    # Get test directory
    fida_dir = fs.utils.get_fidasim_dir()
    test_dir = os.path.join(fida_dir, 'test')
    result_dir = './test_nc_results'
    
    # Create result directory if it doesn't exist
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    
    # Set up inputs with NC weight functions enabled
    inputs = {
        "device": "test", 
        "shot": 1, 
        "time": 1.0,
        "runid": "test_nc",
        "result_dir": result_dir,
        "tables_file": os.path.join(fida_dir, 'tables/atomic_tables.h5'),
        
        # Beam parameters
        "einj": 72.5,  # keV
        "pinj": 1.7,   # MW
        "ab": 2.01410178e0,  # Deuterium mass
        "current_fractions": [0.7, 0.2, 0.1],
        
        # Grid parameters
        "nx": 50, "ny": 60, "nz": 70,
        "xmin": -50.0, "xmax": 50.0,
        "ymin": -230.0, "ymax": -110.0,
        "zmin": -70.0, "zmax": 70.0,
        "alpha": 0.0, "beta": 0.0, "gamma": 0.0,
        "origin": [0.0, 0.0, 0.0],
        
        # Spectral parameters
        "lambdamin": 647.0, "lambdamax": 667.0, "nlambda": 20,
        
        # Monte Carlo particles
        "n_nbi": 50000,
        "n_fida": 0,  # Skip FIDA for this test
        "n_npa": 0,   # Skip NPA for this test
        
        # Calculation flags
        "calc_neutron": 1,      # Basic neutron rate
        "calc_neut_spec": 1,    # Neutron spectroscopy
        "calc_nc_wght": 2,      # NC weights with emissivity
        "calc_fida": 0,
        "calc_npa": 0,
        "calc_brems": 0,
        "calc_bes": 0,
        
        # NC weight function parameters
        "ne_nc": 30,      # Energy bins
        "np_nc": 20,      # Pitch bins
        "emax_nc": 100.0, # Maximum energy [keV]
        
        # Other settings
        "verbose": 1,
        "seed": 42
    }
    
    # Create test beam
    from run_tests import test_beam
    nbi = test_beam()
    
    # Create NC geometry
    from run_tests import test_nc
    nc = test_nc()
    
    # Create grid
    grid = fs.utils.rz_grid(100.0, 240.0, 70, -100.0, 100.0, 100)
    
    # Read equilibrium
    equil, rho, btipsign = fs.utils.read_geqdsk(
        os.path.join(test_dir, 'g000001.01000'), 
        grid, ccw_phi=True, exp_Bp=0
    )
    
    # Read fast-ion distribution
    fbm = fs.utils.read_nubeam(
        os.path.join(test_dir, 'test_fi_1.cdf'), 
        grid, btipsign=btipsign
    )
    
    # Read plasma profiles
    from run_tests import test_profiles
    plasma = test_profiles(
        os.path.join(test_dir, 'test_profiles.cdf'),
        grid, rho
    )
    
    # Adjust thermal ion density
    plasma['deni'] = plasma['deni'] - fbm['denf'].reshape(1, grid['nr'], grid['nz'])
    plasma['deni'] = np.where(plasma['deni'] > 0.0, plasma['deni'], 0.0).astype('float64')
    
    print("\n=== Running FIDASIM with NC Weight Functions ===")
    print(f"NC Channels: {nc['nchan']}")
    print(f"Energy bins: {inputs['ne_nc']}")
    print(f"Pitch bins: {inputs['np_nc']}")
    print(f"Max energy: {inputs['emax_nc']} keV")
    
    # Run FIDASIM
    fs.prefida(inputs, grid, nbi, plasma, equil, fbm, nc=nc)
    
    # Read and analyze results
    nc_weights_file = os.path.join(result_dir, 'test_nc_nc_weights.h5')
    
    if os.path.exists(nc_weights_file):
        print(f"\n=== Reading NC Weight Functions from {nc_weights_file} ===")
        
        # Use the new read_nc_weights function
        nc_data = fs.utils.read_nc_weights(nc_weights_file)
        
        print(f"Loaded NC weights with shape: {nc_data['weight'].shape}")
        print(f"Energy range: {nc_data['energy'][0]:.1f} - {nc_data['energy'][-1]:.1f} keV")
        print(f"Pitch range: {nc_data['pitch'][0]:.2f} - {nc_data['pitch'][-1]:.2f}")
        
        # Calculate some statistics
        for ch in range(min(3, nc['nchan'])):  # Show first 3 channels
            total_weight = np.sum(nc_data['weight'][:, :, ch])
            peak_energy_idx = np.sum(nc_data['weight'][:, :, ch], axis=1).argmax()
            peak_energy = nc_data['energy'][peak_energy_idx]
            
            print(f"\nChannel {ch}:")
            print(f"  Total weight: {total_weight:.3e}")
            print(f"  Peak energy: {peak_energy:.1f} keV")
            
            if 'flux' in nc_data:
                total_flux = np.sum(nc_data['flux'][:, ch])
                print(f"  Total flux: {total_flux:.3e} neutrons/s")
        
        # Create a simple plot
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # Plot weight function for first channel
        ch = 0
        im = axes[0].pcolormesh(
            nc_data['energy'], nc_data['pitch'], 
            nc_data['weight'][:, :, ch].T,
            shading='auto', cmap='viridis'
        )
        axes[0].set_xlabel('Energy [keV]')
        axes[0].set_ylabel('Pitch')
        axes[0].set_title(f'NC Weight Function - Channel {ch}')
        plt.colorbar(im, ax=axes[0])
        
        # Plot energy spectrum
        if 'flux' in nc_data:
            for ch in range(min(3, nc['nchan'])):
                axes[1].plot(nc_data['energy'], nc_data['flux'][:, ch], 
                           label=f'Channel {ch}', linewidth=2)
            axes[1].set_xlabel('Energy [keV]')
            axes[1].set_ylabel('Flux [neutrons/(s·keV)]')
            axes[1].set_title('Neutron Energy Spectrum')
            axes[1].legend()
            axes[1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(os.path.join(result_dir, 'nc_weights_test.png'), dpi=150)
        print(f"\nPlot saved to {os.path.join(result_dir, 'nc_weights_test.png')}")
        plt.show()
        
    else:
        print(f"\nWarning: NC weights file not found at {nc_weights_file}")
    
    # Also check neutron rate file
    neutron_file = os.path.join(result_dir, 'test_nc_neutrons.h5')
    if os.path.exists(neutron_file):
        import h5py
        with h5py.File(neutron_file, 'r') as f:
            if 'rate' in f:
                rate = f['rate'][()]
                print(f"\nTotal neutron rate: {rate:.3e} neutrons/s")
            if 'flux' in f:
                flux = f['flux'][:]
                print(f"NC flux array shape: {flux.shape}")
                for ch in range(min(3, flux.shape[1])):
                    print(f"  Channel {ch} total flux: {np.sum(flux[:, ch]):.3e}")
    
    return nc_data

if __name__ == '__main__':
    print("Testing Neutron Collimator Weight Functions")
    print("=" * 50)
    
    # Run the test
    nc_data = test_nc_weights()
    
    print("\n" + "=" * 50)
    print("Test completed successfully!")
    print("\nYou can now use the plot_weights script to visualize the results:")
    print("  plot_weights test_nc_results/test_nc_nc_weights.h5")
    print("  plot_weights test_nc_results/test_nc_nc_weights.h5 -a  # All channels")
    print("  plot_weights test_nc_results/test_nc_nc_weights.h5 -ch 0 1 2  # Compare channels")
