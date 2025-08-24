#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Unit tests for fidasim.preprocessing module
"""

import sys
import os
import numpy as np
import pytest
from unittest.mock import patch, MagicMock, mock_open
import tempfile
import h5py

# Add FIDASIM lib/python to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../lib/python'))

import fidasim.preprocessing as prep
import fidasim.utils as utils


class TestCheckDictSchema:
    """Test check_dict_schema function"""

    def test_valid_dict(self):
        """Test with valid dictionary"""
        dic = {
            'a': 0,
            'b': np.array([1.0, 2.0]),
            'c': "example"
        }

        schema = {
            'a': {'dims': 0, 'type': [int]},
            'b': {'dims': [2], 'type': [float, np.float64]},
            'c': {'dims': 0, 'type': [str]}
        }

        err = prep.check_dict_schema(schema, dic, desc="Test dict")
        assert err == False

    def test_missing_key(self, capsys):
        """Test with missing required key"""
        dic = {'a': 0}
        schema = {
            'a': {'dims': 0, 'type': [int]},
            'b': {'dims': [2], 'type': [float, np.float64]}
        }

        err = prep.check_dict_schema(schema, dic, desc="Test dict")
        assert err == True

        captured = capsys.readouterr()
        assert 'is missing from' in captured.out

    def test_wrong_type(self, capsys):
        """Test with wrong data type"""
        dic = {'a': 1.5}  # float instead of int
        schema = {'a': {'dims': 0, 'type': [int]}}

        err = prep.check_dict_schema(schema, dic, desc="Test dict")
        assert err == True

        captured = capsys.readouterr()
        assert 'has the wrong type' in captured.out

    def test_wrong_shape(self, capsys):
        """Test with wrong array shape"""
        dic = {'a': np.array([1, 2, 3])}  # Shape (3,) instead of (2,)
        schema = {'a': {'dims': [2], 'type': [np.int64, int]}}

        err = prep.check_dict_schema(schema, dic, desc="Test dict")
        assert err == True

        captured = capsys.readouterr()
        assert 'has the wrong shape' in captured.out

    def test_extra_key(self, capsys):
        """Test with extra key (should note but not error)"""
        dic = {'a': 0, 'extra': 'value'}
        schema = {'a': {'dims': 0, 'type': [int]}}

        err = prep.check_dict_schema(schema, dic, desc="Test dict")
        assert err == False  # Extra keys don't cause errors

        captured = capsys.readouterr()
        assert 'Extra variable' in captured.out

    def test_nan_detection(self, capsys):
        """Test NaN detection in arrays"""
        dic = {'a': np.array([1.0, np.nan, 3.0])}
        schema = {'a': {'dims': [3], 'type': [np.float64, float]}}

        err = prep.check_dict_schema(schema, dic, desc="Test dict")
        assert err == True

        captured = capsys.readouterr()
        assert 'NaN or Infinity detected' in captured.out

    def test_inf_detection(self, capsys):
        """Test Infinity detection in arrays"""
        dic = {'a': np.array([1.0, np.inf, 3.0])}
        schema = {'a': {'dims': [3], 'type': [np.float64, float]}}

        err = prep.check_dict_schema(schema, dic, desc="Test dict")
        assert err == True

        captured = capsys.readouterr()
        assert 'NaN or Infinity detected' in captured.out


class TestCheckInputs:
    """Test check_inputs function"""

    @patch('fidasim.preprocessing.os.path.isfile')
    @patch('fidasim.preprocessing.os.path.isdir')
    def test_check_inputs_basic(self, mock_isdir, mock_isfile):
        """Test basic input checking"""
        mock_isfile.return_value = True
        mock_isdir.return_value = True

        inputs = {
            'comment': 'Test run',
            'device': 'TEST',
            'shot': 1,
            'time': 1.0,
            'runid': 'test001',
            'result_dir': '/tmp/results',
            'tables_file': '/tmp/tables.h5',

            # Beam parameters
            'ab': 2.014,
            'pinj': 1.7,
            'einj': 75.0,
            'current_fractions': np.array([0.7, 0.2, 0.1]),

            # Grid parameters
            'nx': 50,
            'ny': 60,
            'nz': 70,
            'xmin': -50.0,
            'xmax': 50.0,
            'ymin': -100.0,
            'ymax': 100.0,
            'zmin': -50.0,
            'zmax': 50.0,
            'origin': np.array([0.0, 0.0, 0.0]),
            'alpha': 0.0,
            'beta': 0.0,
            'gamma': 0.0,

            # Wavelength parameters
            'lambdamin': 650.0,
            'lambdamax': 660.0,
            'nlambda': 1000,

            # Monte Carlo settings
            'n_nbi': 50000,
            'n_halo': 500000,
            'n_dcx': 500000,
            'n_birth': 10000,
            'n_fida': 5000000,
            'n_pfida': 50000000,
            'n_npa': 5000000,
            'n_pnpa': 50000000,

            # Weight function settings
            'ne_wght': 50,
            'np_wght': 50,
            'nphi_wght': 100,
            'emax_wght': 100.0,
            'nlambda_wght': 1000,
            'lambdamin_wght': 650.0,
            'lambdamax_wght': 660.0,

            # Calculation flags
            'calc_spec': 1,
            'calc_beam': 1,
            'calc_brems': 1,
            'calc_fida': 1,
            'calc_npa': 1,
            'calc_neutron': 0,
            'calc_cfpd': 0,
            'calc_birth': 1,
            'calc_fida_wght': 0,
            'calc_npa_wght': 0,
            'calc_res': 0,
            'calc_bes': 0,
            'calc_dcx': 0,
            'calc_halo': 0,
            'calc_cold': 0,
            'calc_pfida': 0,
            'calc_pnpa': 0,

            # Distribution function settings
            'dist_type': 1,
            'impurity_charge': 6
        }

        # Test that valid inputs pass
        result = prep.check_inputs(inputs, use_abs_path=False)
        assert result is not None
        assert result['comment'] == 'Test run'

    def test_check_inputs_missing_required(self):
        """Test with missing required fields"""
        inputs = {'comment': 'Incomplete'}

        with pytest.raises(SystemExit):
            prep.check_inputs(inputs)


class TestCheckPlasma:
    """Test check_plasma function"""

    def test_check_plasma_valid(self):
        """Test valid plasma dictionary"""
        nr, nz = 10, 20
        plasma = {
            'time': 1.0,
            'data_source': 'test_file.cdf',
            'mask': np.ones((nr, nz), dtype=np.int64),
            'te': np.ones((nr, nz)) * 2.0,  # 2 keV
            'ti': np.ones((nr, nz)) * 2.5,  # 2.5 keV
            'dene': np.ones((nr, nz)) * 1e13,  # 1e19 m^-3
            'deni': np.ones((1, nr, nz)) * 1e13,
            'denimp': np.ones((nr, nz)) * 1e11,
            'denn': np.ones((3, nr, nz)) * 1e8,
            'zeff': np.ones((nr, nz)) * 1.5,
            'vr': np.zeros((nr, nz)),
            'vt': np.zeros((nr, nz)),
            'vz': np.zeros((nr, nz)),
            'species_mass': np.array([2.014]),
            'impurity_charge': 6,
            'nthermal': 1
        }

        result = prep.check_plasma(plasma, (nr, nz))
        assert result is not None
        assert result['time'] == 1.0

    def test_check_plasma_wrong_shape(self):
        """Test plasma with wrong array shapes"""
        nr, nz = 10, 20
        plasma = {
            'time': 1.0,
            'data_source': 'test_file.cdf',
            'mask': np.ones((nr, nz), dtype=np.int64),
            'te': np.ones((nr+1, nz)),  # Wrong shape!
            'ti': np.ones((nr, nz)),
            'dene': np.ones((nr, nz)),
            'deni': np.ones((1, nr, nz)),
            'denimp': np.ones((nr, nz)),
            'denn': np.ones((3, nr, nz)),
            'zeff': np.ones((nr, nz)),
            'vr': np.zeros((nr, nz)),
            'vt': np.zeros((nr, nz)),
            'vz': np.zeros((nr, nz)),
            'species_mass': np.array([2.014]),
            'impurity_charge': 6,
            'nthermal': 1
        }

        with pytest.raises(SystemExit):
            prep.check_plasma(plasma, (nr, nz))


class TestCheckFields:
    """Test check_fields function"""

    def test_check_fields_valid(self):
        """Test valid fields dictionary"""
        nr, nz = 10, 20
        fields = {
            'time': 1.0,
            'data_source': 'test_eqdsk.geq',
            'mask': np.ones((nr, nz), dtype=np.int64),
            'br': np.zeros((nr, nz)),
            'bt': np.ones((nr, nz)) * 2.0,  # 2 Tesla
            'bz': np.zeros((nr, nz)),
            'er': np.zeros((nr, nz)),
            'et': np.zeros((nr, nz)),
            'ez': np.zeros((nr, nz))
        }

        result = prep.check_fields(fields, (nr, nz))
        assert result is not None
        assert result['time'] == 1.0

    def test_check_fields_missing_key(self):
        """Test fields with missing required key"""
        nr, nz = 10, 20
        fields = {
            'time': 1.0,
            'data_source': 'test_eqdsk.geq',
            'mask': np.ones((nr, nz), dtype=np.int64),
            # Missing br, bt, bz, er, et, ez
        }

        with pytest.raises(SystemExit):
            prep.check_fields(fields, (nr, nz))


class TestCheckDistribution:
    """Test check_distribution function"""

    def test_check_fbm_valid(self):
        """Test valid fast-ion distribution"""
        nr, nz = 10, 20
        nenergy, npitch = 15, 25

        fbm = {
            'time': 1.0,
            'data_source': 'test_dist.cdf',
            'type': 1,  # Fast-ion distribution
            'nenergy': nenergy,
            'npitch': npitch,
            'energy': np.logspace(0, 2, nenergy),  # 1-100 keV
            'pitch': np.linspace(-1, 1, npitch),
            'denf': np.random.rand(nr, nz) * 1e12,
            'f': np.random.rand(nenergy, npitch, nr, nz) * 1e10
        }

        result = prep.check_distribution(fbm, (nr, nz))
        assert result is not None
        assert result['type'] == 1

    def test_check_fbm_missing_f(self):
        """Test fast-ion distribution missing f array"""
        nr, nz = 10, 20
        nenergy, npitch = 15, 25

        fbm = {
            'time': 1.0,
            'data_source': 'test_dist.cdf',
            'type': 1,
            'nenergy': nenergy,
            'npitch': npitch,
            'energy': np.logspace(0, 2, nenergy),
            'pitch': np.linspace(-1, 1, npitch),
            'denf': np.random.rand(nr, nz) * 1e12
            # Missing 'f' array
        }

        with pytest.raises(SystemExit):
            prep.check_distribution(fbm, (nr, nz))

    def test_check_fbm_type_2(self):
        """Test Guiding center Monte Carlo distribution"""
        nr, nz = 10, 20
        npart = 100000

        fbm = {
            'time': 1.0,
            'data_source': 'test_mc.cdf',
            'type': 2,  # Guiding center MC
            'nparticle': npart,
            'nclass': 1,
            'r': np.random.uniform(100, 200, npart),
            'z': np.random.uniform(-50, 50, npart),
            'phi_enter': np.random.uniform(0, 2*np.pi, npart),
            'energy': np.random.uniform(10, 100, npart),
            'pitch': np.random.uniform(-1, 1, npart),
            'weight': np.ones(npart),
            'class': np.ones(npart, dtype=np.int32)
        }

        result = prep.check_distribution(fbm, (nr, nz))
        assert result is not None
        assert result['type'] == 2


class TestCheckNBI:
    """Test check_nbi function"""

    def test_check_nbi_valid(self):
        """Test valid NBI dictionary"""
        nbi = {
            'name': 'test_beam',
            'data_source': 'test_geometry.dat',
            'shape': 1,  # Rectangular
            'src': np.array([0.0, -200.0, 0.0]),
            'axis': np.array([0.0, 1.0, 0.0]),
            'widy': 10.0,
            'widz': 20.0,
            'divy': np.array([0.01, 0.01, 0.01]),
            'divz': np.array([0.02, 0.02, 0.02]),
            'focy': 1000.0,
            'focz': 1000.0,
            'naperture': 1,
            'ashape': np.array([1]),
            'awidy': np.array([8.0]),
            'awidz': np.array([15.0]),
            'aoffy': np.array([0.0]),
            'aoffz': np.array([0.0]),
            'adist': np.array([150.0])
        }

        result = prep.check_nbi(nbi)
        assert result is not None
        assert result['name'] == 'test_beam'

    def test_check_nbi_wrong_shape(self):
        """Test NBI with wrong array dimensions"""
        nbi = {
            'name': 'test_beam',
            'data_source': 'test_geometry.dat',
            'shape': 1,
            'src': np.array([0.0, -200.0]),  # Wrong: should be 3D
            'axis': np.array([0.0, 1.0, 0.0]),
            'widy': 10.0,
            'widz': 20.0,
            'divy': np.array([0.01, 0.01, 0.01]),
            'divz': np.array([0.02, 0.02, 0.02]),
            'focy': 1000.0,
            'focz': 1000.0,
            'naperture': 0
        }

        with pytest.raises(SystemExit):
            prep.check_nbi(nbi)


class TestCheckSpec:
    """Test check_spec function"""

    def test_check_spec_valid(self):
        """Test valid spectroscopy chord dictionary"""
        nchan = 3
        spec = {
            'system': 'SPECTRAL',
            'data_source': 'test_chords.dat',
            'nchan': nchan,
            'id': np.array([b'ch1', b'ch2', b'ch3']),
            'lens': np.random.rand(3, nchan),
            'axis': np.random.rand(3, nchan),
            'sigma_pi': np.ones(nchan),
            'spot_size': np.ones(nchan) * 0.5,
            'radius': np.linspace(150, 200, nchan)
        }

        result = prep.check_spec(spec)
        assert result is not None
        assert result['nchan'] == nchan

    def test_check_spec_missing_field(self):
        """Test spec with missing required field"""
        spec = {
            'system': 'SPECTRAL',
            'data_source': 'test_chords.dat',
            'nchan': 3,
            # Missing required fields
        }

        with pytest.raises(SystemExit):
            prep.check_spec(spec)


class TestPrefidaIntegration:
    """Integration tests for prefida function"""

    @patch('fidasim.preprocessing.os.path.isfile')
    @patch('fidasim.preprocessing.os.path.isdir')
    @patch('fidasim.preprocessing.h5py.File')
    def test_prefida_minimal(self, mock_h5file, mock_isdir, mock_isfile):
        """Test prefida with minimal valid inputs"""
        mock_isfile.return_value = True
        mock_isdir.return_value = True
        mock_h5file.return_value.__enter__.return_value = MagicMock()

        # Create minimal valid inputs
        nr, nz = 10, 20

        inputs = {
            'comment': 'Test run',
            'device': 'TEST',
            'shot': 1,
            'time': 1.0,
            'runid': 'test001',
            'result_dir': '/tmp/results',
            'tables_file': '/tmp/tables.h5',
            'ab': 2.014,
            'pinj': 1.7,
            'einj': 75.0,
            'current_fractions': np.array([0.7, 0.2, 0.1]),
            'nx': 50, 'ny': 60, 'nz': 70,
            'xmin': -50.0, 'xmax': 50.0,
            'ymin': -100.0, 'ymax': 100.0,
            'zmin': -50.0, 'zmax': 50.0,
            'origin': np.array([0.0, 0.0, 0.0]),
            'alpha': 0.0, 'beta': 0.0, 'gamma': 0.0,
            'lambdamin': 650.0, 'lambdamax': 660.0, 'nlambda': 1000,
            'n_nbi': 50000, 'n_halo': 500000, 'n_dcx': 500000,
            'n_birth': 10000, 'n_fida': 5000000, 'n_pfida': 50000000,
            'n_npa': 5000000, 'n_pnpa': 50000000,
            'ne_wght': 50, 'np_wght': 50, 'nphi_wght': 100,
            'emax_wght': 100.0, 'nlambda_wght': 1000,
            'lambdamin_wght': 650.0, 'lambdamax_wght': 660.0,
            'calc_spec': 1, 'calc_beam': 1, 'calc_brems': 1,
            'calc_fida': 1, 'calc_npa': 1, 'calc_neutron': 0,
            'calc_cfpd': 0, 'calc_birth': 1, 'calc_fida_wght': 0,
            'calc_npa_wght': 0, 'calc_res': 0, 'calc_bes': 0,
            'calc_dcx': 0, 'calc_halo': 0, 'calc_cold': 0,
            'calc_pfida': 0, 'calc_pnpa': 0,
            'dist_type': 1, 'impurity_charge': 6
        }

        grid = {
            'nr': nr, 'nz': nz, 'nphi': 1,
            'r': np.linspace(100, 200, nr*nz).reshape(nr, nz),
            'z': np.linspace(-50, 50, nr*nz).reshape(nr, nz),
            'phi': np.array([0.0]),
            'r2d': np.linspace(100, 200, nr*nz).reshape(nr, nz),
            'z2d': np.linspace(-50, 50, nr*nz).reshape(nr, nz)
        }

        nbi = {
            'name': 'test_beam',
            'data_source': 'test',
            'shape': 1,
            'src': np.array([0.0, -200.0, 0.0]),
            'axis': np.array([0.0, 1.0, 0.0]),
            'widy': 10.0, 'widz': 20.0,
            'divy': np.array([0.01, 0.01, 0.01]),
            'divz': np.array([0.02, 0.02, 0.02]),
            'focy': 1000.0, 'focz': 1000.0,
            'naperture': 0
        }

        plasma = {
            'time': 1.0,
            'data_source': 'test',
            'mask': np.ones((nr, nz), dtype=np.int64),
            'te': np.ones((nr, nz)) * 2.0,
            'ti': np.ones((nr, nz)) * 2.5,
            'dene': np.ones((nr, nz)) * 1e13,
            'deni': np.ones((1, nr, nz)) * 1e13,
            'denimp': np.ones((nr, nz)) * 1e11,
            'denn': np.ones((3, nr, nz)) * 1e8,
            'zeff': np.ones((nr, nz)) * 1.5,
            'vr': np.zeros((nr, nz)),
            'vt': np.zeros((nr, nz)),
            'vz': np.zeros((nr, nz)),
            'species_mass': np.array([2.014]),
            'impurity_charge': 6,
            'nthermal': 1
        }

        fields = {
            'time': 1.0,
            'data_source': 'test',
            'mask': np.ones((nr, nz), dtype=np.int64),
            'br': np.zeros((nr, nz)),
            'bt': np.ones((nr, nz)) * 2.0,
            'bz': np.zeros((nr, nz)),
            'er': np.zeros((nr, nz)),
            'et': np.zeros((nr, nz)),
            'ez': np.zeros((nr, nz))
        }

        fbm = {
            'time': 1.0,
            'data_source': 'test',
            'type': 1,
            'nenergy': 15,
            'npitch': 25,
            'energy': np.logspace(0, 2, 15),
            'pitch': np.linspace(-1, 1, 25),
            'denf': np.zeros((nr, nz)),
            'f': np.zeros((15, 25, nr, nz))
        }

        # Should not raise any exceptions
        prep.prefida(inputs, grid, nbi, plasma, fields, fbm, use_abs_path=False)
