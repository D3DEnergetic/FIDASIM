#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Unit tests for fidasim.utils module
"""

import sys
import os
import numpy as np
import pytest
from unittest.mock import patch, MagicMock
import tempfile
import h5py

# Add FIDASIM lib/python to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../lib/python'))

import fidasim.utils as utils


class TestCOCOS:
    """Test COCOS class functionality"""

    def test_cocos_defaults(self):
        """Test default COCOS initialization"""
        cc = utils.COCOS()
        assert cc.cocos == 5
        assert cc.exp_Bp == 0
        assert cc.sigma_Bp == 1
        assert cc.sigma_RphZ == 1
        assert cc.sigma_rhothph == -1
        assert cc.sign_q == -1
        assert cc.sign_pprime == -1

    def test_cocos_index_1(self):
        """Test COCOS index 1"""
        cc = utils.COCOS(1)
        assert cc.cocos == 1
        assert cc.exp_Bp == 0
        assert cc.sigma_Bp == 1
        assert cc.sigma_RphZ == 1
        assert cc.sigma_rhothph == 1
        assert cc.sign_q == 1
        assert cc.sign_pprime == -1

    def test_cocos_index_11(self):
        """Test COCOS index 11 (exp_Bp = 1)"""
        cc = utils.COCOS(11)
        assert cc.cocos == 11
        assert cc.exp_Bp == 1
        assert cc.sigma_Bp == 1
        assert cc.sigma_RphZ == 1
        assert cc.sigma_rhothph == 1
        assert cc.sign_q == 1
        assert cc.sign_pprime == -1

    def test_cocos_even_index(self):
        """Test even COCOS index (sigma_RphZ = -1)"""
        cc = utils.COCOS(2)
        assert cc.cocos == 2
        assert cc.sigma_RphZ == -1


class TestGetFidasimDir:
    """Test get_fidasim_dir function"""

    def test_get_fidasim_dir(self):
        """Test that get_fidasim_dir returns a valid path"""
        fida_dir = utils.get_fidasim_dir()
        assert isinstance(fida_dir, str)
        assert os.path.exists(fida_dir)
        # Check for expected FIDASIM structure
        assert os.path.exists(os.path.join(fida_dir, 'lib'))
        assert os.path.exists(os.path.join(fida_dir, 'src'))


class TestGetVersion:
    """Test get_version function"""

    @patch('subprocess.Popen')
    @patch('subprocess.check_output')
    def test_get_version_with_git(self, mock_check_output, mock_popen):
        """Test version retrieval with git available"""
        # Mock git command check
        mock_proc = MagicMock()
        mock_proc.communicate.return_value = (b'/usr/bin/git\n', b'')
        mock_popen.return_value = mock_proc

        # Mock git describe output
        mock_check_output.return_value = b'v1.0.0-dirty\n'

        with patch('os.path.isfile', return_value=True):
            with patch('os.path.isdir', return_value=True):
                version = utils.get_version('/fake/path')
                assert version == b'v1.0.0-dirty'

    def test_get_version_fallback(self):
        """Test version retrieval fallback to VERSION file"""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create VERSION file
            version_file = os.path.join(tmpdir, 'VERSION')
            with open(version_file, 'w') as f:
                f.write('v2.0.0')

            # Mock git not available
            with patch('os.path.isfile', return_value=False):
                version = utils.get_version(tmpdir)
                assert version == 'v2.0.0'


class TestCoordinateTransforms:
    """Test coordinate transformation functions"""

    def test_tb_zyx_identity(self):
        """Test Tait-Bryan angles with zero rotation"""
        basis, inv_basis = utils.tb_zyx(0, 0, 0)

        # Should get identity matrix
        np.testing.assert_array_almost_equal(basis, np.eye(3))
        np.testing.assert_array_almost_equal(inv_basis, np.eye(3))

        # Test that basis * inv_basis = identity
        identity = np.dot(basis, inv_basis)
        np.testing.assert_array_almost_equal(identity, np.eye(3))

    def test_tb_zyx_90_deg_z(self):
        """Test 90 degree rotation about z-axis"""
        basis, inv_basis = utils.tb_zyx(np.pi/2, 0, 0)

        # Check orthogonality
        identity = np.dot(basis, inv_basis)
        np.testing.assert_array_almost_equal(identity, np.eye(3))

    def test_uvw_to_xyz_identity(self):
        """Test uvw_to_xyz with identity transformation"""
        uvw = np.array([1.0, 2.0, 3.0])
        xyz = utils.uvw_to_xyz(0, 0, 0, uvw)
        np.testing.assert_array_almost_equal(xyz, uvw)

    def test_xyz_to_uvw_identity(self):
        """Test xyz_to_uvw with identity transformation"""
        xyz = np.array([1.0, 2.0, 3.0])
        uvw = utils.xyz_to_uvw(0, 0, 0, xyz)
        np.testing.assert_array_almost_equal(uvw, xyz)

    def test_coordinate_transform_inverse(self):
        """Test that uvw_to_xyz and xyz_to_uvw are inverses"""
        alpha, beta, gamma = 0.1, 0.2, 0.3
        origin = np.array([10.0, 20.0, 30.0])
        xyz_original = np.array([1.0, 2.0, 3.0])

        # Transform xyz -> uvw -> xyz
        uvw = utils.xyz_to_uvw(alpha, beta, gamma, xyz_original, origin)
        xyz_back = utils.uvw_to_xyz(alpha, beta, gamma, uvw, origin)

        np.testing.assert_array_almost_equal(xyz_original, xyz_back)


class TestLineBasis:
    """Test line_basis function"""

    def test_line_basis_x_axis(self):
        """Test line basis along x-axis"""
        r0 = np.array([0.0, 0.0, 0.0])
        v0 = np.array([1.0, 0.0, 0.0])

        basis, inv_basis = utils.line_basis(r0, v0)

        # First column should be normalized v0
        np.testing.assert_array_almost_equal(basis[:, 0], v0)

        # Basis should be orthonormal
        identity = np.dot(basis, inv_basis)
        np.testing.assert_array_almost_equal(identity, np.eye(3))

    def test_line_basis_diagonal(self):
        """Test line basis with diagonal direction"""
        r0 = np.array([0.0, 0.0, 0.0])
        v0 = np.array([1.0, 1.0, 1.0])

        basis, inv_basis = utils.line_basis(r0, v0)

        # Check orthonormality
        identity = np.dot(basis, inv_basis)
        np.testing.assert_array_almost_equal(identity, np.eye(3))


class TestRZGrid:
    """Test rz_grid function"""

    def test_rz_grid_basic(self):
        """Test basic RZ grid creation"""
        grid = utils.rz_grid(100.0, 200.0, 10, -50.0, 50.0, 20)

        assert grid['nr'] == 10
        assert grid['nz'] == 20
        assert grid['nphi'] == 1
        assert grid['r'].shape == (10, 20)
        assert grid['z'].shape == (10, 20)
        assert grid['phi'].shape == (1,)

        # Check ranges
        assert np.isclose(grid['r'].min(), 100.0, rtol=1e-5)
        assert np.isclose(grid['r'].max(), 200.0, rtol=1e-5)
        assert np.isclose(grid['z'].min(), -50.0, rtol=1e-5)
        assert np.isclose(grid['z'].max(), 50.0, rtol=1e-5)

    def test_rz_grid_with_phi(self):
        """Test RZ grid with phi dimension"""
        grid = utils.rz_grid(100.0, 200.0, 5, -50.0, 50.0, 10, 
                            phimin=0.0, phimax=2*np.pi, nphi=4)

        assert grid['nphi'] == 4
        assert grid['phi'].shape == (4,)
        assert grid['r3d'].shape == (5, 10, 4)
        assert grid['z3d'].shape == (5, 10, 4)
        assert grid['phi3d'].shape == (5, 10, 4)


class TestBeamGrid:
    """Test beam_grid function"""

    def test_beam_grid_basic(self):
        """Test basic beam grid creation"""
        # Create simple NBI dict
        nbi = {
            'src': np.array([0.0, -200.0, 0.0]),
            'axis': np.array([0.0, 1.0, 0.0])
        }

        rstart = 100.0

        grid = utils.beam_grid(nbi, rstart, 
                              nx=10, ny=20, nz=15,
                              length=100.0, width=50.0, height=30.0)

        assert grid['nx'] == 10
        assert grid['ny'] == 20
        assert grid['nz'] == 15
        assert 'xc' in grid
        assert 'yc' in grid
        assert 'zc' in grid
        assert 'origin' in grid
        assert 'basis' in grid
        assert 'inv_basis' in grid


class TestColoredOutput:
    """Test colored output functions"""

    def test_colored(self):
        """Test colored text function"""
        text = utils.colored("Test", "green")
        assert "Test" in text
        assert "\033[" in text  # Check for ANSI color code

    def test_info(self, capsys):
        """Test info message"""
        utils.info("Test info")
        captured = capsys.readouterr()
        assert "INFO" in captured.out
        assert "Test info" in captured.out

    def test_warn(self, capsys):
        """Test warning message"""
        utils.warn("Test warning")
        captured = capsys.readouterr()
        assert "WARNING" in captured.out
        assert "Test warning" in captured.out

    def test_error_no_halt(self, capsys):
        """Test error message without halting"""
        utils.error("Test error", halt=False)
        captured = capsys.readouterr()
        assert "ERROR" in captured.out
        assert "Test error" in captured.out

    def test_success(self, capsys):
        """Test success message"""
        utils.success("Test success")
        captured = capsys.readouterr()
        assert "SUCCESS" in captured.out
        assert "Test success" in captured.out


class TestWriteData:
    """Test write_data function for HDF5"""

    def test_write_data_basic(self):
        """Test basic HDF5 data writing"""
        with tempfile.NamedTemporaryFile(suffix='.h5', delete=False) as tmp:
            tmp_name = tmp.name

        try:
            with h5py.File(tmp_name, 'w') as f:
                data = {
                    'scalar': 1.0,
                    'array': np.array([1, 2, 3]),
                    'matrix': np.array([[1, 2], [3, 4]]),
                    'string': 'test'
                }

                desc = {
                    'scalar': 'A scalar value',
                    'array': 'A 1D array',
                    'matrix': 'A 2D matrix',
                    'string': 'A test string'
                }

                units = {
                    'scalar': 'm',
                    'array': 'cm',
                    'matrix': 'keV'
                }

                utils.write_data(f, data, desc, units, name='test_group')

                # Verify data was written
                assert 'test_group' in f
                assert 'scalar' in f['test_group']
                assert 'array' in f['test_group']
                assert 'matrix' in f['test_group']
                assert 'string' in f['test_group']

                # Check attributes
                assert f['test_group/scalar'].attrs['description'] == b'A scalar value'
                assert f['test_group/scalar'].attrs['units'] == b'm'

        finally:
            os.unlink(tmp_name)


class TestReadNCDF:
    """Test read_ncdf function"""

    def test_read_ncdf_mock(self):
        """Test reading NetCDF file with mock"""
        with patch('fidasim.utils.netcdf.netcdf_file') as mock_netcdf:
            # Create mock NetCDF file
            mock_file = MagicMock()
            mock_file.variables = {
                'time': MagicMock(data=np.array([0.0, 1.0, 2.0])),
                'data': MagicMock(data=np.array([[1, 2], [3, 4], [5, 6]]))
            }
            mock_netcdf.return_value.__enter__.return_value = mock_file

            result = utils.read_ncdf('fake_file.nc')

            assert 'time' in result
            assert 'data' in result
            np.testing.assert_array_equal(result['time'], [0.0, 1.0, 2.0])


class TestAABBIntersect:
    """Test axis-aligned bounding box intersection"""

    def test_aabb_intersect_hit(self):
        """Test AABB intersection with hit"""
        # Ray from origin along x-axis
        r0 = np.array([0.0, 0.0, 0.0])
        d0 = np.array([1.0, 0.0, 0.0])

        # Box centered at (5, 0, 0) with size 2x2x2
        rc = np.array([5.0, 0.0, 0.0])
        dr = np.array([1.0, 1.0, 1.0])

        tmin, tmax = utils.aabb_intersect(rc, dr, r0, d0)

        # Ray should hit the box
        assert tmin < tmax
        assert tmin > 0  # Box is in front of ray origin

    def test_aabb_intersect_miss(self):
        """Test AABB intersection with miss"""
        # Ray from origin along x-axis
        r0 = np.array([0.0, 0.0, 0.0])
        d0 = np.array([1.0, 0.0, 0.0])

        # Box off to the side
        rc = np.array([5.0, 5.0, 0.0])
        dr = np.array([1.0, 1.0, 1.0])

        tmin, tmax = utils.aabb_intersect(rc, dr, r0, d0)

        # Ray should miss the box
        assert tmin > tmax  # Indicates no intersection
