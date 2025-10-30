"""
Visualization service for handling FIDASIM output data with xarray
"""

import os
import json
import numpy as np
import xarray as xr
import h5py
from typing import Dict, List, Optional, Any, Union, Tuple
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

class VisualizationService:
    """Service for handling data visualization using xarray"""

    def __init__(self):
        self.cache = {}  # Simple in-memory cache for datasets

    def load_dataset(self, file_path: str) -> xr.Dataset:
        """
        Load a FIDASIM output file as an xarray Dataset

        Args:
            file_path: Path to the HDF5/NetCDF file

        Returns:
            xarray.Dataset: The loaded dataset
        """
        if file_path in self.cache:
            return self.cache[file_path]

        try:
            # Try loading as NetCDF first (xarray native)
            ds = xr.open_dataset(file_path, engine='netcdf4')
        except:
            # Fall back to HDF5 with custom loading
            ds = self._load_hdf5_as_xarray(file_path)

        # Cache the dataset
        self.cache[file_path] = ds
        return ds

    def _load_hdf5_as_xarray(self, file_path: str) -> xr.Dataset:
        """
        Load HDF5 file as xarray Dataset with proper dimension handling
        """
        data_vars = {}
        coords = {}
        attrs = {}

        with h5py.File(file_path, 'r') as f:
            # Get global attributes
            for attr_name in f.attrs:
                attrs[attr_name] = f.attrs[attr_name]

            # Process datasets and groups
            for key in f.keys():
                item = f[key]

                if isinstance(item, h5py.Dataset):
                    # Handle scalar values
                    if item.shape == ():
                        attrs[key] = item[()]
                    else:
                        # Get data and attributes
                        data = item[:]
                        var_attrs = {k: v for k, v in item.attrs.items()}

                        # Check if this is a coordinate variable
                        # (1D array that could serve as a dimension)
                        if len(item.shape) == 1:
                            # Check for common coordinate names
                            coord_names = ['x', 'y', 'z', 'r', 'theta', 'phi',
                                         'lambda', 'time', 'energy', 'level',
                                         'radius', 'channel']
                            if any(name in key.lower() for name in coord_names):
                                coords[key] = (key, data, var_attrs)
                            else:
                                # Create a data variable with auto-generated dims
                                data_vars[key] = (['dim_0'], data, var_attrs)
                        else:
                            # Multi-dimensional array
                            # Try to get dimension labels if available
                            dim_labels = None
                            if 'DIMENSION_LABELS' in item.attrs:
                                dim_labels = [str(label) for label in item.attrs['DIMENSION_LABELS']]

                            if dim_labels:
                                dims = dim_labels
                            else:
                                # Generate dimension names
                                dims = [f'{key}_dim_{i}' for i in range(len(item.shape))]

                            data_vars[key] = (dims, data, var_attrs)

                            # Create coordinate arrays if they don't exist
                            for i, (dim, size) in enumerate(zip(dims, item.shape)):
                                if dim not in coords:
                                    coords[dim] = (dim, np.arange(size), {'units': 'index'})

                elif isinstance(item, h5py.Group):
                    # Process group recursively
                    group_data = self._process_hdf5_group(item, f'{key}/')
                    data_vars.update(group_data['data_vars'])
                    coords.update(group_data['coords'])

        # Create the dataset
        ds = xr.Dataset(data_vars=data_vars, coords=coords, attrs=attrs)
        return ds

    def _process_hdf5_group(self, group: h5py.Group, prefix: str = '') -> Dict:
        """Process HDF5 group recursively"""
        data_vars = {}
        coords = {}

        for key in group.keys():
            item = group[key]
            full_key = f'{prefix}{key}'

            if isinstance(item, h5py.Dataset):
                data = item[:]
                var_attrs = {k: v for k, v in item.attrs.items()}

                if len(item.shape) == 1:
                    # Could be a coordinate
                    if 'level' in key or 'channel' in key or any(
                        name in key.lower() for name in ['x', 'y', 'z', 'r']
                    ):
                        coords[full_key] = (full_key, data, var_attrs)
                    else:
                        data_vars[full_key] = ([f'{full_key}_dim'], data, var_attrs)
                else:
                    # Multi-dimensional data
                    dims = [f'{full_key}_dim_{i}' for i in range(len(item.shape))]

                    # Check for dimension labels
                    if 'DIMENSION_LABELS' in item.attrs:
                        try:
                            labels = [str(label) for label in item.attrs['DIMENSION_LABELS']]
                            if len(labels) == len(dims):
                                dims = labels
                        except:
                            pass

                    data_vars[full_key] = (dims, data, var_attrs)

                    # Auto-create coordinates
                    for i, (dim, size) in enumerate(zip(dims, item.shape)):
                        if dim not in coords:
                            coords[dim] = (dim, np.arange(size), {'units': 'index'})

            elif isinstance(item, h5py.Group):
                # Recursive processing
                sub_group = self._process_hdf5_group(item, f'{full_key}/')
                data_vars.update(sub_group['data_vars'])
                coords.update(sub_group['coords'])

        return {'data_vars': data_vars, 'coords': coords}

    def get_dataset_info(self, file_path: str) -> Dict[str, Any]:
        """
        Get metadata about a dataset

        Returns:
            Dictionary containing dataset information
        """
        ds = self.load_dataset(file_path)

        info = {
            'dimensions': {name: int(size) for name, size in ds.dims.items()},
            'coordinates': {
                name: {
                    'dims': list(coord.dims),
                    'shape': list(coord.shape),
                    'dtype': str(coord.dtype),
                    'attrs': dict(coord.attrs)
                }
                for name, coord in ds.coords.items()
            },
            'data_variables': {
                name: {
                    'dims': list(var.dims),
                    'shape': list(var.shape),
                    'dtype': str(var.dtype),
                    'attrs': dict(var.attrs)
                }
                for name, var in ds.data_vars.items()
            },
            'attributes': dict(ds.attrs)
        }

        return info

    def get_variable_data(
        self,
        file_path: str,
        variable: str,
        slices: Optional[Dict[str, Union[int, slice]]] = None
    ) -> Dict[str, Any]:
        """
        Get data for a specific variable with optional slicing

        Args:
            file_path: Path to the dataset file
            variable: Name of the variable to extract
            slices: Dictionary of dimension slices

        Returns:
            Dictionary containing the data and metadata
        """
        ds = self.load_dataset(file_path)

        if variable not in ds.data_vars:
            raise ValueError(f"Variable '{variable}' not found in dataset")

        var = ds[variable]

        # Apply slicing if provided
        if slices:
            slice_dict = {}
            for dim, slice_val in slices.items():
                if isinstance(slice_val, dict):
                    # Handle slice object passed as dict
                    start = slice_val.get('start')
                    stop = slice_val.get('stop')
                    step = slice_val.get('step')
                    slice_dict[dim] = slice(start, stop, step)
                elif isinstance(slice_val, (int, np.integer)):
                    slice_dict[dim] = slice_val
                else:
                    slice_dict[dim] = slice_val

            var = var.isel(slice_dict)

        # Convert to numpy array for JSON serialization
        data = var.values

        # Handle complex data types
        if np.iscomplexobj(data):
            data = {
                'real': np.real(data).tolist(),
                'imag': np.imag(data).tolist()
            }
        else:
            # Convert to list for JSON serialization
            data = data.tolist()

        # Get coordinate values for the sliced dimensions
        coords = {}
        for coord_name in var.coords:
            coord_data = var.coords[coord_name].values
            if np.iscomplexobj(coord_data):
                coords[coord_name] = {
                    'real': np.real(coord_data).tolist(),
                    'imag': np.imag(coord_data).tolist()
                }
            else:
                coords[coord_name] = coord_data.tolist()

        return {
            'variable': variable,
            'dims': list(var.dims),
            'shape': list(var.shape),
            'data': data,
            'coords': coords,
            'attrs': dict(var.attrs)
        }

    def get_slice_for_plot(
        self,
        file_path: str,
        variable: str,
        x_dim: str,
        y_dim: Optional[str] = None,
        fixed_coords: Optional[Dict[str, int]] = None,
        integrate_dims: Optional[List[str]] = None
    ) -> Dict[str, Any]:
        """
        Get a data slice formatted for plotting with optional integration

        Args:
            file_path: Path to the dataset
            variable: Variable name to plot
            x_dim: X-axis dimension
            y_dim: Y-axis dimension (for 2D plots)
            fixed_coords: Fixed coordinates for other dimensions (for slicing)
            integrate_dims: Dimensions to integrate over instead of slicing

        Returns:
            Dictionary with plot-ready data
        """
        ds = self.load_dataset(file_path)

        if variable not in ds.data_vars:
            raise ValueError(f"Variable '{variable}' not found in dataset")

        var = ds[variable]

        # Apply integration first
        if integrate_dims:
            for dim in integrate_dims:
                if dim in var.dims:
                    # Integrate over the dimension (sum for discrete data)
                    var = var.sum(dim=dim)
                    logger.info(f"Integrated over dimension: {dim}")

        # Apply fixed coordinates (slicing)
        if fixed_coords:
            # Only apply slicing for dimensions that weren't integrated
            coords_to_apply = {k: v for k, v in fixed_coords.items()
                             if k not in (integrate_dims or [])}
            if coords_to_apply:
                var = var.isel(coords_to_apply)

        # Prepare plot data based on dimensionality
        if y_dim is None:
            # 1D plot
            if x_dim not in var.dims:
                raise ValueError(f"Dimension '{x_dim}' not found in variable")

            # Ensure we have a 1D array
            if len(var.dims) > 1:
                # Need to fix all other dimensions
                other_dims = [d for d in var.dims if d != x_dim]
                if not all(d in fixed_coords for d in other_dims):
                    raise ValueError(f"Must specify fixed coordinates for dimensions: {other_dims}")

            x_data = var.coords[x_dim].values
            y_data = var.values

            return {
                'type': '1D',
                'x': x_data.tolist(),
                'y': y_data.tolist(),
                'x_label': x_dim,
                'y_label': variable,
                'x_units': var.coords[x_dim].attrs.get('units', ''),
                'y_units': var.attrs.get('units', '')
            }
        else:
            # 2D plot (heatmap)
            if x_dim not in var.dims or y_dim not in var.dims:
                raise ValueError(f"Dimensions '{x_dim}' or '{y_dim}' not found in variable")

            # Fix all other dimensions
            other_dims = [d for d in var.dims if d not in [x_dim, y_dim]]
            if other_dims and not all(d in fixed_coords for d in other_dims):
                raise ValueError(f"Must specify fixed coordinates for dimensions: {other_dims}")

            # Transpose to ensure correct orientation
            var = var.transpose(y_dim, x_dim, ...)

            x_data = var.coords[x_dim].values
            y_data = var.coords[y_dim].values
            z_data = var.values

            return {
                'type': '2D',
                'x': x_data.tolist(),
                'y': y_data.tolist(),
                'z': z_data.tolist(),
                'x_label': x_dim,
                'y_label': y_dim,
                'z_label': variable,
                'x_units': var.coords[x_dim].attrs.get('units', ''),
                'y_units': var.coords[y_dim].attrs.get('units', ''),
                'z_units': var.attrs.get('units', '')
            }

    def list_available_files(self, core_job_id: str) -> List[Dict[str, str]]:
        """
        List available output files for a job

        Args:
            core_job_id: Core job ID to list files for

        Returns:
            List of available files with metadata
        """
        # Define the output directory based on core_job_id
        base_dir = Path('/data') / f'job_{core_job_id}'

        if not base_dir.exists():
            return []

        files = []
        # Look for HDF5 and NetCDF files
        for pattern in ['*.h5', '*.hdf5', '*.nc', '*.cdf']:
            for file_path in base_dir.glob(pattern):
                files.append({
                    'name': file_path.name,
                    'path': str(file_path),
                    'size': file_path.stat().st_size,
                    'modified': file_path.stat().st_mtime
                })

        return files

    def clear_cache(self, file_path: Optional[str] = None):
        """Clear cached datasets"""
        if file_path:
            self.cache.pop(file_path, None)
        else:
            self.cache.clear()

# Singleton instance
visualization_service = VisualizationService()