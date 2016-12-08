#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


def write_data(h5_obj, dic, desc, units, name=''):
    #+##`write_data(h5_obj, dic, desc, units)`
    #+ Write h5 datasets with attributes 'description' and 'units'
    #+###Arguments
    #+     **h5_obj**: An h5 file or group object from h5py
    #+
    #+     **dic**: Dict of data to save as h5 datasets
    #+
    #+     **desc**: Dict with same keys as dic describing each item in dic
    #+
    #+     **units**: Dict with same keys as dic providing units of data in dic, doesn't have to be all keys of dic.
    #+
    #+###Keyword Arguments
    #+     **name**: Name/description of dic for clarity in raising errors
    #+
    #+###Example Usage
    #+```python
    #+>>>import h5py
    #+>>>
    #+>>> uvw_to_xyz(h5_obj, dic, desc, units)
    for key in dic:
        # Transpose data to match expected by Fortran and historically provided by IDL
        if isinstance(dic[key], np.ndarray):
            if dic[key].ndim == 2:
                dic[key] = dic[key].T
            elif dic[key].ndim != 1:
                raise ValueError('Dict {}, key {}, has shape {}, need fix.'.format(name, key, dic[key].shape))

        # Create dataset
        ds = h5_obj.create_dataset(key, data = dic[key])

        # Add descrption attribute
        ds.attrs['description'] = desc[key]

        # Add units attribute (if present)
        if key in units:
            ds.attrs['units'] = units[key]

