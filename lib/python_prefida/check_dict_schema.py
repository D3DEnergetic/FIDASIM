#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from lib.python_prefida.info import info
from lib.python_prefida.error import error
import numpy as np


def check_dict_schema(schema, dic, desc=None):
    #+#check_dict_schema
    #+ Check dict `dec` is formatted according to `schema`
    #+***
    #+##Input Arguments
    #+     **schema**: dict schema
    #+
    #+     **dic**: dict to check
    #+
    #+##Output Arguments
    #+     **err**: error code
    #+
    #+##Keyword Arguments
    #+     **desc**: description of dict `dic`
    #+
    #+##Example usage
    #+```python
    #+>>> dic = {'a':0, 'b':[1.d0,2.d0], 'c':"example"}
    #+>>> schema = {'a':{'dims':0,'type':[int]}, 'b':{'dims':[2],'type':[float, np.float64]}, 'c':{'dims':0,'type':[str]}  }
    #+
    #+>>> err = check_dict_schema(schema, dic, desc="Example dict")
    #+>>> print(err)
    #+    False
    #+```
    if desc is None:
        desc = 'dict'

    err = False
    schema_keys = list(schema.keys())
    dic_keys = list(dic.keys())

    # Note extra variables
    for key in dic_keys:
        if key not in schema_keys:
            info('Extra variable "{}" found in "{}"'.format(key, desc))

    for key in schema_keys:
        # Note missing data
        if key not in dic_keys:
            error('"{}" is missing from "{}"'.format(key, desc))
            err = True
        else:
            # Check type
            if (schema[key]['dims'] == 0):
                if not isinstance(dic[key], tuple(schema[key]['type'])):
                    error('"{}" has the wrong type of {}. Expected {}'.format(key, type(dic[key]), schema[key]['type']))
                    err = True
            elif dic[key].dtype.type not in schema[key]['type']:
                error('"{}" has the wrong type of {}. Expected {}'.format(key, dic[key].dtype.type, schema[key]['type']))
                err = True

            # Check for NaNs or Inf
            if (not isinstance(dic[key], (str, dict, float, int))) and (str not in schema[key]['type']):
                if (dic[key][np.isnan(dic[key])].size > 0) or (dic[key][np.isinf(dic[key])].size > 0):
                    error('NaN or Infinity detected in "{}"'.format(key))
                    err = True

                # Check shape
                if not np.array_equal(dic[key].shape, schema[key]['dims']):
                    error('"{}" has the wrong shape of {}. Expected ({})'.format(key, dic[key].shape, schema[key]['dims']))
                    print('ndim({}) = {}'.format(key, dic[key].ndim))
                    err = True

            # Check shape
            if isinstance(dic[key], (str, int, float)):
                if (schema[key]['dims'] != 0) and (schema[key]['dims'] != [0]):
                    error('"{}" has the wrong shape. Expected ({})'.format(key, schema[key]['dims']))
                    err = True

    return err
