#!/usr/bin/env python
# -*- coding: utf-8 -*-

from lib.python_prefida.info import info
from lib.python_prefida.error import error
import numpy as np
import sys


def check_dict_schema(schema, dic, desc=None):
    """
    #check_struct_schema
     Check structure `dic` is formatted according to `schema`
    ***
    ##Input Arguments
         **schema**: structure schema

         **dic**: structure to check

    ##Output Arguments
         **err**: error code

    ##Keyword Arguments
         **desc**: description of structure `dic`

    ##Example usage
    ```idl
    IDL> dic = {a:0, b:[1.d0,2.d0], c:"example"}
    IDL> schema = {a:{dims:0,type:"INT"}, b:{dims:[2],type:"DOUBLE"}, c:{dims:0,type:"STRING"}  }

    IDL> check_struct_schema, schema, dic, err, desc="Example structure"
    IDL> print, err
        0
    ```
    """
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
#            if not isinstance(dic[key], schema[key]['type']):
#                error('"{}" has the wrong type of {}. Expected {}'.format(key, type(dic[key]), schema[key]['type']))
#                print('type({}) = {}'.format(key, type(dic[key])))
#                err = True
#            print(key)
            if (schema[key]['dims'] == 0):  # or (schema[key]['dims'] ==[0]):
#                print(schema[key]['type'], type(schema[key]['type']))
                if not isinstance(dic[key], schema[key]['type']):
                    error('"{}" has the wrong type of {}. Expected {}'.format(key, type(dic[key]), schema[key]['type']))
#                    print('type({}) = {}'.format(key, type(dic[key])))
                    err = True
            elif dic[key].dtype != schema[key]['type']:
                error('"{}" has the wrong type of {}. Expected {}'.format(key, dic[key].dtype, schema[key]['type']))
#                print('type({}) = {}'.format(key, type(dic[key])))
                err = True

            # Check for NaNs or Inf
#            if (not isinstance(dic[key], str)) and (not isinstance(dic[key], dict)):
            if not isinstance(dic[key], (str, dict, float, int)):
#                print(key, type(dic[key]))
                if (dic[key][np.isnan(dic[key])].size > 0) or (dic[key][np.isinf(dic[key])].size > 0):
#                    print(np.isnan(dic[key]))
#                    print(np.isinf(dic[key]))
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
#                    print('ndim({}) = {}'.format(key, dic[key].ndim))
                    err = True
#            elif dic[key].ndim != schema[key]['dims']:
#                error('"{}" has the wrong dimensions. Expected ({})'.format(key, schema[key]['dims']))
#                print('ndim({}) = {}'.format(key, dic[key].ndim))
#                err = True

    return err
