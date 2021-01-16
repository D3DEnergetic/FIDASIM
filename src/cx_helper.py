from cffi import FFI
import numpy as np
ffi = FFI()
lib = ffi.dlopen("/Users/lmorton/Code/FIDASIM_lib/src/cx_helper.so")
ffi.cdef("void cx_matrix(double arg1[], int arg2[], double arg3[]);")
ffi.cdef("void cx_rates(double arg1[], double arg2[], double arg3[], long arg4[], long arg5[], double arg6[]);")
converter = {'float64':'double*',
            'int32':'int*',
            'int64':'long*'}
def ref(x):
    """
    Convert array to pointer suitable for consumption by C/fortran
    """
    xdtype = x.dtype
    ctype = converter[xdtype.name]
    return ffi.cast(ctype,x.__array_interface__['data'][0])

def scalar(x,ctype):
    """
    Array-wrap non-arrays, especially scalars
    """
    return np.array(x,dtype=ctype,order='F')

def unscalar(x):
    """
    Unwrap single-element array
    """
    return x.item()

def scref(x,ctype):
    return ref(scalar(x,ctype))

def cx_matrix(velocity):
    nlevels = 6
    arg3 = np.empty((nlevels,nlevels),order='F',dtype= 'float64')
    lib.cx_matrix(scref(velocity,'float64'),
                 scref(nlevels,'int32'),
                 ref(arg3),)
    return arg3

def cx_rates(neutral_density,neutral_velocity,ion_velocities):
    space_dimensions = 3
    nlevels = 6
    assert neutral_density.shape == (nlevels,)
    assert neutral_velocity.shape == (space_dimensions,)
    assert ion_velocities.shape[0] == space_dimensions
    assert ion_velocities.ndim == 2
    num_ions = ion_velocities.shape[1]
    rates = np.empty((nlevels,num_ions),dtype='float64')
    lib.cx_rates(ref(neutral_density),
                ref(neutral_velocity),
                ref(ion_velocities),
                scref(num_ions,'int64'),
                scref(nlevels,'int64'),
                ref(rates))
    return rates