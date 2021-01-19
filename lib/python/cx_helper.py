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
    assert neutral_density.dtype.name == 'float64'
    assert neutral_velocity.dtype.name == 'float64'
    assert neutral_velocity.shape == (space_dimensions,)
    
    assert ion_velocities.dtype.name == 'float64'
    assert ion_velocities.ndim == 2
    assert ion_velocities.shape[0] == space_dimensions
    assert ion_velocities.flags['F_CONTIGUOUS']
    
    
    num_ions = ion_velocities.shape[1]
    num_ions_arr = np.array([num_ions],dtype='int64')
    rates = np.zeros((nlevels,num_ions),dtype='float64',order='F')
    lib.cx_rates(ref(neutral_density),
                ref(neutral_velocity),
                ref(np.asfortranarray(ion_velocities)),
                ref(num_ions_arr),
                ref(np.array([nlevels],dtype='int64')),
                ref(rates))
    return rates

from matplotlib.pyplot import loglog,legend,xlabel,ylabel
def test_rates():
    nlevels = 6
    space_dims = 3
    velocity_resolution = 50
    neutral_density = np.zeros((nlevels))
    neutral_density[0] = 1.0
    neutral_velocity = np.zeros((3))
    ion_velocities = np.zeros((space_dims,velocity_resolution),dtype='float64',order='F')
    v_range = np.logspace(6,9,velocity_resolution)
    ion_velocities[0,:]=v_range
    rates = cx_rates(neutral_density,neutral_velocity,ion_velocities)
    for ii in range(nlevels):
        loglog(v_range,rates[ii],label=(ii+1))
    legend(title='n=')
    xlabel('Velocity [cm/s]')
    ylabel('CX rate [cm^-3*s^-1]')
    return rates
