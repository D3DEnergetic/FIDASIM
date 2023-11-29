!+ This file contains the main routines for FIDASIM {!../VERSION!}
module libfida
!+ Main FIDASIM library
USE ISO_C_BINDING
USE H5LT !! High level HDF5 Interface
USE HDF5 !! Base HDF5
USE hdf5_utils !! Additional HDF5 routines
USE eigensystem, ONLY : eigen, linsolve
USE utilities
#ifdef _MPI
USE mpi_utils
#endif

implicit none

character(30) :: version = ''
    !+ FIDASIM version number
integer, dimension(8) :: time_start
    !+ Start time
integer, parameter, private   :: Int32   = 4
    !+ Defines a 32 bit integer
integer, parameter, private   :: Int64   = 8
    !+ Defines a 64 bit integer
integer, parameter, private   :: Float32 = 4
    !+ Defines a 32 bit floating point real
integer, parameter, private   :: Float64 = 8
    !+ Defines a 64 bit floating point real
integer, parameter :: charlim = 150
    !+ Defines character limit for files and directories

character(charlim) :: namelist_file
    !+ Input namelist file
integer, parameter :: nbif_type  = 1
    !+ Identifier for full energy NBI neutral interaction
integer, parameter :: nbih_type  = 2
    !+ Identifier for half energy NBI neutral interaction
integer, parameter :: nbit_type  = 3
    !+ Identifier for third energy NBI neutral interaction
integer, parameter :: dcx_type   = 4
    !+ Identifier for dcx neutral interaction
integer, parameter :: halo_type  = 5
    !+ Identifier for halo neutral interaction
integer, parameter :: fida_type  = 6
    !+ Identifier for fida neutral interaction
integer, parameter :: brems_type = 7
    !+ Identifier for bremsstrahlung interaction. Acts as dummy type
integer, parameter :: ntypes     = 7
    !+ Number of different types of neutrals

!! Physical units
real(Float64), parameter :: e_amu = 5.48579909070d-4
    !+ Atomic mass of an electron [amu]
real(Float64), parameter :: H1_amu = 1.007276466879d0
    !+ Atomic mass of Hydrogen-1 (protium) [amu]
real(Float64), parameter :: H2_amu = 2.013553212745d0
    !+ Atomic mass of Hydrogen-2 (deuterium) [amu]
real(Float64), parameter :: H3_amu = 3.01550071632d0
    !+ Atomic mass of Hydrogen-3 (tritium) [amu]
real(Float64), parameter :: He3_amu = 3.01602931914d0
    !+ Atomic mass of Helium-3 [amu]
real(Float64), parameter :: He4_amu = 4.00260325413d0
    !+ Atomic mass of Helium-4 [amu]
real(Float64), parameter :: B10_amu = 10.81d0
    !+ Atomic mass of Boron [amu]
real(Float64), parameter :: C12_amu = 12.011d0
    !+ Atomic mass of Carbon [amu]
real(Float64), parameter :: mass_u    = 1.660539040d-27
    !+ Atomic mass unit [kg]
real(Float64), parameter :: e0        = 1.60217733d-19
    !+ Electron charge [C]
real(Float64), parameter :: pi        = 3.14159265358979323846264d0
    !+ Pi
real(Float64), parameter :: c0        = 2.99792458d+08
    !+ Speed of light [m/s]
real(Float64), parameter :: h_planck  = 4.135667516d-15
    !+ Planck's constant [eV*s]
real(Float64), parameter :: v2_to_E_per_amu = mass_u/(2.*e0*1.d3)*1.d-4
    !+ \(cm^2/s^2\) to keV conversion factor
real(Float64), parameter :: log_10 = log(10.d0)
    !+ Natural log of 10.0
real(Float64), parameter :: a_0 = 5.29177210903d-11
    !+ bohr radius in [m]

real(Float64) :: line_lambda0(3)   = [0.00000d0, 0.00000d0, 0.00000d0]
    !+ H/D/T emission lines [nm]
integer :: n_stark = 0
    !+ Number of Stark lines
real(Float64), dimension(:), allocatable :: stark_wavel
    !+ Stark wavelengths [nm*m/V]
real(Float64), dimension(:), allocatable :: stark_intens
    !+ Stark Intensities
integer, dimension(:), allocatable :: stark_pi
    !+ Pi line indicators
integer, dimension(:), allocatable :: stark_sigma
    !+ Sigma line indicators
integer :: initial_state = 0
    !+ intial state energy level of transition
integer :: final_state = 0
    !+ final state energy level of transition

!!Numerical Settings
integer, parameter :: nlevs=6
    !+ Number of atomic energy levels
integer :: nbi_outside = 0
    !+ Keeps track of how many beam neutrals do not hit the [[libfida:beam_grid]]
integer, parameter :: reservoir_size = 100
    !+ Size of Neutral Particle Reservoir

!!Loop Parallization Settings
integer :: istart = 1
    !+ Starting loop counter (1 if OpenMP, processor number if MPI)
integer :: istep = 1
    !+ Loop step size (1 if OpenMP, number of processes if MPI)

!! Plasma Composition
integer, parameter :: max_species = 3
    !+ Maximum number of thermal isotopes
real(Float64), parameter :: supported_masses(max_species) = [H1_amu, H2_amu, H3_amu]
    !+ Supported Hydrogen isotopes masses [amu]

integer :: n_thermal = 1
    !+ Number of thermal hydrogen species/isotopes
real(Float64) :: thermal_mass(max_species) = 0.d0
    !+ Thermal ion species mass [amu]
real(Float64) :: thermal_lambda0(max_species) = 0.d0
    !+ Reference wavelengths for thermal species/isotopes [nm]
integer :: impurity_charge = 1
    !+ Charge of main impurity (boron=5, carbon=6,...)
real(Float64) :: impurity_mass = 0.d0
    !+ Impurity species mass [amu]
real(Float64) :: beam_mass = 0.d0
    !+ Beam species mass [amu]
real(Float64) :: beam_lambda0 = 0.d0
    !+ Reference wavelength for beam species

!integer :: n_fast = 1
!    !+ Number of Fast hydrogen species/isotopes
!real(Float64) :: fast_mass(max_species) = 0.d0
!    !+ Fast ion species mass [amu]
!real(Float64) :: fast_lambda0(max_species) = 0.d0
!    !+ Reference wavelength for fast species

type BeamGrid
    !+ Defines a 3D grid for neutral beam calculations
    integer(Int32) :: nx
        !+ Number of cells in the x direction
    integer(Int32) :: ny
        !+ Number of cells in the y direction
    integer(Int32) :: nz
        !+ Number of cells in the z direction
    real(Float64)  :: xmin
        !+ Minimum x value
    real(Float64)  :: xmax
        !+ Maximum x value
    real(Float64)  :: ymin
        !+ Minimum y value
    real(Float64)  :: ymax
        !+ Maximum y value
    real(Float64)  :: zmin
        !+ Minimum z value
    real(Float64)  :: zmax
        !+ Maximum z value
    real(Float64)  :: alpha
        !+ Tait-Bryan angle for a rotation about z [radians]
    real(Float64)  :: beta
        !+ Tait-Bryan angle for a rotation about y' [radians]
    real(Float64)  :: gamma
        !+ Tait-Bryan angle for a rotation about x" [radians]
    real(Float64)  :: drmin
        !+ Minimum cell spacing: `min(dx,dy,dz)`
    real(Float64)  :: dv
        !+ Cell volume [\(cm^3\)]
    real(Float64)  :: volume
        !+ Grid volume [\(cm^3\)]
    integer(Int32) :: ntrack
        !+ Maximum number of cell for particle tracking
    integer(Int32) :: ngrid
        !+ Number of cells
    integer(Int32), dimension(3)  :: dims
        !+ Dimensions of beam grid
    real(Float64), dimension(3)   :: origin
        !+ Origin of beam grid in machine coordinates
    real(Float64), dimension(3)   :: center
        !+ Center of beam grid in beam coordinates
    real(Float64), dimension(3)   :: dr
        !+ Cell spacings [dx, dy, dz]
    real(Float64), dimension(3)   :: lwh
        !+ Grid [length(x), width(y), height(z)]
    real(Float64), dimension(3,3) :: basis
        !+Beam grid basis for converting from beam coordinates(xyz)
        !+to machine coordinates(uvw): (\uvw = B*xyz + origin\)
    real(Float64), dimension(3,3) :: inv_basis
        !+Inverse basis for reverse transformation: (\xyz = B^{-1}*(uvw - origin)\)
    real(Float64), dimension(:), allocatable :: xc
        !+ x positions of cell centers
    real(Float64), dimension(:), allocatable :: yc
        !+ y positions of cell centers
    real(Float64), dimension(:), allocatable :: zc
        !+ z positions of cell centers
end type BeamGrid

type InterpolationGrid
    !+ Defines a 3D R-Z-phi grid for interpolating plasma parameters and fields
    integer(Int32) :: nr
        !+ Number of Radii
    integer(Int32) :: nz
        !+ Number of Z values
    integer(Int32) :: nphi
        !+ Number of phi values
    real(Float64)  :: dr
        !+ Radial spacing [cm]
    real(Float64)  :: dz
        !+ Vertical spacing [cm]
    real(Float64)  :: dphi
        !+ Angular spacing [rad]
    real(Float64)  :: da
        !+ Grid element area [\(cm^2\)]
    real(Float64)  :: dv
        !+ dr*dz*dphi [\(rad*cm^2\)]
    integer(Int32) :: dims(3)
        !+ Dimension of the interpolation grid
    real(Float64), dimension(:),   allocatable :: r
        !+ Radii values [cm]
    real(Float64), dimension(:),   allocatable :: z
        !+ Z values [cm]
    real(Float64), dimension(:),   allocatable :: phi
        !+ Angular values [rad]
    integer(Int32) :: ntrack
        !+ Maximum number of cells for particle tracking
    integer(Int32) :: ngrid
        !+ Number of cells
end type InterpolationGrid

type Profiles
    !+ Torodial symmetric plasma parameters at a given R-Z
    real(Float64) :: dene = 0.d0
        !+ Electron density [\(cm^{-3}\)]
    real(Float64) :: deni(max_species) = 0.d0
        !+ Ion density [\(cm^{-3}\)]
    real(Float64) :: denimp = 0.d0
        !+ Impurity density [\(cm^{-3}\)]
    real(Float64) :: denf = 0.d0
        !+ Fast-ion density [\(cm^{-3}\)]
    real(Float64) :: te = 0.d0
        !+ Electron temperature [kev]
    real(Float64) :: ti = 0.d0
        !+ Ion temperature [kev]
    real(Float64) :: zeff = 0.d0
        !+ Effective Nuclear Charge
    real(Float64) :: vr = 0.d0
        !+ Plasma rotation in radial direction
    real(Float64) :: vt = 0.d0
        !+ Plasma rotation in torodial/phi direction
    real(Float64) :: vz = 0.d0
        !+ Plasma rotation in z direction
    real(Float64) :: denn(nlevs,max_species) = 0.d0
        !+ Cold neutral densities [\(cm^{-3}\)]
end type Profiles

type, extends( Profiles ) :: LocalProfiles
    !+ Plasma parameters at given position
    logical :: in_plasma = .False.
        !+ Indicates whether plasma parameters are valid/known
    integer :: coords = 0
        !+ Indicates coordinate system of vectors. Beam grid (0), machine (1) and cylindrical (2)
    real(Float64), dimension(3) :: pos = 0.d0
        !+ Position in beam grid coordinates
    real(Float64), dimension(3) :: uvw = 0.d0
        !+ Position in machine coordinates
    real(Float64), dimension(3) :: vrot = 0.d0
        !+ Plasma rotation in beam grid coordinates
    real(Float64), dimension(3) :: vrot_uvw = 0.d0
        !+ Plasma rotation in machine coordinates
    type(InterpolCoeffs3D) :: b
        !+ Cylindrical Interpolation Coefficients and indicies for interpolation at `pos`
end type LocalProfiles

type EMFields
    !+ Torodial symmetric electro-magnetic fields at given R-Z
    real(Float64) :: br = 0.d0
        !+ Radial magnetic field [T]
    real(Float64) :: bt = 0.d0
        !+ Torodial magnetic field [T]
    real(Float64) :: bz = 0.d0
        !+ Vertical magnetic field [T]
    real(Float64) :: er = 0.d0
        !+ Radial electric field [V/m]
    real(Float64) :: et = 0.d0
        !+ Torodial electric field [V/m]
    real(Float64) :: ez = 0.d0
        !+ Vertical electric field [V/m]
    real(Float64) :: dbr_dr = 0.d0
        !+ Radial derivative of the radial magnetic field [T/m]
    real(Float64) :: dbr_dphi = 0.d0
        !+ Angular derivative of the radial magnetic field [T/m]
    real(Float64) :: dbr_dz = 0.d0
        !+ Vertical derivative of the radial magnetic field [T/m]
    real(Float64) :: dbt_dr = 0.d0
        !+ Radial derivative of the torodial magnetic field [T/m]
    real(Float64) :: dbt_dphi = 0.d0
        !+ Angular derivative of the torodial magnetic field [T/m]
    real(Float64) :: dbt_dz = 0.d0
        !+ Vertical derivative of the torodial magnetic field [T/m]
    real(Float64) :: dbz_dr = 0.d0
        !+ Radial derivative of the radial magnetic field [T/m]
    real(Float64) :: dbz_dphi = 0.d0
        !+ Angular derivative of the radial magnetic field [T/m]
    real(Float64) :: dbz_dz = 0.d0
        !+ Vertical derivative of the vertical magnetic field [T/m]
end type EMFields

type, extends( EMFields ) :: LocalEMFields
    !+ Electro-magnetic fields at given position
    logical       :: in_plasma = .False.
        !+ Indicates whether fields are valid/known
    integer :: coords = 0
        !+ Indicates coordinate system of vectors. Beam grid (0), machine (1) and cylindrical (2)
    real(Float64) :: b_abs = 0.d0
        !+ Magnitude of magnetic field
    real(Float64) :: e_abs = 0.d0
        !+ Magnitude of electrin field
    real(Float64), dimension(3) :: pos = 0.d0
        !+ Position in beam grid coordinates
    real(Float64), dimension(3) :: uvw = 0.d0
        !+ Position in machine coordinates
    real(Float64), dimension(3) :: b_norm = 0.d0
        !+ Direction of magnetic field in beam grid coordinates
    real(Float64), dimension(3) :: a_norm = 0.d0
        !+ Vector perpendicular to `b_norm` and `c_norm`
    real(Float64), dimension(3) :: c_norm = 0.d0
        !+ Vector perpendicular to `b_norm` and `a_norm`
    real(Float64), dimension(3) :: e_norm = 0.d0
        !+ Direction of electric field in beam grid coordinates
    type(InterpolCoeffs3D) :: b
        !+ Cylindrical Interpolation Coefficients and indicies for interpolation at `pos`
end type LocalEMFields

type Equilibrium
    !+MHD Equilbrium
    type(EMFields), dimension(:,:,:), allocatable :: fields
        !+ Electro-magnetic fields at points defined in [[libfida:inter_grid]]
    type(Profiles), dimension(:,:,:), allocatable :: plasma
        !+ Plasma parameters at points defined in [[libfida:inter_grid]]
    real(Float64), dimension(:,:,:), allocatable  :: mask
        !+ Indicates whether fields and plasma are well-defined at points defined in [[libfida:inter_grid]]
end type Equilibrium

type FastIonDistribution
    !+ Defines a Guiding Center Fast-ion Distribution Function: F(E,p,R,Z,Phi)
    real(Float64) :: A = H2_amu
        !+ Atomic Mass
    integer(Int32) :: nenergy
        !+ Number of energies
    integer(Int32) :: npitch
        !+ Number of pitches
    integer(Int32) :: nr
        !+ Number of radii
    integer(Int32) :: nz
        !+ Number of z values
    integer(Int32) :: nphi
        !+ Number of phi values
    real(Float64)  :: dE
        !+ Energy spacing [keV]
    real(Float64)  :: dp
        !+ Pitch spacing
    real(Float64)  :: dr
        !+ Radial spacing [cm]
    real(Float64)  :: dz
        !+ Z spacing [cm]
    real(Float64)  :: dphi
        !+ Angular spacing [rad]
    real(Float64)  :: emin
        !+ Minimum energy [keV]
    real(Float64)  :: emax
        !+ Maximum energy [keV]
    real(Float64)  :: e_range
        !+ Energy interval length [keV]
    real(Float64)  :: pmin
        !+ Minimum pitch
    real(Float64)  :: pmax
        !+ Maximum pitch
    real(Float64)  :: p_range
        !+ Pitch interval length
    real(Float64)  :: rmin
        !+ Minimum radius [cm]
    real(Float64)  :: rmax
        !+ Maximum radius [cm]
    real(Float64)  :: r_range
        !+ Radius interval length [cm]
    real(Float64)  :: zmin
        !+ Minimum Z [cm]
    real(Float64)  :: zmax
        !+ Maximum Z [cm]
    real(Float64)  :: z_range
        !+ Z interval length [cm]
    real(Float64)  :: phimin
        !+ Minimum Phi [rad]
    real(Float64)  :: phimax
        !+ Maximum Phi [rad]
    real(Float64)  :: phi_range
        !+ Phi interval length [rad]
    real(Float64)  :: n_tot = 0.d0
        !+ Total Number of fast-ions
    real(Float64), dimension(:), allocatable       :: energy
        !+ Energy values [keV]
    real(Float64), dimension(:), allocatable       :: pitch
        !+ Pitch w.r.t. the magnetic field
    real(Float64), dimension(:), allocatable       :: r
        !+ Radius [cm]
    real(Float64), dimension(:), allocatable       :: z
        !+ Z [cm]
    real(Float64), dimension(:), allocatable       :: phi
        !+ Angles [rad]
    real(Float64), dimension(:,:,:), allocatable     :: denf
        !+ Fast-ion density defined on the [[libfida:inter_grid]]: denf(R,Z,Phi)
    real(Float64), dimension(:,:,:,:,:), allocatable :: f
        !+ Fast-ion distribution function defined on the [[libfida:inter_grid]]: F(E,p,R,Z,Phi)
end type FastIonDistribution

type FastIon
    !+ Defines a fast-ion
    real(Float64)  :: A = H2_amu
        !+ Fast-ion atomic mass
    logical        :: beam_grid_cross_grid = .False.
        !+ Indicates whether the fast-ion crosses the [[libfida:beam_grid]]
    real(Float64)  :: r = 0.d0
        !+ Radial position of fast-ion [cm]
    real(Float64)  :: phi = 0.d0
        !+ Angular position of fast-ion [rad]
    real(Float64)  :: z = 0.d0
        !+ Vertical position of fast-ion [cm]
    real(Float64)  :: beam_grid_phi_enter = 0.d0
        !+ Torodial/phi position where fast-ion enters the [[libfida:beam_grid]] [radians]
    real(Float64)  :: delta_phi = 2*pi
        !+ Angle subtended by the [[libfida:beam_grid]] at (r,z)
    real(Float64)  :: energy = 0.d0
        !+ Energy [keV]
    real(Float64)  :: pitch = 0.d0
        !+ Pitch w.r.t. the magnetic field
    real(Float64)  :: vabs = 0.d0
        !+ Speed [cm/s]
    real(Float64)  :: vr = 0.d0
        !+ Radial velocity [cm/s]
    real(Float64)  :: vt = 0.d0
        !+ Torodial velocity [cm/s]
    real(Float64)  :: vz = 0.d0
        !+ Z velocity [cm/s]
    real(Float64)  :: weight = 0.d0
        !+ Particle weight: How many fast-ions does particle represent.
    integer(Int32) :: class = 0
        !+ Orbit class id
end type FastIon

type FastIonParticles
    !+ Collection of fast-ion particles
    integer(Int32) :: nparticle = 0
        !+ Number of particles
    integer(Int32) :: nclass = 1
        !+ Number of orbit classes
    logical :: axisym = .True.
        !+ Indicates whether distribution function is axisymmetric
    type(FastIon), dimension(:), allocatable :: fast_ion
        !+ Fast-ion particles
end type FastIonParticles

type NeutralBeam
    !+ Defines a neutral beam with +x defined to be into the plasma
    character(25) :: name = ''
        !+ Beam name
    integer       :: shape
        !+ Beam source shape 1="rectangular", 2="circular"
    real(Float64) :: widy
        !+ Half width of source in y direction
    real(Float64) :: widz
        !+ Half height of source in z direction
    real(Float64) :: focy
        !+ Focal length in y direction
    real(Float64) :: focz
        !+ Focal length in z direction
    real(Float64) :: einj
        !+ NBI voltage  [kV]
    real(Float64) :: pinj
        !+ NBI power    [MW]
    real(Float64) :: vinj
        !+ NBI velocity [cm/s]
    real(Float64) :: alpha
        !+ Z rotation not same as [[libfida:beam_grid]] alpha
    real(Float64) :: beta
        !+ Tilt rotation not same as [[libfida:beam_grid]] beta
    real(Float64), dimension(3)   :: divy
        !+ Energy dependent divergence in y direction
    real(Float64), dimension(3)   :: divz
        !+ Energy dependent divergence in z direction
    real(Float64), dimension(3)   :: current_fractions
        !+ Fractions of full, half, and third energy neutrals
    real(Float64), dimension(3)   :: src
        !+ Position of source in beam grid coordinates [cm]
    real(Float64), dimension(3)   :: axis
        !+ Beam centerline
    integer :: naperture
        !+ Number of beam apertures
    integer, dimension(:), allocatable       :: ashape
        !+ Aperture shape 1="rectangular", 2="circular"
    real(Float64), dimension(:), allocatable :: awidy
        !+ Half width of the aperture(s) in y direction
    real(Float64), dimension(:), allocatable :: awidz
        !+ Half height of the aperture(s) in z direction
    real(Float64), dimension(:), allocatable :: aoffy
        !+ Horizontal (y) offset of the aperture(s) relative to the beam centerline [cm]
    real(Float64), dimension(:), allocatable :: aoffz
        !+ Vertical (z) offset of the aperture(s) relative to the beam centerline [cm]
    real(Float64), dimension(:), allocatable :: adist
        !+ Distance from the center of the beam source grid to the aperture(s) plane [cm]
    real(Float64), dimension(3,3) :: basis
        !+ Beam basis for converting from centerline coordinates to beam grid coordinates
    real(Float64), dimension(3,3) :: inv_basis
        !+ Inverse basis for reverse transfomation
end type NeutralBeam

type AtomicCrossSection
    !+ Defines a n/m-resolved atomic cross section table
    integer       :: nenergy = 1
        !+ Number of beam energies
    real(Float64) :: logemin = 0.d0
        !+ Log-10 minimum energy
    real(Float64) :: logemax = 0.d0
        !+ Log-10 maximum energy
    integer       :: n_max = nlevs
        !+ Number of initial atomic energy levels
    integer       :: m_max = nlevs
        !+ Number of final atomic energy levels
    real(Float64) :: dlogE = 0.d0
        !+ Log-10 energy spacing
    real(Float64) :: minlog_cross
        !+ Log-10 minimum cross section
    real(Float64), dimension(:,:,:), allocatable :: log_cross
        !+ Log-10 cross sections
end type AtomicCrossSection

type AtomicRates
    !+ Defines a n/m-resolved atomic cross section table
    integer       :: nenergy = 1
        !+ Number of beam energies/amu
    real(Float64) :: logemin = 0.d0
        !+ Log-10 minimum energy/amu
    real(Float64) :: logemax = 0.d0
        !+ Log-10 maximum energy/amu
    integer       :: ntemp = 1
        !+ Number of target temperatures
    real(Float64) :: logtmin = 0.d0
        !+ Log-10 minimum temperature
    real(Float64) :: logtmax = 0.d0
        !+ Log-10 maximum temperature
    integer       :: n_max = nlevs
        !+ Number of initial atomic energy levels
    integer       :: m_max = nlevs
        !+ Number of final atomic energy levels
    real(Float64) :: dlogE = 0.d0
        !+ Log-10 energy spacing
    real(Float64) :: dlogT = 0.d0
        !+ Log-10 temperature spacing
    real(Float64) :: minlog_rate = 0.d0
        !+ Log-10 minimum reaction rate
    real(Float64), dimension(:,:,:,:), allocatable :: log_rate
        !+ Log-10 beam-target rates
end type AtomicRates

type AtomicTransitions
    !+ Defines an atomic table for populating and de-populating reaction rates
    integer       :: nenergy = 1
        !+ Number of beam energies/amu
    real(Float64) :: logemin = 0.d0
        !+ Log-10 minimum energy/amu
    real(Float64) :: logemax = 0.d0
        !+ Log-10 maximum energy/amu
    integer       :: ntemp = 1
        !+ Number of target temperatures
    real(Float64) :: logtmin = 0.d0
        !+ Log-10 minimum temperature
    real(Float64) :: logtmax = 0.d0
        !+ Log-10 maximum temperature
    integer       :: n_max = nlevs
        !+ Number of initial atomic energy levels
    integer       :: m_max = nlevs
        !+ Number of final atomic energy levels
    real(Float64) :: dlogE = 0.d0
        !+ Log-10 energy spacing
    real(Float64) :: dlogT = 0.d0
        !+ Log-10 temperature spacing
    real(Float64) :: minlog_pop = 0.d0
        !+ Log-10 minimum reaction rates for populating transistions
    real(Float64) :: minlog_depop = 0.d0
        !+ Log-10 minimum reaction rates for de-populating transistions
    real(Float64), dimension(:,:,:,:), allocatable :: log_pop
        !+ Log-10 reaction rates for populating transistions
    real(Float64), dimension(:,:,:), allocatable   :: log_depop
        !+ Log-10 reaction rates for de-populating transistions
end type AtomicTransitions

type NuclearRates
    !+ Nuclear reaction rates
    integer       :: nbranch = 1
        !+ Number of reaction branches
    integer       :: nenergy = 1
        !+ Number of beam energies
    real(Float64) :: logemin = 0.d0
        !+ Log-10 minimum energy
    real(Float64) :: logemax = 0.d0
        !+ Log-10 maximum energy
    integer       :: ntemp = 1
        !+ Number of target temperatures
    real(Float64) :: logtmin = 0.d0
        !+ Log-10 minimum temperature
    real(Float64) :: logtmax = 0.d0
        !+ Log-10 maximum temperature
    real(Float64) :: dlogE = 0.d0
        !+ Log-10 energy spacing
    real(Float64) :: dlogT = 0.d0
        !+ Log-10 temperature spacing
    real(Float64) :: minlog_rate = 0.d0
        !+ Log-10 minimum reaction rate
    real(Float64), dimension(2) :: bt_amu = 0.d0
        !+ Isotope mass of beam and thermal ions respectively [amu]
    real(Float64), dimension(:,:,:), allocatable :: log_rate
        !+ Log-10 reaction rates: log_rate(energy, temperature, branch)
end type NuclearRates

type AtomicTables
    !+ Atomic tables for various types of interactions
    type(AtomicCrossSection) :: H_H_cx_cross
        !+ Hydrogen-Hydrogen charge exchange n/m-resolved cross sections
    type(AtomicRates)        :: H_H_cx_rate
        !+ Hydrogen-Hydrogen charge exchange n/m-resolved beam-target rates
    type(AtomicTransitions)  :: H_H
        !+ Hydrogen-Hydrogen atomic transitions
    type(AtomicTransitions)  :: H_e
        !+ Hydrogen-Electron atomic transitions
    type(AtomicTransitions)  :: H_Aq
        !+ Hydrogen-Impurity atomic transitions
    real(Float64), dimension(nlevs,nlevs) :: einstein
        !+ Einstein coefficients for spontaneous emission
    type(NuclearRates)       :: D_D
        !+ Deuterium-Deuterium reaction rates
end type AtomicTables

type LineOfSight
    !+ Defines a line of sight
    real(Float64) :: sigma_pi = 1.d0
        !+ Ratio of sigma to pi line intensity
    real(Float64) :: spot_size = 0.d0
        !+ Radius of spot size [cm]
    real(Float64), dimension(3) :: lens = 0.d0
        !+ Lens location in beam grid coordinates
    real(Float64), dimension(3) :: axis = 0.d0
        !+ Optical axis in beam grid coordinates
    real(Float64), dimension(3) :: lens_uvw = 0.d0
        !+ Lens location in machine coordinates
    real(Float64), dimension(3) :: axis_uvw = 0.d0
        !+ Optical axis in machine coordinates
end type LineOfSight

type LOSElement
    !+ Defines a element of a line of sight and cell intersection
    integer :: id
        !+ Line of sight index
    real(Float64) :: length
        !+ Length of crossing
end type LOSElement

type LOSInters
    !+ Defines the channels that intersect a cell
    integer :: nchan = 0
        !+ Number of channels that intersect
    type(LOSElement), dimension(:), allocatable :: los_elem
        !+ Array of crossing
end type LOSInters

type SpectralChords
    !+ Defines an spectral diagnostic system
    integer :: nchan = 0
        !+ Number of channels
    integer :: ncell = 0
        !+ Number of beam_grid cells with intersections
    integer :: cyl_ncell = 0
        !+ Number of pass_grid cells with intersections
    type(LineOfSight), dimension(:), allocatable :: los
        !+ Line of sight array
    real(Float64), dimension(:), allocatable     :: radius
        !+ Radius of each line of sight
    type(LOSInters), dimension(:,:,:), allocatable :: inter
        !+ Array of LOS intersections with [[libfida:beam_grid]]
    type(LOSInters), dimension(:,:,:), allocatable :: cyl_inter
        !+ Array of LOS intersections with [[libfida:pass_grid]]
    integer, dimension(:), allocatable :: cell
        !+ Linear indices of beam_grid that have intersections
    integer, dimension(:), allocatable :: cyl_cell
        !+ Linear indices of pass_grid that have intersections
end type SpectralChords

type BoundedPlane
    !+ Defines a plane with a circular or rectangular boundary
    integer       :: shape    = 0
        !+ Boundary shape 1="Rectangular", 2="circular"
    real(Float64) :: hh       = 0.d0
        !+ Half height of boundary [cm]
    real(Float64) :: hw       = 0.d0
        !+ Half width of boundary [cm]
    real(Float64), dimension(3)   :: origin   = 0.d0
        !+ Origin of plane in machine coordinates
    real(Float64), dimension(3,3) :: basis    = 0.d0
        !+ Basis vectors basis(:,1) = u_1 is plane normal
    real(Float64), dimension(3,3) :: inv_basis= 0.d0
        !+ Inverse basis
end type BoundedPlane

type NPADetector
    !+ Defines a NPA detector
    type(BoundedPlane) :: detector
        !+ Detecting plane of NPA detector
    type(BoundedPlane) :: aperture
        !+ Aperture plane of NPA detector
end type NPADetector

type NPAProbability
    !+ Type to contain the probability of hitting a NPA detector
    real(Float64) :: p = 0.d0
        !+ Hit probability
    real(Float64) :: pitch = -2.d0
        !+ Pitch
    real(Float64), dimension(3) :: eff_rd = 0.d0
        !+ Effective position of detector
    real(Float64), dimension(3) :: dir = 0.d0
        !+ Trajectory direction
end type NPAProbability

type NPAChords
    !+ Defines a NPA system
    integer :: nchan = 0
         !+ Number of channels
    type(NPADetector), dimension(:), allocatable          :: det
         !+ NPA detector array
    real(Float64), dimension(:), allocatable              :: radius
         !+ Radius [cm]
    logical, dimension(:,:,:), allocatable                :: hit
         !+ Indicates whether a particle can hit any NPA detector from a grid cell: hit(x,y,z)
    type(NPAProbability), dimension(:,:,:,:), allocatable :: phit
         !+ Probability of hitting a detector from a grid cell: phit(x,y,z,chan)
end type NPAChords

type NPAParticle
    !+ Defines a NPA particle
    integer       :: detector = 0
        !+ Detector NPA particle hit
    integer(Int32) :: class = 0
        !+ Orbit class id
    real(Float64) :: A
        !+ Atomic Mass
    real(Float64) :: xi = 0.d0
        !+ Initial x position
    real(Float64) :: yi = 0.d0
        !+ Initial y position
    real(Float64) :: zi = 0.d0
        !+ Initial z position
    real(Float64) :: xf = 0.d0
        !+ Final x position
    real(Float64) :: yf = 0.d0
        !+ Final y position
    real(Float64) :: zf = 0.d0
        !+ Final z position
    real(Float64) :: weight = 0.d0
        !+ NPA particle weight
    real(Float64) :: energy = 0.d0
        !+ Birth Energy [keV]
    real(Float64) :: pitch = 0.d0
        !+ Birth Pitch
end type NPAParticle

type NPAResults
    !+ MC NPA result structure
    integer(Int32) :: nchan = 0
        !+ Number of NPA channels
    integer(Int32) :: npart = 0
        !+ Number of particles that hit a detector
    integer(Int32) :: nmax = 1000000
        !+ Maximum allowed number of particles grows if necessary
    integer(Int32) :: nenergy = 122
        !+ Number of energy values
    type(NPAParticle), dimension(:), allocatable :: part
        !+ Array of NPA particles
    real(Float64), dimension(:), allocatable     :: energy
        !+ Energy array [keV]
    real(Float64), dimension(:,:,:), allocatable :: flux
        !+ Neutral particle flux: flux(energy,chan, orbit_type) [neutrals/(s*dE)]
end type NPAResults

type BirthParticle
    !+ Defines a Birth particle
    integer :: neut_type = 0
        !+ Birth type (1=Full, 2=Half, 3=Third)
    integer(Int32), dimension(3) :: ind = 0
        !+ Initial [[libfida:beam_grid]] indices
    real(Float64), dimension(3) :: ri = 0.d0
        !+ Initial position in beam grid coordinates [cm]
    real(Float64), dimension(3) :: vi = 0.d0
        !+ Initial velocity in beam grid coordinates [cm/s]
    real(Float64), dimension(3) :: ri_gc = 0.d0
        !+ Initial guiding-center position in beam grid coordinates [cm]
    real(Float64) :: weight = 0.d0
        !+ NPA particle weight [fast-ions/s]
    real(Float64) :: energy = 0.d0
        !+ Birth Energy [keV]
    real(Float64) :: pitch = 0.d0
        !+ Birth Pitch w.r.t. the magnetic field
end type BirthParticle

type BirthProfile
    !+ Birth profile structure
    integer :: cnt = 1
        !+ Particle counter
    type(BirthParticle), dimension(:), allocatable :: part
        !+ Array of birth particles
    real(Float64), dimension(:,:,:,:), allocatable :: dens
        !+ Birth density: dens(neutral_type,x,y,z) [fast-ions/(s*cm^3)]
end type BirthProfile

type Spectra
    !+ Spectra storage structure
    real(Float64), dimension(:,:), allocatable   :: brems
        !+ Bremsstruhlung: brems(lambda,chan)
    real(Float64), dimension(:,:,:), allocatable   :: full
        !+ Full energy beam emission stark components: full(n_stark,lambda,chan)
    real(Float64), dimension(:,:,:), allocatable   :: half
        !+ Half energy beam emission stark components: half(n_stark,lambda,chan)
    real(Float64), dimension(:,:,:), allocatable   :: third
        !+ Third energy beam emission stark components: thirdstokes(n_stark,lambda,chan)
    real(Float64), dimension(:,:,:,:), allocatable   :: fullstokes
        !+ Full energy beam emission stark components: fullstokes(n_stark,4,lambda,chan)
    real(Float64), dimension(:,:,:,:), allocatable   :: halfstokes
        !+ Half energy beam emission stark components: halfstokes(n_stark,4,lambda,chan)
    real(Float64), dimension(:,:,:,:), allocatable   :: thirdstokes
        !+ Third energy beam emission stark components: thirdstokes(n_stark,4,lambda,chan)

    real(Float64), dimension(:,:,:,:), allocatable :: dcx
        !+ Direct CX emission stark components: dcx(n_stark,lambda,chan,species)
    real(Float64), dimension(:,:,:,:,:), allocatable :: dcxstokes
        !+ Direct CX emission stark components: dcxstokes(n_stark,4,lambda,chan,species)
    real(Float64), dimension(:,:,:,:), allocatable :: halo
        !+ Thermal halo emission stark components: halo(n_stark,lambda,chan,species)
    real(Float64), dimension(:,:,:,:,:), allocatable :: halostokes
        !+ Thermal halo emission stark components: halo(n_stark,4,lambda,chan,species)
    real(Float64), dimension(:,:,:,:), allocatable :: cold
        !+ Cold D-alpha emission stark components: cold(n_stark,lambda,chan,species)
    real(Float64), dimension(:,:,:,:,:), allocatable :: coldstokes
        !+ Cold D-alpha emission stark components: coldstokes(n_stark,4,lambda,chan,species)
    real(Float64), dimension(:,:,:,:), allocatable :: fida
        !+ Active FIDA emission stark components: fida(n_stark,lambda,chan,orbit_type)
    real(Float64), dimension(:,:,:,:,:), allocatable :: fidastokes
        !+ Active FIDA emission stark components: fidastokes(n_stark,4,lambda,chan,orbit_type)
    real(Float64), dimension(:,:,:,:), allocatable :: pfida
        !+ Passive FIDA emission stark components: pfida(n_stark,lambda,chan,orbit_type)
    real(Float64), dimension(:,:,:,:,:), allocatable :: pfidastokes
        !+ Passive FIDA emission stark components: pfidastokes(n_stark,4,lambda,chan,orbit_type)
end type Spectra

type NeutronRate
    !+ Neutron storage structure
    real(Float64), dimension(:), allocatable :: rate
        !+ Neutron rate: rate(orbit_type) [neutrons/sec]
    real(Float64), dimension(:,:,:,:,:), allocatable :: weight
        !+ Neutron rate weight: weight(E,p,R,Z,Phi)
    real(Float64), dimension(:,:,:), allocatable :: emis
        !+ Neutron emissivity: emis(R,Z,Phi)
end type NeutronRate

type CFPDRate
    !+ Charged Fusion Product (CFPD) storage structure
    real(Float64), dimension(:,:,:,:), allocatable :: weight
        !+ CFPD rate weight: weight(Ep,Ch,E,p)
    real(Float64), dimension(:,:), allocatable         :: flux
        !+ CFPD flux: flux(E3,chan) [kHz]
    real(Float64), dimension(:,:), allocatable         :: prob
        !+ CFPD flux: probability_gyro(E3,chan) [unity]
    real(Float64), dimension(:,:), allocatable         :: gam
        !+ CFPD flux: gyro(E3,chan) [rad]
end type CFPDRate

type NeutralParticle
    real(Float64) :: w = 1.d0
        !+ Neutral particle weight
    real(Float64), dimension(3)  :: v = 0.d0
        !+ Neutral particle position/velocity
end type NeutralParticle

type NeutralParticleReservoir
    !+ Neutral Particle Reservoir
    !+###References
    !+ Li, Kim-Hung (4 December 1994). "Reservoir-Sampling Algorithms of Time
    !+ Complexity O(n(1+log(N/n)))". ACM Transactions on Mathematical Software. 20
    !+ (4): 481–493. doi:10.1145/198429.198435
    integer(Int32) :: n = 0
        !+ Number of neutral particles seen by the reservoir
    integer(Int32) :: i = 0
        !+ Number of the next particle that will be stored in the reservoir
    integer(Int32) :: k = 0
        !+ Size of the reservoir
    real(Float64)  :: W = 0.d0
        !+ Sampling weight
    type(NeutralParticle), dimension(reservoir_size) :: R
        !+ Neutral Particle Reservoir
end type NeutralParticleReservoir

type NeutralPopulation
    real(Float64), dimension(:,:,:,:), allocatable :: dens
        !+ Neutral density: dens(lev,x,y,z)
    type(NeutralParticleReservoir), dimension(:,:,:), allocatable :: res
        !+ Neutral particle reservoir
end type NeutralPopulation

type Neutrals
    !+ Neutral density structure
    type(NeutralPopulation) :: full
        !+ Full energy neutral population
    type(NeutralPopulation) :: half
        !+ Half energy neutral population
    type(NeutralPopulation) :: third
        !+ Third energy neutral population
    type(NeutralPopulation) :: dcx
        !+ DCX neutral population
    type(NeutralPopulation) :: halo
        !+ Halo neutral population
!TODO: Implement Following
!    type(NeutralPopulation) :: fida
!        !+ FIDA neutral population
!    type(NeutralPopulation) :: fida_halo
!        !+ FIDA-Halo neutral population
end type Neutrals

type SpatialSpectra
    !+ Spatial Storage of photon neutral birth
    type(NeutralParticleReservoir), dimension(:), allocatable :: full
    type(NeutralParticleReservoir), dimension(:), allocatable :: half
    type(NeutralParticleReservoir), dimension(:), allocatable :: third
    type(NeutralParticleReservoir), dimension(:), allocatable :: dcx
    type(NeutralParticleReservoir), dimension(:), allocatable :: halo
    type(NeutralParticleReservoir), dimension(:), allocatable :: fida
    type(NeutralParticleReservoir), dimension(:), allocatable :: pfida
end type SpatialSpectra

type FIDAWeights
    !+ FIDA weights structure
    real(Float64), dimension(:,:,:), allocatable   :: mean_f
        !+ Estimate of mean fast-ion distribution function "seen" by LOS: mean_f(E,p,chan)
    real(Float64), dimension(:,:,:,:), allocatable :: weight
        !+ FIDA weight function: weight(lambda,E,p,chan)
end type FIDAWeights

type NPAWeights
    !+ NPA weights structure
    real(Float64), dimension(:,:,:,:,:), allocatable :: attenuation
        !+ Attenuation fraction: attenuation(E,x,y,z,chan)
    real(Float64), dimension(:,:,:,:,:), allocatable :: cx
        !+ Charge Exchange reaction rates: cx(E,x,y,z,chan)
    real(Float64), dimension(:,:,:,:), allocatable   :: emissivity
        !+ Emissivity: emissivity(x,y,z,chan) [neutrals/(s*dV)]
    real(Float64), dimension(:,:,:), allocatable     :: weight
        !+ NPA weight function: weight(E,p,chan) [neutrals/(s*fast-ion*dE*dP)]
    real(Float64), dimension(:,:), allocatable       :: flux
        !+ Neutral particle flux: flux(E,chan) [neutrals/(s*dE)]
end type NPAWeights

type SimulationInputs
    !+ Simulation settings structure
    integer(Int32) :: shot_number
        !+ Shot Number
    real(Float64)  :: time
        !+ Shot time [s]
    character(charlim) :: runid = ''
        !+ FIDASIM run ID
    character(charlim) :: result_dir = ''
        !+ Result directory
    character(charlim) :: tables_file = ''
        !+ Atomic tables file
    character(charlim) :: geometry_file = ''
        !+ FIDASIM input file containing geometric quantities
    character(charlim) :: equilibrium_file = ''
        !+ FIDASIM input file containing the plasma parameters and fields
    character(charlim) :: distribution_file = ''
        !+ FIDASIM input file containing the fast-ion distribution
    character(charlim) :: neutrals_file = ''
        !+ FIDASIM output/input file containing beam neutral density.
        !+ Used when [[SimulationInputs:load_neutrals]] is set.

    !! Random Number Generator Settings
    integer :: seed
        !+ Random number generator seed

    !! Monte Carlo settings
    integer(Int64) :: n_fida
        !+ Number of Active FIDA mc markers
    integer(Int64) :: n_pfida
        !+ Number of Passive FIDA mc markers
    integer(Int64) :: n_npa
        !+ Number of Passiv NPA mc markers
    integer(Int64) :: n_pnpa
        !+ Number of Passive NPA mc markers
    integer(Int64) :: n_nbi
        !+ Number of neutral beam mc markers
    integer(Int64) :: n_dcx
        !+ Number of direct charge exchange (DCX) mc markers
    integer(Int64) :: n_halo
        !+ Number of halo mc markers
    integer(Int64) :: n_birth
        !+ Number of birth particles per [[SimulationInputs:n_nbi]]

    !! Simulation switches
    integer(Int32) :: calc_spec
        !+ Calculate spectra: 0 = off, 1=on
    integer(Int32) :: calc_beam
        !+ Calculate beam densities: 0 = off, 1=on
    integer(Int32) :: calc_nbi_dens
        !+ Calculate neutral beam density: 0 = off, 1=on
    integer(Int32) :: calc_dcx_dens
        !+ Calculate Direct Charge Exchange (DCX) density: 0 = off, 1=on
    integer(Int32) :: calc_halo_dens
        !+ Calculate Thermal Halo density: 0 = off, 1=on
    integer(Int32) :: calc_brems
        !+ Calculate bremmstruhlung: 0 = off, 1=on
    integer(Int32) :: calc_bes
        !+ Calculate NBI: 0 = off, 1=on
    integer(Int32) :: calc_dcx
        !+ Calculate DCX: 0 = off, 1=on
    integer(Int32) :: calc_halo
        !+ Calculate Halo: 0 = off, 1=on
    integer(Int32) :: calc_cold
        !+ Calculate Cold D-alpha: 0 = off, 1=on
    integer(Int32) :: calc_fida
        !+ Calculate Active FIDA: 0 = off, 1=on
    integer(Int32) :: calc_pfida
        !+ Calculate Passive FIDA: 0 = off, 1=on
    integer(Int32) :: tot_spectra
        !+ Total number of spectral switches on
    integer(Int32) :: write_neutrals
        !+ Write neutrals to file: 0 = off, 1=on
    integer(Int32) :: load_neutrals
        !+ Load neutrals from file: 0 = off, 1=on
    integer(Int32) :: output_neutral_reservoir
        !+ Output neutral reservoir: 0 = off, 1=on
    integer(Int32) :: calc_npa
        !+ Calculate Active NPA: 0 = off, 1=on, 2=on++
    integer(Int32) :: calc_pnpa
        !+ Calculate Passive NPA: 0 = off, 1=on, 2=on++
    integer(Int32) :: calc_fida_wght
        !+ Calculate FIDA weight: 0 = off, 1=on, 2=on++
    integer(Int32) :: calc_npa_wght
        !+ Calculate NPA weights: 0 = off, 1=on, 2=on++
    integer(Int32) :: calc_birth
        !+ Calculate birth profile: 0 = off, 1=on
    integer(Int32) :: calc_neutron
        !+ Calculate neutron flux: 0 = off, 1=on, 2=on++
    integer(Int32) :: calc_cfpd
        !+ Calculate Charged Fusion Product flux: 0 = off, 1=on
    integer(Int32) :: calc_res
        !+ Calculate spatial resolution of LOS
    integer(Int32) :: flr
        !+ FLR correction: 0=off, 1=1st order(vxb/omega), 2=2nd order correction
    integer(Int32) :: split
        !+ Split signals by fast ion class: 0=off, 1=on
    integer(Int32) :: verbose
        !+ Verbosity: <0 = off++, 0 = off, 1=on, 2=on++
    integer(Int32) :: stark_components
        !+ Output spectral stark components : 0=off, 1=on

    !! Distribution settings
    integer(Int32) :: dist_type
        !+ Type of fast-ion distribution

    !! Spectrum parameters
    integer(Int32) :: nlambda
        !+ Number of wavelength to calculate
    real(Float64)  :: dlambda
        !+ Wavelength spacing [nm]
    real(Float64)  :: lambdamin
        !+ Minimum wavelength [nm]
    real(Float64)  :: lambdamax
        !+ Maximum wavelength [nm]

    !! Weight function settings
    integer(Int32) :: ne_wght
        !+ Number of energies in weight functions
    integer(Int32) :: np_wght
        !+ Number of pitches in weight functions
    integer(Int32) :: nphi_wght
        !+ Number of gyro-angles to average over in weight functions
    integer(Int32) :: nlambda_wght
        !+ Number of wavelength to calculate in weight functions
    real(Float64)  :: emax_wght
        !+ Maximum energy in weight functions [keV]
    real(Float64)  :: lambdamin_wght
        !+ Minimum wavelength in weight functions [nm]
    real(Float64)  :: lambdamax_wght
        !+ Maximum wavelength in weight functions [nm]
        
    !! Adaptive time step settings
    integer(Int32) :: adaptive
        !+ Simulation switch for adaptive time step, 0:split off, 1:dene, 2:avg(denn(1,:)), 3:denf, 4:avg(deni), 5:denimp, 6:te, 7:ti
    integer(Int32) :: max_cell_splits
        !+ Maximum number of times a cell can be split
    real(Float64)  :: split_tol
        !+ Tolerance level for splitting cells
end type SimulationInputs

type ParticleTrack
    !+ Stores properties seen when traveling through a 3D grid
    real(Float64) :: time = 0.d0
        !+ Time/distance/... in cell
    real(Float64) :: flux = 0.d0
        !+ Flux/density/... in cell
    integer(Int32), dimension(3) :: ind = 0
        !+ Indices of cell
    real(Float64), dimension(3)  :: pos = 0.d0
        !+ Midpoint of track in cell [cm]
end type ParticleTrack

type GyroSurface
    !+ Surface containing the fast-ion velocity vectors for all values of the
    !+ gyro-angle. It takes the form of a hyperboloid
    !+ \((x(\gamma,t) = \alpha \sqrt{1-\rm{pitch}^2}(cos(\gamma + \pi/2) - \omega_i t sin(\gamma + \pi/2)) \)
    !+ \((y(\gamma,t) = \alpha \sqrt{1-\rm{pitch}^2}(sin(\gamma + \pi/2) + \omega_i t cos(\gamma + \pi/2)) \)
    !+ \((z(\gamma,t) = \alpha \omega_i \rm{pitch} t\)
    !+ where \(\gamma\) is the gyro-angle, \(\omega_i\) is the ion
    !+ gyro-frequency and \(\alpha = V/\omega_i \)
    real(Float64) :: v = 0.d0
        !+ Particle speed
    real(Float64) :: omega = 0.d0
        !+ Ion gyro-frequency
    real(Float64), dimension(3)   :: axes = 0.d0
        !+ Semi-axes of the hyperboloid, i.e. a, b, c coefficients
    real(Float64), dimension(3)   :: center = 0.d0
        !+ Center of the gyrosurface
    real(Float64), dimension(3,3) :: A = 0.d0
        !+ Coefficients of quartic surface i.e. `basis*diagm(1/a^2,1/b^2,1/c^2)*basis'`
    real(Float64), dimension(3,3) :: basis = 0.d0
        !+ Basis of coordinate system of gyrosurface
end type GyroSurface

type CFPDTable
    !+ Defines Charged Fusion Product orbit trajectories
    integer :: nenergy = 0
        !+ Number of cfpd energies
    integer :: nrays = 0
        !+ Number of "rays"
    integer :: nsteps = 0
        !+ Number of total steps
    integer :: nchan = 0
        !+ Number of channels
    real(Float64) :: dl
        !+ Step length [cm]
    real(Float64) :: dE
        !+ CFPD energy bin width [keV]
    real(Float64), dimension(:), allocatable :: earray
        !+ Energies of cfpd orbits [keV]: earray(E3)
    real(Float64), dimension(:,:,:), allocatable :: nactual
        !+ Number of spatial steps: nactual(E3,ray,channel)
    real(Float64), dimension(:,:,:), allocatable :: daomega
        !+ Differntial area times solid angle [cm^2]: daomega(E3,ray,channel)
    real(Float64), dimension(:,:,:,:,:), allocatable :: sightline
        !+ Velocity [cm/s] and position [cm] in cylindrical (R,Phi,Z) coordinates:
        !+ sightline(E3,:,step,ray,channel)
end type CFPDTable

interface assignment(=)
    !+ Allows for assigning [[Profiles]],[[LocalProfiles]],
    !+ [[EMFields]],[[LocalEMFields]],[[FastIon]], [[NPAParticle]], and [[BirthParticle]]
    module procedure pp_assign, lpp_assign, plp_assign, lplp_assign, &
                     ff_assign, lff_assign, flf_assign, lflf_assign, &
                     fast_ion_assign,npa_part_assign,birth_part_assign, &
                     neutral_part_assign, npr_assign, np_assign, oo_assign
end interface

interface operator(+)
    !+ Allows for adding [[Profiles]],[[LocalProfiles]],
    !+ [[EMFields]], and [[LocalEMFields]]
    module procedure pp_add,lplp_add,ff_add,lflf_add
end interface

interface operator(-)
    !+ Allows for subtracting [[Profiles]],[[LocalProfiles]],
    !+ [[EMFields]], and [[LocalEMFields]]
    module procedure pp_subtract,lplp_subtract,ff_subtract,lflf_subtract
end interface

interface operator(*)
    !+ Allows for multiplying [[Profiles]],[[LocalProfiles]],
    !+ [[EMFields]], and [[LocalEMFields]] by scalars
    module procedure sp_multiply, ps_multiply, lps_multiply, slp_multiply, &
                     sf_multiply, fs_multiply, lfs_multiply, slf_multiply
end interface

interface operator(/)
    !+ Allows for dividing [[Profiles]],[[LocalProfiles]],
    !+ [[EMFields]], and [[LocalEMFields]] by scalars
    module procedure ps_divide, lps_divide, fs_divide, lfs_divide
end interface

!! definition of the structures:
type(BeamGrid), save            :: beam_grid
    !+ Variable containing beam grid definition
type(InterpolationGrid), save   :: inter_grid
    !+ Variable containing interpolation grid definition
type(InterpolationGrid), save   :: pass_grid
    !+ Variable containing passive neutral grid definition
type(FastIonDistribution), save :: fbm
    !+ Variable containing the fast-ion distribution function
type(FastIonParticles), save    :: particles
    !+ Variable containing a MC fast-ion distribution
type(Equilibrium), save         :: equil
    !+ Variable containing the plasma parameters and fields
type(NeutralBeam), save         :: nbi
    !+ Variable containing the neutral beam geometry and settings
type(AtomicTables), save        :: tables
    !+ Variable containing the atomic tables
type(NPAResults), save          :: npa
    !+ Variable for storing the calculated active NPA results
type(NPAResults), save          :: pnpa
    !+ Variable for storing the calculated passive NPA results
type(SpectralChords), save      :: spec_chords
    !+ Variable containing the spectral system definition
type(NPAChords), save           :: npa_chords
    !+ Variable containing the NPA system definition
type(SimulationInputs), save    :: inputs
    !+ Variable containing the simulation inputs
type(BirthProfile), save        :: birth
    !+ Variable for storing the calculated birth profile
type(Neutrals), save            :: neut
    !+ Variable for storing the calculated beam density
type(Spectra), save             :: spec
    !+ Variable for storing the calculated spectra
type(NeutronRate), save         :: neutron
    !+ Variable for storing the neutron rate
type(CFPDRate), save            :: cfpd
    !+ Variable for storing the Charged Fusion Product rate
type(FIDAWeights), save         :: fweight
    !+ Variable for storing the calculated FIDA weights
type(NPAWeights), save          :: nweight
    !+ Variable for storing the calculated NPA weights
type(CFPDTable), save           :: ctable
    !+ Variable for storing the calculated Charged Fusion Product orbits
type(SpatialSpectra), save      :: spatres
    !+ Variable for storing birth neutral for spatial resolution

contains

subroutine print_banner()
    !+ Prints FIDASIM banner
    write(*,'(a)') "   ____ ____ ___   ___    ____ ____ __  ___"
    write(*,'(a)') "  / __//  _// _ \ / _ |  / __//  _//  |/  /"
    write(*,'(a)') " / _/ _/ / / // // __ | _\ \ _/ / / /|_/ / "
    write(*,'(a)') "/_/  /___//____//_/ |_|/___//___//_/  /_/  "
    write(*,'(a)') "                                           "
    if(version.ne."") then
        write(*,'(a,a)') "Version: ",trim(version)
    endif
    write(*,'(a)') ""
    write(*,'(a)') "FIDASIM is released as open source code under the MIT Licence."
    write(*,'(a)') "For more information visit http://d3denergetic.github.io/FIDASIM/"
    write(*,'(a)') ""

#ifdef _DEBUG
    write(*,'(a)') "########################### ATTENTION ###########################"
    write(*,'(a)') "# Running in debug mode. All optimizations have been turned off #"
    write(*,'(a)') "#################################################################"
    write(*,'(a)') ""
#endif

#ifdef _PROF
    write(*,'(a)') "########################### ATTENTION ###########################"
    write(*,'(a)') "#                   Running in profiling mode                   #"
    write(*,'(a)') "#################################################################"
    write(*,'(a)') ""
#endif

#ifndef _OMP
#ifndef _MPI
    write(*,'(a)') "########################### ATTENTION ###########################"
    write(*,'(a)') "#              OpenMP threading has been disabled               #"
    write(*,'(a)') "#################################################################"
    write(*,'(a)') ""
#endif
#endif

#ifndef _MPI
#ifndef _OMP
    write(*,'(a)') "########################### ATTENTION ###########################"
    write(*,'(a)') "#                     MPI has been disabled                     #"
    write(*,'(a)') "#################################################################"
    write(*,'(a)') ""
#endif
#endif

end subroutine print_banner

!============================================================================
!---------------------------Operator Overloading-----------------------------
!============================================================================
subroutine fast_ion_assign(p1, p2)
    !+ Defines how to assign [[FastIon]] types to eachother
    type(FastIon), intent(in)  :: p2
    type(FastIon), intent(out) :: p1

    p1%beam_grid_cross_grid = p2%beam_grid_cross_grid
    p1%r                    = p2%r
    p1%z                    = p2%z
    p1%phi                  = p2%phi
    p1%beam_grid_phi_enter  = p2%beam_grid_phi_enter
    p1%delta_phi            = p2%delta_phi
    p1%energy               = p2%energy
    p1%pitch                = p2%pitch
    p1%vabs                 = p2%vabs
    p1%vr                   = p2%vr
    p1%vt                   = p2%vt
    p1%vz                   = p2%vz
    p1%weight               = p2%weight
    p1%class                = p2%class

end subroutine fast_ion_assign

subroutine npa_part_assign(p1, p2)
    !+ Defines how to assign [[NPAParticle]] types to eachother
    type(NPAParticle), intent(in)  :: p2
    type(NPAParticle), intent(out) :: p1

    p1%xi = p2%xi
    p1%yi = p2%yi
    p1%zi = p2%zi
    p1%xf = p2%xf
    p1%yf = p2%yf
    p1%zf = p2%zf
    p1%weight = p2%weight
    p1%energy = p2%energy
    p1%pitch = p2%pitch
    p1%detector = p2%detector

end subroutine npa_part_assign

subroutine birth_part_assign(p1, p2)
    !+ Defines how to assign [[BirthParticle]] types to eachother
    type(BirthParticle), intent(in)  :: p2
    type(BirthParticle), intent(out) :: p1

    p1%neut_type = p2%neut_type
    p1%ind       = p2%ind
    p1%ri        = p2%ri
    p1%vi        = p2%vi
    p1%ri_gc     = p2%ri_gc
    p1%weight    = p2%weight
    p1%energy    = p2%energy
    p1%pitch     = p2%pitch

end subroutine birth_part_assign

subroutine neutral_part_assign(p1, p2)
    !+ Defines how to assign [[NeutralParticle]] types to eachother
    type(NeutralParticle), intent(in)  :: p2
    type(NeutralParticle), intent(out) :: p1

    p1%w = p2%w
    p1%v = p2%v

end subroutine neutral_part_assign

subroutine npr_assign(r1,r2)
    !+ Defines how to assign [[NeutralParticleReservoir]] to eachother
    type(NeutralParticleReservoir), intent(in) :: r2
    type(NeutralParticleReservoir), intent(inout) :: r1

    if(r2%k.eq.0) return
    call init_reservoir(r1) !r2%k)
    r1%n = r2%n
    r1%k = r2%k
    r1%i = r2%i
    r1%W = r2%W
    r1%R = r2%R

end subroutine npr_assign

subroutine np_assign(n1, n2)
    !+ Defines how to assign [[NeutralPopulation]] to eachother
    type(NeutralPopulation), intent(in) :: n2
    type(NeutralPopulation), intent(inout) :: n1

    integer :: i, j, k, ic, ind(3), nx, ny, nz

    if(.not.allocated(n1%dens)) then
        nx = size(n2%dens,2); ny = size(n2%dens,3); nz = size(n2%dens,4)
        allocate(n1%dens(nlevs,nx,ny,nz), source=n2%dens)
    else
        n1%dens = n2%dens
    endif

    do ic=1,beam_grid%ngrid
        call ind2sub(beam_grid%dims,ic,ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)
        if (n2%res(i,j,k)%n.ne.0) then
            call npr_assign(n1%res(i,j,k), n2%res(i,j,k))
        endif
    enddo

end subroutine np_assign

pure subroutine pp_assign(p1, p2)
    !+ Defines how to assign [[Profiles]] types to eachother
    type(Profiles), intent(in)    :: p2
    type(Profiles), intent(inout) :: p1

    p1%dene   = p2%dene
    p1%ti     = p2%ti
    p1%te     = p2%te
    p1%deni   = p2%deni
    p1%denf   = p2%denf
    p1%denimp = p2%denimp
    p1%zeff   = p2%zeff
    p1%vr     = p2%vr
    p1%vt     = p2%vt
    p1%vz     = p2%vz
    p1%denn   = p2%denn

end subroutine pp_assign

pure subroutine lpp_assign(p1, p2)
    !+ Defines how to assign a [[Profiles]] type to a [[LocalProfiles]] type
    type(Profiles), intent(in)         :: p2
    type(LocalProfiles), intent(inout) :: p1

    p1%dene   = p2%dene
    p1%ti     = p2%ti
    p1%te     = p2%te
    p1%deni   = p2%deni
    p1%denf   = p2%denf
    p1%denimp = p2%denimp
    p1%zeff   = p2%zeff
    p1%vr     = p2%vr
    p1%vt     = p2%vt
    p1%vz     = p2%vz
    p1%denn   = p2%denn

end subroutine lpp_assign

pure subroutine plp_assign(p1, p2)
    !+ Defines how to assign a [[LocalProfiles]] type to a [[Profiles]] type
    type(LocalProfiles), intent(in) :: p2
    type(Profiles), intent(inout)   :: p1

    p1%dene   = p2%dene
    p1%ti     = p2%ti
    p1%te     = p2%te
    p1%deni   = p2%deni
    p1%denf   = p2%denf
    p1%denimp = p2%denimp
    p1%zeff   = p2%zeff
    p1%vr     = p2%vr
    p1%vt     = p2%vt
    p1%vz     = p2%vz
    p1%denn   = p2%denn

end subroutine plp_assign

pure subroutine lplp_assign(p1, p2)
    !+ Defines how to assign [[LocalProfiles]] types to eachother
    type(LocalProfiles), intent(in)    :: p2
    type(LocalProfiles), intent(inout) :: p1

    p1%pos    = p2%pos
    p1%uvw    = p2%uvw
    p1%dene   = p2%dene
    p1%ti     = p2%ti
    p1%te     = p2%te
    p1%deni   = p2%deni
    p1%denf   = p2%denf
    p1%denimp = p2%denimp
    p1%zeff   = p2%zeff
    p1%vr     = p2%vr
    p1%vt     = p2%vt
    p1%vz     = p2%vz
    p1%denn   = p2%denn
    p1%vrot   = p2%vrot

end subroutine lplp_assign

pure subroutine ff_assign(p1, p2)
    !+ Defines how to assign [[EMFields]] types to eachother
    type(EMFields), intent(in)    :: p2
    type(EMFields), intent(inout) :: p1

    p1%br   = p2%br
    p1%bt   = p2%bt
    p1%bz   = p2%bz
    p1%er   = p2%er
    p1%et   = p2%et
    p1%ez   = p2%ez

    p1%dbr_dr   = p2%dbr_dr
    p1%dbr_dz   = p2%dbr_dz
    p1%dbt_dr   = p2%dbt_dr
    p1%dbt_dz   = p2%dbt_dz
    p1%dbz_dr   = p2%dbz_dr
    p1%dbz_dz   = p2%dbz_dz

end subroutine ff_assign

pure subroutine lff_assign(p1, p2)
    !+ Defines how to assign a [[EMFields]] type to a [[LocalEMFields]] type
    type(EMFields), intent(in)         :: p2
    type(LocalEMFields), intent(inout) :: p1

    p1%br   = p2%br
    p1%bt   = p2%bt
    p1%bz   = p2%bz
    p1%er   = p2%er
    p1%et   = p2%et
    p1%ez   = p2%ez

    p1%dbr_dr   = p2%dbr_dr
    p1%dbr_dz   = p2%dbr_dz
    p1%dbt_dr   = p2%dbt_dr
    p1%dbt_dz   = p2%dbt_dz
    p1%dbz_dr   = p2%dbz_dr
    p1%dbz_dz   = p2%dbz_dz

end subroutine lff_assign

pure subroutine flf_assign(p1, p2)
    !+ Defines how to assign a [[LocalEMFields]] type to a [[EMFields]] type
    type(LocalEMFields), intent(in) :: p2
    type(EMFields), intent(inout)   :: p1

    p1%br   = p2%br
    p1%bt   = p2%bt
    p1%bz   = p2%bz
    p1%er   = p2%er
    p1%et   = p2%et
    p1%ez   = p2%ez

    p1%dbr_dr   = p2%dbr_dr
    p1%dbr_dz   = p2%dbr_dz
    p1%dbt_dr   = p2%dbt_dr
    p1%dbt_dz   = p2%dbt_dz
    p1%dbz_dr   = p2%dbz_dr
    p1%dbz_dz   = p2%dbz_dz

end subroutine flf_assign

pure subroutine lflf_assign(p1, p2)
    !+ Defines how to assign [[LocalEMFields]] types to eachother
    type(LocalEMFields), intent(in)    :: p2
    type(LocalEMFields), intent(inout) :: p1

    p1%pos  = p2%pos
    p1%uvw  = p2%uvw
    p1%br   = p2%br
    p1%bt   = p2%bt
    p1%bz   = p2%bz
    p1%er   = p2%er
    p1%et   = p2%et
    p1%ez   = p2%ez
    p1%b_abs = p2%b_abs
    p1%e_abs = p2%e_abs
    p1%a_norm = p2%a_norm
    p1%b_norm = p2%b_norm
    p1%c_norm = p2%c_norm
    p1%e_norm = p2%e_norm

    p1%dbr_dr   = p2%dbr_dr
    p1%dbr_dz   = p2%dbr_dz
    p1%dbt_dr   = p2%dbt_dr
    p1%dbt_dz   = p2%dbt_dz
    p1%dbz_dr   = p2%dbz_dr
    p1%dbz_dz   = p2%dbz_dz

end subroutine lflf_assign

elemental function pp_add(p1, p2) result (p3)
    !+ Defines how to add two [[Profiles]] types
    type(Profiles), intent(in) :: p1,p2
    type(Profiles)             :: p3

    p3%dene   = p1%dene   + p2%dene
    p3%ti     = p1%ti     + p2%ti
    p3%te     = p1%te     + p2%te
    p3%deni   = p1%deni   + p2%deni
    p3%denf   = p1%denf   + p2%denf
    p3%denimp = p1%denimp + p2%denimp
    p3%zeff   = p1%zeff   + p2%zeff
    p3%vr     = p1%vr     + p2%vr
    p3%vt     = p1%vt     + p2%vt
    p3%vz     = p1%vz     + p2%vz
    p3%denn   = p1%denn   + p2%denn

end function pp_add

elemental function pp_subtract(p1, p2) result (p3)
    !+ Defines how to subtract two [[Profiles]] types
    type(Profiles), intent(in) :: p1,p2
    type(Profiles)             :: p3

    p3%dene   = p1%dene   - p2%dene
    p3%ti     = p1%ti     - p2%ti
    p3%te     = p1%te     - p2%te
    p3%deni   = p1%deni   - p2%deni
    p3%denf   = p1%denf   - p2%denf
    p3%denimp = p1%denimp - p2%denimp
    p3%zeff   = p1%zeff   - p2%zeff
    p3%vr     = p1%vr     - p2%vr
    p3%vt     = p1%vt     - p2%vt
    p3%vz     = p1%vz     - p2%vz
    p3%denn   = p1%denn   - p2%denn

end function pp_subtract

elemental function lplp_add(p1, p2) result (p3)
    !+ Defines how to add two [[LocalProfiles]] types
    type(LocalProfiles), intent(in) :: p1,p2
    type(LocalProfiles)             :: p3

    p3%pos    = p1%pos    + p2%pos
    p3%uvw    = p1%uvw    + p2%uvw
    p3%dene   = p1%dene   + p2%dene
    p3%ti     = p1%ti     + p2%ti
    p3%te     = p1%te     + p2%te
    p3%deni   = p1%deni   + p2%deni
    p3%denf   = p1%denf   + p2%denf
    p3%denimp = p1%denimp + p2%denimp
    p3%zeff   = p1%zeff   + p2%zeff
    p3%vr     = p1%vr     + p2%vr
    p3%vt     = p1%vt     + p2%vt
    p3%vz     = p1%vz     + p2%vz
    p3%denn   = p1%denn   + p2%denn
    p3%vrot   = p1%vrot   + p2%vrot

end function lplp_add

elemental function lplp_subtract(p1, p2) result (p3)
    !+ Defines how to subtract two [[LocalProfiles]] types
    type(LocalProfiles), intent(in) :: p1,p2
    type(LocalProfiles)             :: p3

    p3%pos    = p1%pos    - p2%pos
    p3%uvw    = p1%uvw    - p2%uvw
    p3%dene   = p1%dene   - p2%dene
    p3%ti     = p1%ti     - p2%ti
    p3%te     = p1%te     - p2%te
    p3%deni   = p1%deni   - p2%deni
    p3%denf   = p1%denf   - p2%denf
    p3%denimp = p1%denimp - p2%denimp
    p3%zeff   = p1%zeff   - p2%zeff
    p3%vr     = p1%vr     - p2%vr
    p3%vt     = p1%vt     - p2%vt
    p3%vz     = p1%vz     - p2%vz
    p3%denn   = p1%denn   - p2%denn
    p3%vrot   = p1%vrot   - p2%vrot

end function lplp_subtract

elemental function ps_multiply(p1, real_scalar) result (p3)
    !+ Defines how to multiply [[Profiles]] types by a scalar
    type(Profiles), intent(in) :: p1
    real(Float64), intent(in)  :: real_scalar
    type(Profiles)             :: p3

    p3%dene   = p1%dene   * real_scalar
    p3%ti     = p1%ti     * real_scalar
    p3%te     = p1%te     * real_scalar
    p3%deni   = p1%deni   * real_scalar
    p3%denf   = p1%denf   * real_scalar
    p3%denimp = p1%denimp * real_scalar
    p3%zeff   = p1%zeff   * real_scalar
    p3%vr     = p1%vr     * real_scalar
    p3%vt     = p1%vt     * real_scalar
    p3%vz     = p1%vz     * real_scalar
    p3%denn   = p1%denn   * real_scalar

end function ps_multiply

elemental function sp_multiply(real_scalar, p1) result (p3)
    !+ Defines how to multiply [[Profiles]] types by a scalar
    type(Profiles), intent(in) :: p1
    real(Float64), intent(in)  :: real_scalar
    type(Profiles)             :: p3

    p3 = p1*real_scalar

end function sp_multiply

elemental function ps_divide(p1, real_scalar) result (p3)
    !+ Defines how to divide [[Profiles]] types by a scalar
    type(Profiles), intent(in) :: p1
    real(Float64), intent(in)  :: real_scalar
    type(Profiles)             :: p3

    p3 = p1*(1.d0/real_scalar)

end function ps_divide

elemental function lps_multiply(p1, real_scalar) result (p3)
    !+ Defines how to multiply [[LocalProfiles]] types by a scalar
    type(LocalProfiles), intent(in) :: p1
    real(Float64), intent(in)       :: real_scalar
    type(LocalProfiles)             :: p3

    p3%pos    = p1%pos    * real_scalar
    p3%uvw    = p1%uvw    * real_scalar
    p3%dene   = p1%dene   * real_scalar
    p3%ti     = p1%ti     * real_scalar
    p3%te     = p1%te     * real_scalar
    p3%deni   = p1%deni   * real_scalar
    p3%denf   = p1%denf   * real_scalar
    p3%denimp = p1%denimp * real_scalar
    p3%zeff   = p1%zeff   * real_scalar
    p3%vr     = p1%vr     * real_scalar
    p3%vt     = p1%vt     * real_scalar
    p3%vz     = p1%vz     * real_scalar
    p3%denn   = p1%denn   * real_scalar
    p3%vrot   = p1%vrot   * real_scalar

end function lps_multiply

elemental function slp_multiply(real_scalar, p1) result (p3)
    !+ Defines how to multiply [[LocalProfiles]] types by a scalar
    type(LocalProfiles), intent(in) :: p1
    real(Float64), intent(in)       :: real_scalar
    type(LocalProfiles)             :: p3

    p3 = p1*real_scalar

end function slp_multiply

elemental function lps_divide(p1, real_scalar) result (p3)
    !+ Defines how to divide [[LocalProfiles]] types by a scalar
    type(LocalProfiles), intent(in) :: p1
    real(Float64), intent(in)       :: real_scalar
    type(LocalProfiles)             :: p3

    p3 = p1*(1.d0/real_scalar)

end function lps_divide

elemental function ff_add(p1, p2) result (p3)
    !+ Defines how to add two [[EMFields]] types
    type(EMFields), intent(in) :: p1,p2
    type(EMFields)             :: p3

    p3%br   = p1%br   + p2%br
    p3%bt   = p1%bt   + p2%bt
    p3%bz   = p1%bz   + p2%bz
    p3%er   = p1%er   + p2%er
    p3%et   = p1%et   + p2%et
    p3%ez   = p1%ez   + p2%ez

    p3%dbr_dr   = p1%dbr_dr + p2%dbr_dr
    p3%dbr_dz   = p1%dbr_dz + p2%dbr_dz
    p3%dbt_dr   = p1%dbt_dr + p2%dbt_dr
    p3%dbt_dz   = p1%dbt_dz + p2%dbt_dz
    p3%dbz_dr   = p1%dbz_dr + p2%dbz_dr
    p3%dbz_dz   = p1%dbz_dz + p2%dbz_dz

end function ff_add

elemental function ff_subtract(p1, p2) result (p3)
    !+ Defines how to subtract two [[EMFields]] types
    type(EMFields), intent(in) :: p1,p2
    type(EMFields)             :: p3

    p3%br   = p1%br   - p2%br
    p3%bt   = p1%bt   - p2%bt
    p3%bz   = p1%bz   - p2%bz
    p3%er   = p1%er   - p2%er
    p3%et   = p1%et   - p2%et
    p3%ez   = p1%ez   - p2%ez

    p3%dbr_dr   = p1%dbr_dr - p2%dbr_dr
    p3%dbr_dz   = p1%dbr_dz - p2%dbr_dz
    p3%dbt_dr   = p1%dbt_dr - p2%dbt_dr
    p3%dbt_dz   = p1%dbt_dz - p2%dbt_dz
    p3%dbz_dr   = p1%dbz_dr - p2%dbz_dr
    p3%dbz_dz   = p1%dbz_dz - p2%dbz_dz

end function ff_subtract

elemental function fs_multiply(p1, real_scalar) result (p3)
    !+ Defines how to multiply [[EMFields]] types by a scalar
    type(EMFields), intent(in) :: p1
    real(Float64), intent(in)  :: real_scalar
    type(EMFields)             :: p3

    p3%br   = p1%br   * real_scalar
    p3%bt   = p1%bt   * real_scalar
    p3%bz   = p1%bz   * real_scalar
    p3%er   = p1%er   * real_scalar
    p3%et   = p1%et   * real_scalar
    p3%ez   = p1%ez   * real_scalar

    p3%dbr_dr   = p1%dbr_dr * real_scalar
    p3%dbr_dz   = p1%dbr_dz * real_scalar
    p3%dbt_dr   = p1%dbt_dr * real_scalar
    p3%dbt_dz   = p1%dbt_dz * real_scalar
    p3%dbz_dr   = p1%dbz_dr * real_scalar
    p3%dbz_dz   = p1%dbz_dz * real_scalar

end function fs_multiply

elemental function sf_multiply(real_scalar, p1) result (p3)
    !+ Defines how to multiply [[EMFields]] types by a scalar
    type(EMFields), intent(in) :: p1
    real(Float64), intent(in)  :: real_scalar
    type(EMFields)             :: p3

    p3 = p1*real_scalar

end function sf_multiply

elemental function fs_divide(p1, real_scalar) result (p3)
    !+ Defines how to divide [[EMFields]] types by a scalar
    type(EMFields), intent(in) :: p1
    real(Float64), intent(in)  :: real_scalar
    type(EMFields)             :: p3

    p3 = p1*(1.d0/real_scalar)

end function fs_divide

elemental function lflf_add(p1, p2) result (p3)
    !+ Defines how to add two [[LocalEMFields]] types
    type(LocalEMFields), intent(in) :: p1,p2
    type(LocalEMFields)             :: p3

    real(Float64), dimension(3) :: bfield,efield

    p3%pos    = p1%pos    + p2%pos
    p3%uvw    = p1%uvw    + p2%uvw
    p3%br     = p1%br     + p2%br
    p3%bt     = p1%bt     + p2%bt
    p3%bz     = p1%bz     + p2%bz
    p3%er     = p1%er     + p2%er
    p3%et     = p1%et     + p2%et
    p3%ez     = p1%ez     + p2%ez

    bfield = p1%b_abs*p1%b_norm + p2%b_abs*p2%b_norm
    p3%b_abs = norm2(bfield)
    if(p3%b_abs.gt.0.d0) then
        p3%b_norm = bfield/p3%b_abs
        call calc_perp_vectors(p3%b_norm,p3%a_norm,p3%c_norm)
    endif

    efield = p1%e_abs*p1%e_norm + p2%e_abs*p2%e_norm
    p3%e_abs = norm2(efield)
    if(p3%e_abs.gt.0.d0) p3%e_norm = efield/p3%e_abs

    p3%dbr_dr   = p1%dbr_dr + p2%dbr_dr
    p3%dbr_dz   = p1%dbr_dz + p2%dbr_dz
    p3%dbt_dr   = p1%dbt_dr + p2%dbt_dr
    p3%dbt_dz   = p1%dbt_dz + p2%dbt_dz
    p3%dbz_dr   = p1%dbz_dr + p2%dbz_dr
    p3%dbz_dz   = p1%dbz_dz + p2%dbz_dz

end function lflf_add

elemental function lflf_subtract(p1, p2) result (p3)
    !+ Defines how to subtract two [[LocalEMFields]] types
    type(LocalEMFields), intent(in) :: p1,p2
    type(LocalEMFields)             :: p3

    real(Float64), dimension(3) :: bfield,efield

    p3%pos    = p1%pos    - p2%pos
    p3%uvw    = p1%uvw    - p2%uvw
    p3%br     = p1%br     - p2%br
    p3%bt     = p1%bt     - p2%bt
    p3%bz     = p1%bz     - p2%bz
    p3%er     = p1%er     - p2%er
    p3%et     = p1%et     - p2%et
    p3%ez     = p1%ez     - p2%ez

    bfield = p1%b_abs*p1%b_norm - p2%b_abs*p2%b_norm
    p3%b_abs = norm2(bfield)
    if(p3%b_abs.gt.0.d0) then
        p3%b_norm = bfield/p3%b_abs
        call calc_perp_vectors(p3%b_norm,p3%a_norm,p3%c_norm)
    endif

    efield = p1%e_abs*p1%e_norm - p2%e_abs*p2%e_norm
    p3%e_abs = norm2(efield)
    if(p3%e_abs.gt.0.d0) p3%e_norm = efield/p3%e_abs

    p3%dbr_dr   = p1%dbr_dr - p2%dbr_dr
    p3%dbr_dz   = p1%dbr_dz - p2%dbr_dz
    p3%dbt_dr   = p1%dbt_dr - p2%dbt_dr
    p3%dbt_dz   = p1%dbt_dz - p2%dbt_dz
    p3%dbz_dr   = p1%dbz_dr - p2%dbz_dr
    p3%dbz_dz   = p1%dbz_dz - p2%dbz_dz

end function lflf_subtract

elemental function lfs_multiply(p1, real_scalar) result (p3)
    !+ Defines how to multiply [[LocalEMFields]] types by a scalar
    type(LocalEMFields), intent(in) :: p1
    real(Float64), intent(in)       :: real_scalar
    type(LocalEMFields)             :: p3

    p3%pos  = p1%pos  * real_scalar
    p3%uvw  = p1%uvw  * real_scalar
    p3%br   = p1%br   * real_scalar
    p3%bt   = p1%bt   * real_scalar
    p3%bz   = p1%bz   * real_scalar
    p3%er   = p1%er   * real_scalar
    p3%et   = p1%et   * real_scalar
    p3%ez   = p1%ez   * real_scalar
    p3%b_abs  = p1%b_abs  * real_scalar
    p3%e_abs  = p1%e_abs  * real_scalar
    p3%a_norm = p1%a_norm
    p3%b_norm = p1%b_norm
    p3%c_norm = p1%c_norm
    p3%e_norm = p1%e_norm

    p3%dbr_dr   = p1%dbr_dr * real_scalar
    p3%dbr_dz   = p1%dbr_dz * real_scalar
    p3%dbt_dr   = p1%dbt_dr * real_scalar
    p3%dbt_dz   = p1%dbt_dz * real_scalar
    p3%dbz_dr   = p1%dbz_dr * real_scalar
    p3%dbz_dz   = p1%dbz_dz * real_scalar

end function lfs_multiply

elemental function slf_multiply(real_scalar, p1) result (p3)
    !+ Defines how to multiply [[LocalEMFields]] types by a scalar
    type(LocalEMFields), intent(in) :: p1
    real(Float64), intent(in)       :: real_scalar
    type(LocalEMFields)             :: p3

    p3 = p1*real_scalar

end function slf_multiply

elemental function lfs_divide(p1, real_scalar) result (p3)
    !+ Defines how to divide [[LocalEMFields]] types by a scalar
    type(LocalEMFields), intent(in) :: p1
    real(Float64), intent(in)       :: real_scalar
    type(LocalEMFields)             :: p3

    p3 = p1*(1.d0/real_scalar)

end function lfs_divide

pure subroutine oo_assign(p1, p2)
    !+ Defines how to assign [[CFPDTable]] types to eachother
    type(CFPDTable), intent(in)  :: p2
    type(CFPDTable), intent(out) :: p1

    p1%nenergy   = p2%nenergy
    p1%nrays     = p2%nrays
    p1%nsteps    = p2%nsteps
    p1%nchan     = p2%nchan
    p1%dl        = p2%dl
    p1%dE        = p2%dE
    p1%earray    = p2%earray
    p1%nactual   = p2%nactual
    p1%daomega   = p2%daomega
    p1%sightline = p2%sightline

end subroutine oo_assign

!============================================================================
!-------------------------------I/O Routines---------------------------------
!============================================================================
subroutine read_inputs
    !+ Reads input namelist file and stores the results into [[libfida:inputs]],
    !+ [[libfida:nbi]], and [[libfida:beam_grid]]
    character(charlim) :: runid,result_dir, tables_file
    character(charlim) :: distribution_file, equilibrium_file
    character(charlim) :: geometry_file, neutrals_file
    integer            :: pathlen, calc_neutron, seed, calc_cfpd
    integer            :: calc_brems, calc_dcx, calc_halo, calc_cold, calc_bes
    integer            :: calc_fida, calc_pfida, calc_npa, calc_pnpa
    integer            :: calc_birth,calc_fida_wght,calc_npa_wght, calc_res
    integer            :: write_neutrals,load_neutrals,verbose,flr,split,stark_components
    integer            :: output_neutral_reservoir
    integer(Int64)     :: n_fida,n_pfida,n_npa,n_pnpa,n_nbi,n_halo,n_dcx,n_birth
    integer(Int32)     :: shot,nlambda,ne_wght,np_wght,nphi_wght,nlambda_wght
    integer(Int32)     :: adaptive, max_cell_splits
    real(Float64)      :: time,lambdamin,lambdamax,emax_wght
    real(Float64)      :: lambdamin_wght,lambdamax_wght
    real(Float64)      :: ab,pinj,einj,current_fractions(3)
    integer(Int32)     :: nx,ny,nz
    real(Float64)      :: xmin,xmax,ymin,ymax,zmin,zmax
    real(Float64)      :: alpha,beta,gamma,origin(3)
    real(Float64)      :: split_tol
    logical            :: exis, error

    NAMELIST /fidasim_inputs/ result_dir, tables_file, distribution_file, &
        geometry_file, equilibrium_file, neutrals_file, shot, time, runid, &
        calc_brems, calc_dcx,calc_halo, calc_cold, calc_fida, calc_bes,&
        calc_pfida, calc_npa, calc_pnpa,calc_birth, calc_res, seed, flr, split, &
        calc_fida_wght, calc_npa_wght, write_neutrals, load_neutrals, verbose, stark_components, &
        calc_neutron, calc_cfpd, n_fida, n_pfida, n_npa, n_pnpa, n_nbi, n_halo, n_dcx, n_birth, &
        ab, pinj, einj, current_fractions, output_neutral_reservoir, &
        nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, &
        origin, alpha, beta, gamma, &
        ne_wght, np_wght, nphi_wght, &
        nlambda, lambdamin,lambdamax,emax_wght, &
        nlambda_wght,lambdamin_wght,lambdamax_wght, &
        adaptive, max_cell_splits, split_tol

    inquire(file=namelist_file,exist=exis)
    if(.not.exis) then
        write(*,'(a,a)') 'READ_INPUTS: Input file does not exist: ', trim(namelist_file)
        stop
    endif

    ! variables that are not changed not be auto-initalized
    ! provide reasonable defaults here
    result_dir="."
    tables_file="."
    distribution_file="."
    geometry_file ="."
    equilibrium_file="."
    neutrals_file="."
    shot=0
    time=0
    runid="0"
    seed = -1
    calc_brems=0
    calc_bes=0
    calc_dcx=0
    calc_halo=0
    calc_cold=0
    calc_fida=0
    calc_pfida=0
    calc_npa=0
    calc_pnpa=0
    calc_birth=0
    calc_res=0
    flr=2
    stark_components=0
    split=1
    calc_fida_wght=0
    calc_npa_wght=0
    write_neutrals=1
    load_neutrals=0
    output_neutral_reservoir=1
    verbose=0
    calc_neutron=0
    calc_cfpd=0
    n_fida=0
    n_pfida=0
    n_npa=0
    n_pnpa=0
    n_nbi=0
    n_halo=0
    n_dcx=0
    n_birth=0
    ab=0
    pinj=0
    einj=0
    current_fractions=0
    nx=0
    ny=0
    nz=0
    xmin=0
    xmax=0
    ymin=0
    ymax=0
    zmin=0
    zmax=0
    origin=0
    alpha=0
    beta=0
    gamma=0
    ne_wght=0
    np_wght=0
    nphi_wght=0
    nlambda=0
    lambdamin=0
    lambdamax=0
    emax_wght=0
    nlambda_wght=0
    lambdamin_wght=0
    lambdamax_wght=0
    adaptive=0
    max_cell_splits=1
    split_tol=0

    open(13,file=namelist_file)
    read(13,NML=fidasim_inputs)
    close(13)

    !!General Information
    inputs%shot_number=shot
    inputs%time=time
    inputs%runid=runid
    inputs%result_dir=result_dir

    !!Input Files
    inputs%tables_file=tables_file
    inputs%geometry_file=geometry_file
    inputs%equilibrium_file=equilibrium_file
    inputs%distribution_file=distribution_file
    inputs%neutrals_file=neutrals_file

    !! RNG seed
    inputs%seed = seed
    if(inputs%seed.lt.0) inputs%seed = rng_seed()

    !!Simulation Switches
    if((calc_brems+calc_bes+calc_dcx+calc_halo+&
        calc_cold+calc_fida+calc_pfida).gt.0) then
        inputs%calc_spec=1
        inputs%tot_spectra=calc_brems+calc_bes+calc_dcx+calc_halo+&
                           calc_cold+calc_fida+calc_pfida
    else
        inputs%calc_spec=0
        inputs%tot_spectra=0
    endif

    inputs%calc_beam = 0
    if((calc_bes+calc_birth+calc_dcx+&
        calc_halo+calc_fida+calc_npa+&
        calc_fida_wght+calc_npa_wght).gt.0) then
        inputs%calc_nbi_dens=1
        inputs%calc_beam=1
    else
        inputs%calc_nbi_dens=0
    endif

    if((calc_dcx+calc_halo+calc_fida+calc_npa+&
        calc_fida_wght+calc_npa_wght).gt.0) then
        inputs%calc_dcx_dens=1
        inputs%calc_beam=1
    else
        inputs%calc_dcx_dens=0
    endif

    if((calc_halo+calc_fida+calc_npa+&
        calc_fida_wght+calc_npa_wght).gt.0) then
        inputs%calc_halo_dens=1
        inputs%calc_beam=1
    else
        inputs%calc_halo_dens=0
    endif

    inputs%calc_brems=calc_brems
    inputs%calc_bes=calc_bes
    inputs%calc_dcx=calc_dcx
    inputs%calc_halo=calc_halo
    inputs%calc_cold=calc_cold
    inputs%calc_fida=calc_fida
    inputs%calc_pfida=calc_pfida
    inputs%calc_npa=calc_npa
    inputs%calc_pnpa=calc_pnpa
    inputs%calc_birth=calc_birth
    inputs%calc_fida_wght=calc_fida_wght
    inputs%calc_npa_wght=calc_npa_wght
    inputs%calc_neutron=calc_neutron
    inputs%calc_cfpd=calc_cfpd
    inputs%calc_res = calc_res

    !! Misc. Settings
    inputs%write_neutrals=write_neutrals
    inputs%load_neutrals=load_neutrals
    inputs%output_neutral_reservoir=output_neutral_reservoir

    inputs%verbose=verbose
    inputs%flr = flr
    inputs%stark_components = stark_components
    inputs%split = split

    !!Monte Carlo Settings
    inputs%n_fida=max(10,n_fida)
    inputs%n_pfida=max(10,n_pfida)
    inputs%n_npa=max(10,n_npa)
    inputs%n_pnpa=max(10,n_pnpa)
    inputs%n_nbi=max(10,n_nbi)
    inputs%n_halo=max(10,n_halo)
    inputs%n_dcx=max(10,n_dcx)
    inputs%n_birth= max(1,nint(n_birth/real(n_nbi)))

    !!Neutral Beam Settings
    beam_mass=ab
    nbi%current_fractions=current_fractions
    nbi%einj=einj
    nbi%pinj=pinj

    !!Weight Function Settings
    inputs%ne_wght=ne_wght
    inputs%np_wght=np_wght
    inputs%nphi_wght=nphi_wght
    inputs%emax_wght=emax_wght
    inputs%nlambda_wght = nlambda_wght
    inputs%lambdamin_wght=lambdamin_wght
    inputs%lambdamax_wght=lambdamax_wght

    !!Wavelength Grid Settings
    inputs%nlambda=nlambda
    inputs%lambdamin=lambdamin
    inputs%lambdamax=lambdamax
    inputs%dlambda=(inputs%lambdamax-inputs%lambdamin)/inputs%nlambda

    !!Adaptive Time Step Settings
    inputs%adaptive=adaptive
    if(inputs%adaptive.eq.0) then
        inputs%max_cell_splits=1
        inputs%split_tol=0.0
    else
        inputs%max_cell_splits=max_cell_splits
        inputs%split_tol=split_tol
    endif

    !!Beam Grid Settings
    beam_grid%nx=nx
    beam_grid%ny=ny
    beam_grid%nz=nz
    beam_grid%xmin=xmin
    beam_grid%xmax=xmax
    beam_grid%ymin=ymin
    beam_grid%ymax=ymax
    beam_grid%zmin=zmin
    beam_grid%zmax=zmax
    beam_grid%alpha=alpha
    beam_grid%beta=beta
    beam_grid%gamma=gamma
    beam_grid%origin=origin

#ifdef _MPI
    if(my_rank().ne.0) inputs%verbose=0
#endif

    if(inputs%verbose.ge.1) then
        write(*,'(a)') "---- Shot settings ----"
        write(*,'(T2,"Shot: ",i8)') inputs%shot_number
        write(*,'(T2,"Time: ",i4," [ms]")') int(inputs%time*1.d3)
        write(*,'(T2,"Runid: ",a)') trim(adjustl(inputs%runid))
        write(*,*) ''
        write(*,'(a)') "---- Input files ----"
    endif

    error = .False.

    inquire(file=inputs%tables_file,exist=exis)
    if(exis) then
        if(inputs%verbose.ge.1) then
            write(*,'(T2,"Tables file: ",a)') trim(inputs%tables_file)
        endif
    else
        if(inputs%verbose.ge.0) then
            write(*,'(a,a)') 'READ_INPUTS: Tables file does not exist: ', &
                             trim(inputs%tables_file)
        endif
        error = .True.
    endif

    inquire(file=inputs%geometry_file,exist=exis)
    if(exis) then
        if(inputs%verbose.ge.1) then
            write(*,'(T2,"Geometry file: ",a)') trim(inputs%geometry_file)
        endif
    else
        if(inputs%verbose.ge.0) then
            write(*,'(a,a)') 'READ_INPUTS: Geometry file does not exist: ', &
                             trim(inputs%geometry_file)
        endif
        error = .True.
    endif

    inquire(file=inputs%equilibrium_file,exist=exis)
    if(exis) then
        if(inputs%verbose.ge.1) then
            write(*,'(T2,"Equilibrium file: ",a)') trim(inputs%equilibrium_file)
        endif
    else
        if(inputs%verbose.ge.0) then
            write(*,'(a,a)') 'READ_INPUTS: Equilibrium file does not exist: ', &
                              trim(inputs%equilibrium_file)
        endif
        error = .True.
    endif

    inquire(file=inputs%distribution_file,exist=exis)
    if(exis) then
        if(inputs%verbose.ge.1) then
            write(*,'(T2,"Distribution file: ",a)') trim(inputs%distribution_file)
        endif
    else
        if(inputs%verbose.ge.0) then
            write(*,'(a,a)') 'READ_INPUTS: Distribution file does not exist: ', &
                             trim(inputs%distribution_file)
        endif
        error = .True.
    endif

    pathlen = len_trim(inputs%result_dir)+len_trim(inputs%runid) + 20
    !+20 for suffixes and seperators e.g. /, _npa.h5, ...
    if(pathlen.gt.charlim) then
        if(inputs%verbose.ge.0) then
            write(*,'(a,i3,a,i3)') 'READ_INPUTS: Result directory path + runID use too many characters: ', &
                                   pathlen-20,'>', charlim-20
        endif
        error = .True.
    endif

    if(inputs%verbose.ge.1) then
        write(*,*) ''
    endif

    if(inputs%adaptive.ge.0.and.inputs%adaptive.le.7) then
        write(*,'(a)') "---- Adaptive time step settings ----"
        write(*,'(T2,"Adaptive: ",i2)') inputs%adaptive
        if(inputs%adaptive.gt.0) then
            if(inputs%max_cell_splits.gt.1) then
                write(*,'(T2,"Max cell splits: ",i4)') inputs%max_cell_splits
            else
                write(*,'(a,i4)') 'READ_INPUTS: max_cell_splits must be greater than 1 if adaptive time stepping is on: ', &
                                  inputs%max_cell_splits
                error=.True.
            endif
            if(inputs%split_tol.gt.0.0.and.inputs%split_tol.lt.1.0) then
                write(*, '(T2,"Split tolerance: ",1f9.6)') inputs%split_tol
            else
                write(*,'(a,f2.6)') 'READ_INPUTS: split_tol must be positive if adaptive time stepping is on: ', &
                                inputs%split_tol
                error=.True.
            endif
        else
            write(*,'(a)') ' Adaptive time stepping is off'
        endif
        write(*,*) ''
    else
        write(*,'(a,i2)') 'READ_INPUTS: Invalid adaptive switch setting, must be within [0,7]: ', &
                          inputs%adaptive
        error=.True.
    endif

    if(error) then
        stop
    endif

    ! Identify the transition from lambdamax and lambdamin
    call identify_transition(n_stark, stark_pi, stark_sigma, &
                                stark_intens, stark_wavel, line_lambda0)

end subroutine read_inputs

subroutine identify_transition(n_stark, stark_pi, stark_sigma, &
                                stark_intens, stark_wavel, line_lambda0)
    ! Determines the type of transition from user defined inputs

    integer, intent(out)                                     :: n_stark
    integer, dimension(:), allocatable, intent(out)          :: stark_pi, stark_sigma
    real(FLoat64), dimension(:), allocatable, intent(out)    :: stark_intens, stark_wavel
    real(Float64), dimension(3), intent(out)                 :: line_lambda0

    if (inputs%lambdamin > 620.0d0 .and. inputs%lambdamax < 680.0d0 ) then
        ! Assigns stark variables to balmer alpha transition
        n_stark = 15
        initial_state = 3
        final_state = 2

        line_lambda0 = [ 656.28d0, 656.104d0, 656.045d0 ]
         !"Tritium Diagnostics by Balmer-alpha emission" CH Skinner 1993
        allocate(stark_intens(n_stark))
        stark_intens = &
            [ 1.000d0, 18.00d0, 16.00d0, 1681.d0, 2304.d0, &
              729.0d0, 1936.d0, 5490.d0, 1936.d0, 729.0d0, &
              2304.d0, 1681.d0, 16.00d0, 18.00d0, 1.000d0  ]

        allocate(stark_wavel(n_stark))
        stark_wavel = &
            [-2.20200d-07,-1.65200d-07,-1.37700d-07,-1.10200d-07, &
             -8.26400d-08,-5.51000d-08,-2.75600d-08, 0.00000d0,   &
              2.75700d-08, 5.51500d-08, 8.27400d-08, 1.10300d-07, &
              1.38000d-07, 1.65600d-07, 2.20900d-07               ]

        allocate(stark_pi(n_stark))
        stark_pi = &
            [1, 0 , 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1]

        allocate(stark_sigma(n_stark))
        stark_sigma = 1 - stark_pi
        if(inputs%verbose.ge.1) then
           write(*, '(a)') "---- Atomic Transition ----"
           write(*, '(T2, a)') "Balmer-alpha: 3 -> 2"
           write(*, '(T2, "H/D/T wavelengths: [", f8.3, ",", f8.3, ",", f8.3, "] nm")') line_lambda0
           write(*,*) ''
        endif


    else if (inputs%lambdamin > 103.0d0 .and. inputs%lambdamax < 136.0d0) then
        ! Assigns stark varibales to Lyman alpha transition
        n_stark = 3
        initial_state = 2
        final_state = 1

        line_lambda0 = [ 121.57d0, 121.53d0, 121.52d0 ]
        ! "Ultra-sensitive Detection of Hydrogen Isotopes by Lyman-alpha RIS"
        ! Yasuhiro Miyake et al 2012
        allocate(stark_intens(n_stark))
        stark_intens = [ 1.000d0, 2.000d0, 1.000d0 ]

        allocate(stark_wavel(n_stark))
        stark_wavel = [ -8.25800d-8, 0.00000d0, 8.25800d-8 ]

        allocate(stark_pi(n_stark))
        stark_pi = [1, 0, 1]

        allocate(stark_sigma(n_stark))
        stark_sigma = 1 - stark_pi

        if(inputs%verbose.ge.1) then
           write(*, '(a)') "---- Atomic Transition ----"
           write(*, '(T2, a)') "Lyman-alpha: 2 -> 1"
           write(*, '(T2, "H/D/T wavelengths: [", f8.3, ",", f8.3, ",", f8.3, "] nm")') line_lambda0
           write(*,*) ''
        endif

   endif

end subroutine identify_transition

subroutine make_beam_grid
    !+ Makes [[libfida:beam_grid] from user defined inputs
    integer(Int32) :: i, j, k, n
    real(Float64) :: dx, dy, dz, ri(3)
    logical :: inp

    allocate(beam_grid%xc(beam_grid%nx),  &
             beam_grid%yc(beam_grid%ny),  &
             beam_grid%zc(beam_grid%nz))

    dx = (beam_grid%xmax - beam_grid%xmin)/beam_grid%nx
    dy = (beam_grid%ymax - beam_grid%ymin)/beam_grid%ny
    dz = (beam_grid%zmax - beam_grid%zmin)/beam_grid%nz

    do i=1, beam_grid%nx
        beam_grid%xc(i) = beam_grid%xmin + (i-0.5)*dx
    enddo
    do i=1, beam_grid%ny
        beam_grid%yc(i) = beam_grid%ymin + (i-0.5)*dy
    enddo
    do i=1, beam_grid%nz
        beam_grid%zc(i) = beam_grid%zmin + (i-0.5)*dz
    enddo

    beam_grid%dr(1) = abs(beam_grid%xc(2)-beam_grid%xc(1))
    beam_grid%dr(2) = abs(beam_grid%yc(2)-beam_grid%yc(1))
    beam_grid%dr(3) = abs(beam_grid%zc(2)-beam_grid%zc(1))

    beam_grid%lwh(1) = abs(beam_grid%xc(beam_grid%nx) - beam_grid%xc(1)) + beam_grid%dr(1)
    beam_grid%lwh(2) = abs(beam_grid%yc(beam_grid%ny) - beam_grid%yc(1)) + beam_grid%dr(2)
    beam_grid%lwh(3) = abs(beam_grid%zc(beam_grid%nz) - beam_grid%zc(1)) + beam_grid%dr(3)

    beam_grid%volume = beam_grid%lwh(1)*beam_grid%lwh(2)*beam_grid%lwh(3)

    beam_grid%center(1) = (minval(beam_grid%xc) - 0.5*beam_grid%dr(1)) + 0.5*beam_grid%lwh(1)
    beam_grid%center(2) = (minval(beam_grid%yc) - 0.5*beam_grid%dr(2)) + 0.5*beam_grid%lwh(2)
    beam_grid%center(3) = (minval(beam_grid%zc) - 0.5*beam_grid%dr(3)) + 0.5*beam_grid%lwh(3)

    beam_grid%drmin  = minval(beam_grid%dr)
    beam_grid%dv     = beam_grid%dr(1)*beam_grid%dr(2)*beam_grid%dr(3)
    beam_grid%ntrack = (beam_grid%nx+beam_grid%ny+beam_grid%nz)*inputs%max_cell_splits
    beam_grid%ngrid  = beam_grid%nx*beam_grid%ny*beam_grid%nz

    beam_grid%dims(1) = beam_grid%nx
    beam_grid%dims(2) = beam_grid%ny
    beam_grid%dims(3) = beam_grid%nz

    call tb_zyx(beam_grid%alpha,beam_grid%beta,beam_grid%gamma, &
                beam_grid%basis, beam_grid%inv_basis)

    !! Check if beam grid is in the plasma
    n = 0
    do k=1,beam_grid%nz
        do j=1,beam_grid%ny
            do i=1,beam_grid%nx
                ri = [beam_grid%xc(i),beam_grid%yc(j), beam_grid%zc(k)]
                call in_plasma(ri, inp)
                if(inp) n = n + 1
            enddo
        enddo
    enddo

    if(inputs%verbose.ge.1) then
        write(*,'(a)') "---- Beam grid settings ----"
        write(*,'(T2,"Nx: ", i3)') beam_grid%nx
        write(*,'(T2,"Ny: ", i3)') beam_grid%ny
        write(*,'(T2,"Nz: ", i3)') beam_grid%nz
        write(*,'(T2,"dV: ", f5.2," [cm^3]")') beam_grid%dv
        write(*,'(T2,"alpha: ",f5.2," [rad]")') beam_grid%alpha
        write(*,'(T2,"beta:  ",f5.2," [rad]")') beam_grid%beta
        write(*,'(T2,"gamma: ",f5.2," [rad]")') beam_grid%gamma
        write(*,'(T2,"origin: [",f7.2,",",f7.2,",",f7.2,"] [cm]")') beam_grid%origin
        write(*,'(T2,"Number of cells in plasma: ",i8)') n
        write(*,*) ''
    endif

    if(n.le.(0.1*beam_grid%ngrid)) then
        if(inputs%verbose.ge.0) then
            write(*,'(a)') "MAKE_BEAM_GRID: Beam grid definition is poorly defined. &
                            &Less than 10% of the beam grid cells fall within the plasma."
        endif
        stop
    endif

end subroutine make_beam_grid

subroutine make_passive_grid
    !+ Makes [[libfida:pass_grid] from user defined inputs
    real(Float64), dimension(3,spec_chords%nchan+npa_chords%nchan) :: r0_arr, v0_arr
    real(Float64), dimension(2,spec_chords%nchan+npa_chords%nchan) :: xy_enter, xy_exit
    logical, dimension(8+2*(spec_chords%nchan+npa_chords%nchan))   :: yle, ygt
    logical, dimension(spec_chords%nchan+npa_chords%nchan)         :: skip
    real(Float64), dimension(:), allocatable  :: xarr, yarr
    real(Float64), dimension(3,8) :: vertices_xyz, vertices_uvw
    real(Float64), dimension(2,3) :: extrema
    real(Float64), dimension(3,3) :: xyz_axis
    real(Float64), dimension(3)   :: r0, vi
    real(Float64), dimension(8)   :: xarr_beam_grid, yarr_beam_grid
    real(Float64) :: xmin, ymin, xmax, ymax, zmin, zmax, max_length
    real(Float64) :: dlength = 3.0 !cm
    integer :: i, iin, iout, dim_le, dim_gt, dim
    logical :: inp, phi_pos

    !! Get beam grid boundaries
    ! Convert vertices_xyz to vertices_uvw
    ! Note: vertices_xyz has the following coordinate definitions:
    ! 111 = 1, 222 = 2, 112 = 3, 211 = 4, 121 = 5, 212 = 6, 122 = 7, 221 = 8
    xmin = beam_grid%xmin ; xmax = beam_grid%xmax
    ymin = beam_grid%ymin ; ymax = beam_grid%ymax
    zmin = beam_grid%zmin ; zmax = beam_grid%zmax

    !Initialize minimum and maximum vertices
    vertices_xyz(1,1) = xmin   ; vertices_xyz(2,1) = ymin   ; vertices_xyz(3,1) = zmin
    vertices_xyz(1,2) = xmax   ; vertices_xyz(2,2) = ymax   ; vertices_xyz(3,2) = zmax

    !Initialize
    vertices_xyz(:,3) = vertices_xyz(:,1)
    vertices_xyz(:,4) = vertices_xyz(:,1)
    vertices_xyz(:,5) = vertices_xyz(:,1)
    !Update
    vertices_xyz(3,3) = zmax
    vertices_xyz(1,4) = xmax
    vertices_xyz(2,5) = ymax
    !Initialize
    vertices_xyz(:,7) = vertices_xyz(:,2)
    vertices_xyz(:,8) = vertices_xyz(:,2)
    vertices_xyz(:,6) = vertices_xyz(:,2)
    !Update
    vertices_xyz(1,7) = xmin
    vertices_xyz(3,8) = zmin
    vertices_xyz(2,6) = ymin

    do i=1, 8
        call xyz_to_uvw(vertices_xyz(:,i),vertices_uvw(:,i))
        xarr_beam_grid(i) = vertices_uvw(1,i)
        yarr_beam_grid(i) = vertices_uvw(2,i)
    enddo

    !! Next consider passive diagnostic extrema relative to the plasma
    if((inputs%calc_pfida.gt.0).and.(inputs%calc_pnpa.gt.0)) then
        do i=1,(spec_chords%nchan)
            r0_arr(:,i) = spec_chords%los(i)%lens_uvw
            v0_arr(:,i) = spec_chords%los(i)%axis_uvw
        enddo
        do i=1, npa_chords%nchan
            call xyz_to_uvw(npa_chords%det(i)%detector%origin, r0_arr(:,i+spec_chords%nchan))
            xyz_axis = npa_chords%det(i)%detector%basis
            v0_arr(:,i+spec_chords%nchan) = matmul(beam_grid%basis, xyz_axis(:,3))
        enddo
    else if(inputs%calc_pfida.gt.0) then
        do i=1, spec_chords%nchan
            r0_arr(:,i) = spec_chords%los(i)%lens_uvw
            v0_arr(:,i) = spec_chords%los(i)%axis_uvw
        enddo
    else !pnpa>=1 case
        do i=1, npa_chords%nchan
            call xyz_to_uvw(npa_chords%det(i)%detector%origin, r0_arr(:,i))
            xyz_axis = npa_chords%det(i)%detector%basis
            v0_arr(:,i) = matmul(beam_grid%basis, xyz_axis(:,3))
        enddo
    endif

    call get_plasma_extrema(r0_arr,v0_arr,extrema,xarr_beam_grid,yarr_beam_grid)

    !! Store the passive neutral grid
    pass_grid%dr = inter_grid%dr
    pass_grid%dz = inter_grid%dz
    pass_grid%nr = inter_grid%nr
    pass_grid%nz = inter_grid%nz
    allocate(pass_grid%r(pass_grid%nr), pass_grid%z(pass_grid%nz))
    pass_grid%r = inter_grid%r
    pass_grid%z = inter_grid%z
    pass_grid%da = pass_grid%dr*pass_grid%dz

    pass_grid%dphi = 2*pi/100 !TODO: make this user input
    pass_grid%nphi = int(ceiling((extrema(2,3)-extrema(1,3))/pass_grid%dphi))

    allocate(pass_grid%phi(pass_grid%nphi))
    do i=1, pass_grid%nphi
        pass_grid%phi(i) = extrema(1,3) + (i-1)*pass_grid%dphi
    enddo

    pass_grid%dv = pass_grid%dr*pass_grid%dphi*pass_grid%dz
    pass_grid%dims = [pass_grid%nr, pass_grid%nz, pass_grid%nphi]

    pass_grid%ntrack = (pass_grid%nr+pass_grid%nz+pass_grid%nphi)*inputs%max_cell_splits
    pass_grid%ngrid  = pass_grid%nr*pass_grid%nz*pass_grid%nphi

    if(inputs%verbose.ge.1) then
        write(*,'(a)') "---- Passive grid settings ----"
        write(*,'(T2,"Nr: ", i3)') pass_grid%nr
        write(*,'(T2,"Nz: ", i3)') pass_grid%nz
        write(*,'(T2,"Nphi: ", i3)') pass_grid%nphi
        write(*,'(T2,"R  range = [",f6.2,",",f6.2,"]")') &
              pass_grid%r(1),pass_grid%r(pass_grid%nr)
        write(*,'(T2,"Z  range = [",f7.2,",",f6.2,"]")') &
              pass_grid%z(1),pass_grid%z(pass_grid%nz)
        write(*,'(T2,"Phi  range = [",f5.2,",",f5.2,"]")') &
              pass_grid%phi(1),pass_grid%phi(pass_grid%nphi)
        write(*,'(T2,"dA: ", f5.2," [cm^3]")') pass_grid%da
        write(*,*) ''
    endif

end subroutine make_passive_grid

subroutine make_diagnostic_grids
    !+ Makes [[libfida:pass_grid] from user defined inputs, and stores the quantities in
    !+ [[libfida:spec_chords]] and [[libfida:npa_chords]]
    real(Float64), dimension(:,:,:), allocatable :: dlength
    type(LOSElement), dimension(:), allocatable :: los_elem
    type(ParticleTrack), dimension(:), allocatable :: tracks
    real(Float64) :: r0(3), v0(3), r_enter(3), r_exit(3), r0_cyl(3), ri(3)
    real(Float64), dimension(3,3) :: basis
    real(Float64), dimension(2) :: randomu
    real(Float64) :: theta, sqrt_rho, dl, length
    character(len=20) :: system = ''
    real(Float64), dimension(:), allocatable :: probs
    real(Float64), dimension(:,:), allocatable :: eff_rds
    real(Float64), parameter :: inv_4pi = (4.d0*pi)**(-1)
    real(Float64), dimension(3) :: eff_rd, rd, rd_d, r0_d
    real(Float64), dimension(3,3) :: inv_basis
    real(Float64), dimension(:),allocatable :: xd, yd
    type(LocalEMFields) :: fields
    real(Float64) :: total_prob, hh, hw, dprob, dx, dy, r, pitch
    integer :: ichan,k,id,ix,iy,d_index,nd,ind_d(2)
    integer :: i, j, ic, nc, ntrack, ind(3), ii, jj, kk
    integer :: error

    if((inputs%calc_pfida+inputs%calc_pnpa+inputs%calc_cold).gt.0) then
        if(inter_grid%nphi.gt.1) then
            pass_grid = inter_grid
        else
            call make_passive_grid()
        endif
    endif

    !! Spectral line-of-sight passive neutral grid intersection calculations
    allocate(spec_chords%cyl_inter(pass_grid%nr,pass_grid%nz,pass_grid%nphi))
    if((inputs%calc_pfida+inputs%calc_cold).gt.0) then
        allocate(tracks(pass_grid%ntrack))
        allocate(dlength(pass_grid%nr, &
                         pass_grid%nz, &
                         pass_grid%nphi) )
        pass_grid_chan_loop: do i=1,spec_chords%nchan
            r0 = spec_chords%los(i)%lens_uvw
            v0 = spec_chords%los(i)%axis_uvw
            v0 = v0/norm2(v0)
            call line_basis(r0,v0,basis)

            call grid_intersect(r0,v0,length,r_enter,r_exit,passive=.True.)
            if(length.le.0.d0) then
                if(inputs%verbose.ge.1) then
                    WRITE(*,'("Channel ",i5," missed the passive neutral grid or starts inside the plasma")') i
                endif
                cycle pass_grid_chan_loop
            endif

            if(spec_chords%los(i)%spot_size.le.0.d0) then
                nc = 1
            else
                nc = 100
            endif

            dlength = 0.d0
            !$OMP PARALLEL DO schedule(guided) private(ic,randomu,sqrt_rho,theta,r0, &
            !$OMP& length, r_enter, r_exit, j, tracks, ntrack, ind)
            do ic=1,nc
                ! Uniformally sample within spot size
                call randu(randomu)
                sqrt_rho = sqrt(randomu(1))
                theta = 2*pi*randomu(2)
                r0(1) = 0.d0
                r0(2) = spec_chords%los(i)%spot_size*sqrt_rho*cos(theta)
                r0(3) = spec_chords%los(i)%spot_size*sqrt_rho*sin(theta)
                r0 = matmul(basis,r0) + spec_chords%los(i)%lens_uvw

                call grid_intersect(r0,v0,length,r_enter,r_exit,passive=.True.)
                call track_cylindrical(r_enter, v0, tracks, ntrack)
                pass_grid_track_loop: do j=1, ntrack
                    ind = tracks(j)%ind
                    !inds can repeat so add rather than assign
                    !$OMP ATOMIC UPDATE
                    dlength(ind(1),ind(2),ind(3)) = &
                    dlength(ind(1),ind(2),ind(3)) + tracks(j)%time/real(nc) !time == distance
                    !$OMP END ATOMIC
                enddo pass_grid_track_loop
            enddo
            !$OMP END PARALLEL DO
            do kk=1,pass_grid%nphi
                do jj=1,pass_grid%nz
                    rloop: do ii=1, pass_grid%nr
                        if(dlength(ii,jj,kk).ne.0.d0) then
                            dl = dlength(ii,jj,kk)
                            nc = spec_chords%cyl_inter(ii,jj,kk)%nchan + 1
                            if(nc.eq.1) then
                                allocate(spec_chords%cyl_inter(ii,jj,kk)%los_elem(nc))
                                spec_chords%cyl_inter(ii,jj,kk)%los_elem(nc) = LOSElement(i, dl)
                            else
                                allocate(los_elem(nc))
                                los_elem(1:(nc-1)) = spec_chords%cyl_inter(ii,jj,kk)%los_elem
                                los_elem(nc) = LOSElement(i, dl)
                                deallocate(spec_chords%cyl_inter(ii,jj,kk)%los_elem)
                                call move_alloc(los_elem, spec_chords%cyl_inter(ii,jj,kk)%los_elem)
                            endif
                            spec_chords%cyl_inter(ii,jj,kk)%nchan = nc
                        endif
                    enddo rloop
                enddo
            enddo
        enddo pass_grid_chan_loop

        spec_chords%cyl_ncell = count(spec_chords%cyl_inter%nchan.gt.0)
        allocate(spec_chords%cyl_cell(spec_chords%cyl_ncell))

        nc = 0
        do ic=1,pass_grid%ngrid
            call ind2sub(pass_grid%dims,ic,ind)
            ii = ind(1) ; jj = ind(2) ; kk = ind(3)
            if(spec_chords%cyl_inter(ii,jj,kk)%nchan.gt.0) then
                nc = nc + 1
                spec_chords%cyl_cell(nc) = ic
            endif
        enddo
        deallocate(dlength, tracks)
    endif

    !! Spectral line-of-sight beam grid intersection calculations
    if((inputs%tot_spectra+inputs%calc_fida_wght-inputs%calc_pfida).gt.0) then
        allocate(dlength(beam_grid%nx, &
                         beam_grid%ny, &
                         beam_grid%nz) )
        allocate(tracks(beam_grid%ntrack))
        spec_chan_loop: do i=1,spec_chords%nchan
            r0 = spec_chords%los(i)%lens
            v0 = spec_chords%los(i)%axis
            v0 = v0/norm2(v0)
            call line_basis(r0,v0,basis)

            call grid_intersect(r0,v0,length,r_enter,r_exit)
            if(length.le.0.d0) then
                if(inputs%verbose.ge.1) then
                    WRITE(*,'("Channel ",i5," missed the beam grid")') i
                endif
                cycle spec_chan_loop
            endif

            if(spec_chords%los(i)%spot_size.le.0.d0) then
                nc = 1
            else
                nc = 100
            endif

            dlength = 0.d0
            !$OMP PARALLEL DO schedule(guided) private(ic,randomu,sqrt_rho,theta,r0, &
            !$OMP& length, r_enter, r_exit, j, tracks, ntrack, ind)
            do ic=1,nc
                ! Uniformally sample within spot size
                call randu(randomu)
                sqrt_rho = sqrt(randomu(1))
                theta = 2*pi*randomu(2)
                r0(1) = 0.d0
                r0(2) = spec_chords%los(i)%spot_size*sqrt_rho*cos(theta)
                r0(3) = spec_chords%los(i)%spot_size*sqrt_rho*sin(theta)
                r0 = matmul(basis,r0) + spec_chords%los(i)%lens

                call grid_intersect(r0, v0, length, r_enter, r_exit)
                call track(r_enter, v0, tracks, ntrack)
                track_loop: do j=1, ntrack
                    ind = tracks(j)%ind
                    !inds can repeat so add rather than assign
                    !$OMP ATOMIC UPDATE
                    dlength(ind(1),ind(2),ind(3)) = &
                    dlength(ind(1),ind(2),ind(3)) + tracks(j)%time/real(nc) !time == distance
                    !$OMP END ATOMIC
                enddo track_loop
            enddo
            !$OMP END PARALLEL DO
            do kk=1,beam_grid%nz
                do jj=1,beam_grid%ny
                    xloop: do ii=1, beam_grid%nx
                        if(dlength(ii,jj,kk).ne.0.d0) then
                            dl = dlength(ii,jj,kk)
                            nc = spec_chords%inter(ii,jj,kk)%nchan + 1
                            if(nc.eq.1) then
                                allocate(spec_chords%inter(ii,jj,kk)%los_elem(nc))
                                spec_chords%inter(ii,jj,kk)%los_elem(nc) = LOSElement(i, dl)
                            else
                                allocate(los_elem(nc))
                                los_elem(1:(nc-1)) = spec_chords%inter(ii,jj,kk)%los_elem
                                los_elem(nc) = LOSElement(i, dl)
                                deallocate(spec_chords%inter(ii,jj,kk)%los_elem)
                                call move_alloc(los_elem, spec_chords%inter(ii,jj,kk)%los_elem)
                            endif
                            spec_chords%inter(ii,jj,kk)%nchan = nc
                        endif
                    enddo xloop
                enddo
            enddo
        enddo spec_chan_loop

        spec_chords%ncell = count(spec_chords%inter%nchan.gt.0)
        allocate(spec_chords%cell(spec_chords%ncell))

        nc = 0
        do ic=1,beam_grid%ngrid
            call ind2sub(beam_grid%dims,ic,ind)
            ii = ind(1) ; jj = ind(2) ; kk = ind(3)
            if(spec_chords%inter(ii,jj,kk)%nchan.gt.0) then
                nc = nc + 1
                spec_chords%cell(nc) = ic
            endif
        enddo
    endif

    !! NPA probability calculations
    allocate(xd(50),yd(50))
    allocate(probs(beam_grid%ngrid))
    allocate(eff_rds(3,beam_grid%ngrid))
    allocate(npa_chords%phit(beam_grid%nx, &
                             beam_grid%ny, &
                             beam_grid%nz, &
                             npa_chords%nchan) )
    allocate(npa_chords%hit(beam_grid%nx, &
                            beam_grid%ny, &
                            beam_grid%nz) )
    npa_chords%hit = .False.

    npa_chan_loop: do ichan=1,npa_chords%nchan
        v0 = npa_chords%det(ichan)%aperture%origin - npa_chords%det(ichan)%detector%origin
        v0 = v0/norm2(v0)
        call grid_intersect(npa_chords%det(ichan)%detector%origin,v0,length,r0,r0_d)
        if(length.le.0.0) then
            if(inputs%verbose.ge.0) then
                WRITE(*,'("Channel ",i3," centerline missed the beam grid")') ichan
            endif
        endif

        if(inputs%calc_npa_wght.ge.1) then
            hw = npa_chords%det(ichan)%detector%hw
            hh = npa_chords%det(ichan)%detector%hh
            nd = size(xd)
            do i=1,nd
                xd(i) = -hw + 2*hw*(i-0.5)/real(nd)
                yd(i) = -hh + 2*hh*(i-0.5)/real(nd)
            enddo
            dx = abs(xd(2) - xd(1))
            dy = abs(yd(2) - yd(1))
            basis = npa_chords%det(ichan)%detector%basis
            inv_basis = npa_chords%det(ichan)%detector%inv_basis
            eff_rds = 0.d0
            probs = 0.d0
            ! For each grid point find the probability of hitting the detector given an isotropic source
            !$OMP PARALLEL DO schedule(guided) private(ic,i,j,k,ix,iy,total_prob,eff_rd,r0,r0_d, &
            !$OMP& rd_d,rd,d_index,v0,dprob,r,fields,id,ind_d,ind)
            do ic=istart,beam_grid%ngrid,istep
                call ind2sub(beam_grid%dims,ic,ind)
                i = ind(1) ; j = ind(2) ; k = ind(3)
                r0 = [beam_grid%xc(i),beam_grid%yc(j),beam_grid%zc(k)]
                r0_d = matmul(inv_basis,r0-npa_chords%det(ichan)%detector%origin)
                do id = 1, nd*nd
                    call ind2sub([nd,nd],id,ind_d)
                    ix = ind_d(1) ; iy = ind_d(2)
                    rd_d = [xd(ix),yd(iy),0.d0]
                    rd = matmul(basis,rd_d) + npa_chords%det(ichan)%detector%origin
                    v0 = rd - r0
                    d_index = 0
                    call hit_npa_detector(r0,v0,d_index,det=ichan)
                    if(d_index.ne.0) then
                        r = norm2(rd_d - r0_d)**2
                        dprob = (dx*dy) * inv_4pi * r0_d(3)/(r*sqrt(r))
                        eff_rds(:,ic) = eff_rds(:,ic) + dprob*rd
                        probs(ic) = probs(ic) + dprob
                    endif
                enddo
            enddo
            !$OMP END PARALLEL DO
#ifdef _MPI
            call parallel_sum(eff_rds)
            call parallel_sum(probs)
#endif
            do ic = 1, beam_grid%ngrid
                if(probs(ic).gt.0.0) then
                    call ind2sub(beam_grid%dims, ic, ind)
                    i = ind(1) ; j = ind(2) ; k = ind(3)
                    r0 = [beam_grid%xc(i),beam_grid%yc(j),beam_grid%zc(k)]
                    eff_rd = eff_rds(:,ic)/probs(ic)
                    call get_fields(fields,pos=r0)
                    v0 = (eff_rd - r0)/norm2(eff_rd - r0)
                    npa_chords%phit(i,j,k,ichan)%pitch = dot_product(fields%b_norm,v0)
                    npa_chords%phit(i,j,k,ichan)%p = probs(ic)
                    npa_chords%phit(i,j,k,ichan)%dir = v0
                    npa_chords%phit(i,j,k,ichan)%eff_rd = eff_rd
                    npa_chords%hit(i,j,k) = .True.
                endif
            enddo
            total_prob = sum(probs)
            if(total_prob.le.0.d0) then
                if(inputs%verbose.ge.0) then
                    WRITE(*,'("Channel ",i3," missed the beam grid")') ichan
                endif
                cycle npa_chan_loop
            endif
        endif
    enddo npa_chan_loop
    deallocate(probs,eff_rds,xd,yd)

end subroutine make_diagnostic_grids

subroutine read_beam
    !+ Reads neutral beam geometry and stores the quantities in [[libfida:nbi]]
    integer(HID_T) :: fid, gid
    integer(HSIZE_T), dimension(1) :: dims

    real(Float64), dimension(3) ::uvw_src, uvw_axis, pos
    real(Float64) :: dis
    logical :: path_valid
    integer :: error

    !!Initialize HDF5 interface
    call h5open_f(error)

    !!Open HDF5 file
    call h5fopen_f(inputs%geometry_file, H5F_ACC_RDONLY_F, fid, error)

    !!Check if SPEC group exists
    call h5ltpath_valid_f(fid, "/nbi", .True., path_valid, error)
    if(.not.path_valid) then
        if(inputs%verbose.ge.0) then
            write(*,'(a)') 'READ_BEAM: NBI geometry is not in the geometry file'
        endif
        stop
    endif

    !!Open NBI group
    call h5gopen_f(fid, "/nbi", gid, error)

    !!Read in beam definitions
    call h5ltread_dataset_string_f(gid, "/nbi/name",nbi%name, error)
    dims(1) = 3
    call h5ltread_dataset_double_f(gid, "/nbi/src", uvw_src, dims, error)
    call h5ltread_dataset_double_f(gid, "/nbi/axis", uvw_axis, dims, error)
    call h5ltread_dataset_double_f(gid, "/nbi/divy", nbi%divy, dims, error)
    call h5ltread_dataset_double_f(gid, "/nbi/divz", nbi%divz, dims, error)
    call h5ltread_dataset_int_scalar_f(gid, "/nbi/shape", nbi%shape, error)
    call h5ltread_dataset_double_scalar_f(gid, "/nbi/focy", nbi%focy, error)
    call h5ltread_dataset_double_scalar_f(gid, "/nbi/focz", nbi%focz, error)
    call h5ltread_dataset_double_scalar_f(gid, "/nbi/widy", nbi%widy, error)
    call h5ltread_dataset_double_scalar_f(gid, "/nbi/widz", nbi%widz, error)

    !!Read in aperture definitions
    !! Check for naperture for compatibility with old runs
    call h5ltpath_valid_f(gid, "/nbi/naperture", .True., path_valid, error)
    if(path_valid) then
        call h5ltread_dataset_int_scalar_f(gid,"/nbi/naperture",nbi%naperture, error)
    else
        nbi%naperture = 0
    endif
    if(nbi%naperture.gt.0) then
        allocate(nbi%ashape(nbi%naperture), nbi%adist(nbi%naperture), &
                 nbi%awidy(nbi%naperture), nbi%awidz(nbi%naperture),  &
                 nbi%aoffy(nbi%naperture), nbi%aoffz(nbi%naperture)   )

        dims(1) = nbi%naperture
        call h5ltread_dataset_int_f(gid, "/nbi/ashape", nbi%ashape, dims, error)
        call h5ltread_dataset_double_f(gid, "/nbi/awidy", nbi%awidy, dims, error)
        call h5ltread_dataset_double_f(gid, "/nbi/awidz", nbi%awidz, dims, error)
        call h5ltread_dataset_double_f(gid, "/nbi/aoffy", nbi%aoffy, dims, error)
        call h5ltread_dataset_double_f(gid, "/nbi/aoffz", nbi%aoffz, dims, error)
        call h5ltread_dataset_double_f(gid, "/nbi/adist", nbi%adist, dims, error)
    endif

    !!Close NBI group
    call h5gclose_f(gid, error)

    !!Close file id
    call h5fclose_f(fid, error)

    !!Close HDF5 interface
    call h5close_f(error)

    !!Convert to beam grid coordinates
    call uvw_to_xyz(uvw_src,nbi%src)
    nbi%axis = matmul(beam_grid%inv_basis,uvw_axis)

    nbi%vinj=sqrt(2.d0*nbi%einj*1.d3 *e0/(beam_mass*mass_u))*1.d2 !! [cm/s]
    pos = nbi%src + 200.0*nbi%axis
    dis = sqrt(sum((pos - nbi%src)**2))
    nbi%beta = asin((nbi%src(3) - pos(3))/dis)
    nbi%alpha = atan2(pos(2)-nbi%src(2),pos(1)-nbi%src(1))

    call tb_zyx(nbi%alpha,nbi%beta,0.d0,nbi%basis,nbi%inv_basis)

    if(inputs%verbose.ge.1) then
        write(*,'(a)') '---- Neutral beam settings ----'
        write(*,'(T2,"Beam: ",a)') nbi%name
        write(*,'(T2,"Power:   ",f5.2," [MW]")') nbi%pinj
        write(*,'(T2,"Voltage: ",f6.2," [keV]")') nbi%einj
        write(*,'(T2,"Current fractions: [",f5.2,",",f5.2,",",f5.2,"]")') nbi%current_fractions
        write(*,*) ''
    endif

end subroutine read_beam

subroutine read_chords
    !+ Reads the spectral geometry and stores the quantities in [[libfida:spec_chords]]
    integer(HID_T) :: fid, gid
    integer(HSIZE_T), dimension(2) :: dims
    logical :: path_valid

    real(Float64), dimension(:,:), allocatable :: lenses
    real(Float64), dimension(:,:), allocatable :: axes
    real(Float64) :: xyz_lens(3), xyz_axis(3)
    character(charlim) :: system = ''

    integer :: i
    integer :: error

    if(inputs%verbose.ge.1) then
        write(*,'(a)') '---- FIDA/BES settings ----'
    endif

    !!Initialize HDF5 interface
    call h5open_f(error)

    !!Open HDF5 file
    call h5fopen_f(inputs%geometry_file, H5F_ACC_RDONLY_F, fid, error)

    !!Check if SPEC group exists
    call h5ltpath_valid_f(fid, "/spec", .True., path_valid, error)
    if(.not.path_valid) then
        if(inputs%verbose.ge.1) then
            write(*,'(a)') 'FIDA/BES geometry is not in the geometry file'
            write(*,'(a)') 'Continuing without spectral diagnostics'
        endif
        inputs%calc_spec = 0
        inputs%tot_spectra=0
        inputs%calc_fida = 0
        inputs%calc_pfida = 0
        inputs%calc_bes = 0
        inputs%calc_dcx = 0
        inputs%calc_halo = 0
        inputs%calc_cold = 0
        inputs%calc_brems = 0
        inputs%calc_fida_wght = 0
        call h5fclose_f(fid, error)
        call h5close_f(error)
        return
    endif

    !!Open SPEC group
    call h5gopen_f(fid, "/spec", gid, error)

    call h5ltread_dataset_string_f(gid, "/spec/system", system, error)
    call h5ltread_dataset_int_scalar_f(gid, "/spec/nchan", spec_chords%nchan, error)

    allocate(lenses(3, spec_chords%nchan))
    allocate(axes(3, spec_chords%nchan))
    allocate(spec_chords%los(spec_chords%nchan))
    allocate(spec_chords%radius(spec_chords%nchan))

    dims = [3,spec_chords%nchan]
    call h5ltread_dataset_double_f(gid, "/spec/lens", lenses, dims, error)
    call h5ltread_dataset_double_f(gid, "/spec/axis", axes, dims, error)
    call h5ltread_dataset_double_f(gid, "/spec/spot_size", spec_chords%los%spot_size, dims(2:2), error)
    call h5ltread_dataset_double_f(gid, "/spec/sigma_pi", spec_chords%los%sigma_pi, dims(2:2), error)
    call h5ltread_dataset_double_f(gid, "/spec/radius", spec_chords%radius, dims(2:2), error)

    !!Close SPEC group
    call h5gclose_f(gid, error)

    !!Close file id
    call h5fclose_f(fid, error)

    !!Close HDF5 interface
    call h5close_f(error)

    chan_loop: do i=1,spec_chords%nchan
        call uvw_to_xyz(lenses(:,i), xyz_lens)
        xyz_axis = matmul(beam_grid%inv_basis, axes(:,i))
        spec_chords%los(i)%lens = xyz_lens
        spec_chords%los(i)%axis = xyz_axis
        spec_chords%los(i)%lens_uvw = lenses(:,i)
        spec_chords%los(i)%axis_uvw = axes(:,i)
    enddo chan_loop

    if(inputs%verbose.ge.1) then
        write(*,'(T2,"FIDA/BES System: ",a)') trim(adjustl(system))
        write(*,'(T2,"Number of channels: ",i5)') spec_chords%nchan
        write(*,*) ''
    endif

    deallocate(lenses,axes)

end subroutine read_chords

subroutine read_npa
    !+ Reads the NPA geometry and stores the quantities in [[libfida:npa_chords]]
    integer(HID_T) :: fid, gid
    integer(HSIZE_T), dimension(2) :: dims
    logical :: path_valid

    real(Float64), dimension(:,:), allocatable :: a_tedge,a_redge,a_cent
    real(Float64), dimension(:,:), allocatable :: d_tedge,d_redge,d_cent
    character(len=20) :: system = ''

    real(Float64), dimension(3) :: xyz_a_tedge,xyz_a_redge,xyz_a_cent
    real(Float64), dimension(3) :: xyz_d_tedge,xyz_d_redge,xyz_d_cent
    real(Float64), dimension(3,3) :: basis, inv_basis
    real(Float64) :: hh, hw
    integer :: ichan
    integer :: error

    !!Initialize HDF5 interface
    call h5open_f(error)

    !!Open HDF5 file
    call h5fopen_f(inputs%geometry_file, H5F_ACC_RDWR_F, fid, error)

    !!Check if NPA group exists
    call h5ltpath_valid_f(fid, "/npa", .True., path_valid, error)
    if(.not.path_valid) then
        if(inputs%verbose.ge.0) then
            write(*,'(a)') 'NPA geometry is not in the geometry file'
            write(*,'(a)') 'Continuing without NPA diagnostics'
        endif
        inputs%calc_npa = 0
        inputs%calc_pnpa = 0
        inputs%calc_npa_wght = 0
        call h5fclose_f(fid, error)
        call h5close_f(error)
        return
    endif

    !!Open NPA group
    call h5gopen_f(fid, "/npa", gid, error)

    call h5ltread_dataset_string_f(gid, "/npa/system", system, error)
    call h5ltread_dataset_int_scalar_f(gid, "/npa/nchan", npa_chords%nchan, error)

    if(inputs%verbose.ge.1) then
        write(*,'(a)') "---- NPA settings ----"
        write(*,'(T2,"NPA System: ", a)') trim(adjustl(system))
        write(*,'(T2,"Number of channels: ",i3)') npa_chords%nchan
    endif

    allocate(a_tedge(3, npa_chords%nchan))
    allocate(a_redge(3, npa_chords%nchan))
    allocate(a_cent(3,  npa_chords%nchan))
    allocate(d_tedge(3, npa_chords%nchan))
    allocate(d_redge(3, npa_chords%nchan))
    allocate(d_cent(3,  npa_chords%nchan))
    allocate(npa_chords%radius(npa_chords%nchan))
    allocate(npa_chords%det(npa_chords%nchan))

    dims = [3,spec_chords%nchan]
    call h5ltread_dataset_double_f(gid,"/npa/radius", npa_chords%radius, dims(2:2), error)
    call h5ltread_dataset_int_f(gid, "/npa/a_shape", npa_chords%det%aperture%shape, dims(2:2), error)
    call h5ltread_dataset_double_f(gid, "/npa/a_tedge", a_tedge, dims, error)
    call h5ltread_dataset_double_f(gid, "/npa/a_redge", a_redge, dims, error)
    call h5ltread_dataset_double_f(gid, "/npa/a_cent",  a_cent, dims, error)

    call h5ltread_dataset_int_f(gid, "/npa/d_shape", npa_chords%det%detector%shape, dims(2:2), error)
    call h5ltread_dataset_double_f(gid, "/npa/d_tedge", d_tedge, dims, error)
    call h5ltread_dataset_double_f(gid, "/npa/d_redge", d_redge, dims, error)
    call h5ltread_dataset_double_f(gid, "/npa/d_cent",  d_cent, dims, error)

    !!Close NPA group
    call h5gclose_f(gid, error)

    !!Close file id
    call h5fclose_f(fid, error)

    !!Close HDF5 interface
    call h5close_f(error)

    chan_loop: do ichan=1,npa_chords%nchan
        ! Convert to beam grid coordinates
        call uvw_to_xyz(a_cent(:,ichan), xyz_a_cent)
        call uvw_to_xyz(a_redge(:,ichan),xyz_a_redge)
        call uvw_to_xyz(a_tedge(:,ichan),xyz_a_tedge)
        call uvw_to_xyz(d_cent(:,ichan), xyz_d_cent)
        call uvw_to_xyz(d_redge(:,ichan),xyz_d_redge)
        call uvw_to_xyz(d_tedge(:,ichan),xyz_d_tedge)

        ! Define detector/aperture hh/hw
        npa_chords%det(ichan)%detector%hw = norm2(xyz_d_redge - xyz_d_cent)
        npa_chords%det(ichan)%aperture%hw = norm2(xyz_a_redge - xyz_a_cent)

        npa_chords%det(ichan)%detector%hh = norm2(xyz_d_tedge - xyz_d_cent)
        npa_chords%det(ichan)%aperture%hh = norm2(xyz_a_tedge - xyz_a_cent)

        ! Define detector/aperture origin
        npa_chords%det(ichan)%detector%origin = xyz_d_cent
        npa_chords%det(ichan)%aperture%origin = xyz_a_cent

        ! Define detector/aperture basis
        call plane_basis(xyz_d_cent, xyz_d_redge, xyz_d_tedge, &
             npa_chords%det(ichan)%detector%basis, &
             npa_chords%det(ichan)%detector%inv_basis)
        call plane_basis(xyz_a_cent, xyz_a_redge, xyz_a_tedge, &
             npa_chords%det(ichan)%aperture%basis, &
             npa_chords%det(ichan)%aperture%inv_basis)
    enddo chan_loop

    if(inputs%verbose.ge.1) write(*,'(50X,a)') ""

    deallocate(a_cent,a_redge,a_tedge)
    deallocate(d_cent,d_redge,d_tedge)

end subroutine read_npa

subroutine read_cfpd
    !+ Reads the CFPD geometry and stores the quantities in [[libfida:ctable]]
    integer(HID_T) :: fid, gid
    integer(HSIZE_T), dimension(3) :: dims3
    integer(HSIZE_T), dimension(5) :: dims5
    logical :: path_valid
    real(Float64),dimension(3) :: uvwi,uvwf,rpzi,rpzf

    character(len=20) :: system = ''

    integer :: error

    !!Initialize HDF5 interface
    call h5open_f(error)

    !!Open HDF5 file
    call h5fopen_f(inputs%geometry_file, H5F_ACC_RDWR_F, fid, error)

    if(inputs%verbose.ge.1) write(*,'(a)') "---- CFPD settings ----"
    !!Check if CFPD group exists
    call h5ltpath_valid_f(fid, "/cfpd", .True., path_valid, error)
    if(.not.path_valid) then
        if(inputs%verbose.ge.0) then
            write(*,'(a)') 'CFPD geometry is not in the geometry file'
            write(*,'(a)') 'Continuing without CFPD diagnostics'
            write(*,*) ''
        endif
        inputs%calc_cfpd = 0
        call h5fclose_f(fid, error)
        call h5close_f(error)
        return
    endif

    !!Open CFPD group
    call h5gopen_f(fid, "/cfpd", gid, error)

    call h5ltread_dataset_string_f(gid, "/cfpd/system", system, error)
    call h5ltread_dataset_int_scalar_f(gid, "/cfpd/nchan", ctable%nchan, error)

    if(inputs%verbose.ge.1) then
        write(*,'(T2,"CFPD System: ", a)') trim(adjustl(system))
        write(*,'(T2,"Number of channels: ",i3)') ctable%nchan
    endif

    call h5ltread_dataset_int_scalar_f(gid,"/cfpd/nenergy", ctable%nenergy, error)
    call h5ltread_dataset_int_scalar_f(gid,"/cfpd/nrays", ctable%nrays, error)
    call h5ltread_dataset_int_scalar_f(gid,"/cfpd/nsteps", ctable%nsteps, error)

    allocate(ctable%earray(ctable%nenergy))
    allocate(ctable%nactual(ctable%nenergy, ctable%nrays, ctable%nchan))
    allocate(ctable%daomega(ctable%nenergy, ctable%nrays, ctable%nchan))
    allocate(ctable%sightline(ctable%nenergy, 6, ctable%nsteps, ctable%nrays, ctable%nchan))

    dims3 = [ctable%nenergy, ctable%nrays, ctable%nchan]
    dims5 = [ctable%nenergy, 6, ctable%nsteps, ctable%nrays, ctable%nchan]
    call h5ltread_dataset_double_f(gid, "/cfpd/earray", ctable%earray, dims3(1:1), error)
    call h5ltread_dataset_double_f(gid, "/cfpd/nactual", ctable%nactual, dims3, error)
    call h5ltread_dataset_double_f(gid, "/cfpd/daomega", ctable%daomega, dims3, error)
    call h5ltread_dataset_double_f(gid, "/cfpd/sightline", ctable%sightline, dims5, error)

    !!Close CFPD group
    call h5gclose_f(gid, error)

    !!Close file id
    call h5fclose_f(fid, error)

    !!Close HDF5 interface
    call h5close_f(error)

    rpzi(1) = ctable%sightline(1,4,1,1,1)
    rpzi(2) = ctable%sightline(1,5,1,1,1)
    rpzi(3) = ctable%sightline(1,6,1,1,1)
    rpzf(1) = ctable%sightline(1,4,2,1,1)
    rpzf(2) = ctable%sightline(1,5,2,1,1)
    rpzf(3) = ctable%sightline(1,6,2,1,1)

    uvwi(1) = rpzi(1)*cos(rpzi(2))
    uvwi(2) = rpzi(1)*sin(rpzi(2))
    uvwi(3) = rpzi(3)
    uvwf(1) = rpzf(1)*cos(rpzf(2))
    uvwf(2) = rpzf(1)*sin(rpzf(2))
    uvwf(3) = rpzf(3)
    ctable%dl = norm2(uvwf-uvwi)

    ctable%dE = ctable%earray(2)-ctable%earray(1)

    if(inputs%verbose.ge.1) write(*,'(50X,a)') ""

end subroutine read_cfpd

subroutine read_plasma
    !+ Reads in Plasma composition e.g. the thermal, fast-ion, and impurity species
    integer(HID_T) :: fid, gid
    integer(HSIZE_T), dimension(1) :: dims
    integer :: error,i,w

    !!Initialize HDF5 interface
    call h5open_f(error)

    !!Open HDF5 file
    call h5fopen_f(inputs%equilibrium_file, H5F_ACC_RDONLY_F, fid, error)

    !!Open PLASMA group
    call h5gopen_f(fid, "/plasma", gid, error)

    !! Read in number of thermal species
    call h5ltread_dataset_int_scalar_f(gid, "/plasma/nthermal", n_thermal, error)
    if (n_thermal.gt.max_species) then
        write(*,'(a,i1,">",i1)') 'READ_PLASMA: Number of ion species exceeds maximum: ',n_thermal, max_species
        stop
    endif

    !! Read in species mass and charge
    dims(1) = n_thermal
    call h5ltread_dataset_double_f(gid, "/plasma/species_mass", thermal_mass(1:n_thermal), dims, error)
    call h5ltread_dataset_int_scalar_f(gid,"/plasma/impurity_charge", impurity_charge, error)

    !!Close PLASMA group
    call h5gclose_f(gid, error)

    !!Close file
    call h5fclose_f(fid, error)

    !!Close HDF5 interface
    call h5close_f(error)

    !! Correct misspecified masses
    w = minloc(abs(supported_masses - beam_mass),1)
    if(abs(beam_mass - supported_masses(w)).gt.0.5) then
        write(*,'("READ_PLASMA: Got unsupported beam ion mass: ",f6.3," amu")') beam_mass
        write(*,'(a)') 'Only Hydrogen, Deuterium, and Tritium are supported'
        stop
    else
        beam_mass = supported_masses(w)
        beam_lambda0 = line_lambda0(w)
    endif

    do i=1, n_thermal
        w = minloc(abs(supported_masses-thermal_mass(i)),1)
        if(abs(thermal_mass(i)-supported_masses(w)).gt.0.5) then
            write(*,'("READ_PLASMA: Got unsupported thermal-ion mass: ",f6.3," amu")') thermal_mass(i)
            write(*,'(a)') 'Only Hydrogen, Deuterium, and Tritium are supported'
            stop
        else
            thermal_mass(i) = supported_masses(w)
            thermal_lambda0(i) = line_lambda0(w)
        endif
    enddo

    select case (impurity_charge)
        case (2)
            impurity_mass = He4_amu
        case (5)
            impurity_mass = B10_amu
        case (6)
            impurity_mass = C12_amu
        case DEFAULT
            impurity_mass = 2.d0*impurity_charge*H1_amu
    end select

    if(inputs%verbose.ge.1) then
        write(*,'(a)') '---- Plasma Composition ----'
        write(*,'(T2,"Thermal-ion mass: ",f6.3," [amu]")') thermal_mass(1:n_thermal)
        write(*,'(T2,"Impurity mass: ",f6.3," [amu]")') impurity_mass
        write(*,'(T2,"Beam-ion mass: ",f6.3," [amu]")') beam_mass
        write(*,'(T2,"Fast-ion mass: ",f6.3," [amu]")') beam_mass
        write(*,*) ''
    endif

end subroutine

subroutine read_equilibrium
    !+ Reads in the interpolation grid, plasma parameters, and fields
    !+ and stores the quantities in [[libfida:inter_grid]] and [[libfida:equil]]
    integer(HID_T) :: fid, gid
    integer(HSIZE_T), dimension(3) :: dims
    integer(HSIZE_T), dimension(4) :: dims4

    integer :: impc, ic, ir, iz, iphi, it, ind(3), i
    type(LocalProfiles) :: plasma
    real(Float64) :: photons, smix(max_species)
    real(Float64), dimension(nlevs) :: rates, denn, rates_avg
    real(Float64), dimension(3) :: vi
    integer :: error
    integer :: n = 50
    logical :: path_valid

    integer, dimension(:,:,:), allocatable :: p_mask, f_mask
    real(Float64), dimension(:,:,:), allocatable :: denn3d
    real(Float64), dimension(:,:,:,:), allocatable :: deni

    !!Initialize HDF5 interface
    call h5open_f(error)

    !!Open HDF5 file
    call h5fopen_f(inputs%equilibrium_file, H5F_ACC_RDONLY_F, fid, error)

    !!Open PLASMA group
    call h5gopen_f(fid, "/plasma", gid, error)

    !!Read in interpolation grid
    call h5ltread_dataset_int_scalar_f(gid, "/plasma/nr", inter_grid%nr, error)
    call h5ltread_dataset_int_scalar_f(gid, "/plasma/nz", inter_grid%nz, error)
    call h5ltpath_valid_f(gid, "/plasma/nphi", .True., path_valid, error)
    if(path_valid) then
        call h5ltread_dataset_int_scalar_f(gid, "/plasma/nphi", inter_grid%nphi, error)
    else
        inter_grid%nphi=1
    endif

    inter_grid%dims = [inter_grid%nr, inter_grid%nz, inter_grid%nphi]

    allocate(inter_grid%r(inter_grid%nr),inter_grid%z(inter_grid%nz),inter_grid%phi(inter_grid%nphi))
    allocate(p_mask(inter_grid%nr,inter_grid%nz,inter_grid%nphi))
    allocate(f_mask(inter_grid%nr,inter_grid%nz,inter_grid%nphi))
    allocate(denn3d(inter_grid%nr,inter_grid%nz,inter_grid%nphi))

    dims = [inter_grid%nr, inter_grid%nz, inter_grid%nphi]

    call h5ltread_dataset_double_f(gid, "/plasma/r", inter_grid%r, dims(1:1), error)
    call h5ltread_dataset_double_f(gid, "/plasma/z", inter_grid%z, dims(2:2), error)
    if(path_valid) then
        call h5ltread_dataset_double_f(gid, "/plasma/phi", inter_grid%phi, dims(3:3), error)
    else
        inter_grid%phi=0.d0
    endif

    inter_grid%dr = abs(inter_grid%r(2)-inter_grid%r(1))
    inter_grid%dz = abs(inter_grid%z(2)-inter_grid%z(1))
    inter_grid%da = inter_grid%dr*inter_grid%dz
    if (inter_grid%nphi .eq. 1) then
        inter_grid%dphi = 2*pi
    else
        inter_grid%dphi = abs(inter_grid%phi(2)-inter_grid%phi(1))
    endif
    inter_grid%dv = inter_grid%dr*inter_grid%dphi*inter_grid%dz

    inter_grid%ntrack = inter_grid%nr+inter_grid%nz+inter_grid%nphi
    inter_grid%ngrid  = inter_grid%nr*inter_grid%nz*inter_grid%nphi

    if(inputs%verbose.ge.1) then
        write(*,'(a)') '---- Interpolation grid settings ----'
        write(*,'(T2,"Nr: ",i3)') inter_grid%nr
        write(*,'(T2,"Nz: ",i3)') inter_grid%nz
        if(inter_grid%nphi.gt.1) then
            write(*,'(T2,"Nphi: ",i3)') inter_grid%nphi
        endif
        write(*,'(T2,"dA: ", f5.2," [cm^2]")') inter_grid%da
        write(*,*) ''
    endif

    !!Read in plasma parameters
    allocate(equil%plasma(inter_grid%nr,inter_grid%nz,inter_grid%nphi))

    call h5ltread_dataset_double_f(gid, "/plasma/dene", equil%plasma%dene, dims, error)
    call h5ltread_dataset_double_f(gid, "/plasma/denimp", equil%plasma%denimp, dims, error)
    call h5ltread_dataset_double_f(gid, "/plasma/te", equil%plasma%te, dims, error)
    call h5ltread_dataset_double_f(gid, "/plasma/ti", equil%plasma%ti, dims, error)
    call h5ltread_dataset_double_f(gid, "/plasma/zeff", equil%plasma%zeff, dims, error)
    call h5ltread_dataset_double_f(gid, "/plasma/vr", equil%plasma%vr, dims, error)
    call h5ltread_dataset_double_f(gid, "/plasma/vt", equil%plasma%vt, dims, error)
    call h5ltread_dataset_double_f(gid, "/plasma/vz", equil%plasma%vz, dims, error)
    call h5ltread_dataset_int_f(gid, "/plasma/mask", p_mask, dims,error)

    ! Read in ion densities
    dims4(1) = n_thermal
    dims4(2:4) = dims

    allocate(deni(n_thermal,inter_grid%nr, inter_grid%nz, inter_grid%nphi))
    call h5ltread_dataset_double_f(gid, "/plasma/deni", deni, dims4, error)
    do i=1,n_thermal
        equil%plasma(:,:,:)%deni(i) = deni(i,:,:,:)
    enddo
    deallocate(deni)

    impc = impurity_charge

    where(equil%plasma%zeff.lt.1.0)
        equil%plasma%zeff = 1
    endwhere

    where(equil%plasma%zeff.gt.impc)
        equil%plasma%zeff = impc
    endwhere

    where(equil%plasma%dene.lt.0.0)
        equil%plasma%dene = 0.0
    endwhere

    where(equil%plasma%denimp.lt.0.0)
        equil%plasma%denimp = 0.0
    endwhere

    where(equil%plasma%te.lt.0.0)
        equil%plasma%te = 0.0
    endwhere

    where(equil%plasma%ti.lt.0.0)
        equil%plasma%ti = 0.0
    endwhere

    do i=1,n_thermal
        where(equil%plasma%deni(i).lt.0.0)
            equil%plasma%deni(i) = 0.0
        endwhere
    enddo


    call h5ltpath_valid_f(fid, "/plasma/denn", .True., path_valid, error)
    if(path_valid) then
        call h5ltread_dataset_double_f(gid, "/plasma/denn", denn3d, dims, error)
        where(denn3d.lt.0.0)
            denn3d = 0.0
        endwhere
        loop_over_cells: do ic=1, inter_grid%nr*inter_grid%nz*inter_grid%nphi
            call ind2sub(inter_grid%dims,ic,ind)
            ir = ind(1) ; iz = ind(2) ; iphi = ind(3)
            if(p_mask(ir,iz,iphi).lt.0.5) cycle loop_over_cells
            if(denn3d(ir,iz,iphi).le.0.0) cycle loop_over_cells
            plasma = equil%plasma(ir,iz,iphi)
            plasma%vrot = [plasma%vr, plasma%vt, plasma%vz]
            plasma%in_plasma = .True.
            smix = plasma%deni/sum(plasma%deni)
            do i=1,n_thermal
                rates_avg = 0.0
                do it=1,n
                    rates = 0.0
                    rates(1) = 1.d19
                    call mc_halo(plasma, thermal_mass(i), vi)
                    call colrad(plasma, thermal_mass(i), vi, 1.0d-7, rates, denn, photons)
                    rates_avg = rates_avg + rates/n
                enddo
                if(sum(rates_avg).le.0.0) cycle loop_over_cells
                equil%plasma(ir,iz,iphi)%denn(:,i) = smix(i)*denn3d(ir,iz,iphi)*(rates_avg)/sum(rates_avg)
            enddo
        enddo loop_over_cells
    else
        if((inputs%calc_pnpa + inputs%calc_pfida + inputs%calc_cold).gt.0) then
            if(inputs%verbose.ge.0) then
                write(*,'(a)') "READ_EQUILIBRIUM: Cold neutral density was not provided"
                write(*,'(a)') "Continuing without passive calculations"
            endif
        endif
        inputs%calc_pnpa = 0
        inputs%calc_pfida = 0
        inputs%calc_cold = 0
    endif

    !!Close PLASMA group
    call h5gclose_f(gid, error)

    !!Open FIELDS group
    call h5gopen_f(fid, "/fields", gid, error)

    allocate(equil%fields(inter_grid%nr,inter_grid%nz,inter_grid%nphi))

    !!Read in electromagnetic fields
    call h5ltread_dataset_double_f(gid, "/fields/br", equil%fields%br, dims, error)
    call h5ltread_dataset_double_f(gid, "/fields/bt", equil%fields%bt, dims, error)
    call h5ltread_dataset_double_f(gid, "/fields/bz", equil%fields%bz, dims, error)
    call h5ltread_dataset_double_f(gid, "/fields/er", equil%fields%er, dims, error)
    call h5ltread_dataset_double_f(gid, "/fields/et", equil%fields%et, dims, error)
    call h5ltread_dataset_double_f(gid, "/fields/ez", equil%fields%ez, dims, error)
    call h5ltread_dataset_int_f(gid, "/fields/mask", f_mask, dims,error)

    !!Calculate B field derivatives
    call deriv(inter_grid%r, inter_grid%z, inter_grid%phi, equil%fields%br, &
        equil%fields%dbr_dr, equil%fields%dbr_dz, equil%fields%dbr_dphi)
    call deriv(inter_grid%r, inter_grid%z, inter_grid%phi, equil%fields%bt, &
        equil%fields%dbt_dr, equil%fields%dbt_dz, equil%fields%dbt_dphi)
    call deriv(inter_grid%r, inter_grid%z, inter_grid%phi, equil%fields%bz, &
        equil%fields%dbz_dr, equil%fields%dbz_dz, equil%fields%dbz_dphi)

    !!Close FIELDS group
    call h5gclose_f(gid, error)

    !!Close file
    call h5fclose_f(fid, error)

    !!Close HDF5 interface
    call h5close_f(error)

    allocate(equil%mask(inter_grid%nr,inter_grid%nz,inter_grid%nphi))
    equil%mask = 0.d0
    where ((p_mask.eq.1).and.(f_mask.eq.1)) equil%mask = 1.d0
    if (sum(equil%mask).le.0.d0) then
        if(inputs%verbose.ge.0) then
            write(*,'(a)') "READ_EQUILIBRIUM: Plasma and/or fields are not well defined anywhere"
        endif
        stop
    endif

end subroutine read_equilibrium

subroutine read_f(fid, error)
    !+ Reads in the fast-ion distribution function and stores the quantities in [[libfida:fbm]]
    integer(HID_T), intent(inout) :: fid
        !+ HDF5 file ID
    integer, intent(out)          :: error
        !+ Error code

    integer(HSIZE_T), dimension(5) :: dims
    real(Float64) :: deni_tot
    integer :: ir,is
    logical :: path_valid

    if(inputs%verbose.ge.1) then
        write(*,'(a)') '---- Fast-ion distribution settings ----'
    endif

    call h5ltread_dataset_int_scalar_f(fid,"/nenergy", fbm%nenergy, error)
    call h5ltread_dataset_int_scalar_f(fid,"/npitch", fbm%npitch, error)
    call h5ltread_dataset_int_scalar_f(fid,"/nr", fbm%nr, error)
    call h5ltread_dataset_int_scalar_f(fid,"/nz", fbm%nz, error)
    call h5ltpath_valid_f(fid, "/nphi", .True., path_valid, error)
    if(path_valid) then
        call h5ltread_dataset_int_scalar_f(fid,"/nphi", fbm%nphi, error)
    else
        fbm%nphi=1
    endif

    if(((fbm%nr.ne.inter_grid%nr).or.(fbm%nz.ne.inter_grid%nz)).or.(fbm%nphi.ne.inter_grid%nphi)) then
        if(inputs%verbose.ge.0) then
            write(*,'(a)') "READ_F: Distribution file has incompatable grid dimensions"
        endif
        stop
    endif

    if (fbm%nphi .eq. 1) then
        allocate(fbm%energy(fbm%nenergy), fbm%pitch(fbm%npitch), fbm%r(fbm%nr), fbm%z(fbm%nz), fbm%phi(1))
        allocate(fbm%denf(fbm%nr, fbm%nz,1))
        allocate(fbm%f(fbm%nenergy, fbm%npitch, fbm%nr, fbm%nz,1))
    else
        allocate(fbm%energy(fbm%nenergy), fbm%pitch(fbm%npitch), fbm%r(fbm%nr), fbm%z(fbm%nz), fbm%phi(fbm%nphi))
        allocate(fbm%denf(fbm%nr, fbm%nz, fbm%nphi))
        allocate(fbm%f(fbm%nenergy, fbm%npitch, fbm%nr, fbm%nz, fbm%nphi))
    endif

    if (fbm%nphi .eq. 1) then
        dims = [fbm%nenergy, fbm%npitch, fbm%nr, fbm%nz, 1]
        call h5ltread_dataset_double_f(fid, "/denf",fbm%denf, dims(3:4), error)
        call h5ltread_dataset_double_f(fid, "/f", fbm%f, dims(1:4), error)
    else
        dims = [fbm%nenergy, fbm%npitch, fbm%nr, fbm%nz, fbm%nphi]
        call h5ltread_dataset_double_f(fid, "/denf",fbm%denf, dims(3:5), error)
        call h5ltread_dataset_double_f(fid, "/f", fbm%f, dims(1:5), error)
    endif
    call h5ltread_dataset_double_f(fid, "/energy", fbm%energy, dims(1:1), error)
    call h5ltread_dataset_double_f(fid, "/pitch", fbm%pitch, dims(2:2), error)
    call h5ltread_dataset_double_f(fid, "/r", fbm%r, dims(3:3), error)
    call h5ltread_dataset_double_f(fid, "/z", fbm%z, dims(4:4), error)
    call h5ltpath_valid_f(fid, "/phi", .True., path_valid, error)
    if(path_valid) then
        call h5ltread_dataset_double_f(fid, "/phi", fbm%phi, dims(5:5), error)
    else
        fbm%phi=0.d0
    endif
    !call h5ltread_dataset_double_scalar_f(fid,"/A",fbm%A, error)
    fbm%A = beam_mass

    equil%plasma%denf = fbm%denf

    fbm%dE = abs(fbm%energy(2) - fbm%energy(1))
    fbm%dp = abs(fbm%pitch(2) - fbm%pitch(1))
    fbm%dr = abs(fbm%r(2) - fbm%r(1))
    fbm%dz = abs(fbm%z(2) - fbm%z(1))
    if (fbm%nphi .eq. 1) then
        fbm%dphi = 2*pi
    else
        fbm%dphi = abs(fbm%phi(2)-fbm%phi(1))
    endif

    fbm%emin = minval(fbm%energy,1)
    fbm%emax = maxval(fbm%energy,1)
    fbm%e_range = fbm%emax - fbm%emin
    fbm%pmin = minval(fbm%pitch,1)
    fbm%pmax = maxval(fbm%pitch,1)
    fbm%p_range = fbm%pmax - fbm%pmin
    fbm%rmin = minval(fbm%r,1)
    fbm%rmax = maxval(fbm%r,1)
    fbm%r_range = fbm%rmax - fbm%rmin
    fbm%zmin = minval(fbm%z,1)
    fbm%zmax = maxval(fbm%z,1)
    fbm%z_range = fbm%zmax - fbm%zmin
    fbm%phimin = minval(fbm%phi,1)
    fbm%phimax = maxval(fbm%phi,1)
    fbm%phi_range = fbm%phimax - fbm%phimin

    deni_tot = 0.0
    do ir=1,fbm%nr
        fbm%n_tot = fbm%n_tot + fbm%dphi*fbm%dr*fbm%dz*sum(fbm%denf(ir,:,:))*fbm%r(ir)
        do is=1,n_thermal
            deni_tot = deni_tot + fbm%dphi*fbm%dr*fbm%dz*sum(equil%plasma(ir,:,:)%deni(is))*fbm%r(ir)
        enddo
    enddo

    if(fbm%n_tot.ge.deni_tot) then
        if(inputs%verbose.ge.0) then
            write(*,'(a," (",ES10.3," >=",ES10.3,")")') &
                "READ_F: The total of number of fast ions exceeded the total number of thermal ions.", &
                 fbm%n_tot, deni_tot
            write(*,'(a)') "This is usually caused by zeff being incorrect."
        endif
        stop
    endif

    if(inputs%verbose.ge.1) then
        if(fbm%nphi.gt.1) then
            write(*,'(T2,"Distribution type: ",a)') "Non-axisymmetric Fast-ion Density Function F(energy,pitch,R,Z,Phi)"
        else
            write(*,'(T2,"Distribution type: ",a)') "Axisymmetric Fast-ion Density Function F(energy,pitch,R,Z)"
        endif
        write(*,'(T2,"Nenergy = ",i3)') fbm%nenergy
        write(*,'(T2,"Npitch  = ",i3)') fbm%npitch
        write(*,'(T2,"Nr  = ",i3)') fbm%nr
        write(*,'(T2,"Nz  = ",i3)') fbm%nz
        if(fbm%nphi.gt.1) then
            write(*,'(T2,"Nphi  = ",i3)') fbm%nphi
        endif
        write(*,'(T2,"Energy range = [",f5.2,",",f6.2,"]")') fbm%emin,fbm%emax
        write(*,'(T2,"Pitch  range = [",f5.2,",",f5.2,"]")') fbm%pmin,fbm%pmax
        write(*,'(T2,"R  range = [",f6.2,",",f6.2,"]")') fbm%rmin,fbm%rmax
        write(*,'(T2,"Z  range = [",f7.2,",",f6.2,"]")') fbm%zmin,fbm%zmax
        if(fbm%nphi.gt.1) then
            write(*,'(T2,"Phi  range = [",f5.2,",",f5.2,"]")') fbm%phimin,fbm%phimax
        endif
        write(*,'(T2,"Ntotal = ",ES10.3)') fbm%n_tot
        write(*,*) ''
    endif

end subroutine read_f

subroutine read_mc(fid, error)
    !+ Reads in a MC particle fast-ion distribution and puts them in [[libfida:particles]]
    integer(HID_T), intent(inout) :: fid
        !+ HDF5 file ID
    integer, intent(out)          :: error
        !+ Error code

    integer(HSIZE_T), dimension(1) :: dims
    integer(Int32) :: i,j,ii,ir,iz,iphi,nphi
    real(Float64) :: phi,beam_grid_phi_enter,beam_grid_phi_exit,delta_phi,xp,yp,zp
    real(Float64), dimension(3) :: uvw,xyz,ri,vi,e1_xyz,e2_xyz,C_xyz,dum
    real(Float64), dimension(:), allocatable :: weight
    type(LocalEMFields) :: fields
    integer :: cnt,num
    logical :: inp,path_valid
    character(len=50) :: dist_type_name = ''

    if(inputs%verbose.ge.1) then
        write(*,'(a)') '---- Fast-ion distribution settings ----'
    endif

    call h5ltread_dataset_int_scalar_f(fid, "/nparticle", particles%nparticle, error)

    !!ALLOCATE SPACE
    allocate(particles%fast_ion(particles%nparticle))
    allocate(weight(particles%nparticle))

    dims(1) = particles%nparticle
    !call h5ltread_dataset_double_f(fid, "/A", particles%fast_ion%A, dims, error)
    particles%fast_ion%A = beam_mass
    call h5ltread_dataset_double_f(fid, "/r", particles%fast_ion%r, dims, error)
    call h5ltread_dataset_double_f(fid, "/z", particles%fast_ion%z, dims, error)
    call h5ltpath_valid_f(fid, "/phi", .True., path_valid, error)
    if(path_valid) then
        call h5ltread_dataset_double_f(fid, "/phi", particles%fast_ion%phi, dims, error)
        particles%axisym = .False.
    endif
    call h5ltread_dataset_int_f(fid, "/class", particles%fast_ion%class, dims, error)

    if(any(particles%fast_ion%class.le.0)) then
        if(inputs%verbose.ge.0) then
            write(*,'(a)') 'READ_MC: Orbit class ID must be greater than 0'
        endif
        stop
    endif

    if(inputs%split.ge.1) then
        call h5ltread_dataset_int_scalar_f(fid, "/nclass", particles%nclass, error)
        if(any(particles%fast_ion%class.gt.particles%nclass)) then
            if(inputs%verbose.ge.0) then
                write(*,'(a)') 'READ_MC: Orbit class ID greater than the number of classes'
            endif
            stop
        endif
    endif

    if(inputs%dist_type.eq.2) then
        if(particles%axisym) then
            dist_type_name = "Axisymmetric Guiding Center Monte Carlo"
        else
            dist_type_name = "Non-axisymmetric Guiding Center Monte Carlo"
        endif
        call h5ltread_dataset_double_f(fid, "/energy", particles%fast_ion%energy, dims, error)
        call h5ltread_dataset_double_f(fid, "/pitch", particles%fast_ion%pitch, dims, error)
        particles%fast_ion%vabs  = sqrt(particles%fast_ion%energy/(v2_to_E_per_amu*particles%fast_ion%A))
    else
        if(particles%axisym) then
            dist_type_name = "Axisymmetric Full Orbit Monte Carlo"
        else
            dist_type_name = "Non-axisymmetric Full Orbit Monte Carlo"
        endif
        call h5ltread_dataset_double_f(fid, "/vr", particles%fast_ion%vr, dims, error)
        call h5ltread_dataset_double_f(fid, "/vt", particles%fast_ion%vt, dims, error)
        call h5ltread_dataset_double_f(fid, "/vz", particles%fast_ion%vz, dims, error)
        particles%fast_ion%vabs = sqrt(particles%fast_ion%vr**2 + &
                                       particles%fast_ion%vt**2 + &
                                       particles%fast_ion%vz**2)
        particles%fast_ion%energy = v2_to_E_per_amu*particles%fast_ion%A*particles%fast_ion%vabs**2
    endif

    call h5ltread_dataset_double_f(fid, "/weight", weight, dims, error)

    cnt=0
    e1_xyz = matmul(beam_grid%inv_basis,[1.0,0.0,0.0])
    e2_xyz = matmul(beam_grid%inv_basis,[0.0,1.0,0.0])
    !$OMP PARALLEL DO schedule(guided) private(i,ii,j,ir,iz,iphi,fields,uvw,phi,ri,vi, &
    !$OMP& delta_phi,beam_grid_phi_enter,beam_grid_phi_exit,C_xyz,xyz,xp,yp,zp,dum,inp)
    particle_loop: do i=1,particles%nparticle
        if(inputs%verbose.ge.2) then
            WRITE(*,'(f7.2,"% completed",a,$)') cnt/real(particles%nparticle)*100,char(13)
        endif

        if(particles%axisym) then
            uvw = [particles%fast_ion(i)%r, 0.d0 , particles%fast_ion(i)%z]
        else
            xp = particles%fast_ion(i)%r * cos(particles%fast_ion(i)%phi)
            yp = particles%fast_ion(i)%r * sin(particles%fast_ion(i)%phi)
            zp = particles%fast_ion(i)%z
            uvw = [xp,yp,zp]
        endif
        call in_plasma(uvw,inp,input_coords=1)
        if(.not.inp) cycle particle_loop

        if(particles%axisym) then
            beam_grid_phi_enter = 0.0
            beam_grid_phi_exit = 0.0
            dum = [0.d0, 0.d0, particles%fast_ion(i)%z]
            call uvw_to_xyz(dum, C_xyz)
            call circle_grid_intersect(C_xyz,e1_xyz,e2_xyz,particles%fast_ion(i)%r,beam_grid_phi_enter,beam_grid_phi_exit)
            delta_phi = beam_grid_phi_exit-beam_grid_phi_enter
            if(delta_phi.gt.0) then
                particles%fast_ion(i)%beam_grid_cross_grid = .True.
            else
                particles%fast_ion(i)%beam_grid_cross_grid = .False.
                delta_phi = 2*pi
            endif
            particles%fast_ion(i)%beam_grid_phi_enter = beam_grid_phi_enter
        else
            delta_phi = 2*pi
            call uvw_to_xyz(uvw,xyz)
            particles%fast_ion(i)%beam_grid_cross_grid = in_grid(xyz)
            particles%fast_ion(i)%beam_grid_phi_enter = particles%fast_ion(i)%phi
        endif
        particles%fast_ion(i)%delta_phi = delta_phi
        particles%fast_ion(i)%weight = weight(i)

        ir = minloc(abs(inter_grid%r - particles%fast_ion(i)%r),1)
        iphi = minloc(abs(inter_grid%phi - particles%fast_ion(i)%phi),1)
        iz = minloc(abs(inter_grid%z - particles%fast_ion(i)%z),1)

        !$OMP ATOMIC UPDATE
        equil%plasma(ir,iz,iphi)%denf = equil%plasma(ir,iz,iphi)%denf + weight(i) / &
                                   (particles%fast_ion(i)%r*inter_grid%dv)
        !$OMP END ATOMIC
        cnt=cnt+1
    enddo particle_loop
    !$OMP END PARALLEL DO

    num = count(particles%fast_ion%beam_grid_cross_grid)
    if(num.le.0) then
        if(inputs%verbose.ge.0) then
            write(*,'(a)') 'READ_MC: No mc particles in beam grid'
        endif
        stop
    endif

    if(inputs%verbose.ge.1) then
        write(*,'(T2,"Distribution type: ",a)') dist_type_name
        write(*,'(T2,"Number of mc particles: ",i10)') particles%nparticle
        write(*,'(T2,"Number of orbit classes: ",i6)') particles%nclass
        write(*,*) ''
    endif

end subroutine read_mc

subroutine read_distribution
    !+ Reads in the fast-ion distribution
    integer(HID_T) :: fid
    integer :: error

    !!Initialize HDF5 interface
    call h5open_f(error)

    !!Open HDF5 file
    call h5fopen_f(inputs%distribution_file, H5F_ACC_RDONLY_F, fid, error)

    !!Get distribution type
    call h5ltread_dataset_int_scalar_f(fid, "/type", inputs%dist_type, error)

    if(inputs%dist_type.eq.1) then
        call read_f(fid, error)
    else !2 or 3
        call read_mc(fid, error)
    endif

    !!Close file
    call h5fclose_f(fid, error)

    !!Close HDF5 interface
    call h5close_f(error)

end subroutine read_distribution

subroutine quasineutrality_check
    !+ Checks whether quasi-neutrality is satisfied to some tol
    real(Float64) :: tol = 0.01
    real(Float64), dimension(:,:,:), allocatable :: quasi, deni
    integer :: i, quasi_cnt

    allocate(deni(inter_grid%nr,inter_grid%nz,inter_grid%nphi))
    deni = 0.d0
    do i=1,n_thermal
            deni = deni + equil%plasma%deni(i)
    enddo

    !! Quasi-neutrality check
    allocate(quasi(inter_grid%nr,inter_grid%nz,inter_grid%nphi))
    where(equil%mask.ge.0.5)
        quasi = equil%plasma%dene - (deni + impurity_charge*equil%plasma%denimp + equil%plasma%denf)
    end where

    quasi_cnt = count(abs(quasi).gt.(tol*equil%plasma%dene))
    if(inputs%verbose.ge.1) then
        write(*,'(a)') '---- Quasi-neutrality Check ----'
        write(*,'(T2, a, f5.2, "%")') 'Tolerance: ',tol*100
        write(*,'(T2, a, f6.2, "%")') 'Percent Failure: ', &
                                  100.0*quasi_cnt/size(quasi)
        write(*,*) ''
    endif
    deallocate(deni,quasi)

end subroutine quasineutrality_check

subroutine read_atomic_cross(fid, grp, cross)
    !+ Reads in a cross section table from file
    !+ and puts it into a [[AtomicCrossSection]] type
    integer(HID_T), intent(in)              :: fid
        !+ HDF5 file ID
    character(len=*), intent(in)            :: grp
        !+ HDF5 group to read from
    type(AtomicCrossSection), intent(inout) :: cross
        !+ Atomic cross section

    integer(HSIZE_T), dimension(3) :: dim3
    real(Float64) :: emin, emax, rmin
    integer :: i, n_max, m_max, error
    real(Float64), dimension(:,:,:), allocatable :: dummy3
    logical :: path_valid

    call h5ltpath_valid_f(fid, grp, .True., path_valid, error)
    if(.not.path_valid) then
        if(inputs%verbose.ge.0) then
            write(*,'(a,a)') 'READ_ATOMIC_CROSS: Unknown atomic interaction: ', trim(grp)
        endif
        stop
    endif

    call h5ltread_dataset_int_scalar_f(fid, grp//"/nenergy", cross%nenergy, error)
    call h5ltread_dataset_int_scalar_f(fid, grp//"/n_max", n_max, error)
    call h5ltread_dataset_int_scalar_f(fid, grp//"/m_max", m_max, error)
    call h5ltread_dataset_double_scalar_f(fid,grp//"/emin", emin, error)
    call h5ltread_dataset_double_scalar_f(fid,grp//"/emax", emax, error)
    call h5ltread_dataset_double_scalar_f(fid,grp//"/dlogE", cross%dlogE, error)
    cross%logemin = log10(emin)
    cross%logemax = log10(emax)

    allocate(dummy3(n_max, m_max, cross%nenergy))

    allocate(cross%log_cross(cross%m_max,cross%n_max, cross%nenergy))

    dim3 = [n_max, m_max, cross%nenergy]
    call h5ltread_dataset_double_f(fid,grp//"/cx", dummy3, dim3, error)
    rmin = minval(dummy3,dummy3.gt.0.d0)
    where (dummy3.le.0.0)
        dummy3 = 0.9*rmin
    end where
    cross%minlog_cross = log10(rmin)
    do i=1, cross%nenergy
        cross%log_cross(:,:,i) = log10(transpose(dummy3(1:nlevs,1:nlevs,i)))
    enddo
    deallocate(dummy3)

end subroutine read_atomic_cross

subroutine read_atomic_rate(fid, grp, rates)
    !+ Reads in a atomic rate table from file
    !+ and puts it into a [[AtomicRates]] type
    integer(HID_T), intent(in)              :: fid
        !+ HDF5 file ID
    character(len=*), intent(in)            :: grp
        !+ HDF5 group to read from
    type(AtomicRates), intent(inout)        :: rates
        !+ Atomic reaction rates

    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(3) :: dim3
    integer(HSIZE_T), dimension(4) :: dim4
    integer(HSIZE_T), dimension(5) :: dim5
    logical :: path_valid
    integer :: i, j, n, n_max, m_max, error
    integer :: drank
    real(Float64) :: emin,emax,tmin,tmax,rmin
    real(Float64), dimension(:,:), allocatable :: dummy2
    real(Float64), dimension(:,:,:), allocatable :: dummy3
    real(Float64), dimension(:,:,:,:), allocatable :: dummy4
    real(Float64), dimension(:,:,:,:,:), allocatable :: dummy5

    call h5ltpath_valid_f(fid, grp, .True., path_valid, error)
    if(.not.path_valid) then
        if(inputs%verbose.ge.0) then
            write(*,'(a,a)') 'READ_ATOMIC_RATE: Unknown atomic interaction: ', trim(grp)
        endif
        stop
    endif

    call h5ltread_dataset_int_scalar_f(fid, grp//"/n_max", n_max, error)
    call h5ltread_dataset_int_scalar_f(fid, grp//"/m_max", m_max, error)
    call h5ltread_dataset_int_scalar_f(fid, grp//"/nenergy", rates%nenergy, error)
    call h5ltread_dataset_double_scalar_f(fid, grp//"/emin", emin, error)
    call h5ltread_dataset_double_scalar_f(fid, grp//"/emax", emax, error)
    call h5ltread_dataset_double_scalar_f(fid, grp//"/dlogE", rates%dlogE, error)
    call h5ltread_dataset_int_scalar_f(fid, grp//"/ntemp", rates%ntemp, error)
    call h5ltread_dataset_double_scalar_f(fid, grp//"/tmin", tmin, error)
    call h5ltread_dataset_double_scalar_f(fid, grp//"/tmax", tmax, error)
    call h5ltread_dataset_double_scalar_f(fid, grp//"/dlogT", rates%dlogT, error)
    rates%logemin = log10(emin)
    rates%logemax = log10(emax)
    rates%logtmin = log10(tmin)
    rates%logtmax = log10(tmax)


    allocate(rates%log_rate(&
                    rates%m_max, &
                    rates%n_max, &
                    rates%nenergy, &
                    rates%ntemp))
    rates%log_rate = 0.d0

    !!Read CX
    call h5ltpath_valid_f(fid, grp//"/cx", .True., path_valid, error)
    if(path_valid) then
        call h5ltget_dataset_ndims_f(fid, grp//"/cx", drank, error)
        if(drank.eq.4) then
            allocate(dummy4(n_max, m_max, &
                           rates%nenergy, &
                           rates%ntemp))
            dim4 = [n_max, m_max, rates%nenergy, rates%ntemp]
            call h5ltread_dataset_double_f(fid, grp//"/cx", dummy4, dim4, error)

            do j=1,rates%ntemp
                do i=1,rates%nenergy
                        rates%log_rate(:,:,i,j) = transpose(dummy4(1:nlevs,1:nlevs,i,j))
                        rates%log_rate(:,:,i,j) = transpose(dummy4(1:nlevs,1:nlevs,i,j))
                enddo
            enddo
            deallocate(dummy4)
        else
            if(inputs%verbose.ge.0) then
                write(*,'(a,a)') 'READ_ATOMIC_RATE: Unsupported atomic interaction: ', trim(grp)
            endif
            stop
        endif
    endif

    rmin = minval(rates%log_rate, rates%log_rate.gt.0.d0)
    where (rates%log_rate.le.0.d0)
        rates%log_rate = 0.9*rmin
    end where
    rates%minlog_rate = log10(rmin)
    rates%log_rate = log10(rates%log_rate)

end subroutine read_atomic_rate

subroutine read_atomic_transitions(fid, grp,rates)
    !+ Reads in a atomic transitions table from file
    !+ and puts it into a [[AtomicTransitions]] type
    integer(HID_T), intent(in)              :: fid
        !+ HDF5 file ID
    character(len=*), intent(in)            :: grp
        !+ HDF5 group to read from
    type(AtomicTransitions), intent(inout)  :: rates
        !+ Atomic transitions

    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(3) :: dim3
    integer(HSIZE_T), dimension(4) :: dim4
    integer(HSIZE_T), dimension(5) :: dim5
    logical :: path_valid
    integer :: i, j, n, n_max, m_max, error
    integer :: drank
    real(Float64) :: emin,emax,tmin,tmax,rmin
    real(Float64), dimension(:,:), allocatable :: dummy2
    real(Float64), dimension(:,:,:), allocatable :: dummy3
    real(Float64), dimension(:,:,:,:), allocatable :: dummy4
    real(Float64), dimension(:,:,:,:,:), allocatable :: dummy5

    call h5ltpath_valid_f(fid, grp, .True., path_valid, error)
    if(.not.path_valid) then
        if(inputs%verbose.ge.0) then
            write(*,'(a,a)') 'READ_ATOMIC_TRANSITIONS: Unknown atomic interaction: ', trim(grp)
        endif
        stop
    endif

    call h5ltread_dataset_int_scalar_f(fid, grp//"/n_max", n_max, error)
    call h5ltread_dataset_int_scalar_f(fid, grp//"/m_max", m_max, error)
    call h5ltread_dataset_int_scalar_f(fid, grp//"/nenergy", rates%nenergy, error)
    call h5ltread_dataset_double_scalar_f(fid, grp//"/emin", emin, error)
    call h5ltread_dataset_double_scalar_f(fid, grp//"/emax", emax, error)
    call h5ltread_dataset_double_scalar_f(fid, grp//"/dlogE", rates%dlogE, error)
    call h5ltread_dataset_int_scalar_f(fid, grp//"/ntemp", rates%ntemp, error)
    call h5ltread_dataset_double_scalar_f(fid, grp//"/tmin", tmin, error)
    call h5ltread_dataset_double_scalar_f(fid, grp//"/tmax", tmax, error)
    call h5ltread_dataset_double_scalar_f(fid, grp//"/dlogT", rates%dlogT, error)
    rates%logemin = log10(emin)
    rates%logemax = log10(emax)
    rates%logtmin = log10(tmin)
    rates%logtmax = log10(tmax)

    allocate(rates%log_pop(&
                    rates%m_max, &
                    rates%n_max, &
                    rates%nenergy, &
                    rates%ntemp))
    allocate(rates%log_depop(&
                    rates%n_max, &
                    rates%nenergy, &
                    rates%ntemp))
    rates%log_pop = 0.d0
    rates%log_depop = 0.d0

    !!Read CX
    call h5ltpath_valid_f(fid, grp//"/cx", .True., path_valid, error)
    if(path_valid) then
        call h5ltget_dataset_ndims_f(fid, grp//"/cx", drank, error)
        if(drank.eq.3) then
            allocate(dummy3(n_max, &
                           rates%nenergy, &
                           rates%ntemp))
            dim3 = [n_max, rates%nenergy, rates%ntemp]
            call h5ltread_dataset_double_f(fid, grp//"/cx", dummy3, dim3, error)
            do j=1,rates%ntemp
                do i=1,rates%nenergy
                    do n=1,rates%n_max
                        rates%log_depop(n,i,j) = dummy3(n,i,j)
                    enddo
                enddo
            enddo
            deallocate(dummy3)
        endif
        if(drank.eq.4) then
            allocate(dummy4(n_max, m_max, &
                           rates%nenergy, &
                           rates%ntemp))
            dim4 = [n_max, m_max, rates%nenergy, rates%ntemp]
            call h5ltread_dataset_double_f(fid, grp//"/cx", dummy4, dim4, error)
            do j=1,rates%ntemp
                do i=1,rates%nenergy
                    do n=1,rates%n_max
                        rates%log_depop(n,i,j) = sum(dummy4(n,:,i,j))
                    enddo
                enddo
            enddo
            deallocate(dummy4)
        endif
    endif

    !!Read ionization
    call h5ltpath_valid_f(fid, grp//"/ionization", .True., path_valid, error)
    if(path_valid) then
        allocate(dummy3(n_max, &
                       rates%nenergy, &
                       rates%ntemp))
        dim3 = [n_max, rates%nenergy, rates%ntemp]
        call h5ltread_dataset_double_f(fid, grp//"/ionization", dummy3, dim3, error)
        do j=1,rates%ntemp
            do i=1,rates%nenergy
                do n=1,rates%n_max
                    rates%log_depop(n,i,j) = rates%log_depop(n,i,j) + &
                                               dummy3(n,i,j)
                enddo
            enddo
        enddo
        deallocate(dummy3)
    endif

    !!Read excitation
    call h5ltpath_valid_f(fid, grp//"/excitation", .True., path_valid, error)
    if(path_valid) then
        allocate(dummy4(n_max, m_max,&
                       rates%nenergy, &
                       rates%ntemp))
        dim4 = [n_max, m_max, rates%nenergy, rates%ntemp]
        call h5ltread_dataset_double_f(fid, grp//"/excitation", dummy4, dim4, error)
        do j=1,rates%ntemp
            do i=1,rates%nenergy
                rates%log_pop(:,:,i,j) = transpose(dummy4(1:nlevs,1:nlevs,i,j))
                do n=1,rates%n_max
                    rates%log_depop(n,i,j) = rates%log_depop(n,i,j) + &
                                               sum(dummy4(n,:,i,j))
                enddo
            enddo
        enddo
        deallocate(dummy4)
    endif

    rmin = minval(rates%log_depop, rates%log_depop.gt.0.d0)
    where (rates%log_depop.le.0.d0)
        rates%log_depop = 0.9*rmin
    end where
    rates%minlog_depop = log10(rmin)
    rates%log_depop = log10(rates%log_depop)

    rmin = minval(rates%log_pop, rates%log_pop.gt.0.d0)
    where (rates%log_pop.le.0.d0)
        rates%log_pop = 0.9*rmin
    end where
    rates%minlog_pop = log10(rmin)
    rates%log_pop = log10(rates%log_pop)

end subroutine read_atomic_transitions

subroutine read_nuclear_rates(fid, grp, rates)
    !+ Reads in a nuclear reaction rates table from file
    !+ and puts it into a [[NuclearRates]] type
    integer(HID_T), intent(in)              :: fid
        !+ HDF5 file ID
    character(len=*), intent(in)            :: grp
        !+ HDF5 group to read from
    type(NuclearRates), intent(inout)       :: rates
        !+ Atomic reaction rates

    integer(HSIZE_T), dimension(1) :: dim1
    integer(HSIZE_T), dimension(3) :: dim3
    logical :: path_valid, err
    integer :: i, j, error
    real(Float64) :: emin, emax, tmin, tmax, rmin

    err = .False.
    call h5ltpath_valid_f(fid, grp, .True., path_valid, error)
    if(.not.path_valid) then
        if(inputs%verbose.ge.0) then
            write(*,'(a,a)') 'READ_NUCLEAR_RATES: Unknown nuclear interaction: ', trim(grp)
            write(*,'(a)') 'Continuing without neutron calculation'
        endif
        inputs%calc_neutron=0
        inputs%calc_cfpd=0
        return
    endif

    dim1 = [2]
    call h5ltread_dataset_double_f(fid, grp//"/bt_amu", rates%bt_amu, dim1, error)

    if(abs(beam_mass-rates%bt_amu(1)).gt.0.2) then
        if(inputs%verbose.ge.0) then
            write(*,'(a,f6.3,a,f6.3,a)') 'READ_NUCLEAR_RATES: Unexpected beam species mass. Expected ',&
                rates%bt_amu(1),' amu got ', beam_mass, ' amu'
        endif
        err = .True.
    endif

    if(abs(thermal_mass(1)-rates%bt_amu(2)).gt.0.2) then
        if(inputs%verbose.ge.0) then
            write(*,'(a,f6.3,a,f6.3,a)') 'READ_NUCLEAR_RATES: Unexpected thermal species mass. Expected ',&
                 rates%bt_amu(2),' amu got ', thermal_mass(1), ' amu'
        endif
     !!!err = .True.
    endif

    if(err) then
        if(inputs%verbose.ge.0) then
            write(*,'(a)') 'Continuing without neutron calculation'
        endif
        inputs%calc_neutron=0
        inputs%calc_cfpd=0
        return
    endif

    call h5ltread_dataset_int_scalar_f(fid, grp//"/nbranch", rates%nbranch, error)
    call h5ltread_dataset_int_scalar_f(fid, grp//"/nenergy", rates%nenergy, error)
    call h5ltread_dataset_double_scalar_f(fid, grp//"/emin", emin, error)
    call h5ltread_dataset_double_scalar_f(fid, grp//"/emax", emax, error)
    call h5ltread_dataset_double_scalar_f(fid, grp//"/dlogE", rates%dlogE, error)
    call h5ltread_dataset_int_scalar_f(fid, grp//"/ntemp", rates%ntemp, error)
    call h5ltread_dataset_double_scalar_f(fid, grp//"/tmin", tmin, error)
    call h5ltread_dataset_double_scalar_f(fid, grp//"/tmax", tmax, error)
    call h5ltread_dataset_double_scalar_f(fid, grp//"/dlogT", rates%dlogT, error)
    rates%logemin = log10(emin)
    rates%logemax = log10(emax)
    rates%logtmin = log10(tmin)
    rates%logtmax = log10(tmax)

    allocate(rates%log_rate(rates%nenergy, &
                            rates%ntemp,   &
                            rates%nbranch))

    dim3 = [rates%nenergy, rates%ntemp, rates%nbranch]
    call h5ltread_dataset_double_f(fid, grp//"/fusion", rates%log_rate, dim3, error)

    rmin = minval(rates%log_rate, rates%log_rate.gt.0.d0)
    where (rates%log_rate.le.0.d0)
        rates%log_rate = 0.9*rmin
    end where
    rates%minlog_rate = log10(rmin)
    rates%log_rate = log10(rates%log_rate)

end subroutine read_nuclear_rates

subroutine read_tables
    !+ Reads in atomic tables from file and stores them in [[libfida:tables]]
    integer(HID_T) :: fid
    integer(HSIZE_T), dimension(2) :: dim2
    integer :: error

    integer :: n_max, m_max
    character(len=4) :: impname
    real(Float64), dimension(2) :: b_amu
    real(Float64), dimension(:,:), allocatable :: dummy2

    !!Initialize HDF5 interface
    call h5open_f(error)

    !!Open HDF5 file
    call h5fopen_f(inputs%tables_file, H5F_ACC_RDONLY_F, fid, error)

    !!Read Hydrogen-Hydrogen CX Cross Sections
    call read_atomic_cross(fid,"/cross/H_H",tables%H_H_cx_cross)

    !!Read Hydrogen-Hydrogen CX Rates
    call read_atomic_rate(fid,"/rates/H_H", tables%H_H_cx_rate)

    !!Read Hydrogen-Hydrogen Transitions
    call read_atomic_transitions(fid,"/rates/H_H", tables%H_H)

    !!Read Hydrogen-Electron Transitions
    call read_atomic_transitions(fid,"/rates/H_e", tables%H_e)

    !!Read Hydrogen-Impurity Transitions
    impname = ''
    select case (impurity_charge)
        case (2)
            impname = "He2"
        case (5)
            impname = "B5"
        case (6)
            impname = "C6"
        case DEFAULT
            write(impname,'("A",i1)') impurity_charge
    end select

    call read_atomic_transitions(fid,"/rates/H_"//trim(adjustl(impname)), tables%H_Aq)

    !!Read Einstein coefficients
    call h5ltread_dataset_int_scalar_f(fid,"/rates/spontaneous/n_max", n_max, error)
    call h5ltread_dataset_int_scalar_f(fid,"/rates/spontaneous/m_max", m_max, error)
    allocate(dummy2(n_max,m_max))
    dim2 = [n_max, m_max]
    call h5ltread_dataset_double_f(fid,"/rates/spontaneous/einstein",dummy2, dim2, error)
    tables%einstein(:,:) = transpose(dummy2(1:nlevs,1:nlevs))
    deallocate(dummy2)

    !!Read nuclear Deuterium-Deuterium rates
    if(inputs%calc_neutron.ge.1.or.inputs%calc_cfpd.ge.1) then
        call read_nuclear_rates(fid, "/rates/D_D", tables%D_D)
    endif

    !!Close file
    call h5fclose_f(fid, error)

    !!Close HDF5 interface
    call h5close_f(error)

    if(inputs%verbose.ge.1) then
        write(*,'(a)') "---- Atomic tables settings ----"
        write(*,'(T2,"Maximum n/m: ",i2)') nlevs
        write(*,*) ''
    endif

end subroutine read_tables

subroutine write_beam_grid(id, error)
    !+ Write [[libfida:beam_grid]] to an HDF5 file
    integer(HID_T), intent(inout) :: id
        !+ HDF5 file ID
    integer, intent(out)          :: error
        !+ Error code

    integer(HID_T) :: gid
    integer(HSIZE_T), dimension(3) :: dims
    real(Float64), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz) :: u_grid, v_grid, w_grid
    real(Float64) :: xyz(3),uvw(3)
    integer :: i,j,k

    !Create uvw grids
    do k=1, beam_grid%nz
        do j=1, beam_grid%ny
            do i=1, beam_grid%nx
                xyz = [beam_grid%xc(i), &
                       beam_grid%yc(j), &
                       beam_grid%zc(k)]
                call xyz_to_uvw(xyz,uvw)
                u_grid(i,j,k) = uvw(1)
                v_grid(i,j,k) = uvw(2)
                w_grid(i,j,k) = uvw(3)
           enddo
        enddo
    enddo

    !Create grid group
    call h5gcreate_f(id, "grid", gid, error)

    !Write variables
    dims(1) = 1
    call h5ltmake_dataset_int_f(gid,"nx", 0, dims(1:1), [beam_grid%nx], error)
    call h5ltmake_dataset_int_f(gid,"ny", 0, dims(1:1), [beam_grid%ny], error)
    call h5ltmake_dataset_int_f(gid,"nz", 0, dims(1:1), [beam_grid%nz], error)

    dims = [beam_grid%nx, beam_grid%ny, beam_grid%nz]
    call h5ltmake_compressed_dataset_double_f(gid,"x", 1, dims(1:1), beam_grid%xc, error)
    call h5ltmake_compressed_dataset_double_f(gid,"y", 1, dims(2:2), beam_grid%yc, error)
    call h5ltmake_compressed_dataset_double_f(gid,"z", 1, dims(3:3), beam_grid%zc, error)

    call h5ltmake_compressed_dataset_double_f(gid,"x_grid", 3, dims, u_grid, error)
    call h5ltmake_compressed_dataset_double_f(gid,"y_grid", 3, dims, v_grid, error)
    call h5ltmake_compressed_dataset_double_f(gid,"z_grid", 3, dims, w_grid, error)

    !Write attributes
    call h5ltset_attribute_string_f(gid,"nx","description", &
         "Number of cells in the X direction", error)
    call h5ltset_attribute_string_f(gid,"ny","description", &
         "Number of cells in the Y direction", error)
    call h5ltset_attribute_string_f(gid,"nz","description", &
         "Number of cells in the Z direction", error)
    call h5ltset_attribute_string_f(gid,"x","description", &
         "X value of cell center in beam grid coordinates", error)
    call h5ltset_attribute_string_f(gid,"x","units", "cm", error)
    call h5ltset_attribute_string_f(gid,"y","description", &
         "Y value of cell center in beam grid coordinates", error)
    call h5ltset_attribute_string_f(gid,"y","units", "cm", error)
    call h5ltset_attribute_string_f(gid,"z","description", &
         "Z value of cell center in beam grid coordinates", error)
    call h5ltset_attribute_string_f(gid,"z","units", "cm", error)

    call h5ltset_attribute_string_f(gid,"x_grid","description", &
         "X value of cell center in machine coordinates: x_grid(x,y,z)", error)
    call h5ltset_attribute_string_f(gid,"x_grid","units", "cm", error)
    call h5ltset_attribute_string_f(gid,"y_grid","description", &
         "Y value of cell center in machine coordinates: y_grid(x,y,z)", error)
    call h5ltset_attribute_string_f(gid,"y_grid","units", "cm", error)
    call h5ltset_attribute_string_f(gid,"z_grid","description", &
         "Z value of cell center in machine coordinates: z_grid(x,y,z)", error)
    call h5ltset_attribute_string_f(gid,"z_grid","units", "cm", error)

    call h5ltset_attribute_string_f(id,"grid","coordinate_system", &
         "Right-handed cartesian",error)

    !Close grid group
    call h5gclose_f(gid, error)

end subroutine write_beam_grid

subroutine write_birth_profile
    !+ Writes [[libfida:birth]] to a HDF5 file
    integer(HID_T) :: fid
    integer(HSIZE_T), dimension(4) :: dim4
    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(1) :: d
    integer :: error, i, c, npart, start_index, end_index

    character(charlim) :: filename
    type(LocalEMFields) :: fields
    type(BirthParticle) :: part
    real(Float64), dimension(:,:), allocatable :: ri, ri_gc
    real(Float64), dimension(:,:), allocatable :: vi
    real(Float64), dimension(:), allocatable :: energy, pitch, weight
    integer, dimension(:,:), allocatable :: inds
    integer, dimension(:), allocatable :: neut_types
    real(Float64), dimension(3) :: xyz,v_xyz,uvw,v_uvw,r_gyro,uvw_gc
    logical :: do_write

#ifdef _MPI
    integer :: rank
    integer, dimension(:), allocatable :: npart_image
#endif

    npart = birth%cnt-1
#ifdef _MPI
    call parallel_sum(npart)
#endif

    allocate(ri(3,npart))
    allocate(vi(3,npart))
    allocate(ri_gc(3,npart))
    allocate(energy(npart),pitch(npart),weight(npart))
    allocate(inds(3,npart))
    allocate(neut_types(npart))
    ri = 0.d0
    vi = 0.d0
    ri_gc = 0.d0
    energy = 0.d0
    pitch = 0.d0
    weight = 0.d0
    inds = 0
    neut_types = 0

    c = 0
#ifdef _MPI
    rank = my_rank()
    allocate(npart_image(0:num_ranks()-1))
    npart_image(:) = 0
    npart_image(rank) = birth%cnt - 1
    call parallel_sum(npart_image)

    do i=0,rank-1
        c = c +  npart_image(i)
    enddo
    deallocate(npart_image)
#endif
    start_index = 1 + c
    end_index = birth%cnt - 1 + c

    c = 1
    do i=start_index, end_index
        part = birth%part(c)
        ! Convert position to rzphi
        xyz = part%ri
        v_xyz = part%vi

        inds(:,i) = part%ind
        neut_types(i) = part%neut_type
        energy(i) = part%energy
        weight(i) = part%weight

        ! Get guiding center positions
        pitch(i) = part%pitch
        call xyz_to_uvw(part%ri_gc, uvw_gc)
        ri_gc(1,i) = sqrt(uvw_gc(1)*uvw_gc(1) + uvw_gc(2)*uvw_gc(2))
        ri_gc(2,i) = uvw_gc(3)
        ri_gc(3,i) = atan2(uvw_gc(2),uvw_gc(1))

        ! Get position in cylindrical
        call xyz_to_uvw(xyz,uvw)
        ri(1,i) = sqrt(uvw(1)*uvw(1) + uvw(2)*uvw(2))
        ri(2,i) = uvw(3)
        ri(3,i) = atan2(uvw(2),uvw(1))

        ! Convert velocity to cylindrical
        v_uvw = matmul(beam_grid%basis, v_xyz)
        vi(1,i) = v_uvw(1)*cos(ri(3,i)) + v_uvw(2)*sin(ri(3,i))
        vi(2,i) = v_uvw(3)
        vi(3,i) = -v_uvw(1)*sin(ri(3,i)) + v_uvw(2)*cos(ri(3,i))
        c = c + 1
    enddo

    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_birth.h5"

    do_write = .True.
#ifdef _MPI
    call parallel_sum(ri_gc)
    call parallel_sum(energy)
    call parallel_sum(pitch)
    call parallel_sum(weight)
    call parallel_sum(ri)
    call parallel_sum(vi)
    call parallel_sum(inds)
    call parallel_sum(neut_types)
    if(my_rank().ne.0) do_write = .False.
#endif

    if(do_write) then
        !Open HDF5 interface
        call h5open_f(error)

        !Create file overwriting any existing file
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, fid, error)

        !Write variables
        call write_beam_grid(fid, error)
        d(1) = 1
        call h5ltmake_dataset_int_f(fid, "/n_birth", 0, d, [npart], error)
        dim4 = shape(birth%dens)
        call h5ltmake_compressed_dataset_double_f(fid,"/dens", 4, dim4, birth%dens, error)
        dim2 = [3, npart]
        call h5ltmake_compressed_dataset_double_f(fid,"/ri_gc", 2, dim2, ri_gc, error)
        call h5ltmake_compressed_dataset_double_f(fid,"/ri", 2, dim2, ri, error)
        call h5ltmake_compressed_dataset_double_f(fid,"/vi", 2, dim2, vi, error)
        call h5ltmake_compressed_dataset_int_f(fid,"/ind", 2, dim2, inds, error)
        call h5ltmake_compressed_dataset_double_f(fid,"/energy", 1, dim2(2:2), energy, error)
        call h5ltmake_compressed_dataset_double_f(fid,"/pitch", 1, dim2(2:2), pitch, error)
        call h5ltmake_compressed_dataset_double_f(fid,"/weight", 1, dim2(2:2), weight, error)
        call h5ltmake_compressed_dataset_int_f(fid,"/type", 1, dim2(2:2), neut_types, error)

        !Add attributes
        call h5ltset_attribute_string_f(fid, "/n_birth","description", &
             "Number of birth mc particles deposited", error)
        call h5ltset_attribute_string_f(fid, "/dens", "description", &
             "Birth density: dens(beam_component,x,y,z)", error)
        call h5ltset_attribute_string_f(fid, "/dens", "units", &
             "fast-ions/(s*cm^3)", error)
        call h5ltset_attribute_string_f(fid, "/ri_gc", "description", &
             "Fast-ion guiding-center birth position in R-Z-Phi: ri_gc([r,z,phi],particle)", error)
        call h5ltset_attribute_string_f(fid, "/ri_gc", "units", "cm, radians", error)
        call h5ltset_attribute_string_f(fid, "/ri", "description", &
             "Fast-ion birth position in R-Z-Phi: ri([r,z,phi],particle)", error)
        call h5ltset_attribute_string_f(fid, "/ri", "units", "cm, radians", error)
        call h5ltset_attribute_string_f(fid, "/vi", "description", &
             "Fast-ion birth velocity in R-Z-Phi: vi([r,z,phi],particle)", error)
        call h5ltset_attribute_string_f(fid, "/vi", "units", "cm/s", error)
        call h5ltset_attribute_string_f(fid, "/energy", "description", &
             "Fast-ion birth energy: energy(particle)", error)
        call h5ltset_attribute_string_f(fid, "/energy", "units", "keV", error)
        call h5ltset_attribute_string_f(fid, "/weight", "description", &
             "Fast-ion birth weight: weight(particle)", error)
        call h5ltset_attribute_string_f(fid, "/weight", "units", "fast-ions/s", error)
        call h5ltset_attribute_string_f(fid, "/pitch", "description", &
             "Fast-ion birth pitch w.r.t. the magnetic field: pitch(particle)", error)
        call h5ltset_attribute_string_f(fid, "/ind", "description", &
             "Fast-ion birth beam grid indices: ind([i,j,k],particle)", error)
        call h5ltset_attribute_string_f(fid, "/type", "description", &
             "Fast-ion birth type (1=Full, 2=Half, 3=Third)", error)

        call h5ltset_attribute_string_f(fid, "/", "coordinate_system", &
             "Cylindrical (R,Z,Phi)",error)
        call h5ltset_attribute_string_f(fid, "/", "version", version, error)
        call h5ltset_attribute_string_f(fid, "/", "description", &
             "Birth density and particles calculated by FIDASIM", error)

        !!Close file
        call h5fclose_f(fid, error)

        !!Close HDF5 interface
        call h5close_f(error)
    endif

    ! Deallocate arrays since they aren't needed anymore
    deallocate(ri,vi,ri_gc,energy,pitch,neut_types,inds)
    deallocate(birth%dens,birth%part)

    if(inputs%verbose.ge.1) then
        write(*,'(T4,a,a)') 'birth profile written to: ',trim(filename)
    endif

end subroutine write_birth_profile

subroutine write_neutral_population(id, pop, error)
    !+ Writes Neutral Population to HDF5 group
    integer(HID_T), intent(inout)       :: id
        !+ HDF5 group ID
    type(NeutralPopulation), intent(in) :: pop
        !+ Neutral Population to write
    integer, intent(out)                :: error
        !+ Error code

    integer(HID_T) :: gid
    integer(HSIZE_T), dimension(1) :: d
    integer(HSIZE_T), dimension(4) :: dims4
    integer(HSIZE_T), dimension(5) :: dims5

    real(Float64), dimension(:,:,:,:,:), allocatable :: v
    real(Float64), dimension(:,:,:,:), allocatable   :: w
    integer(Int32), dimension(:,:,:), allocatable    :: n

    integer :: ic, ir, i, j, k, ind(3), nx, ny, nz, nk

    nx = beam_grid%nx
    ny = beam_grid%ny
    nz = beam_grid%nz
    dims4 = [nlevs, nx, ny, nz]
    d(1) =1

    call h5ltmake_dataset_int_f(id,"nlevel", 0, d, [nlevs], error)
    call h5ltset_attribute_string_f(id,"nlevel","description", &
         "Number of atomic energy levels", error)

    call h5ltmake_compressed_dataset_double_f(id, "dens", 4, dims4, pop%dens, error)
    call h5ltset_attribute_string_f(id,"dens","units","neutrals*cm^-3",error)
    call h5ltset_attribute_string_f(id,"dens","description", &
         "Neutral density dens(level,x,y,z)", error)

    dims5 = [3, reservoir_size, nx, ny, nz]
    allocate(v(3,reservoir_size, nx, ny, nz))
    allocate(w(reservoir_size, nx, ny, nz))
    allocate(n(nx, ny, nz))
    v = 0.d0; w = 0.d0
    n = 0

    do ic=1,beam_grid%ngrid
        call ind2sub(beam_grid%dims,ic,ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)
        n(i,j,k) = pop%res(i,j,k)%n
        nk = min(pop%res(i,j,k)%k, pop%res(i,j,k)%n)
        w(1:nk,i,j,k) = pop%res(i,j,k)%R(1:nk)%w
        do ir=1, nk
            v(:,ir,i,j,k) = pop%res(i,j,k)%R(ir)%v
        enddo
    enddo

    !Create reservoir group
    if(inputs%output_neutral_reservoir.eq.1) then
       call h5gcreate_f(id, "reservoir", gid, error)
       call h5ltset_attribute_string_f(id,"reservoir","description", &
            "Neutral Particle Reservoir",error)

       call h5ltmake_dataset_int_f(gid,"k", 0, d, [maxval(pop%res%k)], error)
       call h5ltset_attribute_string_f(gid,"k","description", &
            "Reservoir Size", error)

       call h5ltmake_compressed_dataset_int_f(gid, "n", 3, dims4(2:4), n, error)
       call h5ltset_attribute_string_f(gid, "n", "description", &
            "Number of particles in each reservoir", error)

       call h5ltmake_compressed_dataset_double_f(gid, "w", 4, dims5(2:5), w, error, compress=.False.)
       call h5ltset_attribute_string_f(gid, "w", "description", &
            "Neutral Particle Weight", error)

       call h5ltmake_compressed_dataset_double_f(gid, "v", 5, dims5, v, error, compress=.False.)
       call h5ltset_attribute_string_f(gid,"v","units","cm/s",error)
       call h5ltset_attribute_string_f(gid,"v","description", &
            "Neutral Particle velocity in beam grid coordinates v(:,particle,i,j,k)", error)
       !Close reservoir group
       call h5gclose_f(gid, error)
    endif

    deallocate(v,w,n)

end subroutine write_neutral_population

subroutine write_neutrals
    !+ Writes [[libfida:neut]] to a HDF5 file
    integer(HID_T) :: fid, gid
    integer :: error

    integer :: i
    character(charlim) :: filename
    integer(HSIZE_T), dimension(4) :: dims4

    dims4 = [nlevs, beam_grid%nx, beam_grid%ny, beam_grid%nz]

    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_neutrals.h5"

    !Open HDF5 interface
    call h5open_f(error)

    !Create file overwriting any existing file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, fid, error)

    !Write variables
    call write_beam_grid(fid, error)

    if(inputs%calc_nbi_dens.ge.1) then
        call h5gcreate_f(fid, "/full", gid, error)
        call write_neutral_population(gid, neut%full, error)
        call h5gclose_f(gid, error)
        call h5ltset_attribute_string_f(fid,"/full","description", &
             "Full Energy Neutral Population", error)

        call h5gcreate_f(fid, "/half", gid, error)
        call write_neutral_population(gid, neut%half, error)
        call h5gclose_f(gid, error)
        call h5ltset_attribute_string_f(fid,"/half","description", &
             "Half Energy Neutral Population", error)

        call h5gcreate_f(fid, "/third", gid, error)
        call write_neutral_population(gid, neut%third, error)
        call h5gclose_f(gid, error)
        call h5ltset_attribute_string_f(fid,"/third","description", &
             "Third Energy Neutral Population", error)
    endif

    if(inputs%calc_dcx_dens.ge.1) then
        call h5gcreate_f(fid, "/dcx", gid, error)
        call write_neutral_population(gid, neut%dcx, error)
        call h5gclose_f(gid, error)
        call h5ltset_attribute_string_f(fid,"/dcx","description", &
             "Direct Charge Exchange (DCX) Neutral Population", error)
    endif

    if(inputs%calc_halo_dens.ge.1) then
        call h5gcreate_f(fid, "/halo", gid, error)
        call write_neutral_population(gid, neut%halo, error)
        call h5gclose_f(gid, error)
        call h5ltset_attribute_string_f(fid,"/halo","description", &
             "Thermal Halo Neutral Population", error)
    endif

    call h5ltset_attribute_string_f(fid, "/", "version", version, error)
    call h5ltset_attribute_string_f(fid,"/","description", &
         "Neutral Populations calculated by FIDASIM", error)

    !Close file
    call h5fclose_f(fid, error)

    !Close HDF5 interface
    call h5close_f(error)

    if(inputs%verbose.ge.1) then
        write(*,'(T4,a,a)') 'neutral populations written to: ',trim(filename)
    endif

end subroutine write_neutrals

subroutine write_npa
    !+ Writes [[libfida:npa]] to a HDF5 file
    integer(HID_T) :: fid, gid
    integer(HSIZE_T), dimension(3) :: dim3
    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(1) :: d
    integer :: error

    integer, dimension(:), allocatable :: dcount
    real(Float64), dimension(:,:), allocatable :: ri, rf
    real(Float64), dimension(:), allocatable :: weight, energy, pitch
    integer, dimension(:), allocatable :: det, orbit_type
    integer :: i, npart, c, start_index, end_index
    character(charlim) :: filename = ''
    logical :: do_write = .True.

#ifdef _MPI
    integer :: rank
    integer, dimension(:), allocatable :: npart_image

    rank = my_rank()
    allocate(npart_image(0:num_ranks()-1))
#endif

#ifdef _MPI
    if(my_rank().ne.0) do_write = .False.
#endif

    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_npa.h5"

    if(do_write) then
        !Open HDF5 interface
        call h5open_f(error)

        !Create file overwriting any existing file
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, fid, error)
    endif

    !! Active
    allocate(dcount(npa_chords%nchan))
    npart = npa%npart
    if(npart.gt.0) then
        do i=1,npa_chords%nchan
            dcount(i) = count(npa%part%detector.eq.i)
        enddo
    else
        dcount = 0
    endif

#ifdef _MPI
    call parallel_sum(dcount)
    call parallel_sum(npart)
#endif

    c = 0
#ifdef _MPI
    npart_image(:) = 0
    npart_image(rank) = npa%npart
    call parallel_sum(npart_image)

    do i=0,rank-1
        c = c +  npart_image(i)
    enddo
#endif
    start_index = 1 + c
    end_index = npa%npart + c

    if(npart.gt.0) then
        allocate(ri(3,npart),rf(3,npart))
        allocate(weight(npart),energy(npart),pitch(npart))
        allocate(det(npart),orbit_type(npart))
        ri = 0.d0     ; rf = 0.d0
        weight = 0.d0 ; energy = 0.d0
        pitch = 0.d0  ; det = 0
        orbit_type = 0
        c = 1
        do i=start_index, end_index
            ri(1,i) = npa%part(c)%xi
            ri(2,i) = npa%part(c)%yi
            ri(3,i) = npa%part(c)%zi
            rf(1,i) = npa%part(c)%xf
            rf(2,i) = npa%part(c)%yf
            rf(3,i) = npa%part(c)%zf
            weight(i) = npa%part(c)%weight
            energy(i) = npa%part(c)%energy
            pitch(i) = npa%part(c)%pitch
            det(i) = npa%part(c)%detector
            orbit_type(i) = npa%part(c)%class
            c = c + 1
        enddo
#ifdef _MPI
        call parallel_sum(ri)
        call parallel_sum(rf)
        call parallel_sum(weight)
        call parallel_sum(energy)
        call parallel_sum(pitch)
        call parallel_sum(det)
        call parallel_sum(orbit_type)
#endif
    endif

    if(do_write.and.(inputs%calc_npa.ge.1)) then
        !Write Active Flux
        d(1) = 1
        dim2 = [npa%nenergy, npa%nchan]
        dim3 = [npa%nenergy, npa%nchan, particles%nclass]
        if(particles%nclass.gt.1) then
            call h5ltmake_dataset_int_f(fid,"/nclass", 0, d, [particles%nclass], error)
            call h5ltmake_compressed_dataset_double_f(fid,"/flux",3,dim3,npa%flux, error)
            call h5ltset_attribute_string_f(fid,"/flux", "description", &
                 "Active Neutral flux: flux(energy,chan,class)", error)
        else
            call h5ltmake_compressed_dataset_double_f(fid,"/flux",2,dim3(1:2),npa%flux(:,:,1), error)
            call h5ltset_attribute_string_f(fid,"/flux", "description", &
                 "Active Neutral flux: flux(energy,chan)", error)
        endif
        call h5ltset_attribute_string_f(fid,"/flux", "units","neutrals/(s*dE)", error)

        call h5ltmake_dataset_int_f(fid,"/nenergy", 0, d, [npa%nenergy], error)
        call h5ltmake_dataset_int_f(fid,"/nchan", 0, d, [npa%nchan], error)
        call h5ltmake_compressed_dataset_double_f(fid,"/energy",1,dim2(1:1),&
             npa%energy, error)
        call h5ltmake_compressed_dataset_double_f(fid,"/radius",1,dim2(2:2),&
             npa_chords%radius, error)
        call h5ltmake_compressed_dataset_int_f(fid,"/count",1,dim2(2:2), dcount, error)

        !Add attributes
        call h5ltset_attribute_string_f(fid, "/", "version", version, error)
        call h5ltset_attribute_string_f(fid,"/","description", &
             "NPA flux calculated by FIDASIM",error)
        call h5ltset_attribute_string_f(fid,"/nenergy","description",&
             "Number of energy values",error)
        call h5ltset_attribute_string_f(fid,"/nchan","description",&
             "Number of channels",error)
        call h5ltset_attribute_string_f(fid,"/energy","description", &
             "Energy array", error)
        call h5ltset_attribute_string_f(fid,"/energy","units","keV", error)
        call h5ltset_attribute_string_f(fid,"/radius","description", &
             "Detector line of sight radius at midplane or tangency point", error)
        call h5ltset_attribute_string_f(fid,"/radius","units","cm",error)
        call h5ltset_attribute_string_f(fid,"/count","description", &
             "Number of particles that hit the detector: count(chan)", error)


        if((npart.gt.0).and.(inputs%calc_npa.ge.2)) then
            !Create Group
            call h5gcreate_f(fid,"/particles",gid, error)
            call h5ltmake_dataset_int_f(gid, "nparticle", 0, d, [npart], error)
            d(1) = npart
            dim2 = [3, npart]
            call h5ltmake_compressed_dataset_double_f(gid,"ri",2,dim2, ri, error)
            call h5ltmake_compressed_dataset_double_f(gid,"rf",2,dim2, rf, error)
            call h5ltmake_compressed_dataset_double_f(gid,"pitch",1,d, pitch, error)
            call h5ltmake_compressed_dataset_double_f(gid,"energy",1,d, energy, error)
            call h5ltmake_compressed_dataset_double_f(gid,"weight",1,d, weight, error)
            call h5ltmake_compressed_dataset_int_f(gid,"detector",1,d, det, error)
            call h5ltmake_compressed_dataset_int_f(gid,"class",1,d, orbit_type, error)

            !Add attributes
            call h5ltset_attribute_string_f(gid,"nparticle","description", &
                 "Number of particles that hit a detector", error)
            call h5ltset_attribute_string_f(gid,"ri","description", &
                 "Neutral particle's birth position in machine coordinates: ri([x,y,z],particle)", error)
            call h5ltset_attribute_string_f(gid,"ri","units", "cm", error)
            call h5ltset_attribute_string_f(gid,"rf","description", &
                 "Neutral particle's hit position in machine coordinates: rf([x,y,z],particle)", error)
            call h5ltset_attribute_string_f(gid,"rf","units", "cm", error)
            call h5ltset_attribute_string_f(gid,"pitch","description", &
                 "Pitch value of the neutral particle: p = v_parallel/v  w.r.t. the magnetic field", error)
            call h5ltset_attribute_string_f(gid,"energy","description", &
                 "Energy value of the neutral particle", error)
            call h5ltset_attribute_string_f(gid,"energy","units","keV",error)
            call h5ltset_attribute_string_f(gid,"weight","description", &
                 "Neutral particle's contribution to the flux", error)
            call h5ltset_attribute_string_f(gid,"weight","units","neutrals/s",error)
            call h5ltset_attribute_string_f(gid,"detector","description", &
                 "Detector that the neutral particle hit", error)
            call h5ltset_attribute_string_f(gid,"class","description", &
                 "Class of the neutral particle", error)

            call h5ltset_attribute_string_f(fid,"/particles","coordinate_system", &
                 "Right-handed cartesian",error)
            call h5ltset_attribute_string_f(fid,"/particles","description", &
                 "Active NPA Monte Carlo particles",error)

            !Close group
            call h5gclose_f(gid, error)
        endif
    endif

    deallocate(dcount)
    if(npart.gt.0) then
        deallocate(ri,rf)
        deallocate(energy,pitch,weight,det,orbit_type)
    endif

    !! Passive
    allocate(dcount(npa_chords%nchan))
    npart = pnpa%npart
    if(npart.gt.0) then
        do i=1,npa_chords%nchan
            dcount(i) = count(pnpa%part%detector.eq.i)
        enddo
    else
        dcount = 0
    endif

#ifdef _MPI
    call parallel_sum(dcount)
    call parallel_sum(npart)
#endif

    c = 0
#ifdef _MPI
    npart_image(:) = 0
    npart_image(rank) = pnpa%npart
    call parallel_sum(npart_image)

    do i=0,rank-1
        c = c +  npart_image(i)
    enddo
#endif
    start_index = 1 + c
    end_index = pnpa%npart + c

    if(npart.gt.0) then
        allocate(ri(3,npart),rf(3,npart))
        allocate(weight(npart),energy(npart),pitch(npart))
        allocate(det(npart),orbit_type(npart))
        ri = 0.d0     ; rf = 0.d0
        weight = 0.d0 ; energy = 0.d0
        pitch = 0.d0  ; det = 0
        orbit_type = 0
        c = 1
        do i=start_index, end_index
            ri(1,i) = pnpa%part(c)%xi
            ri(2,i) = pnpa%part(c)%yi
            ri(3,i) = pnpa%part(c)%zi
            rf(1,i) = pnpa%part(c)%xf
            rf(2,i) = pnpa%part(c)%yf
            rf(3,i) = pnpa%part(c)%zf
            weight(i) = pnpa%part(c)%weight
            energy(i) = pnpa%part(c)%energy
            pitch(i) = pnpa%part(c)%pitch
            det(i) = pnpa%part(c)%detector
            orbit_type(i) = pnpa%part(c)%class
            c = c + 1
        enddo
#ifdef _MPI
        call parallel_sum(ri)
        call parallel_sum(rf)
        call parallel_sum(weight)
        call parallel_sum(energy)
        call parallel_sum(pitch)
        call parallel_sum(det)
        call parallel_sum(orbit_type)
#endif
    endif

    if(do_write.and.(inputs%calc_pnpa.ge.1)) then
        !Write Passive Flux
        d(1) = 1
        dim2 = [pnpa%nenergy, pnpa%nchan]
        dim3 = [pnpa%nenergy, pnpa%nchan, particles%nclass]
        if(particles%nclass.gt.1) then
            call h5ltmake_compressed_dataset_double_f(fid,"/pflux",3,dim3,pnpa%flux, error)
            call h5ltset_attribute_string_f(fid,"/pflux", "description", &
                 "Passive Neutral flux: pflux(energy,chan,class)", error)
        else
            call h5ltmake_compressed_dataset_double_f(fid,"/pflux",2,dim3(1:2),pnpa%flux(:,:,1), error)
            call h5ltset_attribute_string_f(fid,"/pflux", "description", &
                 "Passive Neutral flux: pflux(energy,chan)", error)
        endif
        call h5ltset_attribute_string_f(fid,"/pflux", "units","neutrals/(s*dE)", error)
        call h5ltmake_compressed_dataset_int_f(fid,"/pcount",1,dim2(2:2), dcount, error)
        call h5ltset_attribute_string_f(fid,"/pcount","description", &
             "Number of passive particles that hit the detector: pcount(chan)", error)

        if(inputs%calc_npa.le.0) then
            call h5ltmake_dataset_int_f(fid,"/nenergy", 0, d, [pnpa%nenergy], error)
            call h5ltmake_dataset_int_f(fid,"/nchan", 0, d, [pnpa%nchan], error)
            call h5ltmake_compressed_dataset_double_f(fid,"/energy",1,dim2(1:1),&
                 pnpa%energy, error)
            call h5ltmake_compressed_dataset_double_f(fid,"/radius",1,dim2(2:2),&
                 npa_chords%radius, error)

            !Add attributes
            call h5ltset_attribute_string_f(fid, "/", "version", version, error)
            call h5ltset_attribute_string_f(fid,"/","description", &
                 "NPA flux calculated by FIDASIM",error)
            call h5ltset_attribute_string_f(fid,"/nenergy","description",&
                 "Number of energy values",error)
            call h5ltset_attribute_string_f(fid,"/nchan","description",&
                 "Number of channels",error)
            call h5ltset_attribute_string_f(fid,"/energy","description", &
                 "Energy array", error)
            call h5ltset_attribute_string_f(fid,"/energy","units","keV", error)
            call h5ltset_attribute_string_f(fid,"/radius","description", &
                 "Detector line of sight radius at midplane or tangency point", error)
            call h5ltset_attribute_string_f(fid,"/radius","units","cm",error)
        endif

        if((npart.gt.0).and.(inputs%calc_pnpa.ge.2)) then
            !Create Group
            call h5gcreate_f(fid,"/passive_particles",gid, error)
            call h5ltmake_dataset_int_f(gid, "nparticle", 0, d, [npart], error)
            d(1) = npart
            dim2 = [3, npart]
            call h5ltmake_compressed_dataset_double_f(gid,"ri",2,dim2, ri, error)
            call h5ltmake_compressed_dataset_double_f(gid,"rf",2,dim2, rf, error)
            call h5ltmake_compressed_dataset_double_f(gid,"pitch",1,d, pitch, error)
            call h5ltmake_compressed_dataset_double_f(gid,"energy",1,d, energy, error)
            call h5ltmake_compressed_dataset_double_f(gid,"weight",1,d, weight, error)
            call h5ltmake_compressed_dataset_int_f(gid,"detector",1,d, det, error)
            call h5ltmake_compressed_dataset_int_f(gid,"class",1,d, orbit_type, error)

            !Add attributes
            call h5ltset_attribute_string_f(gid,"nparticle","description", &
                 "Number of particles that hit a detector", error)
            call h5ltset_attribute_string_f(gid,"ri","description", &
                 "Neutral particle's birth position in machine coordinates: ri([x,y,z],particle)", error)
            call h5ltset_attribute_string_f(gid,"ri","units", "cm", error)
            call h5ltset_attribute_string_f(gid,"rf","description", &
                 "Neutral particle's hit position in machine coordinates: rf([x,y,z],particle)", error)
            call h5ltset_attribute_string_f(gid,"rf","units", "cm", error)
            call h5ltset_attribute_string_f(gid,"pitch","description", &
                 "Pitch value of the neutral particle: p = v_parallel/v  w.r.t. the magnetic field", error)
            call h5ltset_attribute_string_f(gid,"energy","description", &
                 "Energy value of the neutral particle", error)
            call h5ltset_attribute_string_f(gid,"energy","units","keV",error)
            call h5ltset_attribute_string_f(gid,"weight","description", &
                 "Neutral particle's contribution to the flux", error)
            call h5ltset_attribute_string_f(gid,"weight","units","neutrals/s",error)
            call h5ltset_attribute_string_f(gid,"detector","description", &
                 "Detector that the neutral particle hit", error)
            call h5ltset_attribute_string_f(gid,"class","description", &
                 "Class of the neutral particle", error)

            call h5ltset_attribute_string_f(fid,"/passive_particles","coordinate_system", &
                 "Right-handed cartesian",error)
            call h5ltset_attribute_string_f(fid,"/passive_particles","description", &
                 "Passive NPA Monte Carlo particles",error)

            !Close group
            call h5gclose_f(gid, error)
        endif
    endif

    if(do_write) then
        !Close file
        call h5fclose_f(fid, error)

        !Close HDF5 interface
        call h5close_f(error)
    endif

    if(inputs%verbose.ge.1) then
        write(*,'(T4,a,a)') 'NPA data written to: ',trim(filename)
    endif

#ifdef _MPI
    deallocate(npart_image)
#endif

end subroutine write_npa

subroutine write_spectra
    !+ Writes [[libfida:spectra]] to a HDF5 file
    integer(HID_T) :: fid, gid, sid
    integer(HSIZE_T), dimension(4) :: dims
    integer(HSIZE_T), dimension(5) :: dims_stokes
    integer(HSIZE_T), dimension(1) :: d
    integer :: error

    character(charlim) :: filename
    integer :: i, ir, ic
    real(Float64) :: factor
    real(Float64), dimension(:), allocatable :: lambda_arr

    real(Float64), dimension(inputs%nlambda,spec_chords%nchan) :: full, half, third
    real(Float64), dimension(4,inputs%nlambda,spec_chords%nchan) :: fullstokes, halfstokes, thirdstokes
    real(Float64), dimension(inputs%nlambda,spec_chords%nchan,n_thermal) :: dcx, halo, cold
    real(Float64), dimension(4,inputs%nlambda,spec_chords%nchan,n_thermal) :: dcxstokes, halostokes, coldstokes
    real(Float64), dimension(inputs%nlambda,spec_chords%nchan,particles%nclass) :: fida, pfida
    real(Float64), dimension(4,inputs%nlambda,spec_chords%nchan,particles%nclass) :: fidastokes, pfidastokes

    real(Float64), dimension(3,reservoir_size,spec_chords%nchan) :: spat
    real(Float64), dimension(reservoir_size,spec_chords%nchan) :: photons
    real(Float64), dimension(reservoir_size,spec_chords%nchan) :: distance
    real(Float64) :: src(3),axis(3),pos(3)
    integer, dimension(n_stark) :: stark_sign
    stark_sign = +1*stark_sigma - 1*stark_pi

    allocate(lambda_arr(inputs%nlambda))
    do i=1,inputs%nlambda
        lambda_arr(i) = (i-0.5)*inputs%dlambda + inputs%lambdamin
    enddo

    !! convert [Ph/(s*wavel_bin*cm^2*all_directions)] to [Ph/(s*nm*sr*m^2)]!
    factor = 1.d0/(inputs%dlambda)/(4.d0*pi)*1.d4

    !! write to file
    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_spectra.h5"

    !Open HDF5 interface
    call h5open_f(error)

    !Create file overwriting any existing file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, fid, error)
    call h5gcreate_f(fid, "stokes", gid, error)
    call h5ltset_attribute_string_f(fid,"stokes","description", &
            "Stoke Parameters",error)

    !Write variables
    d(1) = 1
    call h5ltmake_dataset_int_f(fid, "/nchan", 0, d, [spec_chords%nchan], error)
    call h5ltmake_dataset_int_f(fid, "/nlambda", 0, d, [inputs%nlambda], error)

    dims(1) = n_stark
    dims(2) = inputs%nlambda
    dims(3) = spec_chords%nchan
    dims(4) = n_thermal
    dims_stokes(1) = n_stark
    dims_stokes(2) = 4
    dims_stokes(3) = inputs%nlambda
    dims_stokes(4) = spec_chords%nchan
    dims_stokes(5) = n_thermal

    if(inputs%stark_components.ge.1) then
        call h5ltmake_dataset_int_f(fid, "nstark", 0, d, [n_stark], error)
        call h5ltmake_compressed_dataset_int_f(fid, "stark_sign", 1, dims(1:1), stark_sign, error)
        call h5ltset_attribute_string_f(fid,"stark_sign", "description", &
         "Stark line indicator: 1=sigma, -1=pi ", error)
    endif

    call h5ltmake_compressed_dataset_double_f(fid, "/lambda", 1, dims(2:2), &
         lambda_arr, error)
    call h5ltmake_compressed_dataset_double_f(fid, "/radius", 1, dims(3:3), &
         spec_chords%radius, error)

    !Add attributes
    call h5ltset_attribute_string_f(fid,"/nchan", "description", &
         "Number of channels", error)
    call h5ltset_attribute_string_f(fid,"/nlambda", "description", &
         "Number of wavelengths", error)
    call h5ltset_attribute_string_f(fid,"/lambda","description", &
         "Wavelength array", error)
    call h5ltset_attribute_string_f(fid,"/lambda","units","nm", error)
    call h5ltset_attribute_string_f(fid,"/radius", "description", &
         "Line of sight radius at midplane or tangency point", error)
    call h5ltset_attribute_string_f(fid,"/radius","units","cm", error)

    if(inputs%calc_brems.ge.1) then
        spec%brems = factor*spec%brems
        !Write variables
        call h5ltmake_compressed_dataset_double_f(fid, "/brems", 2, &
             dims(2:3), spec%brems, error)
        !Add attributes
        call h5ltset_attribute_string_f(fid,"/brems","description", &
             "Visible Bremsstrahlung: brems(lambda,chan)", error)
        call h5ltset_attribute_string_f(fid,"/brems","units",&
             "Ph/(s*nm*sr*m^2)",error )
    endif

    if(inputs%calc_bes.ge.1) then
        spec%full = factor*spec%full
        spec%fullstokes = factor*spec%fullstokes
        spec%half = factor*spec%half
        spec%halfstokes = factor*spec%halfstokes
        spec%third = factor*spec%third
        spec%thirdstokes = factor*spec%thirdstokes
        if (inputs%stark_components.eq.0) then
            full  = sum(spec%full, dim=1)
            half  = sum(spec%half, dim=1)
            third = sum(spec%third, dim=1)
            fullstokes = sum(spec%fullstokes,dim=1)
            halfstokes = sum(spec%halfstokes, dim=1)
            thirdstokes = sum(spec%thirdstokes, dim=1)
            !Write variables
            call h5ltmake_compressed_dataset_double_f(fid, "/full", 2, dims(2:3), &
                 full, error)
            call h5ltmake_compressed_dataset_double_f(fid, "/half", 2, dims(2:3), &
                 half, error)
            call h5ltmake_compressed_dataset_double_f(fid, "/third", 2, dims(2:3),&
                 third, error)

            call h5ltmake_compressed_dataset_double_f(gid, "full", 3, dims_stokes(2:4), &
                 fullstokes, error)
            call h5ltmake_compressed_dataset_double_f(gid, "half", 3, dims_stokes(2:4), &
                 halfstokes, error)
            call h5ltmake_compressed_dataset_double_f(gid, "third", 3, dims_stokes(2:4),&
                 thirdstokes, error)

            !Add attributes
            call h5ltset_attribute_string_f(fid,"/full","description", &
                 "Full energy component of the beam emmision: full(lambda,chan)", error)
            call h5ltset_attribute_string_f(fid,"/full","units","Ph/(s*nm*sr*m^2)",error )
            call h5ltset_attribute_string_f(fid,"/half","description", &
                 "Half energy component of the beam emmision: half(lambda,chan)", error)
            call h5ltset_attribute_string_f(fid,"/half","units","Ph/(s*nm*sr*m^2)",error )
            call h5ltset_attribute_string_f(fid,"/third","description", &
                 "Third energy component of the beam emmision: third(lambda,chan)", error)
            call h5ltset_attribute_string_f(fid,"/third","units","Ph/(s*nm*sr*m^2)",error )

            call h5ltset_attribute_string_f(gid,"full","description", &
                 "Full energy component of the beam emmision stokes parameters: &
                 &full(4,lambda,chan)", error)
            call h5ltset_attribute_string_f(gid,"full","units","Ph/(s*nm*sr*m^2)",error )
            call h5ltset_attribute_string_f(gid,"half","description", &
                 "Half energy component of the beam emmision stokes parameters: &
                 &half(4,lambda,chan)", error)
            call h5ltset_attribute_string_f(gid,"half","units","Ph/(s*nm*sr*m^2)",error )
            call h5ltset_attribute_string_f(gid,"third","description", &
                 "Third energy component of the beam emmision stokes parameters: &
                 &third(4,lambda,chan)", error)
            call h5ltset_attribute_string_f(gid,"third","units","Ph/(s*nm*sr*m^2)",error )
        else
            call h5ltmake_compressed_dataset_double_f(fid, "/full", 3, dims(1:3), &
                 spec%full, error)
            call h5ltmake_compressed_dataset_double_f(fid, "/half", 3, dims(1:3), &
                 spec%half, error)
            call h5ltmake_compressed_dataset_double_f(fid, "/third", 3, dims(1:3),&
                 spec%third, error)

            call h5ltmake_compressed_dataset_double_f(gid, "full", 4, dims_stokes(1:4), &
                spec%fullstokes, error)
            call h5ltmake_compressed_dataset_double_f(gid, "half", 4, dims_stokes(1:4), &
                spec%halfstokes, error)
            call h5ltmake_compressed_dataset_double_f(gid, "third", 4, dims_stokes(1:4),&
                spec%thirdstokes, error)

            call h5ltset_attribute_string_f(fid,"/full","description", &
                 "Full energy component of the beam emmision stark components: full(stark,lambda,chan)", &
                 error)
            call h5ltset_attribute_string_f(fid,"/full","units","Ph/(s*nm*sr*m^2)",error )
            call h5ltset_attribute_string_f(fid,"/half","description", &
                 "Half energy component of the beam emmision stark components: half(stark,lambda,chan)", &
                 error)
            call h5ltset_attribute_string_f(fid,"/half","units","Ph/(s*nm*sr*m^2)",error )
            call h5ltset_attribute_string_f(fid,"/third","description", &
                 "Third energy component of the beam emmision stark components: third(stark,lambda,chan)", &
                 error)
            call h5ltset_attribute_string_f(fid,"/third","units","Ph/(s*nm*sr*m^2)",error )

            call h5ltset_attribute_string_f(gid,"full","description", &
                 "Full energy component of the beam emmision stark components stokes parameters: &
                 &full(stark,4,lambda,chan)", &
                 error)
            call h5ltset_attribute_string_f(gid,"full","units","Ph/(s*nm*sr*m^2)",error )
            call h5ltset_attribute_string_f(gid,"half","description", &
                 "Half energy component of the beam emmision stark components stokes parameters: &
                 &half(stark,4,lambda,chan)", &
                 error)
            call h5ltset_attribute_string_f(gid,"half","units","Ph/(s*nm*sr*m^2)",error )
            call h5ltset_attribute_string_f(gid,"third","description", &
                 "Third energy component of the beam emmision stark components stokes parameters: &
                 &third(stark,4,lambda,chan)", &
                 error)
            call h5ltset_attribute_string_f(gid,"third","units","Ph/(s*nm*sr*m^2)",error )
        endif
    endif

    if(inputs%calc_dcx.ge.1) then
        spec%dcx = factor*spec%dcx
        spec%dcxstokes = factor*spec%dcxstokes
        if (inputs%stark_components.eq.0) then
            dcx = sum(spec%dcx, dim=1)
            dcxstokes = sum(spec%dcxstokes, dim=1)
            call h5ltmake_compressed_dataset_double_f(fid, "/dcx", 3, dims(2:4), &
                 dcx, error)
            call h5ltset_attribute_string_f(fid,"/dcx","description", &
                 "Direct Charge Exchange (DCX) emission: dcx(lambda,chan,species)", error)
            call h5ltset_attribute_string_f(fid,"/dcx","units","Ph/(s*nm*sr*m^2)",error )

            call h5ltmake_compressed_dataset_double_f(gid,"dcx", 4, dims_stokes(2:5),&
                 dcxstokes, error)
            call h5ltset_attribute_string_f(gid,"dcx","description", &
                 "Direct Charge Exchange (DCX) emission stokes parameters: &
                 &dcx(4,lambda,chan,species)", error)
            call h5ltset_attribute_string_f(gid,"dcx","units","Ph/(s*nm*sr*m^2)",error )
        else
            call h5ltmake_compressed_dataset_double_f(fid, "/dcx", 4, dims, &
                 spec%dcx, error)

            call h5ltmake_compressed_dataset_double_f(gid, "dcx", 5, dims_stokes, &
                 spec%dcxstokes, error)

            call h5ltset_attribute_string_f(fid,"/dcx","description", &
                 "Direct Charge Exchange (DCX) emission stark components: dcx(stark,lambda,chan,species)", error)
            call h5ltset_attribute_string_f(fid,"/dcx","units","Ph/(s*nm*sr*m^2)",error )

            call h5ltset_attribute_string_f(gid,"dcx","description", &
                 "Direct Charge Exchange (DCX) emission stark components stokes parameters: &
                 &dcx(stark,4,lambda,chan,species)", error)
            call h5ltset_attribute_string_f(gid,"dcx","units","Ph/(s*nm*sr*m^2)",error )
        endif
    endif

    if(inputs%calc_halo.ge.1) then
        spec%halo = factor*spec%halo
        spec%halostokes = factor*spec%halostokes
        if (inputs%stark_components.eq.0) then
            halo = sum(spec%halo, dim=1)
            halostokes = sum(spec%halostokes, dim=1)
            call h5ltmake_compressed_dataset_double_f(fid, "/halo", 3, dims(2:4), &
                 halo, error)
            call h5ltset_attribute_string_f(fid,"/halo","description", &
                 "Halo component of the beam emmision: halo(lambda,chan,species)", error)
            call h5ltset_attribute_string_f(fid,"/halo","units","Ph/(s*nm*sr*m^2)",error )

            call h5ltmake_compressed_dataset_double_f(gid, "halo", 4, dims_stokes(2:5), &
                 halostokes, error)
            call h5ltset_attribute_string_f(gid,"halo","description", &
                 "Halo component of the beam emmision stokes parameters: halo(4,lambda,chan,species)", error)
            call h5ltset_attribute_string_f(gid,"halo","units","Ph/(s*nm*sr*m^2)",error )
        else
            call h5ltmake_compressed_dataset_double_f(fid, "/halo", 4, dims, &
                 spec%halo, error)
            call h5ltset_attribute_string_f(fid,"/halo","description", &
                 "Halo component of the beam emmision stark components: halo(stark,lambda,chan,species)", error)
            call h5ltset_attribute_string_f(fid,"/halo","units","Ph/(s*nm*sr*m^2)",error )

            call h5ltmake_compressed_dataset_double_f(gid, "halo", 5, dims_stokes(1:5), &
                 spec%halostokes, error)
            call h5ltset_attribute_string_f(gid,"halo","description", &
                 "Halo component of the beam emmision stark components stokes parameters: &
                 &halo(stark,4,lambda,chan,species)", error)
            call h5ltset_attribute_string_f(gid,"halo","units","Ph/(s*nm*sr*m^2)",error )
        endif
    endif

    if((inputs%calc_halo.ge.1) .or. (inputs%calc_dcx.ge.1)) then
       call h5ltmake_compressed_dataset_double_f(fid, "/thermal_mass", 1, dims(4:4), thermal_mass(1:n_thermal), error)
       call h5ltset_attribute_string_f(fid,"/thermal_mass","description", &
            "Mass of each of the thermal ions: thermal_mass(species)", error)
       call h5ltset_attribute_string_f(fid,"/thermal_mass", "units", "amu", error)

       call h5ltmake_compressed_dataset_double_f(fid, "/thermal_lambda0", 1, dims(4:4), thermal_lambda0(1:n_thermal), error)
       call h5ltset_attribute_string_f(fid,"/thermal_lambda0","description", &
            "Rest wavelength of the thermal species lines: thermal_lambda0(species)", error)
       call h5ltset_attribute_string_f(fid,"/thermal_lambda0", "units", "nm", error)
    endif


    if(inputs%calc_cold.ge.1) then
        spec%cold = factor*spec%cold
        spec%coldstokes = factor*spec%coldstokes
        if (inputs%stark_components.eq.0) then
            cold = sum(spec%cold, dim=1)
            coldstokes = sum(spec%coldstokes, dim=1)
            call h5ltmake_compressed_dataset_double_f(fid, "/cold", 3, dims(2:4), &
                 cold, error)
            call h5ltset_attribute_string_f(fid,"/cold","description", &
                 "Cold D-alpha emission: cold(lambda,chan,species)", error)
            call h5ltset_attribute_string_f(fid,"/cold","units","Ph/(s*nm*sr*m^2)",error )

            call h5ltmake_compressed_dataset_double_f(gid, "cold", 4, dims_stokes(2:5), &
                 coldstokes, error)
            call h5ltset_attribute_string_f(gid,"cold","description", &
                 "Cold D-alpha emission stokes parameters: cold(4,lambda,chan,species)", error)
            call h5ltset_attribute_string_f(gid,"cold","units","Ph/(s*nm*sr*m^2)",error )
        else
            call h5ltmake_compressed_dataset_double_f(fid, "/cold", 4, dims, &
                 spec%cold, error)
            call h5ltset_attribute_string_f(fid,"/cold","description", &
                 "Cold D-alpha emission stark components: cold(stark,lambda,chan,species)", &
                 error)
            call h5ltset_attribute_string_f(fid,"/cold","units","Ph/(s*nm*sr*m^2)",error )

            call h5ltmake_compressed_dataset_double_f(gid, "cold", 5, dims_stokes, &
                 spec%coldstokes, error)
            call h5ltset_attribute_string_f(gid,"cold","description", &
                 "Cold D-alpha emission stokes parameters: cold(stark,4,lambda,chan,species)", error)
            call h5ltset_attribute_string_f(gid,"cold","units","Ph/(s*nm*sr*m^2)",error )
        endif
    endif

    dims(4) = particles%nclass
    if(inputs%calc_fida.ge.1) then
        spec%fida = factor*spec%fida
        spec%fidastokes = factor*spec%fidastokes
        if(inputs%stark_components.eq.0) then
            fida = sum(spec%fida, dim=1)
            fidastokes = sum(spec%fidastokes, dim=1)
        endif
        !Write variables
        if(particles%nclass.le.1) then
            if (inputs%stark_components.eq.0) then
                call h5ltmake_compressed_dataset_double_f(fid, "/fida", 2, &
                     dims(2:3), fida(:,:,1), error)

                call h5ltmake_compressed_dataset_double_f(gid, "fida", 3, &
                     dims_stokes(2:4), fidastokes(:,:,:,1), error)

                !Add attributes
                call h5ltset_attribute_string_f(fid,"/fida","description", &
                     "Active Fast-ion D-alpha (FIDA) emmision: fida(lambda,chan)", error)

                call h5ltset_attribute_string_f(gid,"fida","description", &
                     "Active Fast-ion D-alpha (FIDA) emmision stokes parameters: &
                     &fida(4,lambda,chan)", error)
            else
                call h5ltmake_compressed_dataset_double_f(fid, "/fida", 3, &
                     dims(1:3), spec%fida(:,:,:,1), error)

                call h5ltmake_compressed_dataset_double_f(gid, "fida", 4, &
                     dims_stokes(1:4), spec%fidastokes(:,:,:,:,1), error)

                !Add attributes
                call h5ltset_attribute_string_f(fid,"/fida","description", &
                     "Active Fast-ion D-alpha (FIDA) emmision stark components: fida(stark,lambda,chan)", &
                     error)

                call h5ltset_attribute_string_f(gid,"fida","description", &
                     "Active Fast-ion D-alpha (FIDA) emmision stark components stokes parameters: &
                     &fida(stark,4,lambda,chan)", &
                     error)
            endif
        else
            if (inputs%stark_components.eq.0) then
                call h5ltmake_dataset_int_f(fid,"/nclass", 0, d, [particles%nclass], error)
                call h5ltmake_compressed_dataset_double_f(fid, "/fida", 3, &
                     dims(2:4), fida, error)

                call h5ltmake_compressed_dataset_double_f(gid, "fida", 4, &
                     dims_stokes(2:5), fidastokes, error)
                !Add attributes
                call h5ltset_attribute_string_f(fid,"/fida","description", &
                     "Active Fast-ion D-alpha (FIDA) emmision: fida(lambda,chan,class)", error)

                call h5ltset_attribute_string_f(gid,"fida","description", &
                     "Active Fast-ion D-alpha (FIDA) emmision stokes parameters: &
                     &fida(4,lambda,chan,class)", error)
            else
                call h5ltmake_compressed_dataset_double_f(fid, "/fida", 4, &
                     dims, spec%fida, error)

                call h5ltmake_compressed_dataset_double_f(gid, "fida", 5, &
                     dims_stokes, spec%fidastokes, error)

                !Add attributes
                call h5ltset_attribute_string_f(fid,"/fida","description", &
                     "Active Fast-ion D-alpha (FIDA) emmision stark components: fida(stark,lambda,chan,class)", &
                     error)

                call h5ltset_attribute_string_f(gid,"fida","description", &
                     "Active Fast-ion D-alpha (FIDA) emmision stark components stokes parameters: &
                     &fida(stark,4,lambda,chan,class)", &
                     error)
            endif
       endif
       call h5ltset_attribute_string_f(fid,"/fida","units","Ph/(s*nm*sr*m^2)",error )

       call h5ltset_attribute_string_f(gid,"fida","units","Ph/(s*nm*sr*m^2)",error )
    endif

    if(inputs%calc_pfida.ge.1) then
        spec%pfida = factor*spec%pfida
        spec%pfidastokes = factor*spec%pfidastokes
        if (inputs%stark_components.eq.0) then
            pfida = sum(spec%pfida, dim=1)
            pfidastokes = sum(spec%pfidastokes, dim=1)
        endif
        !Write variables
        if(particles%nclass.le.1) then
            if (inputs%stark_components.eq.0) then
                call h5ltmake_compressed_dataset_double_f(fid, "/pfida", 2, &
                     dims(2:3), pfida(:,:,1), error)

                call h5ltmake_compressed_dataset_double_f(gid, "pfida", 3, &
                     dims_stokes(2:4), pfidastokes(:,:,:,1), error)

                !Add attributes
                call h5ltset_attribute_string_f(fid,"/pfida","description", &
                     "Passive Fast-ion D-alpha (p-FIDA) emmision: pfida(lambda,chan)", error)

                call h5ltset_attribute_string_f(gid,"pfida","description", &
                     "Passive Fast-ion D-alpha (p-FIDA) emmision stokes parameters: &
                     &pfida(4,lambda,chan)", error)
            else
                call h5ltmake_compressed_dataset_double_f(fid, "/pfida", 3, &
                     dims(1:3), spec%pfida(:,:,:,1), error)

                call h5ltmake_compressed_dataset_double_f(gid, "pfida", 4, &
                     dims_stokes(1:4), spec%pfidastokes(:,:,:,:,1), error)

                !Add attributes
                call h5ltset_attribute_string_f(fid,"/pfida","description", &
                     "Passive Fast-ion D-alpha (p-FIDA) emmision stark components: pfida(stark,lambda,chan)", &
                     error)

                call h5ltset_attribute_string_f(gid,"pfida","description", &
                     "Passive Fast-ion D-alpha (p-FIDA) emmision stokes parameters: &
                     &pfida(stark,4,lambda,chan)", error)
            endif
        else
            if(inputs%calc_fida.le.0) then
                call h5ltmake_dataset_int_f(fid,"/nclass", 0, d, [particles%nclass], error)
            endif
            if (inputs%stark_components.eq.0) then
                call h5ltmake_compressed_dataset_double_f(fid, "/pfida", 3, &
                     dims(2:4), pfida, error)

                call h5ltmake_compressed_dataset_double_f(gid, "pfida", 4, &
                     dims_stokes(2:5), pfidastokes, error)

                !Add attributes
                call h5ltset_attribute_string_f(fid,"/pfida","description", &
                     "Passive Fast-ion D-alpha (p-FIDA) emmision: pfida(lambda,chan,class)", error)

                call h5ltset_attribute_string_f(gid,"pfida","description", &
                     "Passive Fast-ion D-alpha (p-FIDA) emmision stokes parameters: pfida(4,lambda,chan, class)", error)
            else
                call h5ltmake_compressed_dataset_double_f(fid, "/pfida", 4, &
                     dims, spec%pfida, error)

                call h5ltmake_compressed_dataset_double_f(gid, "pfida", 5, &
                     dims_stokes, spec%pfidastokes, error)

                !Add attributes
                call h5ltset_attribute_string_f(fid,"/pfida","description", &
                     "Passive Fast-ion D-alpha (p-FIDA) emmision stark components: pfida(stark,lambda,chan,class)", &
                     error)

                call h5ltset_attribute_string_f(gid,"pfida","description", &
                     "Passive Fast-ion D-alpha (p-FIDA) emmision stokes parameters: pfida(stark,4,lambda,chan, class)", error)
            endif
       endif
       call h5ltset_attribute_string_f(fid,"/pfida","units","Ph/(s*nm*sr*m^2)",error )

       call h5ltset_attribute_string_f(gid,"pfida","units","Ph/(s*nm*sr*m^2)",error )
    endif
    call h5gclose_f(gid, error)

    if(inputs%calc_res.ge.1) then 
       !Create spatial group
       call h5gcreate_f(fid, "spatial", sid, error)
       call h5ltset_attribute_string_f(fid,"spatial","description", &
            "Spatial Resolution",error)
       dims(1) = 3
       dims(2) = reservoir_size
       dims(3) = spec_chords%nchan
       if(inputs%calc_bes.ge.1) then
           !! Full
           call h5gcreate_f(sid,"full", gid, error)
           call h5ltset_attribute_string_f(sid,"full","description","Full Energy Spatial Resolution",error)
           photons = 0.0
           spat = 0.0
           distance = 0.0
           do ic=1,spec_chords%nchan
               axis = spec_chords%los(ic)%axis
               axis = axis/norm2(axis)
               src = spec_chords%los(ic)%lens
               do ir=1,reservoir_size
                   photons(ir,ic) = spatres%full(ic)%R(ir)%w
                   pos = spatres%full(ic)%R(ir)%v
                   spat(:,ir,ic) = pos
                   distance(ir,ic) = dot_product((pos - src),axis)
               enddo
           enddo
           call h5ltmake_compressed_dataset_double_f(gid, "ri", 3, dims(1:3),&
                spat, error)
           call h5ltset_attribute_string_f(gid,"ri","units","cm",error)
           call h5ltset_attribute_string_f(gid,"ri","description",&
                "Birth position of neutral that produced a photon: ri([x,y,z],sample,channel)",error)
           call h5ltmake_compressed_dataset_double_f(gid, "photons", 2, dims(2:3),&
                photons, error)
           call h5ltset_attribute_string_f(gid,"photons","units","Ph",error)
           call h5ltset_attribute_string_f(gid,"photons","description",&
                "Number of photons produced by neutral: photons(sample,channel)",error)
           call h5ltmake_compressed_dataset_double_f(gid, "distance", 2, dims(2:3),&
                distance, error)
           call h5ltset_attribute_string_f(gid,"distance","units","cm",error)
           call h5ltset_attribute_string_f(gid,"distance","description",&
                "Distance from lens such that lens + distance*axis = point on LOS closest to ri: distance(sample,channel)",error)
           call h5gclose_f(gid, error)

           !! Half 
           call h5gcreate_f(sid,"half", gid, error)
           call h5ltset_attribute_string_f(sid,"half","description","Half Energy Spatial Resolution",error)
           photons = 0.0
           spat = 0.0
           distance = 0.0
           do ic=1,spec_chords%nchan
               axis = spec_chords%los(ic)%axis
               axis = axis/norm2(axis)
               src = spec_chords%los(ic)%lens
               do ir=1,reservoir_size
                   photons(ir,ic) = spatres%half(ic)%R(ir)%w
                   pos = spatres%half(ic)%R(ir)%v
                   spat(:,ir,ic) = pos
                   distance(ir,ic) = dot_product((pos - src),axis)
               enddo
           enddo
           call h5ltmake_compressed_dataset_double_f(gid, "ri", 3, dims(1:3),&
                spat, error)
           call h5ltset_attribute_string_f(gid,"ri","units","cm",error)
           call h5ltset_attribute_string_f(gid,"ri","description",&
                "Birth position of neutral that produced a photon: ri([x,y,z],sample,channel)",error)
           call h5ltmake_compressed_dataset_double_f(gid, "photons", 2, dims(2:3),&
                photons, error)
           call h5ltset_attribute_string_f(gid,"photons","units","Ph",error)
           call h5ltset_attribute_string_f(gid,"photons","description",&
                "Number of photons produced by neutral: photons(sample,channel)",error)
           call h5ltmake_compressed_dataset_double_f(gid, "distance", 2, dims(2:3),&
                distance, error)
           call h5ltset_attribute_string_f(gid,"distance","units","cm",error)
           call h5ltset_attribute_string_f(gid,"distance","description",&
                "Distance from lens such that lens + distance*axis = point on LOS closest to ri: distance(sample,channel)",error)
           call h5gclose_f(gid, error)

           !! Third
           call h5gcreate_f(sid,"third", gid, error)
           call h5ltset_attribute_string_f(sid,"third","description","Third Energy Spatial Resolution",error)
           photons = 0.0
           spat = 0.0
           distance = 0.0
           do ic=1,spec_chords%nchan
               axis = spec_chords%los(ic)%axis
               axis = axis/norm2(axis)
               src = spec_chords%los(ic)%lens
               do ir=1,reservoir_size
                   photons(ir,ic) = spatres%third(ic)%R(ir)%w
                   pos = spatres%third(ic)%R(ir)%v
                   spat(:,ir,ic) = pos
                   distance(ir,ic) = dot_product((pos - src),axis)
               enddo
           enddo
           call h5ltmake_compressed_dataset_double_f(gid, "ri", 3, dims(1:3),&
                spat, error)
           call h5ltset_attribute_string_f(gid,"ri","units","cm",error)
           call h5ltset_attribute_string_f(gid,"ri","description",&
                "Birth position of neutral that produced a photon: ri([x,y,z],sample,channel)",error)
           call h5ltmake_compressed_dataset_double_f(gid, "photons", 2, dims(2:3),&
                photons, error)
           call h5ltset_attribute_string_f(gid,"photons","units","Ph",error)
           call h5ltset_attribute_string_f(gid,"photons","description",&
                "Number of photons produced by neutral: photons(sample,channel)",error)
           call h5ltmake_compressed_dataset_double_f(gid, "distance", 2, dims(2:3),&
                distance, error)
           call h5ltset_attribute_string_f(gid,"distance","units","cm",error)
           call h5ltset_attribute_string_f(gid,"distance","description",&
                "Distance from lens such that lens + distance*axis = point on LOS closest to ri: distance(sample,channel)",error)
           call h5gclose_f(gid, error)
       endif

       if(inputs%calc_dcx.ge.1) then
           call h5gcreate_f(sid,"dcx", gid, error)
           call h5ltset_attribute_string_f(sid,"dcx","description","DCX Spatial Resolution",error)
           photons = 0.0
           spat = 0.0
           distance = 0.0
           do ic=1,spec_chords%nchan
               axis = spec_chords%los(ic)%axis
               axis = axis/norm2(axis)
               src = spec_chords%los(ic)%lens
               do ir=1,reservoir_size
                   photons(ir,ic) = spatres%dcx(ic)%R(ir)%w
                   pos = spatres%dcx(ic)%R(ir)%v
                   spat(:,ir,ic) = pos
                   distance(ir,ic) = dot_product((pos - src),axis)
               enddo
           enddo
           call h5ltmake_compressed_dataset_double_f(gid, "ri", 3, dims(1:3),&
                spat, error)
           call h5ltset_attribute_string_f(gid,"ri","units","cm",error)
           call h5ltset_attribute_string_f(gid,"ri","description",&
                "Birth position of neutral that produced a photon: ri([x,y,z],sample,channel)",error)
           call h5ltmake_compressed_dataset_double_f(gid, "photons", 2, dims(2:3),&
                photons, error)
           call h5ltset_attribute_string_f(gid,"photons","units","Ph",error)
           call h5ltset_attribute_string_f(gid,"photons","description",&
                "Number of photons produced by neutral: photons(sample,channel)",error)
           call h5ltmake_compressed_dataset_double_f(gid, "distance", 2, dims(2:3),&
                distance, error)
           call h5ltset_attribute_string_f(gid,"distance","units","cm",error)
           call h5ltset_attribute_string_f(gid,"distance","description",&
                "Distance from lens such that lens + distance*axis = point on LOS closest to ri: distance(sample,channel)",error)
           call h5gclose_f(gid, error)

       endif

       if(inputs%calc_halo.ge.1) then
           call h5gcreate_f(sid,"halo", gid, error)
           call h5ltset_attribute_string_f(sid,"halo","description","Halo Spatial Resolution",error)
           photons = 0.0
           spat = 0.0
           distance = 0.0
           do ic=1,spec_chords%nchan
               axis = spec_chords%los(ic)%axis
               axis = axis/norm2(axis)
               src = spec_chords%los(ic)%lens
               do ir=1,reservoir_size
                   photons(ir,ic) = spatres%halo(ic)%R(ir)%w
                   pos = spatres%halo(ic)%R(ir)%v
                   spat(:,ir,ic) = pos
                   distance(ir,ic) = dot_product((pos - src),axis)
               enddo
           enddo
           call h5ltmake_compressed_dataset_double_f(gid, "ri", 3, dims(1:3),&
                spat, error)
           call h5ltset_attribute_string_f(gid,"ri","units","cm",error)
           call h5ltset_attribute_string_f(gid,"ri","description",&
                "Birth position of neutral that produced a photon: ri([x,y,z],sample,channel)",error)
           call h5ltmake_compressed_dataset_double_f(gid, "photons", 2, dims(2:3),&
                photons, error)
           call h5ltset_attribute_string_f(gid,"photons","units","Ph",error)
           call h5ltset_attribute_string_f(gid,"photons","description",&
                "Number of photons produced by neutral: photons(sample,channel)",error)
           call h5ltmake_compressed_dataset_double_f(gid, "distance", 2, dims(2:3),&
                distance, error)
           call h5ltset_attribute_string_f(gid,"distance","units","cm",error)
           call h5ltset_attribute_string_f(gid,"distance","description",&
                "Distance from lens such that lens + distance*axis = point on LOS closest to ri: distance(sample,channel)",error)
           call h5gclose_f(gid, error)
       endif

       if(inputs%calc_fida.ge.1) then
           call h5gcreate_f(sid,"fida", gid, error)
           call h5ltset_attribute_string_f(sid,"fida","description","FIDA Spatial Resolution",error)
           photons = 0.0
           spat = 0.0
           distance = 0.0
           do ic=1,spec_chords%nchan
               axis = spec_chords%los(ic)%axis
               axis = axis/norm2(axis)
               src = spec_chords%los(ic)%lens
               do ir=1,reservoir_size
                   photons(ir,ic) = spatres%fida(ic)%R(ir)%w
                   pos = spatres%fida(ic)%R(ir)%v
                   spat(:,ir,ic) = pos
                   distance(ir,ic) = dot_product((pos - src),axis)
               enddo
           enddo
           call h5ltmake_compressed_dataset_double_f(gid, "ri", 3, dims(1:3),&
                spat, error)
           call h5ltset_attribute_string_f(gid,"ri","units","cm",error)
           call h5ltset_attribute_string_f(gid,"ri","description",&
                "Birth position of neutral that produced a photon: ri([x,y,z],sample,channel)",error)
           call h5ltmake_compressed_dataset_double_f(gid, "photons", 2, dims(2:3),&
                photons, error)
           call h5ltset_attribute_string_f(gid,"photons","units","Ph",error)
           call h5ltset_attribute_string_f(gid,"photons","description",&
                "Number of photons produced by neutral: photons(sample,channel)",error)
           call h5ltmake_compressed_dataset_double_f(gid, "distance", 2, dims(2:3),&
                distance, error)
           call h5ltset_attribute_string_f(gid,"distance","units","cm",error)
           call h5ltset_attribute_string_f(gid,"distance","description",&
                "Distance from lens such that lens + distance*axis = point on LOS closest to ri: distance(sample,channel)",error)
           call h5gclose_f(gid, error)
       endif

       if(inputs%calc_pfida.ge.1) then
           call h5gcreate_f(sid,"pfida", gid, error)
           call h5ltset_attribute_string_f(sid,"pfida","description","Passive FIDA Spatial Resolution",error)
           photons = 0.0
           spat = 0.0
           distance = 0.0
           do ic=1,spec_chords%nchan
               axis = spec_chords%los(ic)%axis
               axis = axis/norm2(axis)
               src = spec_chords%los(ic)%lens
               do ir=1,reservoir_size
                   photons(ir,ic) = spatres%fida(ic)%R(ir)%w
                   pos = spatres%fida(ic)%R(ir)%v
                   spat(:,ir,ic) = pos
                   distance(ir,ic) = dot_product((pos - src),axis)
               enddo
           enddo
           call h5ltmake_compressed_dataset_double_f(gid, "ri", 3, dims(1:3),&
                spat, error)
           call h5ltset_attribute_string_f(gid,"ri","units","cm",error)
           call h5ltset_attribute_string_f(gid,"ri","description",&
                "Birth position of neutral that produced a photon: ri([x,y,z],sample,channel)",error)
           call h5ltmake_compressed_dataset_double_f(gid, "photons", 2, dims(2:3),&
                photons, error)
           call h5ltset_attribute_string_f(gid,"photons","units","Ph",error)
           call h5ltset_attribute_string_f(gid,"photons","description",&
                "Number of photons produced by neutral: photons(sample,channel)",error)
           call h5ltmake_compressed_dataset_double_f(gid, "distance", 2, dims(2:3),&
                distance, error)
           call h5ltset_attribute_string_f(gid,"distance","units","cm",error)
           call h5ltset_attribute_string_f(gid,"distance","description",&
                "Distance from lens such that lens + distance*axis = point on LOS closest to ri: distance(sample,channel)",error)
           call h5gclose_f(gid, error)
       endif

       !Close spatial group
       call h5gclose_f(sid, error)
    endif

    call h5ltset_attribute_string_f(fid, "/", "version", version, error)
    call h5ltset_attribute_string_f(fid,"/","description",&
         "Spectra calculated by FIDASIM", error)

    !Close file
    call h5fclose_f(fid, error)

    !Close HDF5 interface
    call h5close_f(error)

    if(inputs%verbose.ge.1) then
        write(*,'(T4,a,a)') 'Spectra written to: ', trim(filename)
    endif

end subroutine write_spectra

subroutine write_neutrons
    !+ Writes [[libfida:neutron]] to a HDF5 file
    integer(HID_T) :: fid
    integer(HSIZE_T), dimension(1) :: dim1
    integer(HSIZE_T), dimension(3) :: dim3
    integer(HSIZE_T), dimension(5) :: dim5
    integer :: error

    character(charlim) :: filename

    !! write to file
    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_neutrons.h5"

    !Open HDF5 interface
    call h5open_f(error)

    !Create file overwriting any existing file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, fid, error)

    !Write variables
    if(particles%nclass.gt.1) then
        dim1(1) = 1
        call h5ltmake_dataset_int_f(fid,"/nclass", 0, dim1, [particles%nclass], error)
        dim1(1) = particles%nclass
        call h5ltmake_compressed_dataset_double_f(fid, "/rate", 1, dim1, neutron%rate, error)
        call h5ltset_attribute_string_f(fid,"/rate","description", &
             "Neutron rate: rate(orbit_class)", error)
    else
        dim1(1) = 1
        call h5ltmake_dataset_double_f(fid, "/rate", 0, dim1, neutron%rate, error)
        call h5ltset_attribute_string_f(fid,"/rate","description", &
             "Neutron rate", error)
    endif
    call h5ltset_attribute_string_f(fid,"/rate","units","neutrons/s",error )

    if((inputs%dist_type.eq.1).and.(inputs%calc_neutron.ge.2)) then
        dim1(1) = 1
        call h5ltmake_dataset_int_f(fid,"/nenergy",0,dim1,[fbm%nenergy], error)
        call h5ltmake_dataset_int_f(fid,"/npitch",0,dim1,[fbm%npitch], error)
        call h5ltmake_dataset_int_f(fid,"/nr",0,dim1,[fbm%nr], error)
        call h5ltmake_dataset_int_f(fid,"/nz",0,dim1,[fbm%nz], error)
        dim5 = shape(neutron%weight)
        call h5ltmake_compressed_dataset_double_f(fid, "/weight", 5, dim5, neutron%weight, error)
        dim3 = shape(neutron%emis)
        call h5ltmake_compressed_dataset_double_f(fid, "/emissivity", 3, dim3, neutron%emis, error)
        call h5ltmake_compressed_dataset_double_f(fid,"/energy", 1, dim5(1:1), fbm%energy, error)
        call h5ltmake_compressed_dataset_double_f(fid,"/pitch", 1, dim5(2:2), fbm%pitch, error)
        call h5ltmake_compressed_dataset_double_f(fid,"/r", 1, dim5(3:3), fbm%r, error)
        call h5ltmake_compressed_dataset_double_f(fid,"/z", 1, dim5(4:4), fbm%z, error)
        call h5ltmake_compressed_dataset_double_f(fid,"/phi", 1, dim5(5:5), fbm%phi, error)

        call h5ltset_attribute_string_f(fid,"/nenergy", "description", &
             "Number of energy values", error)
        call h5ltset_attribute_string_f(fid,"/npitch", "description", &
             "Number of pitch values", error)
        call h5ltset_attribute_string_f(fid,"/nr", "description", &
             "Number of R values", error)
        call h5ltset_attribute_string_f(fid,"/nz", "description", &
             "Number of Z values", error)
        call h5ltset_attribute_string_f(fid,"/weight", "description", &
             "Neutron Weight Function: weight(E,p,R,Z,Phi), rate = sum(f*weight)", error)
        call h5ltset_attribute_string_f(fid,"/weight", "units","neutrons*cm^3*dE*dp/fast-ion*s", error)
        call h5ltset_attribute_string_f(fid,"/emissivity", "description", &
             "Neutron Emissivity: emissivity(R,Z,Phi), rate = sum(emissivity)", error)
        call h5ltset_attribute_string_f(fid,"/emissivity", "units","neutrons*cm^3/fast-ion*s", error)

        call h5ltset_attribute_string_f(fid,"/energy","description", &
             "Energy array", error)
        call h5ltset_attribute_string_f(fid,"/energy", "units","keV", error)
        call h5ltset_attribute_string_f(fid,"/pitch", "description", &
             "Pitch array: p = v_parallel/v  w.r.t. the magnetic field", error)
        call h5ltset_attribute_string_f(fid,"/r","description", &
             "Radius array", error)
        call h5ltset_attribute_string_f(fid,"/r", "units","cm", error)
        call h5ltset_attribute_string_f(fid,"/z","description", &
             "Z array", error)
        call h5ltset_attribute_string_f(fid,"/z", "units","cm", error)
    endif

    call h5ltset_attribute_string_f(fid, "/", "version", version, error)
    call h5ltset_attribute_string_f(fid,"/","description",&
         "Neutron rate calculated by FIDASIM", error)

    !Close file
    call h5fclose_f(fid, error)

    !Close HDF5 interface
    call h5close_f(error)

    if(inputs%verbose.ge.1) then
        write(*,'(T4,a,a)') 'Neutrons written to: ', trim(filename)
    endif

end subroutine write_neutrons

subroutine write_cfpd_weights
    !+ Writes [[libfida:cfpd]] to a HDF5 file
    integer(HID_T) :: fid
    integer(HSIZE_T), dimension(1) :: dim1
    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(4) :: dim4
    integer :: error

    character(charlim) :: filename

    !! write to file
    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_cfpd.h5"

    !Open HDF5 interface
    call h5open_f(error)

    !Create file overwriting any existing file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, fid, error)

    !Write variables
    if(inputs%dist_type.eq.1) then
        dim1(1) = 1
        dim2 = [ctable%nenergy, ctable%nchan]
        dim4 = [ctable%nenergy, ctable%nchan, fbm%nenergy, fbm%npitch]

        call h5ltmake_compressed_dataset_double_f(fid, "/flux", 2, dim2, cfpd%flux, error)
        call h5ltmake_compressed_dataset_double_f(fid, "/prob", 2, dim2, cfpd%prob, error)
        call h5ltmake_compressed_dataset_double_f(fid, "/gam", 2, dim2, cfpd%gam, error)
        call h5ltmake_compressed_dataset_double_f(fid, "/weight", 4, dim4, cfpd%weight, error)

        call h5ltmake_dataset_int_f(fid,"/nenergy",0,dim1,[fbm%nenergy], error)
        call h5ltmake_dataset_int_f(fid,"/npitch",0,dim1,[fbm%npitch], error)
        call h5ltmake_compressed_dataset_double_f(fid,"/energy", 1, dim4(3:3), fbm%energy, error)
        call h5ltmake_compressed_dataset_double_f(fid,"/pitch", 1, dim4(4:4), fbm%pitch, error)
        call h5ltmake_compressed_dataset_double_f(fid,"/earray", 1, dim2(1:1), ctable%earray, error)

        call h5ltset_attribute_string_f(fid,"/flux", "description", &
             "CFPD flux: flux(energy,chan)", error)
        call h5ltset_attribute_string_f(fid,"/flux", "units", &
             "kHz", error)

        call h5ltset_attribute_string_f(fid,"/prob", "description", &
             "CFPD average nonzero probability: prob(energy,chan)", error)
        call h5ltset_attribute_string_f(fid,"/gam", "description", &
             "CFPD average nonzero gyroangle: gam(energy,chan)", error)
        call h5ltset_attribute_string_f(fid,"/gam", "units", "rad", error)
        call h5ltset_attribute_string_f(fid,"/weight", "description", &
             "CFPD Weight Function: weight(Ch,E3,E,p), rate = sum(f*weight)", error)
        call h5ltset_attribute_string_f(fid,"/weight", "units","products*cm^3*dE*dp/fast-ion*s", error)

        call h5ltset_attribute_string_f(fid,"/nenergy", "description", &
             "Number of distribution function energy values", error)
        call h5ltset_attribute_string_f(fid,"/npitch", "description", &
             "Number of distribution function pitch values", error)

        call h5ltset_attribute_string_f(fid,"/energy","description", &
             "Energy array", error)
        call h5ltset_attribute_string_f(fid,"/energy", "units","keV", error)
        call h5ltset_attribute_string_f(fid,"/pitch", "description", &
             "Pitch array: p = v_parallel/v  w.r.t. the magnetic field", error)
        call h5ltset_attribute_string_f(fid,"/earray","description", &
             "E3 energy array", error)
        call h5ltset_attribute_string_f(fid,"/earray", "units","keV", error)
    endif

    call h5ltset_attribute_string_f(fid, "/", "version", version, error)
    call h5ltset_attribute_string_f(fid,"/","description",&
         "CFPD signals calculated by FIDASIM", error)

    !Close file
    call h5fclose_f(fid, error)

    !Close HDF5 interface
    call h5close_f(error)

    if(inputs%verbose.ge.1) then
        write(*,'(T4,a,a)') 'Charged fusion products written to: ', trim(filename)
    endif

end subroutine write_cfpd_weights

subroutine write_fida_weights
    !+ Writes [[libfida:fweight]] to a HDF5 file
    !! HDF5 variables
    integer(HID_T) :: fid
    integer(HSIZE_T), dimension(4) :: dim4
    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(1) :: dim1
    integer :: error

    character(charlim) :: filename
    integer :: i,ie,ip,ic,iwav
    real(Float64), dimension(:),   allocatable :: lambda_arr
    real(Float64), dimension(:),   allocatable :: ebarr,ptcharr
    real(Float64), dimension(:,:), allocatable :: jacobian,e_grid,p_grid
    real(Float64), dimension(:,:), allocatable :: vpa_grid,vpe_grid,fida
    real(Float64) :: dlambda, wtot, dE, dP

    dlambda=(inputs%lambdamax_wght-inputs%lambdamin_wght)/inputs%nlambda_wght
    allocate(lambda_arr(inputs%nlambda_wght))
    do i=1,inputs%nlambda_wght
        lambda_arr(i)=(i-0.5)*dlambda + inputs%lambdamin_wght
    enddo

    !! define arrays
    !! define energy - array
    allocate(ebarr(inputs%ne_wght))
    do i=1,inputs%ne_wght
        ebarr(i)=real(i-0.5)*inputs%emax_wght/real(inputs%ne_wght)
    enddo
    dE = abs(ebarr(2)-ebarr(1))

    !! define pitch - array
    allocate(ptcharr(inputs%np_wght))
    do i=1,inputs%np_wght
        ptcharr(i)=real(i-0.5)*2./real(inputs%np_wght)-1.
    enddo
    dP = abs(ptcharr(2)-ptcharr(1))

    !! define 2d grids
    !! define energy grid
    allocate(e_grid(inputs%ne_wght,inputs%np_wght))
    do i=1,inputs%ne_wght
        e_grid(i,:) = ebarr(i)
    enddo

    !! define pitch grid
    allocate(p_grid(inputs%ne_wght,inputs%np_wght))
    do i=1,inputs%np_wght
        p_grid(:,i) = ptcharr(i)
    enddo

    !! define velocity space grid
    allocate(vpe_grid(inputs%ne_wght,inputs%np_wght)) !! V perpendicular
    allocate(vpa_grid(inputs%ne_wght,inputs%np_wght)) !! V parallel
    vpa_grid = 100*sqrt((((2.0d3)*e0)/(mass_u*beam_mass))*e_grid)*p_grid ! [cm/s]
    vpe_grid = 100*sqrt((((2.0d3)*e0)/(mass_u*beam_mass))*e_grid*(1.0-p_grid**2)) ![cm/s]

    !! define jacobian to convert between E-p to velocity
    allocate(jacobian(inputs%ne_wght,inputs%np_wght))
    jacobian = ((beam_mass*mass_u)/(e0*1.0d3)) *vpe_grid/sqrt(vpa_grid**2 + vpe_grid**2)

    !! normalize mean_f
    do ic=1,spec_chords%nchan
        do ip=1,inputs%np_wght
            do ie=1,inputs%ne_wght
                wtot = sum(fweight%weight(:,ie,ip,ic))
                if((wtot.gt.0.d0)) then
                    fweight%mean_f(ie,ip,ic) = fweight%mean_f(ie,ip,ic)/wtot
                endif
            enddo
        enddo
    enddo

    !! Calculate FIDA estimate
    allocate(fida(inputs%nlambda_wght,spec_chords%nchan))
    do iwav=1,size(fida,1)
        fida(iwav,:) = (dE*dP*1d4)*sum(sum(fweight%mean_f(:,:,:)*fweight%weight(iwav,:,:,:),1),1)
    enddo

    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_fida_weights.h5"

    !Open HDF5 interface
    call h5open_f(error)

    !Create file overwriting any existing file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, fid, error)

    dim1(1) = 1
    dim2 = [inputs%nlambda_wght, spec_chords%nchan]
    dim4 = [inputs%nlambda_wght, inputs%ne_wght, inputs%np_wght, spec_chords%nchan]
    call h5ltmake_dataset_int_f(fid,"/nenergy",0,dim1,[inputs%ne_wght], error)
    call h5ltmake_dataset_int_f(fid,"/npitch",0,dim1,[inputs%np_wght], error)
    call h5ltmake_dataset_int_f(fid,"/nchan",0,dim1,[spec_chords%nchan], error)
    call h5ltmake_compressed_dataset_double_f(fid,"/weight",4,dim4,fweight%weight,error)
    call h5ltmake_compressed_dataset_double_f(fid,"/fida",2,dim2,fida,error)
    call h5ltmake_compressed_dataset_double_f(fid,"/mean_f",3,dim4(2:4),fweight%mean_f,error)
    call h5ltmake_compressed_dataset_double_f(fid,"/lambda",1,dim4(1:1),lambda_arr,error)
    call h5ltmake_compressed_dataset_double_f(fid,"/energy",1,dim4(2:2),ebarr, error)
    call h5ltmake_compressed_dataset_double_f(fid,"/pitch",1,dim4(3:3),ptcharr, error)
    call h5ltmake_compressed_dataset_double_f(fid,"/radius",1,dim4(4:4),spec_chords%radius, error)
    dim2 = [inputs%ne_wght, inputs%np_wght]
    call h5ltmake_compressed_dataset_double_f(fid,"/jacobian",2,dim2, jacobian, error)
    call h5ltmake_compressed_dataset_double_f(fid,"/vpe_grid",2,dim2,vpe_grid, error)
    call h5ltmake_compressed_dataset_double_f(fid,"/vpa_grid",2,dim2,vpa_grid, error)
    call h5ltmake_compressed_dataset_double_f(fid,"/e_grid",2,dim2,e_grid, error)
    call h5ltmake_compressed_dataset_double_f(fid,"/p_grid",2,dim2,p_grid, error)

    !Add attributes
    call h5ltset_attribute_string_f(fid, "/", "version", version, error)
    if(inputs%calc_fida_wght.eq.1) then
        call h5ltset_attribute_string_f(fid,"/", "description", &
             "Line of Sight averaged FIDA E-p space sensitivity/weights " // &
             "and spectra calculated by FIDASIM", error)
    else
        call h5ltset_attribute_string_f(fid,"/", "description", &
             "Full FIDA E-p space sensitivity/weights and spectra calculated " // &
             "by FIDASIM via Monte Carlo method", error)
    endif
    call h5ltset_attribute_string_f(fid,"/weight","description", &
         "E-p space sensivity/weight of FIDA diagnostic: weight(lambda,energy,pitch,chan)", error)
    call h5ltset_attribute_string_f(fid,"/weight","units", &
         "(Ph*cm)/(s*nm*sr*fast-ion*dE*dP)",error)
    call h5ltset_attribute_string_f(fid,"/fida","units", &
         "Ph/(s*nm*sr*m^2)",error )
    call h5ltset_attribute_string_f(fid,"/fida","description", &
         "Estimate of Fast-ion D-alpha (FIDA) emmision calculated by 1e4*weight*mean_f*dEdP: fida(lambda,chan)", error)
    call h5ltset_attribute_string_f(fid,"/mean_f","description", &
         "Estimated mean fast-ion distribution function seen by los: mean_f(energy,pitch,chan)", error)
    call h5ltset_attribute_string_f(fid,"/mean_f","units", &
         "fast-ion/(dE*dP*cm^3)", error)
    call h5ltset_attribute_string_f(fid,"/lambda","description", &
         "Wavelength array", error)
    call h5ltset_attribute_string_f(fid,"/lambda","units","nm", error)
    call h5ltset_attribute_string_f(fid,"/nchan", "description", &
         "Number of channels", error)
    call h5ltset_attribute_string_f(fid,"/nenergy", "description", &
         "Number of energy values", error)
    call h5ltset_attribute_string_f(fid,"/npitch", "description", &
         "Number of pitch value", error)
    call h5ltset_attribute_string_f(fid,"/energy","description", &
         "Energy array", error)
    call h5ltset_attribute_string_f(fid,"/energy", "units","keV", error)
    call h5ltset_attribute_string_f(fid,"/pitch", "description", &
         "Pitch array: p = v_parallel/v  w.r.t. the magnetic field", error)
    call h5ltset_attribute_string_f(fid,"/radius", "description", &
         "Line of sight radius at midplane or tangency point", error)
    call h5ltset_attribute_string_f(fid,"/radius", "units","cm", error)
    call h5ltset_attribute_string_f(fid,"/jacobian","description", &
         "Jacobian used to convert from E-p space to velocity space", error)
    call h5ltset_attribute_string_f(fid,"/jacobian","units", &
         "(dE*dP)/(dvpa*dvpe)", error)
    call h5ltset_attribute_string_f(fid,"/e_grid","description", &
         "2D energy grid", error)
    call h5ltset_attribute_string_f(fid,"/e_grid","units","keV", error)
    call h5ltset_attribute_string_f(fid,"/p_grid","description", &
         "2D pitch grid", error)
    call h5ltset_attribute_string_f(fid,"/vpe_grid","description", &
         "2D perpendicular velocity grid", error)
    call h5ltset_attribute_string_f(fid,"/vpe_grid","units","cm/s", error)
    call h5ltset_attribute_string_f(fid,"/vpa_grid","description", &
         "2D parallel velocity grid", error)
    call h5ltset_attribute_string_f(fid,"/vpa_grid","units","cm/s", error)

    !Close file
    call h5fclose_f(fid, error)

    !Close HDF5 interface
    call h5close_f(error)

    if(inputs%verbose.ge.1) then
        write(*,'(T4,a,a)') 'FIDA weights written to: ', trim(filename)
    endif

end subroutine write_fida_weights

subroutine write_npa_weights
    !+ Writes [[libfida:nweight]] to a HDF5 file
    character(charlim) :: filename
    integer :: i
    real(Float64), dimension(:), allocatable :: ebarr,ptcharr

    !! HDF5 variables
    integer(HID_T) :: fid
    integer(HSIZE_T), dimension(5) :: dim5
    integer(HSIZE_T), dimension(3) :: dim3
    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(1) :: d
    integer :: error

    !! define energy - array
    allocate(ebarr(inputs%ne_wght))
    do i=1,inputs%ne_wght
        ebarr(i)=real(i-0.5)*inputs%emax_wght/real(inputs%ne_wght)
    enddo

    !! define pitch - array
    allocate(ptcharr(inputs%np_wght))
    do i=1,inputs%np_wght
        ptcharr(i)=real(i-0.5)*2./real(inputs%np_wght)-1.
    enddo

    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_npa_weights.h5"

    !Open HDF5 interface
    call h5open_f(error)

    !Create file overwriting any existing file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, fid, error)

    !Write variables
    d(1) = 1
    dim2 = [inputs%ne_wght, npa_chords%nchan]
    dim3 = [inputs%ne_wght, inputs%np_wght, npa_chords%nchan]
    dim5 = [inputs%ne_wght, beam_grid%nx, beam_grid%ny, beam_grid%nz, npa_chords%nchan]
    call h5ltmake_dataset_int_f(fid, "/nchan", 0, d, [npa_chords%nchan], error)
    call h5ltmake_dataset_int_f(fid, "/nenergy", 0, d, [inputs%ne_wght], error)
    call h5ltmake_dataset_int_f(fid, "/npitch", 0, d, [inputs%np_wght], error)
    call h5ltmake_compressed_dataset_double_f(fid, "/radius", 1, &
         dim2(2:2), npa_chords%radius, error)
    call h5ltmake_compressed_dataset_double_f(fid, "/energy", 1, &
         dim2(1:1), ebarr, error)
    call h5ltmake_compressed_dataset_double_f(fid, "/pitch", 1, &
         dim3(2:2), ptcharr, error)
    call h5ltmake_compressed_dataset_double_f(fid, "/flux", 2, &
         dim2, nweight%flux, error)
    call h5ltmake_compressed_dataset_double_f(fid, "/weight", 3, &
         dim3, nweight%weight, error)

    !Add attributes
    call h5ltset_attribute_string_f(fid, "/", "version", version, error)
    call h5ltset_attribute_string_f(fid,"/", "description", &
         "NPA E-p space sensitivity/weights and Flux calculated by FIDASIM", error)
    call h5ltset_attribute_string_f(fid,"/nchan", "description", &
         "Number of channels", error)
    call h5ltset_attribute_string_f(fid,"/nenergy", "description", &
         "Number of energy values", error)
    call h5ltset_attribute_string_f(fid,"/npitch", "description", &
         "Number of pitch value", error)
    call h5ltset_attribute_string_f(fid,"/energy","description", &
         "Energy array", error)
    call h5ltset_attribute_string_f(fid,"/energy", "units","keV", error)
    call h5ltset_attribute_string_f(fid,"/pitch", "description", &
         "Pitch array: p = v_parallel/v  w.r.t. the magnetic field", error)
    call h5ltset_attribute_string_f(fid,"/radius", "description", &
         "Line of sight radius at midplane or tangency point", error)
    call h5ltset_attribute_string_f(fid,"/radius", "units","cm", error)
    call h5ltset_attribute_string_f(fid,"/flux", "description", &
         "Neutral flux: flux(energy,chan)", error)
    call h5ltset_attribute_string_f(fid,"/flux", "units", &
         "neutrals/(s*dE)", error)
    call h5ltset_attribute_string_f(fid,"/weight", "description", &
         "E-p space sensivity/weight of NPA diagnostics: weight(energy,pitch,chan)", error)
    call h5ltset_attribute_string_f(fid,"/weight","units", &
         "neutrals/(s*fast-ion*dE*dP)",error)

    if(inputs%calc_npa_wght.ge.2) then !Write diagnostic variables
        call write_beam_grid(fid, error)
        call h5ltmake_compressed_dataset_double_f(fid, "/emissivity", 4, &
             dim5(2:5), nweight%emissivity, error)
        call h5ltmake_compressed_dataset_double_f(fid, "/attenuation", 5, &
             dim5, nweight%attenuation, error)
        call h5ltmake_compressed_dataset_double_f(fid, "/cx", 5, &
             dim5, nweight%cx, error)
        call h5ltmake_compressed_dataset_double_f(fid, "/phit", 4, &
             dim5(2:5), npa_chords%phit%p, error)

        call h5ltset_attribute_string_f(fid,"/emissivity", "description", &
             "Neutral emissivity: emissivity(x,y,z,chan)", error)
        call h5ltset_attribute_string_f(fid,"/emissivity", "units", &
             "neutrals/(s*dV)", error)
        call h5ltset_attribute_string_f(fid,"/cx", "description", &
             "Charge-exchange rate: cx(energy,x,y,z,chan)", error)
        call h5ltset_attribute_string_f(fid,"/cx", "units", "s^(-1)", error)
        call h5ltset_attribute_string_f(fid,"/attenuation","description", &
             "Attenuation factor i.e. survival probability: attenuation(energy,x,y,z,chan)", error)
        call h5ltset_attribute_string_f(fid,"/phit","description", &
             "Probability of hitting the detector given an isotropic source: phit(x,y,z,chan)", error)
    endif

    !Close file
    call h5fclose_f(fid, error)

    !Close HDF5 interface
    call h5close_f(error)
    if(inputs%verbose.ge.1) then
        write(*,'(T4,a,a)') 'NPA weights written to: ',trim(filename)
    endif

end subroutine write_npa_weights

subroutine read_neutral_population(id, pop, error)
    !+ Reads neutral population from file
    integer(HID_T), intent(in) :: id
        !+ HDF5 id
    type(NeutralPopulation), intent(inout) :: pop
        !+ Neutral Population to populate
    integer, intent(out) :: error
        !+ Error code

    integer(HID_T) :: gid
    integer(HSIZE_T), dimension(4) :: dims4
    integer(HSIZE_T), dimension(5) :: dims5

    integer :: k, nx, ny, nz, ic, ind(3), ii, jj, kk, ir
    real(Float64) :: amu
    integer(Int32), dimension(:,:,:), allocatable :: n
    real(Float64), dimension(:,:,:,:), allocatable :: w
    real(Float64), dimension(:,:,:,:,:), allocatable :: v

    call init_neutral_population(pop)

    nx = beam_grid%nx; ny = beam_grid%ny; nz = beam_grid%nz

    dims4 = [nlevs, nx, ny, nz]
    call h5ltread_dataset_double_f(id, "dens", pop%dens, dims4, error)

    call h5gopen_f(id, "reservoir", gid, error)
    k = reservoir_size
    call h5ltread_dataset_int_scalar_f(gid, "k", k, error)
    if (k.ne.reservoir_size) then
        if(inputs%verbose.ge.0) then
            write(*,'(a,a)') 'READ_NEUTRAL_POPULATION: Unsupported reservoir size: ',k, reservoir_size
        endif
        stop
    endif

    allocate(v(3,k,nx,ny,nz))
    allocate(w(k,nx,ny,nz))
    allocate(n(nx,ny,nz))
    w = 0.d0; v=0.d0; n=0

    dims5 = [3, k, nx, ny, nz]
    call h5ltread_dataset_double_f(gid, "v", v, dims5, error)
    call h5ltread_dataset_double_f(gid, "w", w, dims5(2:5), error)
    call h5ltread_dataset_int_f(gid, "n", n, dims4(2:4), error)

    call h5gclose_f(gid, error)

    do ic=1,beam_grid%ngrid
        call ind2sub(beam_grid%dims,ic,ind)
        ii = ind(1) ; jj = ind(2) ; kk = ind(3)
        call init_reservoir(pop%res(ii,jj,kk))
        pop%res(ii,jj,kk)%n = n(ii,jj,kk)
        do ir=1,min(k,pop%res(ii,jj,kk)%n)
            pop%res(ii,jj,kk)%R(ir) = NeutralParticle(w(ir,ii,jj,kk), v(:,ir,ii,jj,kk))
        enddo
    enddo

    deallocate(v,w,n)

end subroutine read_neutral_population

subroutine read_neutrals
    !+ Reads neutral populations from file and puts it in [[libfida:neut]]
    integer(HID_T) :: fid, gid
    integer(HSIZE_T), dimension(4) :: dims
    integer :: error,nx,ny,nz
    logical :: exis,fatal_error

    if(inputs%verbose.ge.1) then
        write(*,'(a)') '---- loading neutrals ----'
    endif

    inquire(file=inputs%neutrals_file,exist=exis)
    if(exis) then
        if(inputs%verbose.ge.1) then
            write(*,'(T2,"Neutrals file: ",a)') trim(inputs%neutrals_file)
            write(*,*) ''
        endif
    else
        if(inputs%verbose.ge.0) then
            write(*,'(a,a)') 'READ_NEUTRALS: Neutrals file does not exist: ',inputs%neutrals_file
        endif
        stop
    endif

    !Open HDF5 interface
    call h5open_f(error)

    !Create file overwriting any existing file
    call h5fopen_f(inputs%neutrals_file, H5F_ACC_RDONLY_F, fid, error)

    call h5gopen_f(fid, "/grid", gid, error)
    call h5ltread_dataset_int_scalar_f(gid,"nx", nx, error)
    call h5ltread_dataset_int_scalar_f(gid,"ny", ny, error)
    call h5ltread_dataset_int_scalar_f(gid,"nz", nz, error)
    call h5gclose_f(gid, error)

    fatal_error = .False.
    if((nx.ne.beam_grid%nx).or. &
       (ny.ne.beam_grid%ny).or. &
       (nz.ne.beam_grid%nz)) then
        if(inputs%verbose.ge.0) then
            write(*,'(a)') 'READ_NEUTRALS: Neutrals file has incompatable grid dimensions'
        endif
        fatal_error = .True.
    endif

    !Check to make sure the neutrals file has all the needed neutrals
    call h5ltpath_valid_f(fid, "/full", .True., exis, error)
    if((.not.exis).and.(inputs%calc_nbi_dens.ge.1)) then
        if(inputs%verbose.ge.0) then
            write(*,'(a)') 'READ_NEUTRALS: Full energy neutral population is not in the neutrals file'
        endif
        fatal_error = .True.
    endif
    call h5ltpath_valid_f(fid, "/half", .True., exis, error)
    if((.not.exis).and.(inputs%calc_nbi_dens.ge.1)) then
        if(inputs%verbose.ge.0) then
            write(*,'(a)') 'READ_NEUTRALS: Half energy neutral population is not in the neutrals file'
        endif
        fatal_error = .True.
    endif
    call h5ltpath_valid_f(fid, "/third", .True., exis, error)
    if((.not.exis).and.(inputs%calc_nbi_dens.ge.1)) then
        if(inputs%verbose.ge.0) then
            write(*,'(a)') 'READ_NEUTRALS: Third energy neutral population is not in the neutrals file'
        endif
        fatal_error = .True.
    endif
    call h5ltpath_valid_f(fid, "/dcx", .True., exis, error)
    if((.not.exis).and.(inputs%calc_dcx_dens.ge.1)) then
        if(inputs%verbose.ge.0) then
            write(*,'(a)') 'READ_NEUTRALS: Direct Charge Exchange (DCX) neutral population is not in the neutrals file'
        endif
        fatal_error = .True.
    endif
    call h5ltpath_valid_f(fid, "/halo", .True., exis, error)
    if((.not.exis).and.(inputs%calc_halo_dens.ge.1)) then
        if(inputs%verbose.ge.0) then
            write(*,'(a)') 'READ_NEUTRALS: Thermal Halo neutral population is not in the neutrals file'
        endif
        fatal_error = .True.
    endif

    if(fatal_error) stop

    if(inputs%calc_nbi_dens.ge.1) then
        call h5gopen_f(fid, "/full", gid, error)
        call read_neutral_population(gid, neut%full, error)
        call h5gclose_f(gid, error)

        call h5gopen_f(fid, "/half", gid, error)
        call read_neutral_population(gid, neut%half, error)
        call h5gclose_f(gid, error)

        call h5gopen_f(fid, "/third", gid, error)
        call read_neutral_population(gid, neut%third, error)
        call h5gclose_f(gid, error)
    endif
    if(inputs%calc_dcx_dens.ge.1) then
        call h5gopen_f(fid, "/dcx", gid, error)
        call read_neutral_population(gid, neut%dcx, error)
        call h5gclose_f(gid, error)
    endif
    if(inputs%calc_halo_dens.ge.1) then
        call h5gopen_f(fid, "/halo", gid, error)
        call read_neutral_population(gid, neut%halo, error)
        call h5gclose_f(gid, error)
    endif

    !Close file
    call h5fclose_f(fid, error)

    !Close HDF5 interface
    call h5close_f(error)

end subroutine read_neutrals

!=============================================================================
!-----------------------------Geometry Routines-------------------------------
!=============================================================================
function approx_eq(x,y,tol) result(a)
    !+ Inexact equality comparison: `x ~= y` true if `abs(x-y) <= tol` else false
    real(Float64), intent(in) :: x
        !+First value in comparison
    real(Float64), intent(in) :: y
        !+Second value in comparison
    real(Float64), intent(in) :: tol
        !+Equality tolerance

    logical :: a

    a = abs(x-y).le.tol

end function approx_eq

function approx_ge(x,y,tol) result(a)
    !+ Inexact greater than or equal to comparison: `x >~= y`
    real(Float64), intent(in) :: x
        !+First value in comparison
    real(Float64), intent(in) :: y
        !+Second value in comparison
    real(Float64), intent(in) :: tol
        !+Equality tolerance

    logical :: a

    a = (x.gt.y).or.(approx_eq(x,y,tol))

end function approx_ge

function approx_le(x,y,tol) result(a)
    !+ Inexact less then or equal to comparison: `x <~= y`
    real(Float64), intent(in) :: x
        !+First value in comparison
    real(Float64), intent(in) :: y
        !+Second value in comparison
    real(Float64), intent(in) :: tol
        !+Equality tolerance

    logical :: a

    a = (x.lt.y).or.(approx_eq(x,y,tol))

end function approx_le

function cross_product(u, v) result(s)
    !+ Calculates the cross product of two vectors: `u`x`v`
    real(Float64), dimension(3), intent(in) :: u
    real(Float64), dimension(3), intent(in) :: v
    real(Float64), dimension(3)             :: s

    s(1) = u(2)*v(3) - u(3)*v(2)
    s(2) = u(3)*v(1) - u(1)*v(3)
    s(3) = u(1)*v(2) - u(2)*v(1)

end function cross_product

subroutine tb_zyx(alpha, beta, gamma, basis, inv_basis)
    !+ Creates active rotation matrix for z-y'-x" rotation given Tait-Bryan angles
    real(Float64), intent(in)                            :: alpha
        !+ Angle of rotation about z
    real(Float64), intent(in)                            :: beta
        !+ Angle of rotation about y'
    real(Float64), intent(in)                            :: gamma
        !+ Angle of rotation about x"
    real(Float64), dimension(3,3), intent(out)           :: basis
        !+ Rotation matrix/basis for transforming from rotated to non-rotated coordinates
    real(Float64), dimension(3,3), intent(out), optional :: inv_basis
        !+ Inverse basis for reverse transformation

    real(Float64) :: sa, sb, sg, ca, cb, cg

    sa = sin(alpha) ; sb = sin(beta) ; sg = sin(gamma)
    ca = cos(alpha) ; cb = cos(beta) ; cg = cos(gamma)

    basis(1,1) = ca*cb ; basis(1,2) = ca*sb*sg - cg*sa ; basis(1,3) = sa*sg + ca*cg*sb
    basis(2,1) = cb*sa ; basis(2,2) = ca*cg + sa*sb*sg ; basis(2,3) = cg*sa*sb - ca*sg
    basis(3,1) = -sb   ; basis(3,2) = cb*sg            ; basis(3,3) = cb*cg

    if(present(inv_basis)) inv_basis = transpose(basis)

end subroutine tb_zyx

subroutine line_basis(r0, v0, basis, inv_basis)
    !+ Calculates basis from a line with +x in the direction of line
    real(Float64), dimension(3), intent(in)              :: r0
        !+ Starting point of line [cm]
    real(Float64), dimension(3), intent(in)              :: v0
        !+ Direction of line
    real(Float64), dimension(3,3), intent(out)           :: basis
        !+ Basis for transforming from line coordinates to cartesian
    real(Float64), dimension(3,3), intent(out), optional :: inv_basis
        !+ Inverse basis for the reverse transformation cartesian to line

    real(Float64), dimension(3) :: rf
    real(Float64) :: alpha, beta, dis

    rf = r0 + v0
    dis = sqrt(sum((rf - r0)**2))
    beta = asin((r0(3) - rf(3))/dis)
    alpha = atan2(rf(2)-r0(2),rf(1)-r0(1))

    call tb_zyx(alpha,beta,0.d0,basis)

    if(present(inv_basis)) inv_basis = transpose(basis)

end subroutine line_basis

subroutine plane_basis(center, redge, tedge, basis, inv_basis)
    !+ Calculates basis from 3 points on a plane with +z being the plane normal
    real(Float64), dimension(3), intent(in)              :: center
        !+ Plane origin
    real(Float64), dimension(3), intent(in)              :: redge
        !+ Right edge of plane
    real(Float64), dimension(3), intent(in)              :: tedge
        !+ Top edge of plane
    real(Float64), dimension(3,3), intent(out)           :: basis
        !+ Basis for transforming from plane to cartesian coordinates
    real(Float64), dimension(3,3), intent(out), optional :: inv_basis
        !+ Inverse basis for the reverse transformation cartesian to plane

    real(Float64), dimension(3) :: u1,u2,u3

    u1 = (redge - center)
    u1 = u1/norm2(u1)
    u2 = (tedge - center)
    u2 = u2/norm2(u2)
    u3 = cross_product(u1,u2)
    u3 = u3/norm2(u3)

    basis(:,1) = u1
    basis(:,2) = u2
    basis(:,3) = u3

    if(present(inv_basis)) inv_basis = transpose(basis)

end subroutine plane_basis

subroutine line_plane_intersect(l0, l, p0, n, p, t)
    !+ Calculates the intersection of a line and a plane
    real(Float64), dimension(3), intent(in)  :: l0
        !+ Point on line
    real(Float64), dimension(3), intent(in)  :: l
        !+ Ray of line
    real(Float64), dimension(3), intent(in)  :: p0
        !+ Point on plane
    real(Float64), dimension(3), intent(in)  :: n
        !+ Normal vector of plane
    real(Float64), dimension(3), intent(out) :: p
        !+ Line-plane intersect point
    real(Float64), intent(out)               :: t
        !+ "time" to intersect

    real(Float64) :: ldotn

    ldotn = dot_product(l, n)
    if(ldotn.eq.0.0)then
        t = 0.0
    else
        t = dot_product(p0 - l0, n)/ldotn
    endif
    p = l0 + t*l

end subroutine line_plane_intersect

subroutine line_cylinder_intersect(l0, l, p0, p, t)
    !+ Calculates the intersection of a line and a cylinder
    real(Float64), dimension(3), intent(in)  :: l0
        !+ Point on line
    real(Float64), dimension(3), intent(in)  :: l
        !+ Ray of line
    real(Float64), dimension(3), intent(in)  :: p0
        !+ Point on cylinder
    real(Float64), dimension(3), intent(out) :: p
        !+ Line-cylinder intersect point
    real(Float64), intent(out)               :: t
        !+ "time" to intersect

    real(Float64), dimension(2) :: times
    logical, dimension(2) :: mask
    real(Float64) :: r, vx, vy, x0, y0
    real(Float64) :: radicand, npos

    r = sqrt(p0(1) * p0(1) + p0(2) * p0(2))
    x0 = l0(1) ; y0 = l0(2)
    vx = l(1)  ; vy = l(2)

    if((vx.eq.0.d0).and.(vy.eq.0.d0)) then
        t = 0.d0        ! Parallel to a plane tangent to the cylinder
    else
        radicand = r**2 * (vx**2 + vy**2) - (vy * x0 - vx * y0)**2
        if(radicand.lt.0) then
            t = 0.d0    ! Parallel to a plane tangent to the cylinder
        else
            times(1) = (- vx * x0 - vy * y0 - sqrt(radicand)) / (vx**2 + vy**2)
            times(2) = (- vx * x0 - vy * y0 + sqrt(radicand)) / (vx**2 + vy**2)
            mask = times.gt.0
            npos = count(mask)
            if(npos.gt.0) then
                t = minval(times, mask=times.gt.0)
            else
                t = maxval(times, mask=times.le.0)
            endif
        endif
    endif

    p = l0 + l * t

end subroutine line_cylinder_intersect

function in_boundary(bplane, p) result(in_b)
    !+ Indicator function for determining if a point on a plane is within the plane boundary
    type(BoundedPlane), intent(in)          :: bplane
        !+ Plane with boundary
    real(Float64), dimension(3), intent(in) :: p
        !+ Point on plane
    logical :: in_b

    real(Float64), dimension(3) :: pp
    real(Float64) :: hh, hw

    hh = bplane%hh
    hw = bplane%hw
    pp = matmul(bplane%inv_basis, p - bplane%origin)
    in_b = .False.
    SELECT CASE (bplane%shape)
        CASE (1) !Rectangular boundary
            if((abs(pp(1)).le.hw).and. &
               (abs(pp(2)).le.hh)) then
                in_b = .True.
            endif
        CASE (2) !Circular/Ellipsoidal boundary
            if(((hh*pp(1))**2 + (hw*pp(2))**2).le.((hh*hw)**2)) then
                in_b = .True.
            endif
        CASE DEFAULT
            if(inputs%verbose.ge.0) then
                write(*,'("IN_BOUNDARY: Unknown boundary shape: ",i2)') bplane%shape
            endif
            stop
    END SELECT

end function in_boundary

subroutine boundary_edge(bplane, bedge, nb)
    !+ Returns 3 x `nb` array containing points along the BoundedPlane's boundary edge
    type(BoundedPlane), intent(in)             :: bplane
        !+ Bounded plane
    real(Float64), dimension(:,:), intent(out) :: bedge
        !+ Boundary edge points of bounded plane
    integer, intent(out)                       :: nb
        !+ Number of points in boundary edge

    integer :: i
    real(Float64) :: th, dth, x, y
    real(Float64), dimension(4) :: xx, yy

    select case (bplane%shape)
        case (1) !Rectangular boundary
            nb = 4
            if(nb.gt.size(bedge,2)) then
                if(inputs%verbose.ge.0) then
                    write(*,'("BOUNDARY_EDGE: Incompatible boundary edge array : ",i2," > ",i2)') nb, size(bedge,2)
                endif
                stop
            endif
            xx = [-bplane%hw,-bplane%hw,bplane%hw,bplane%hw]
            yy = [-bplane%hh,bplane%hh,bplane%hh,-bplane%hh]
            do i=1,nb
                bedge(:,i) = matmul(bplane%basis,[xx(i),yy(i),0.d0]) + bplane%origin
            enddo
        case (2)
            nb = 50
            if(nb.gt.size(bedge,2)) then
                if(inputs%verbose.ge.0) then
                    write(*,'("BOUNDARY_EDGE: Incompatible boundary edge array : ",i2," > ",i2)') nb, size(bedge,2)
                endif
                stop
            endif
            dth = 2*pi/nb
            do i=1,nb
                th = i*dth
                x = bplane%hw*cos(th)
                y = bplane%hh*sin(th)
                bedge(:,i) = matmul(bplane%basis,[x,y,0.d0]) + bplane%origin
            enddo
        case default
            if(inputs%verbose.ge.0) then
                write(*,'("BOUNDARY_EDGE: Unknown boundary shape: ",i2)') bplane%shape
            endif
            stop
    end select

end subroutine boundary_edge

subroutine gyro_surface(fields, energy, pitch, Ai, gs)
    !+ Calculates the surface of all possible trajectories
    type(LocalEMFields), intent(in) :: fields
        !+ Electromagnetic fields at guiding center
    real(Float64), intent(in)       :: energy
        !+ Energy of particle
    real(Float64), intent(in)       :: pitch
        !+ Particle pitch w.r.t the magnetic field
    real(Float64), intent(in)       :: Ai
        !+ Ion atomic mass
    type(GyroSurface), intent(out)  :: gs
        !+ Gyro-surface

    integer :: i
    real(Float64) :: alpha, vabs, omega
    real(Float64), dimension(3,3) :: s

    vabs  = sqrt(energy/(v2_to_E_per_amu*Ai))
    omega= (fields%b_abs*e0)/(Ai*mass_u)
    alpha = vabs/omega

    gs%omega = omega
    gs%v = vabs
    gs%axes(1) = alpha*sqrt(1-pitch**2)
    gs%axes(2) = alpha*sqrt(1-pitch**2)
    gs%axes(3) = pitch*alpha

    s = 0.d0
    s(1,1) = gs%axes(1)**(-2)
    s(2,2) = gs%axes(2)**(-2)
    s(3,3) = -gs%axes(3)**(-2)

    gs%center = fields%pos

    gs%basis(:,1) = fields%a_norm
    gs%basis(:,2) = fields%c_norm
    gs%basis(:,3) = fields%b_norm

    gs%A = matmul(gs%basis,matmul(s,transpose(gs%basis)))

end subroutine gyro_surface

subroutine line_gyro_surface_intersect(r0, v0, gs, t)
    !+ Calculates the times of intersection of a line and a gyro-surface
    real(Float64), dimension(3), intent(in)  :: r0
        !+ Point on line
    real(Float64), dimension(3), intent(in)  :: v0
        !+ Direction of line
    type(GyroSurface), intent(in)            :: gs
        !+ Gyro-surface
    real(Float64), dimension(2), intent(out) :: t
        !+ "time" to intersect

    real(Float64), dimension(3) :: rr
    real(Float64) :: a, b, c, d, tp, tm

    rr = r0 - gs%center
    a = dot_product(v0, matmul(gs%A,v0))
    b = dot_product(rr, matmul(gs%A,v0)) + dot_product(v0,matmul(gs%A,rr))
    c = dot_product(rr, matmul(gs%A,rr)) - 1.0

    d = b**2 - 4*a*c
    if(d.lt.0.0) then
        t = 0.0
        return
    endif

    t(1) = (-b - sqrt(d))/(2*a)
    t(2) = (-b + sqrt(d))/(2*a)

end subroutine line_gyro_surface_intersect

subroutine gyro_surface_coordinates(gs, p, u)
    !+ Calculates the parametric coordinates, `u`, of point `p` on the gyro_surface
    type(GyroSurface), intent(in)            :: gs
        !+ Gyro_surface
    real(Float64), dimension(3), intent(in)  :: p
        !+ Point on gyro_surface
    real(Float64), dimension(2), intent(out) :: u
        !+ Parametric coordinates (gyro-angle, t)

    real(Float64), dimension(3) :: pp
    real(Float64) :: t, a, b, c, d, thm, thp, dp, dm, th
    integer :: i

    pp = matmul(transpose(gs%basis),p - gs%center)
    t = pp(3)/gs%axes(3)
    a = gs%axes(1) + gs%axes(2)*t
    b = gs%axes(2) - gs%axes(1)*t
    d = pp(1) + pp(2)
    c = max(min(d/sqrt(a**2 + b**2),1.d0),-1.d0)

    thm = -acos(c) + atan2(b,a)
    thp =  acos(c) + atan2(b,a)

    dm = norm2([gs%axes(1)*(cos(thm) - t*sin(thm)), &
                gs%axes(2)*(sin(thm) + t*cos(thm)), &
                gs%axes(3)*t ] - pp)
    dp = norm2([gs%axes(1)*(cos(thp) - t*sin(thp)), &
                gs%axes(2)*(sin(thp) + t*cos(thp)), &
                gs%axes(3)*t ] - pp)

    th = thm - pi/2
    if(dp.le.dm) th = thp - pi/2
    if(th.lt.0.0) th = th + 2*pi
    u = [th, t/gs%omega]

end subroutine gyro_surface_coordinates

subroutine gyro_trajectory(gs, theta, ri, vi)
    !+ Calculate particle trajectory for a given gyro-angle and gyro-surface
    type(GyroSurface), intent(in) :: gs
        !+ Gyro-Surface
    real(Float64), intent(in) :: theta
        !+ Gyro-angle
    real(Float64), dimension(3) :: ri
        !+ Particle position
    real(Float64), dimension(3) :: vi
        !+ Particle Velocity

    real(Float64) :: a,b,c,th
    a = gs%axes(1)
    b = gs%axes(2)
    c = gs%axes(3)
    th = theta + pi/2
    ri = matmul(gs%basis, [a*cos(th), b*sin(th), 0.d0]) + gs%center
    vi = gs%omega*matmul(gs%basis, [-a*sin(th), b*cos(th), c])

end subroutine gyro_trajectory

function in_gyro_surface(gs, p) result(in_gs)
    !+ Indicator function for determining if a point is inside the gyro_surface
    type(GyroSurface), intent(in)           :: gs
        !+ Gyro-surface
    real(Float64), dimension(3), intent(in) :: p
        !+ Point
    logical :: in_gs

    real(Float64), dimension(3) :: pp

    pp = p - gs%center
    in_gs = dot_product(pp, matmul(gs%A, pp)).le.1.d0

end function in_gyro_surface

subroutine gyro_range(b, gs, gyrange, nrange)
    !+ Calculates the range(s) of gyro-angles that would land within a bounded plane
    type(BoundedPlane), intent(in)             :: b
        !+ Bounded Plane
    type(GyroSurface), intent(in)              :: gs
        !+ Gyro-surface
    real(Float64), dimension(2,4), intent(out) :: gyrange
        !+ (theta, dtheta) values
    integer, intent(out) :: nrange
        !+ Number of ranges. `1 <= nrange <= 4`

    integer :: nb, i, j, ninter
    logical :: in_gs
    logical, dimension(8) :: cross = .False.
    real(Float64) :: t_p, th1, th2, dth
    real(Float64), dimension(2) :: u_cur, t_i
    real(Float64), dimension(3) :: rc, p_pre, p_cur, v0, ri
    real(Float64), dimension(2,8) :: u
    real(Float64), dimension(3,50) :: bedge

    nrange = 0
    gyrange = 0.d0

    call line_plane_intersect(gs%center, gs%basis(:,3), b%origin, b%basis(:,3), rc, t_p)
    if(t_p.eq.0.0) return

    call boundary_edge(b, bedge, nb)
    p_pre = bedge(:,1)
    in_gs = in_gyro_surface(gs, p_pre)

    ninter = 0
    u = 0.d0
    boundary_loop: do i=1,nb
        p_cur = bedge(:,modulo(i,nb)+1)
        v0 = p_cur - p_pre
        call line_gyro_surface_intersect(p_pre, v0, gs, t_i)
        do j=1,2
            if((t_i(j).gt.0.0).and.(t_i(j).lt.1.0)) then
                ri = p_pre + t_i(j)*v0
                call gyro_surface_coordinates(gs, ri, u_cur)
                if(u_cur(2).gt.0.0) then
                    in_gs = .not.in_gs
                    ninter = ninter + 1
                    cross(ninter) = in_gs
                    u(:,ninter) = u_cur
                endif
            endif
        enddo
        p_pre = p_cur
    enddo boundary_loop

    if(ninter.eq.0) then
        if(in_boundary(b, rc)) then
            nrange = 1
            gyrange(:,1) = [0.d0,2*pi]
        endif
        return
    endif

    do i=1, ninter
        if(cross(i)) then
            th1 = u(1,i)
            j = modulo(i,ninter) + 1
            th2 = u(1,j)
            dth = th2-th1
            nrange = nrange + 1
            if(dth.gt.0.0) then
                gyrange(:,nrange) = [th1, dth]
            else
                gyrange(:,nrange) = [th2, -dth]
            endif
        endif
        !! OpenMP with multiple threads is duplicating gyro-ranges for some markers
        !! causing double counting and I don't know why.
        !! It should be very unlikely for multiple gyro-ranges to occur so for
        !! now I'm including this cludge to force only one gyro-range when using
        !! OpenMP.
#ifdef _OMP
        if(nrange.eq.1) exit
#endif
    enddo

end subroutine gyro_range

subroutine npa_gyro_range(ichan, gs, gyrange, nrange)
    !+ Calculates range of gyro-angles that would hit the NPA detector
    integer, intent(in) :: ichan
        !+ Index of NPA detector
    type(GyroSurface), intent(in) :: gs
    real(Float64), dimension(2,4), intent(out) :: gyrange
    integer, intent(out) :: nrange

    type(LocalEMFields) :: fields
    integer :: i, j, a_nrange, d_nrange
    real(Float64) :: a0, a, b, c, d
    real(Float64), dimension(2,4) :: a_gyrange, d_gyrange

    nrange = 0
    gyrange = 0.d0

    call gyro_range(npa_chords%det(ichan)%aperture, gs, a_gyrange, a_nrange)
    if(a_nrange.eq.0) return

    call gyro_range(npa_chords%det(ichan)%detector, gs, d_gyrange, d_nrange)
    if(d_nrange.eq.0) return

    if((a_nrange.eq.1).and.approx_eq(a_gyrange(2,1),2*pi,1d-6)) then
        gyrange = d_gyrange
        nrange = d_nrange
        return
    endif

    if((d_nrange.eq.1).and.approx_eq(d_gyrange(2,1),2*pi,1d-6)) then
        gyrange = a_gyrange
        nrange = a_nrange
        return
    endif

    do i=1,a_nrange
        do j=1, d_nrange
            a0 = 0.d0
            if(d_gyrange(1,j).gt.a_gyrange(1,i)) then
                a0 = a_gyrange(1,i)
                a = 0.d0
                b = modulo(a_gyrange(1,i) + a_gyrange(2,i) - a0, 2*pi)
                c = modulo(d_gyrange(1,j) - a0, 2*pi)
                d = modulo(d_gyrange(1,j) + d_gyrange(2,j) - a0, 2*pi)
            else
                a0 = d_gyrange(1,j)
                a = 0.d0
                b = modulo(d_gyrange(1,j) + d_gyrange(2,j) - a0, 2*pi)
                c = modulo(a_gyrange(1,i) - a0, 2*pi)
                d = modulo(a_gyrange(1,i) + a_gyrange(2,i) - a0, 2*pi)
            endif
            if((c.lt.b).or.(d.lt.c)) then
                if(c.lt.d) then
                    nrange = nrange + 1
                    gyrange(:,nrange) = [a0 + c, min(d-c,b-c)]
                else
                    nrange = nrange + 1
                    gyrange(:,nrange) = [a0, d]
                    nrange = nrange + 1
                    gyrange(:,nrange) = [a0+c, b-c]
                endif
            endif
        enddo
    enddo

end subroutine npa_gyro_range

subroutine hit_npa_detector(r0, v0, d_index, rd, det)
    !+ Routine to check if a particle will hit a NPA detector
    real(Float64), dimension(3), intent(in)            :: r0
        !+ Starting point of particle
    real(Float64), dimension(3), intent(in)            :: v0
        !+ Particle velocity
    integer, intent(out)                               :: d_index
        !+ Index of NPA detector. Zero if particle doesn't hit
    real(Float64), dimension(3), intent(out), optional :: rd
        !+ Point where particle hit detector
    integer, intent(in), optional :: det
        !+ Index of NPA detector to check

    real(Float64), dimension(3) :: d, a
    real(Float64) :: t_a,t_d
    integer :: i, s, ndet

    if(present(det)) then
        s = det
        ndet = det
    else
        s = 1
        ndet = npa_chords%nchan
    endif

    d_index = 0
    detector_loop: do i=s,ndet
        !! Find where trajectory crosses detector plane
        call line_plane_intersect(r0,v0,npa_chords%det(i)%detector%origin, &
             npa_chords%det(i)%detector%basis(:,3),d,t_d)

        !! Find where trajectory crosses aperture plane
        call line_plane_intersect(r0,v0,npa_chords%det(i)%aperture%origin, &
             npa_chords%det(i)%aperture%basis(:,3),a,t_a)

        !! If both points are in plane boundaries and the
        !! particle is heading toward the detector then its a hit
        if(in_boundary(npa_chords%det(i)%aperture,a) .and. &
           in_boundary(npa_chords%det(i)%detector,d) .and. &
           (t_d.gt.0.0) ) then
            d_index = i
            exit detector_loop
        endif
    enddo detector_loop
    if(present(rd)) rd = d

end subroutine hit_npa_detector

subroutine xyz_to_uvw(xyz, uvw)
    !+ Convert beam coordinate `xyz` to machine coordinate `uvw`
    real(Float64), dimension(3), intent(in)  :: xyz
    real(Float64), dimension(3), intent(out) :: uvw

    real(Float64), dimension(3) :: origin
    real(Float64), dimension(3,3) :: basis

    origin = beam_grid%origin
    basis = beam_grid%basis
    uvw = matmul(basis,xyz)
    uvw = uvw + origin

end subroutine xyz_to_uvw

subroutine xyz_to_cyl(xyz, cyl)
    !+ Convert beam coordinate `xyz` to cylindrical coordinate `cyl`
    real(Float64), dimension(3), intent(in)  :: xyz
    real(Float64), dimension(3), intent(out) :: cyl

    real(Float64), dimension(3) :: uvw

    call xyz_to_uvw(xyz, uvw)
    call uvw_to_cyl(uvw, cyl)

end subroutine xyz_to_cyl

subroutine uvw_to_xyz(uvw,xyz)
    !+ Convert machine coordinate `uvw` to beam coordinate `xyz`
    real(Float64), dimension(3), intent(in)  :: uvw
    real(Float64), dimension(3), intent(out) :: xyz

    real(Float64), dimension(3) :: origin, uvw_p
    real(Float64), dimension(3,3) :: basis

    origin = beam_grid%origin
    basis = beam_grid%inv_basis
    uvw_p = uvw - origin
    xyz = matmul(basis,uvw_p)

end subroutine uvw_to_xyz

subroutine cyl_to_uvw(cyl, uvw)
    !+ Convert cylindrical coordinate `cyl` to machine coordinate `uvw`
    real(Float64), dimension(3), intent(in)  :: cyl
    real(Float64), dimension(3), intent(out) :: uvw

    uvw(1) = cyl(1) * cos(cyl(3))
    uvw(2) = cyl(1) * sin(cyl(3))
    uvw(3) = cyl(2)

end subroutine cyl_to_uvw

subroutine cyl_to_xyz(cyl, xyz)
    !+ Convert cylindrical coordinate `cyl` to beam coordinate `xyz`
    real(Float64), dimension(3), intent(in)  :: cyl
    real(Float64), dimension(3), intent(out) :: xyz

    real(Float64), dimension(3) :: uvw

    call cyl_to_uvw(cyl, uvw)
    call uvw_to_xyz(uvw, xyz)

end subroutine cyl_to_xyz

subroutine uvw_to_cyl(uvw, cyl)
    !+ Convert machine coordinate `uvw` to cylindrical coordinate `cyl`
    real(Float64), dimension(3), intent(in)  :: uvw
    real(Float64), dimension(3), intent(out) :: cyl

    cyl(1) = sqrt(uvw(1)*uvw(1) + uvw(2)*uvw(2))
    cyl(2) = uvw(3)
    cyl(3) = atan2(uvw(2),uvw(1))

end subroutine uvw_to_cyl

subroutine convert_sightline_to_xyz(ie3, ist, iray, ich, xyz, v3_xyz)
    !+ Convert sightline position and velocity from cylindrical coordinate `rpz` to beam coordinate `xyz`
    integer, intent(in) :: ie3
        !+ CFPD energy index
    integer, intent(in) :: ist
        !+ Step index
    integer, intent(in) :: iray
        !+ "Ray" index
    integer, intent(in) :: ich
        !+ Detector channel
    real(Float64), dimension(3), intent(out) :: xyz
        !+ Sightline position in beam coordinates
    real(Float64), dimension(3), intent(out) :: v3_xyz
        !+ Sightline velocity in beam coordinates

    real(Float64), dimension(3) :: rpz, uvw, v3_rpz, v3_uvw
    real(Float64) :: phi

    !! Position
    rpz(1) = ctable%sightline(ie3,4,ist,iray,ich)
    rpz(2) = ctable%sightline(ie3,5,ist,iray,ich)
    rpz(3) = ctable%sightline(ie3,6,ist,iray,ich)

    phi = rpz(2)
    uvw(1) = rpz(1)*cos(phi)
    uvw(2) = rpz(1)*sin(phi)
    uvw(3) = rpz(3)
    call uvw_to_xyz(uvw, xyz)

    !! Velocity
    v3_rpz(1) = ctable%sightline(ie3,1,ist,iray,ich)
    v3_rpz(2) = ctable%sightline(ie3,2,ist,iray,ich)
    v3_rpz(3) = ctable%sightline(ie3,3,ist,iray,ich)

    v3_uvw(1) = v3_rpz(1)*cos(phi) - v3_rpz(2)*sin(phi)
    v3_uvw(2) = v3_rpz(1)*sin(phi) + v3_rpz(2)*cos(phi)
    v3_uvw(3) = v3_rpz(3)

    v3_xyz = matmul(beam_grid%inv_basis, v3_uvw)

end subroutine convert_sightline_to_xyz

subroutine grid_intersect(r0, v0, length, r_enter, r_exit, center_in, lwh_in, passive)
    !+ Calculates a particles intersection length with the [[libfida:beam_grid]]
    real(Float64), dimension(3), intent(in)           :: r0
        !+ Initial position of particle [cm]
    real(Float64), dimension(3), intent(in)           :: v0
        !+ Velocity of particle [cm/s]
    real(Float64), intent(out)                        :: length
        !+ Intersection length [cm]
    real(Float64), dimension(3), intent(out)          :: r_enter
        !+ Point where particle enters
    real(Float64), dimension(3), intent(out)          :: r_exit
        !+ Point where particle exits
    real(Float64), dimension(3), intent(in), optional :: center_in
        !+ Alternative grid center
    real(Float64), dimension(3), intent(in), optional :: lwh_in
        !+ Alternative grid [length,width,height]
    logical, intent(in), optional                     :: passive
        !+ Calculates a particles intersection length with the [[libfida:pass_grid]]

    real(Float64), dimension(3,6) :: ipnts
    real(Float64), dimension(3) :: ri, vi
    real(Float64), dimension(3) :: center
    real(Float64), dimension(3) :: lwh
    integer, dimension(6) :: side_inter
    integer, dimension(2) :: ind
    integer :: i, j, nunique, ind1, ind2
    real(Float64) :: dlength, max_length
    logical :: ing, pas

    pas = .False.
    if(present(passive)) pas = passive

    if(pas) then
        ing = in_passive_grid(r0)
        if (ing) then
            length = 0.d0
            return
        endif

        ri = r0 ; vi = v0
        dlength = 2.0 !cm
        max_length=0.0

        do while (.not.ing)
            ri = ri + vi*dlength ! move dlength
            ing = in_passive_grid(ri)
            max_length = max_length + dlength
            if(max_length.gt.1d3) then
                length = 0.d0
                return
            endif
        enddo
        r_enter = ri

        do while (ing)
            ri = ri + vi*dlength
            ing =  in_passive_grid(ri)
        enddo
        r_exit = ri

        length = sqrt(sum((r_exit-r_enter)**2))
    else
        if (present(center_in)) then
            center = center_in
        else
            center = beam_grid%center
        endif

        if (present(lwh_in)) then
            lwh = lwh_in
        else
            lwh = beam_grid%lwh
        endif

        side_inter = 0
        ipnts = 0.d0
        do i=1,6
            j = int(ceiling(i/2.0))
            if (j.eq.1) ind = [2,3]
            if (j.eq.2) ind = [1,3]
            if (j.eq.3) ind = [1,2]
            if (abs(v0(j)).gt.0.d0) then
                ipnts(:,i) = r0 + v0*( ( (center(j) + &
                             (mod(i,2) - 0.5)*lwh(j)) - r0(j))/v0(j) )
                if ((abs(ipnts(ind(1),i) - center(ind(1))).le.(0.5*lwh(ind(1)))).and. &
                    (abs(ipnts(ind(2),i) - center(ind(2))).le.(0.5*lwh(ind(2))))) then
                    side_inter(i) = 1
                endif
            endif
        enddo

        length = 0.d0
        r_enter = r0
        r_exit  = r0
        ind1=0
        ind2=0
        if (sum(side_inter).ge.2) then
            ! Find first intersection side
            i=1
            do while (i.le.6)
                if(side_inter(i).eq.1) exit
                i=i+1
            enddo
            ind1=i
            !Find number of unique points
            nunique = 0
            do i=ind1+1,6
                if (side_inter(i).ne.1) cycle
                if (sqrt( sum( ( ipnts(:,i)-ipnts(:,ind1) )**2 ) ).gt.0.001) then
                    ind2=i
                    nunique = 2
                    exit
                endif
            enddo

            if(nunique.eq.2) then
                vi = ipnts(:,ind2) - ipnts(:,ind1)
                if (dot_product(v0,vi).gt.0.0) then
                    r_enter = ipnts(:,ind1)
                    r_exit  = ipnts(:,ind2)
                else
                    r_enter = ipnts(:,ind2)
                    r_exit  = ipnts(:,ind1)
                endif
                length = sqrt(sum((r_exit - r_enter)**2))
            endif
        endif
    endif

end subroutine grid_intersect

function in_grid(xyz) result(ing)
    !+ Determines if a position `pos` is in the [[libfida:beam_grid]]
    real(Float64), dimension(3), intent(in) :: xyz
        !+ Position in beam grid coordinates [cm]
    logical :: ing
        !+ Indicates whether the position is in the beam grid

    real(Float64) :: tol = 1.0d-10

    if((approx_ge(xyz(1),beam_grid%xmin,tol).and.approx_le(xyz(1),beam_grid%xmax,tol)).and. &
       (approx_ge(xyz(2),beam_grid%ymin,tol).and.approx_le(xyz(2),beam_grid%ymax,tol)).and. &
       (approx_ge(xyz(3),beam_grid%zmin,tol).and.approx_le(xyz(3),beam_grid%zmax,tol))) then
        ing = .True.
    else
        ing = .False.
    endif

end function

function in_passive_grid(uvw) result(ing)
    !+ Determines if a position `pos` is in the [[libfida:pass_grid]]
    real(Float64), dimension(3), intent(in) :: uvw
        !+ Position in machine coordinates [cm]
    logical :: ing
        !+ Indicates whether the position is in the passive neutral grid

    real(Float64), dimension(3) :: cyl
    real(Float64) :: phimin, phimax
    real(Float64) :: tol = 1.0d-10

    phimin = minval(pass_grid%phi)
    phimax = maxval(pass_grid%phi)

    call uvw_to_cyl(uvw, cyl)

    if((approx_ge(cyl(1),minval(pass_grid%r),tol).and.approx_le(cyl(1),maxval(pass_grid%r),tol)).and. &
       (approx_ge(cyl(2),minval(pass_grid%z),tol).and.approx_le(cyl(2),maxval(pass_grid%z),tol)).and. &
       ((approx_ge(cyl(3),phimin,tol).and.approx_le(cyl(3),phimax,tol)).or. &
       (approx_ge(modulo(cyl(3),2*pi),phimin,tol).and.approx_le(modulo(cyl(3),2*pi),phimax,tol)))) then
        ing = .True.
    else
        ing = .False.
    endif

end function

subroutine circle_grid_intersect(r0, e1, e2, radius, beam_grid_phi_enter, beam_grid_phi_exit)
    !+ Calculates the intersection arclength of a circle with the [[libfida:beam_grid]]
    real(Float64), dimension(3), intent(in) :: r0
        !+ Position of center enter of the circle in beam grid coordinates [cm]
    real(Float64), dimension(3), intent(in) :: e1
        !+ Unit vector pointing towards (R, 0) (r,phi) position of the circle in beam grid coordinates
    real(Float64), dimension(3), intent(in) :: e2
        !+ Unit vector pointing towards (R, pi/2) (r,phi) position of the circle in beam grid coordinates
    real(Float64), intent(in)               :: radius
        !+ Radius of circle [cm]
    real(Float64), intent(out)              :: beam_grid_phi_enter
        !+ Phi value where the circle entered the [[libfida:beam_grid]] [rad]
    real(Float64), intent(out)              :: beam_grid_phi_exit
        !+ Phi value where the circle exits the [[libfida:beam_grid]] [rad]

    real(Float64), dimension(3) :: i1_p,i1_n,i2_p,i2_n
    real(Float64), dimension(4) :: d
    real(Float64), dimension(6) :: p, gams
    real(Float64), dimension(4,6) :: phi
    logical, dimension(4,6) :: inter
    integer, dimension(6) :: n
    integer :: i
    real(Float64) :: alpha,beta,delta,sinx1,cosx1,sinx2,cosx2,tmp
    real(Float64) :: tol = 1.0d-10
    logical :: r0_ing

    p = [beam_grid%xmin, beam_grid%xmax, &
         beam_grid%ymin, beam_grid%ymax, &
         beam_grid%zmin, beam_grid%zmax ]
    n = [1, 1, 2, 2, 3, 3]
    inter = .False.
    phi = 0.d0

    r0_ing = in_grid(r0)
    do i=1,6
        alpha = e2(n(i))
        beta = e1(n(i))
        if((alpha.eq.0.0).and.(beta.eq.0.0)) cycle
        gams(i) = (p(i) - r0(n(i)))/radius
        delta = alpha**4 + (alpha**2)*(beta**2 - gams(i)**2)
        if (delta.ge.0.0) then
            cosx1 = (gams(i)*beta + sqrt(delta))/(alpha**2 + beta**2)
            if ((cosx1**2).le.1.0) then
                sinx1 = sqrt(1 - cosx1**2)
                i1_p = r0 + radius*cosx1*e1 + radius*sinx1*e2
                i1_n = r0 + radius*cosx1*e1 - radius*sinx1*e2
                if (approx_eq(i1_p(n(i)),p(i),tol).and.in_grid(i1_p)) then
                    inter(1,i) = .True.
                    phi(1,i) = atan2(sinx1,cosx1)
                endif
                if (approx_eq(i1_n(n(i)),p(i),tol).and.in_grid(i1_n)) then
                    inter(2,i) = .True.
                    phi(2,i) = atan2(-sinx1,cosx1)
                endif
            endif

            if(delta.gt.0.0) then
                cosx2 = (gams(i)*beta - sqrt(delta))/(alpha**2 + beta**2)
                if ((cosx2**2).le.1.0) then
                    sinx2 = sqrt(1 - cosx2**2)
                    i2_p = r0 + radius*cosx2*e1 + radius*sinx2*e2
                    i2_n = r0 + radius*cosx2*e1 - radius*sinx2*e2
                    if (approx_eq(i2_p(n(i)),p(i),tol).and.in_grid(i2_p)) then
                        inter(3,i) = .True.
                        phi(3,i) = atan2(sinx2,cosx2)
                    endif
                    if (approx_eq(i2_n(n(i)),p(i),tol).and.in_grid(i2_n)) then
                        inter(4,i) = .True.
                        phi(4,i) = atan2(-sinx2,cosx2)
                    endif
                endif
            endif
        endif
    enddo

    beam_grid_phi_enter = 0.d0
    beam_grid_phi_exit = 0.d0
    if (count(inter).gt.2) then
        write(*,'("CIRCLE_GRID_INTERSECT: Circle intersects grid more than 2 times: ",i2)') count(inter)
        return
    endif

    if(any(inter)) then
        beam_grid_phi_enter = minval(phi,inter)
        beam_grid_phi_exit = maxval(phi,inter)
        if(r0_ing.and.any(count(inter,1).ge.2)) then
            if((beam_grid_phi_exit - beam_grid_phi_enter) .lt. pi) then
                tmp = beam_grid_phi_enter
                beam_grid_phi_enter = beam_grid_phi_exit
                beam_grid_phi_exit = tmp + 2*pi
            endif
        else
            if((beam_grid_phi_exit - beam_grid_phi_enter) .gt. pi) then
                tmp = beam_grid_phi_enter
                beam_grid_phi_enter = beam_grid_phi_exit
                beam_grid_phi_exit = tmp + 2*pi
            endif
        endif
        if(approx_eq(beam_grid_phi_exit-beam_grid_phi_enter,pi,tol).and.r0_ing) then
            beam_grid_phi_enter = 0.0
            beam_grid_phi_exit = 2*pi
        endif
    else
        if(r0_ing) then
            call grid_intersect(r0, e1, tmp, i1_n,i1_p)
            call grid_intersect(r0, e2, tmp, i2_n,i2_p)
            d(1) = norm2(r0 - i1_n)/radius
            d(2) = norm2(r0 - i1_p)/radius
            d(3) = norm2(r0 - i2_n)/radius
            d(4) = norm2(r0 - i2_p)/radius
            if(all(d.ge.1.0)) then
                beam_grid_phi_enter = 0.d0
                beam_grid_phi_exit = 2.d0*pi
            endif
        endif
    endif

end subroutine circle_grid_intersect

subroutine get_indices(pos, ind)
    !+ Find closests [[libfida:beam_grid]] indices `ind` to position `pos`
    real(Float64),  dimension(3), intent(in)  :: pos
        !+ Position [cm]
    integer(Int32), dimension(3), intent(out) :: ind
        !+ Closest indices to position

    real(Float64),  dimension(3) :: mini
    integer(Int32), dimension(3) :: maxind
    integer :: i

    maxind(1) = beam_grid%nx
    maxind(2) = beam_grid%ny
    maxind(3) = beam_grid%nz

    mini(1) = minval(beam_grid%xc) - 0.5*beam_grid%dr(1)
    mini(2) = minval(beam_grid%yc) - 0.5*beam_grid%dr(2)
    mini(3) = minval(beam_grid%zc) - 0.5*beam_grid%dr(3)

    do i=1,3
        ind(i) = floor((pos(i)-mini(i))/beam_grid%dr(i)) + 1
        if (ind(i).gt.maxind(i)) ind(i)=maxind(i)
        if (ind(i).lt.1) ind(i)=1
    enddo

end subroutine get_indices

subroutine get_passive_grid_indices(pos, ind, input_coords)
    !+ Find closest [[libfida:pass_grid]] indices `ind` to position `pos`
    real(Float64),  dimension(3), intent(in)  :: pos
        !+ Position [cm]
    integer(Int32), dimension(3), intent(out) :: ind
        !+ Closest indices to position
    integer, intent(in), optional             :: input_coords
        !+ Indicates coordinate system of `pos`. Beam grid (0), machine (1) and cylindrical (2)

    real(Float64),  dimension(3) :: mini, differentials, loc
    integer(Int32), dimension(3) :: maxind
    integer :: i, ics

    if(present(input_coords)) then
        ics = input_coords
    else
        ics = 2
    endif

    if(ics.eq.0) then
        call xyz_to_cyl(pos, loc)
    endif

    if(ics.eq.1) then
        loc(1) = sqrt(pos(1)*pos(1) + pos(2)*pos(2))
        loc(2) = pos(3)
        loc(3) = atan2(pos(2),pos(1))
    endif

    if(ics.eq.2) then
        loc(1) = pos(1)
        loc(2) = pos(2)
        loc(3) = pos(3)
    endif

    maxind(1) = pass_grid%nr
    maxind(2) = pass_grid%nz
    maxind(3) = pass_grid%nphi

    mini(1) = minval(pass_grid%r)
    mini(2) = minval(pass_grid%z)
    mini(3) = minval(pass_grid%phi)

    differentials(1) = pass_grid%dr
    differentials(2) = pass_grid%dz
    differentials(3) = pass_grid%dphi

    do i=1,3
        ind(i) = floor((loc(i)-mini(i))/differentials(i)) + 1
        if (ind(i).gt.maxind(i)) ind(i)=maxind(i)
        if (ind(i).lt.1) ind(i)=1
    enddo

end subroutine get_passive_grid_indices

subroutine get_interpolation_grid_indices(pos, ind, input_coords)
    !+ Find closest [[libfida:inter_grid]] indices `ind` to position `pos`
    real(Float64),  dimension(3), intent(in)  :: pos
        !+ Position [cm]
    integer(Int32), dimension(3), intent(out) :: ind
        !+ Closest indices to position
    integer, intent(in), optional             :: input_coords
        !+ Indicates coordinate system of `pos`. Beam grid (0), machine (1) and cylindrical (2)

    real(Float64),  dimension(3) :: mini, differentials, loc
    integer(Int32), dimension(3) :: maxind
    integer :: i, ics

    if(present(input_coords)) then
        ics = input_coords
    else
        ics = 2
    endif

    if(ics.eq.1) then
        loc(1) = sqrt(pos(1)*pos(1) + pos(2)*pos(2))
        loc(2) = pos(3)
        loc(3) = atan2(pos(2),pos(1))
    endif

    if(ics.eq.2) then
        loc(1) = pos(1)
        loc(2) = pos(2)
        loc(3) = pos(3)
    endif

    maxind(1) = inter_grid%nr
    maxind(2) = inter_grid%nz
    maxind(3) = inter_grid%nphi

    mini(1) = minval(inter_grid%r)
    mini(2) = minval(inter_grid%z)
    mini(3) = minval(inter_grid%phi)

    differentials(1) = inter_grid%dr
    differentials(2) = inter_grid%dz
    differentials(3) = inter_grid%dphi

    do i=1,3
        ind(i) = floor((loc(i)-mini(i))/differentials(i)) + 1
        if (ind(i).gt.maxind(i)) ind(i)=maxind(i)
        if (ind(i).lt.1) ind(i)=1
    enddo

end subroutine get_interpolation_grid_indices

subroutine get_plasma_extrema(r0, v0, extrema, x0, y0)
    !+ Returns extrema points where line(s) parametrized by `r0` and `v0` intersect the plasma boudnary
    real(Float64), dimension(:,:), intent(in)            :: r0, v0
    !+ Arrays the define line(s) in machine coordinates
    real(Float64), dimension(2,3), intent(out)           :: extrema
    !+ Minimum and maximumm R, Z, and Phi points
    real(Float64), dimension(:), intent(in), optional  :: x0, y0
    !+ Additional x and y points to consider

    real(Float64), dimension(:,:), allocatable :: xy_in, xy_out
    real(Float64), dimension(:,:), allocatable :: cyl_in, cyl_out
    real(Float64), dimension(:), allocatable  :: x, y, xlo, xhi
    real(Float64), dimension(:), allocatable  :: r, z, phi
    logical, dimension(:), allocatable :: skip

    real(Float64), dimension(3) :: ri, vi
    real(Float64) :: max_length, dlength
    integer :: i, iin, iout, ilo, ihi, dlo, dhi, dim, nlines, d
    logical :: inp

    nlines = size(r0(1,:))
    allocate(xy_in(2,nlines), xy_out(2,nlines))
    allocate(cyl_in(3,nlines), cyl_out(3,nlines))
    allocate(skip(nlines))

    dlength = 3.0 !cm
    skip = .False.

    loop_over_channels: do i=1, nlines
        ri = r0(:,i)
        vi = v0(:,i)
        vi = vi/norm2(vi)

        ! Find the position that the los first intersects the plasma
        call in_plasma(ri, inp, input_coords=1)
        max_length=0.0
        do while (.not.inp)
            ri = ri + vi*dlength ! move dlength
            call in_plasma(ri, inp, input_coords=1)
            max_length = max_length + dlength
            if(max_length.gt.1d9) then
                skip(i) = .True. !used below to skip los that do not intersect the plasma
                cycle loop_over_channels
            endif
        enddo

        xy_in(1,i) = ri(1) ; xy_in(2,i) = ri(2)
        call uvw_to_cyl(ri, cyl_in(:,i))

        ! Find the position that the los intersects upon exiting the plasma
        do while (inp)
            ri = ri + vi*dlength
            call in_plasma(ri, inp, input_coords=1)
        enddo

        xy_out(1,i) = ri(1) ; xy_out(2,i) = ri(2)
        call uvw_to_cyl(ri, cyl_out(:,i))

    enddo loop_over_channels

    dim = 2*count(.not.skip) ! 2 for enter and exit

    d = 0
    if (present(x0).and.present(y0)) then
        d = size(x0)
        if (size(x0).ne.size(y0)) then
            if(inputs%verbose.ge.0) then
                write(*,'("GET_PLASMA_EXTREMA: Sizes of X and Y input arrays are not identical")')
                stop
            endif
        endif
    endif

    allocate(x(dim+d), y(dim+d), r(dim), z(dim), phi(dim))

    iin = 1 ; iout = 2
    do i=1, nlines
        if (.not.skip(i)) then
            x(iin) = xy_in(1,i)  ; x(iout) = xy_out(1,i)
            y(iin) = xy_in(2,i)  ; y(iout) = xy_out(2,i)
            r(iin) = cyl_in(1,i) ; r(iout) = cyl_out(1,i)
            z(iin) = cyl_in(2,i) ; z(iout) = cyl_out(2,i)
            iin = iin + 2 ; iout = iout + 2
        endif
    enddo

    extrema(1,1) = minval(r) ; extrema(2,1) = maxval(r)
    extrema(1,2) = minval(z) ; extrema(2,2) = maxval(z)

    ! Append extra x and y points
    if(d.gt.0) then
        x(dim+1:dim+d) = x0
        y(dim+1:dim+d) = y0
    endif

    !! Domain is between 0 and 2 pi if all x points are left of the line x=0
    !! Else domain is between -pi and pi
    dlo = count(y.le.0) ; dhi = count(y.gt.0)
    allocate(xlo(dlo), xhi(dhi))

    ilo = 0 ; ihi = 0
    do i=1, size(x)
        if (y(i).le.0.d10) then
            ilo = ilo + 1
            xlo(ilo) = x(i)
        else
            ihi = ihi + 1
            xhi(ihi) = x(i)
        endif
    enddo

    phi = atan2(y, x)
    if (all(x.le.0) & !none in quadrant 1 and 4
        .or.(all(xlo.le.0).and.(size(x).ne.size(xhi))) & !quadrant 3
        .or.(all(xhi.le.0).and.(size(x).ne.size(xlo)))) then !quadrant 2
        phi = modulo(phi, 2*pi)
    endif

    extrema(1,3) = minval(phi) ; extrema(2,3) = maxval(phi)

end subroutine get_plasma_extrema

subroutine get_position(ind, pos, input_coords)
    !+ Get position `pos` given indices `ind`
    integer(Int32), dimension(3), intent(in) :: ind
        !+ Indices
    real(Float64), dimension(3), intent(out) :: pos
        !+ Position in [[libfida:beam_grid]] coordinates [cm]
    integer, intent(in), optional            :: input_coords
        !+ Indicates coordinate system of `ind`. Beam grid (0) and cylindrical (2)

    real(Float64), dimension(3) :: pos_temp1, pos_temp2
    integer :: ics, ocs

    if(present(input_coords)) then
        ics = input_coords
    else
        ics = 0
    endif

    if(ics.eq.0) then
        pos(1) = beam_grid%xc(ind(1))
        pos(2) = beam_grid%yc(ind(2))
        pos(3) = beam_grid%zc(ind(3))
    endif

    if(ics.eq.2) then
        pos_temp1(1) = inter_grid%r(ind(1))
        pos_temp1(2) = inter_grid%z(ind(2))
        pos_temp1(3) = inter_grid%phi(ind(3))
        call cyl_to_xyz(pos_temp1,pos)
    endif

end subroutine get_position

subroutine track(rin, vin, tracks, ntrack, los_intersect)
    !+ Computes the path of a neutral through the [[libfida:beam_grid]]
    real(Float64), dimension(3), intent(in)          :: rin
        !+ Initial position of particle
    real(Float64), dimension(3), intent(in)          :: vin
        !+ Initial velocity of particle
    type(ParticleTrack), dimension(:), intent(inout) :: tracks
        !+ Array of [[ParticleTrack]] type
    integer(Int32), intent(out)                      :: ntrack
        !+ Number of cells that a particle crosses
    logical, intent(out), optional                   :: los_intersect
        !+ Indicator whether particle intersects a LOS in [[libfida:spec_chords]]

    integer :: cc, i, j, ii, mind, ncross, id, k, adaptive,  max_cell_splits
    integer, dimension(3) :: ind
    logical :: in_plasma1, in_plasma2, in_plasma_tmp, los_inter
    real(Float64) :: dT, dt1, inv_50, dt2, n_cells, split_tol, inv_param, inv_tol, inv_N, inv_dl, param_sum, param1, param2
    real(Float64), dimension(3) :: dt_arr, dr
    real(Float64), dimension(3) :: vn, inv_vn, vp
    real(Float64), dimension(3) :: ri, ri_tmp, ri_cell
    type(LocalEMFields) :: fields
    type(LocalProfiles) :: plasma1, plasma2
    type(LOSInters) :: inter
    real(Float64), dimension(n_stark) :: lambda
    integer, dimension(3) :: sgn
    integer, dimension(3) :: gdims

    vn = vin ; ri = rin ; sgn = 0 ; ntrack = 0
    adaptive = 0; max_cell_splits = 1; split_tol = 0.0

    los_inter=.False.
    if(.not.present(los_intersect)) then
        los_inter = .True. !avoids computation if not needed
    endif

    if(dot_product(vin,vin).eq.0.0) then
        return
    endif

    gdims(1) = beam_grid%nx
    gdims(2) = beam_grid%ny
    gdims(3) = beam_grid%nz

    !! define actual cell
    call get_indices(ri,ind)
    ri_cell = [beam_grid%xc(ind(1)), &
               beam_grid%yc(ind(2)), &
               beam_grid%zc(ind(3))]

    do i=1,3
        if (vn(i).gt.0.0) sgn(i) = 1
        if (vn(i).lt.0.0) sgn(i) =-1
        if (vn(i).eq.0.0) vn(i)  = 1.0d-3
    enddo

    !! Check adaptive switch, adaptive > 0 activates splitting
    adaptive = inputs%adaptive
    max_cell_splits = inputs%max_cell_splits
    if(adaptive.gt.0) then
        split_tol = inputs%split_tol
        if(split_tol.eq.0.0) then
            adaptive = 0
        else
            inv_tol = 1.0/split_tol
        endif
    endif

    dr = beam_grid%dr*sgn
    inv_vn = 1/vn
    inv_50 = 1.0/50.0
    cc=1
    tracks%time = 0.d0
    tracks%flux = 0.d0
    ncross = 0
    call in_plasma(ri,in_plasma1)
    track_loop: do i=1,beam_grid%ntrack
        if(cc.gt.beam_grid%ntrack) exit track_loop
        dt_arr = abs(( (ri_cell + 0.5*dr) - ri)*inv_vn)
        mind = minloc(dt_arr,1)
        dT = dt_arr(mind)
        ri_tmp = ri + dT*vn

        !! Check if velocity intersects LOS and produces wavelength in the right region
        inter = spec_chords%inter(ind(1),ind(2),ind(3))
        if((.not.los_inter).and.(inter%nchan.ne.0))then
            call get_fields(fields, pos=ri_tmp)
            chan_loop: do j=1,inter%nchan
                id = inter%los_elem(j)%id
                vp = ri_tmp - spec_chords%los(id)%lens
                call doppler_stark(vp, vn, fields, beam_lambda0, lambda)
                los_inter = any((lambda.ge.inputs%lambdamin).and.(lambda.le.inputs%lambdamax))
                if(los_inter) exit chan_loop
            enddo chan_loop
        endif

        call in_plasma(ri_tmp,in_plasma2)
        if(in_plasma1.neqv.in_plasma2) then
            dt1 = 0.0
            track_fine: do ii=1,50
                dt1 = dt1 + dT*inv_50
                ri_tmp = ri + vn*dt1
                call in_plasma(ri_tmp,in_plasma_tmp)
                if(in_plasma2.eqv.in_plasma_tmp) exit track_fine
            enddo track_fine
            tracks(cc)%pos = ri + 0.5*dt1*vn
            tracks(cc+1)%pos = ri + 0.5*(dt1 + dT)*vn
            tracks(cc)%time = dt1
            tracks(cc+1)%time = dT - dt1
            tracks(cc)%ind = ind
            tracks(cc+1)%ind = ind
            cc = cc + 2
            ncross = ncross + 1
        elseif(adaptive.eq.0) then
            tracks(cc)%pos = ri + 0.5*dT*vn
            tracks(cc)%time = dT
            tracks(cc)%ind = ind
            cc = cc + 1
        elseif(adaptive.gt.0) then !! If splitting is activated, calculate n_cells based on gradient between entrance and exit points
            dt2 = 0.0
            call get_plasma(plasma1, ri)
            call get_plasma(plasma2, ri_tmp)

            !! Split cell according to gradient of parameter, parameter is selected according to value of adaptive
            select case (adaptive)
                case(1)
                    param1 = plasma1%dene
                    param2 = plasma2%dene
                case(2)
                    param1 = sum(plasma1%denn(1,:))/size(plasma1%denn(1,:))
                    param2 = sum(plasma2%denn(1,:))/size(plasma2%denn(1,:))
                case(3)
                    param1 = plasma1%denf
                    param2 = plasma2%denf
                case(4)
                    param1 = sum(plasma1%deni(:))/size(plasma1%deni(:))
                    param2 = sum(plasma2%deni(:))/size(plasma2%deni(:))
                case(5)
                    param1 = plasma1%denimp
                    param2 = plasma2%denimp
                case(6)
                    param1 = plasma1%te
                    param2 = plasma2%te
                case(7)
                    param1 = plasma1%ti
                    param2 = plasma2%ti
                case default
                    param1 = plasma1%dene
                    param2 = plasma2%dene
            end select

            param_sum = param1 + param2
            if(param_sum.le.0.0) then
                inv_param = 1.0
            else
                inv_param = 2.0/param_sum
            endif
            inv_dl = 1/(dT*norm2(vn))
            n_cells = ceiling(abs(param1 - param2)*inv_param*inv_dl*inv_tol)
            if(n_cells.gt.max_cell_splits) then
                n_cells = max_cell_splits
            elseif(n_cells.lt.1.0) then
                n_cells = 1.0
            endif
            inv_N = 1.0/n_cells
            split_loop2: do k=1,int(n_cells) !! Split time step and position by n_cells
                dt2 = dt2 + dT*inv_N
                tracks(cc)%pos = ri + 0.5*dt2*vn
                tracks(cc)%time = dT*inv_N
                tracks(cc)%ind = ind
                cc = cc + 1
            enddo split_loop2
        else
            tracks(cc)%pos = ri + 0.5*dT*vn
            tracks(cc)%time = dT
            tracks(cc)%ind = ind
            cc = cc + 1
        endif
        in_plasma1 = in_plasma2

        ri = ri + dT*vn
        ind(mind) = ind(mind) + sgn(mind)
        ri_cell(mind) = ri_cell(mind) + dr(mind)

        if (ind(mind).gt.gdims(mind)) exit track_loop
        if (ind(mind).lt.1) exit track_loop
        if (ncross.ge.2) then
            cc = cc - 1 !dont include last segment
            exit track_loop
        endif
    enddo track_loop
    ntrack = cc-1
    if(present(los_intersect)) then
        los_intersect = los_inter
    endif

end subroutine track

subroutine track_cylindrical(rin, vin, tracks, ntrack, los_intersect)
    !+ Computes the path of a neutral through the [[libfida:pass_grid]]
    real(Float64), dimension(3), intent(in)          :: rin
        !+ Initial position of particle
    real(Float64), dimension(3), intent(in)          :: vin
        !+ Initial velocity of particle
    type(ParticleTrack), dimension(:), intent(inout) :: tracks
        !+ Array of [[ParticleTrack]] type
    integer(Int32), intent(out)                      :: ntrack
        !+ Number of cells that a particle crosses
    logical, intent(out), optional                   :: los_intersect
        !+ Indicator whether particle intersects a LOS in [[libfida:spec_chords]]

    real(Float64), dimension(3,3) :: basis
    real(Float64), dimension(3) :: dt_arr, dr !! Need to make this rzphi
    real(Float64), dimension(3) :: vn, vn_cyl, vp
    real(Float64), dimension(3) :: ri, ri_cyl, ri_tmp
    real(Float64), dimension(3) :: p, nz
    real(Float64), dimension(3) :: v_plane_cyl, v_plane
    real(Float64), dimension(3) :: h_plane_cyl, h_plane
    real(Float64), dimension(3) :: arc_cyl, arc
    real(Float64), dimension(3) :: redge, tedge
    real(Float64), dimension(3) :: redge_cyl, tedge_cyl
    integer, dimension(3) :: sgn
    integer, dimension(3) :: gdims
    integer, dimension(1) :: minpos
    integer, dimension(3) :: ind
    real(Float64) :: dT, dt1, inv_50, dt2
    real(Float64) :: s, c, phi
    real(Float64) :: n_cells, split_tol, inv_param, inv_tol, inv_N, inv_dl, param_sum, param1, param2
    integer :: cc, i, j, ii, mind, ncross, id, k, adaptive, max_cell_splits
    type(LocalEMFields) :: fields
    type(LocalProfiles) :: plasma1, plasma2
    type(LOSInters) :: inter
    real(Float64), dimension(n_stark) :: lambda
    logical :: in_plasma1, in_plasma2, in_plasma_tmp, los_inter
    integer :: ir, iz, iphi

    vn = vin ;  ri = rin ; sgn = 0 ; ntrack = 0
    adaptive = 0; max_cell_splits = 1; split_tol = 0.0

    los_inter=.False.
    if(.not.present(los_intersect)) then
        los_inter = .True. !avoids computation if not needed
    endif

    if(dot_product(vn,vn).eq.0.0) then
        return
    endif

    gdims(1) = pass_grid%nr
    gdims(2) = pass_grid%nz
    gdims(3) = pass_grid%nphi

    phi = atan2(rin(2),rin(1))
    s = sin(phi) ; c = cos(phi)
    vn_cyl(1) = c*vn(1) + s*vn(2)
    vn_cyl(3) = -s*vn(1) + c*vn(2)
    vn_cyl(2) = vn(3)

    do i=1,3
        !! sgn is in R-Z-Phi coordinates
        if (vn_cyl(i).gt.0.d0) then
            sgn(i) = 1
        else if (vn_cyl(i).lt.0.d0) then
            sgn(i) =-1
        end if
    enddo

    dr(1) = pass_grid%dr*sgn(1)
    dr(2) = pass_grid%dz*sgn(2)
    dr(3) = pass_grid%dphi*sgn(3)

    !! Define actual cell
    ri_cyl(1) = sqrt(ri(1)*ri(1) + ri(2)*ri(2))
    ri_cyl(2) = ri(3)
    ri_cyl(3) = atan2(ri(2), ri(1))
    call get_passive_grid_indices(ri_cyl,ind)

    arc_cyl(1) = pass_grid%r(ind(1))
    arc_cyl(2) = pass_grid%z(ind(2))
    arc_cyl(3) = pass_grid%phi(ind(3))
    h_plane_cyl = arc_cyl
    v_plane_cyl = arc_cyl

    !! Define surfaces to intersect
    if(sgn(1).gt.0.d0) arc_cyl(1) = pass_grid%r(ind(1)+1)
    if(sgn(2).gt.0.d0) h_plane_cyl(2) = pass_grid%z(ind(2)+1)
    if(sgn(3).gt.0.d0) v_plane_cyl(3) = pass_grid%phi(ind(3)+1)
    ! Special case of the particle being on the surace handled below
    if((sgn(1).lt.0.d0).and.(arc_cyl(1).eq.ri_cyl(1))) then
        arc_cyl(1) = pass_grid%r(ind(1)-1)
        h_plane_cyl(1) = pass_grid%r(ind(1)-1)
        v_plane_cyl(1) = pass_grid%r(ind(1)-1)
    endif
    if((sgn(2).lt.0.d0).and.(h_plane_cyl(2).eq.ri_cyl(2))) then
        arc_cyl(2) = pass_grid%z(ind(2)-1)
        h_plane_cyl(2) = pass_grid%z(ind(2)-1)
        v_plane_cyl(2) = pass_grid%z(ind(2)-1)
    endif
    if((sgn(3).lt.0.d0).and.(v_plane_cyl(3).eq.ri_cyl(3))) then
        arc_cyl(3) = pass_grid%phi(ind(3)-1)
        h_plane_cyl(3) = pass_grid%phi(ind(3)-1)
        v_plane_cyl(3) = pass_grid%phi(ind(3)-1)
    endif
    call cyl_to_uvw(arc_cyl,arc)
    call cyl_to_uvw(h_plane_cyl,h_plane)
    call cyl_to_uvw(v_plane_cyl,v_plane)

    !! Normal vectors
    nz(1) = 0.d0 ; nz(2) = 0.d0 ; nz(3) = 1.d0
    redge_cyl(1) = v_plane_cyl(1) + pass_grid%dr
    redge_cyl(2) = v_plane_cyl(2)
    redge_cyl(3) = v_plane_cyl(3)
    call cyl_to_uvw(redge_cyl,redge)
    tedge_cyl(1) = v_plane_cyl(1)
    tedge_cyl(2) = v_plane_cyl(2) + pass_grid%dz
    tedge_cyl(3) = v_plane_cyl(3)
    call cyl_to_uvw(tedge_cyl,tedge)
    call plane_basis(v_plane, redge, tedge, basis)

    !! Check passive adaptive switch, adaptive > 0 activates splitting
    adaptive = inputs%adaptive
    max_cell_splits = inputs%max_cell_splits
    if(adaptive.gt.0) then
        split_tol = inputs%split_tol
        if(split_tol.eq.0.0) then
            adaptive = 0
        else
            inv_tol = 1.0/split_tol
        endif
    endif

    !! Track the particle
    inv_50 = 1.0/50.0
    cc=1
    tracks%time = 0.d0
    tracks%flux = 0.d0
    ncross = 0
    call in_plasma(ri,in_plasma1,input_coords=1)
    track_loop: do i=1,pass_grid%ntrack
        if(cc.gt.pass_grid%ntrack) exit track_loop

        call line_cylinder_intersect(ri, vn, arc, p, dt_arr(1))
        call line_plane_intersect(ri, vn, h_plane, nz, p, dt_arr(2))
        call line_plane_intersect(ri, vn, v_plane, -basis(:,3), p, dt_arr(3))

        minpos = minloc(dt_arr, mask=dt_arr.gt.0.d0)
        mind = minpos(1)
        dT = dt_arr(mind)
        ri_tmp = ri + dT*vn

        !! Check if velocity intersects LOS and produces wavelength in the right region
        inter = spec_chords%cyl_inter(ind(1),ind(2),ind(3))
        if((.not.los_inter).and.(inter%nchan.ne.0))then
            call get_fields(fields, pos=ri_tmp, input_coords=1)
            chan_loop: do j=1,inter%nchan
                id = inter%los_elem(j)%id
                vp = ri_tmp - spec_chords%los(id)%lens_uvw
                call doppler_stark(vp, vn, fields, beam_lambda0, lambda)
                los_inter = any((lambda.ge.inputs%lambdamin).and.(lambda.le.inputs%lambdamax))
                if(los_inter) exit chan_loop
            enddo chan_loop
        endif

        call in_plasma(ri_tmp,in_plasma2,input_coords=1)
        if(in_plasma1.neqv.in_plasma2) then
            dt1 = 0.0
            track_fine: do ii=1,50
                dt1 = dt1 + dT*inv_50
                ri_tmp = ri + vn*dt1
                call in_plasma(ri_tmp,in_plasma_tmp,input_coords=1)
                if(in_plasma2.eqv.in_plasma_tmp) exit track_fine
            enddo track_fine
            tracks(cc)%pos = ri + 0.5*dt1*vn
            tracks(cc+1)%pos = ri + 0.5*(dt1 + dT)*vn
            tracks(cc)%time = dt1
            tracks(cc+1)%time = dT - dt1
            tracks(cc)%ind = ind
            tracks(cc+1)%ind = ind
            cc = cc + 2
            ncross = ncross + 1
        elseif(adaptive.eq.0) then
            tracks(cc)%pos = ri + 0.5*dT*vn
            tracks(cc)%time = dT
            tracks(cc)%ind = ind
            cc = cc + 1
        elseif(adaptive.gt.0) then !! If splitting is activated, calculate n_cells based on gradient between entrance and exit points
            dt2 = 0.0
            call get_plasma(plasma1, ri, input_coords=1)
            call get_plasma(plasma2, ri_tmp, input_coords=1)

            !! Split cell according to gradient of parameter, parameter is selected according to value of adaptive
            select case (adaptive)
                case(1)
                    param1 = plasma1%dene
                    param2 = plasma2%dene
                case(2)
                    param1 = sum(plasma1%denn(1,:))/size(plasma1%denn(1,:))
                    param2 = sum(plasma2%denn(1,:))/size(plasma2%denn(1,:))
                case(3)
                    param1 = plasma1%denf
                    param2 = plasma2%denf
                case(4)
                    param1 = sum(plasma1%deni(:))/size(plasma1%deni(:))
                    param2 = sum(plasma2%deni(:))/size(plasma2%deni(:))
                case(5)
                    param1 = plasma1%denimp
                    param2 = plasma2%denimp
                case(6)
                    param1 = plasma1%te
                    param2 = plasma2%te
                case(7)
                    param1 = plasma1%ti
                    param2 = plasma2%ti
                case default
                    param1 = plasma1%dene
                    param2 = plasma2%dene
            end select

            param_sum = param1 + param2
            if(param_sum.le.0.0) then
                inv_param = 1.0
            else
                inv_param = 2.0/param_sum
            endif
            inv_dl = 1/(dT*norm2(vn))
            n_cells = ceiling(abs(param1 - param2)*inv_param*inv_dl*inv_tol)
            if(n_cells.gt.max_cell_splits) then
                n_cells = max_cell_splits
            elseif(n_cells.lt.1.0) then
                n_cells = 1.0
            endif
            inv_N = 1.0/n_cells
            split_loop2: do k=1,int(n_cells) !! Split time step and position by n_cells
                dt2 = dt2 + dT*inv_N
                tracks(cc)%pos = ri + 0.5*dt2*vn
                tracks(cc)%time = dT*inv_N
                tracks(cc)%ind = ind
                cc = cc + 1
            enddo split_loop2
        else
            tracks(cc)%pos = ri + 0.5*dT*vn
            tracks(cc)%time = dT
            tracks(cc)%ind = ind
            cc = cc + 1
        endif
        in_plasma1 = in_plasma2

        ri = ri + dT*vn
        ind(mind) = ind(mind) + sgn(mind)

        if (ind(mind).gt.gdims(mind)) exit track_loop
        if (ind(mind).lt.1) exit track_loop
        if (ncross.ge.2) then
            cc = cc - 1 !dont include last segment
            exit track_loop
        endif
        !! Particle advancement and basis update
        arc_cyl(mind) = arc_cyl(mind) + dr(mind)
        h_plane_cyl(mind) = h_plane_cyl(mind) + dr(mind)
        v_plane_cyl(mind) = v_plane_cyl(mind) + dr(mind)
        call cyl_to_uvw(arc_cyl,arc)
        call cyl_to_uvw(h_plane_cyl,h_plane)
        call cyl_to_uvw(v_plane_cyl,v_plane)
        redge_cyl(1) = v_plane_cyl(1) + pass_grid%dr
        redge_cyl(2) = v_plane_cyl(2)
        redge_cyl(3) = v_plane_cyl(3)
        call cyl_to_uvw(redge_cyl,redge)
        tedge_cyl(1) = v_plane_cyl(1)
        tedge_cyl(2) = v_plane_cyl(2) + pass_grid%dz
        tedge_cyl(3) = v_plane_cyl(3)
        call cyl_to_uvw(tedge_cyl,tedge)
        call plane_basis(v_plane, redge, tedge, basis)
    enddo track_loop
    ntrack = cc-1
    if(present(los_intersect)) then
        los_intersect = los_inter
    endif

end subroutine track_cylindrical

!=============================================================================
!-------------------------Profiles and Fields Routines------------------------
!=============================================================================
subroutine in_plasma(xyz, inp, input_coords, coeffs, uvw_out)
    !+ Indicator subroutine to determine if a position is in a region where
    !+ the plasma parameter and fields are valid/known
    real(Float64), dimension(3), intent(in) :: xyz
        !+ Position in beam coordinates
    logical, intent(out)                    :: inp
        !+ Indicates whether plasma parameters and fields are valid/known
    integer, intent(in), optional           :: input_coords
        !+ Indicates coordinate system of xyz. Beam grid (0), machine (1) and cylindrical (2)
    type(InterpolCoeffs3D), intent(out), optional      :: coeffs
        !+ Interpolation coefficients used in calculation
    real(Float64), dimension(3), intent(out), optional :: uvw_out
        !+ Position in machine coordinates

    real(Float64), dimension(3) :: uvw
    type(InterpolCoeffs3D) :: b
    real(Float64) :: R, W, mask
    real(Float64) :: phi
    integer :: i, j, k, k2, err, ics

    err = 1

    if(present(input_coords)) then
        ics = input_coords
    else
        ics = 0
    endif

    if(ics.eq.0) then
        call xyz_to_uvw(xyz,uvw)
    endif
    if(ics.eq.1) then
        uvw = xyz
    endif
    if(ics.eq.2) then
        call cyl_to_uvw(xyz,uvw)
    endif

    R = sqrt(uvw(1)*uvw(1) + uvw(2)*uvw(2))
    W = uvw(3)
    phi = atan2(uvw(2),uvw(1))
    !! Interpolate mask value
    call interpol_coeff(inter_grid%r, inter_grid%z, inter_grid%phi, R, W, phi, b, err)

    inp = .False.
    if(err.eq.0) then
        i = b%i
        j = b%j
        k = b%k
        if(inter_grid%nphi .eq. 1) then
            k2 = min(k+1,inter_grid%nphi)
        else
            k2 = k+1
        endif
        mask = b%b111*equil%mask(i,j,k)     + b%b112*equil%mask(i,j,k2) +   &
               b%b121*equil%mask(i,j+1,k)   + b%b122*equil%mask(i,j+1,k2) + &
               b%b211*equil%mask(i+1,j,k)   + b%b212*equil%mask(i+1,j,k2) + &
               b%b221*equil%mask(i+1,j+1,k) + b%b222*equil%mask(i+1,j+1,k2)

        if((mask.ge.0.5).and.(err.eq.0)) then
            inp = .True.
        endif
    endif


    if(present(coeffs)) coeffs = b
    if(present(uvw_out)) uvw_out = uvw

end subroutine in_plasma

subroutine get_plasma(plasma, pos, ind, input_coords, output_coords)
    !+ Gets plasma parameters at position `pos` or [[libfida:beam_grid]] indices `ind`
    type(LocalProfiles), intent(out)                   :: plasma
        !+ Plasma parameters at `pos`/`ind`
    real(Float64), dimension(3), intent(in), optional  :: pos
        !+ Position in beam grid coordinates
    integer(Int32), dimension(3), intent(in), optional :: ind
        !+ [[libfida:beam_grid]] indices
    integer(Int32), intent(in), optional               :: input_coords
        !+ Indicates coordinate system of inputs. Beam grid (0), machine (1) and cylindrical (2)
    integer(Int32), intent(in), optional               :: output_coords
        !+ Indicates coordinate system of outputs. Beam grid (0), machine (1) and cylindrical (2)

    logical :: inp
    type(InterpolCoeffs3D) :: coeffs
    real(Float64), dimension(3) :: xyz, uvw, cyl, vrot_uvw
    real(Float64) :: phi, s, c
    integer :: i, j, k, k2, ics, ocs

    plasma%in_plasma = .False.

    if(present(input_coords)) then
        ics = input_coords
    else
        ics = 0
    endif
    if(present(output_coords)) then
        ocs = output_coords
    else
        ocs = 0
    endif

    if(present(ind)) then
        if(ics.eq.0) then
            call get_position(ind,xyz)
        endif
        if(ics.eq.2) then
            call get_position(ind,xyz,input_coords=2)
        endif
    endif

    if(present(pos)) then
        if(ics.eq.0) then
            xyz = pos
            call xyz_to_uvw(xyz, uvw)
        endif
        if(ics.eq.1) then
            uvw = pos
            call uvw_to_xyz(uvw, xyz)
        endif
    endif

    call in_plasma(xyz,inp,0,coeffs)
    if(inp) then
        phi = atan2(uvw(2),uvw(1))
        i = coeffs%i
        j = coeffs%j
        k = coeffs%k
        if(inter_grid%nphi .eq. 1) then
            k2 = min(k+1,inter_grid%nphi)
        else
            k2 = k+1
        endif

        plasma = coeffs%b111*equil%plasma(i,j,k)    + coeffs%b121*equil%plasma(i,j+1,k) +   &
                 coeffs%b112*equil%plasma(i,j,k2)   + coeffs%b122*equil%plasma(i,j+1,k2) +  &
                 coeffs%b211*equil%plasma(i+1,j,k)  + coeffs%b221*equil%plasma(i+1,j+1,k) + &
                 coeffs%b212*equil%plasma(i+1,j,k2) + coeffs%b222*equil%plasma(i+1,j+1,k2)

        s = sin(phi) ; c = cos(phi)
        vrot_uvw(1) = plasma%vr*c - plasma%vt*s
        vrot_uvw(2) = plasma%vr*s + plasma%vt*c
        vrot_uvw(3) = plasma%vz
        if(ocs.eq.0) then
            plasma%vrot_uvw = vrot_uvw
            plasma%vrot = matmul(beam_grid%inv_basis,vrot_uvw)
            plasma%pos = xyz
        endif
        if(ocs.eq.1) then
            plasma%vrot = vrot_uvw
            plasma%pos = uvw
        endif
        plasma%uvw = uvw
        plasma%in_plasma = .True.
        plasma%b = coeffs
    endif

end subroutine get_plasma

pure subroutine calc_perp_vectors(b, a, c)
  !+ Calculates normalized vectors that are perpendicular to b
  !+ such that `a` x `c` = `b_norm`
  real(Float64), dimension(3), intent(in)  :: b
  real(Float64), dimension(3), intent(out) :: a
  real(Float64), dimension(3), intent(out) :: c

  real(Float64), dimension(3) :: bnorm

  bnorm=b/norm2(b)

  if (abs(bnorm(3)).eq.1) then
      a=[1.d0,0.d0,0.d0]
      c=[0.d0,1.d0,0.d0]
  else
      if (bnorm(3).eq.0.) then
          a=[0.d0,0.d0,1.d0]
          c=[bnorm(2),-bnorm(1), 0.d0]/sqrt(bnorm(1)**2+bnorm(2)**2)
      else
          a=[bnorm(2),-bnorm(1),0.d0]/sqrt(bnorm(1)**2+bnorm(2)**2)
          c=-[ a(2) , -a(1) , (a(1)*bnorm(2)-a(2)*bnorm(1))/bnorm(3) ]
          c=c/norm2(c)
          if(bnorm(3).lt.0.0) then
              c=-c
          endif
      endif
  endif

end subroutine calc_perp_vectors

subroutine get_fields(fields, pos, ind, input_coords, output_coords)
    !+ Gets electro-magnetic fields at position `pos` or [[libfida:beam_grid]] indices `ind`
    type(LocalEMFields),intent(out)                    :: fields
        !+ Electro-magnetic fields at `pos`/`ind`
    real(Float64), dimension(3), intent(in), optional  :: pos
        !+ Position in beam grid coordinates
    integer(Int32), dimension(3), intent(in), optional :: ind
        !+ [[libfida:beam_grid]] indices
    integer(Int32), intent(in), optional               :: input_coords
        !+ Indicates coordinate system of inputs. Beam grid (0), machine (1) and cylindrical (2)
    integer(Int32), intent(in), optional               :: output_coords
        !+ Indicates coordinate system of outputs. Beam grid (0), machine (1) and cylindrical (2)

    logical :: inp
    real(Float64), dimension(3) :: xyz, uvw
    real(Float64), dimension(3) :: uvw_bfield, uvw_efield
    real(Float64), dimension(3) :: xyz_bfield, xyz_efield
    real(Float64) :: phi, s, c
    type(InterpolCoeffs3D) :: coeffs
    integer :: i, j, k, k2, mc, ocs, ics

    fields%in_plasma = .False.

    if(present(input_coords)) then
        ics = input_coords
    else
        ics = 0
    endif

    if(present(output_coords)) then
        ocs = output_coords
    else
        ocs = 0
    endif

    if(present(ind)) call get_position(ind,xyz)

    if(present(pos)) then
        if(ics.eq.0) then
            xyz = pos
            call xyz_to_uvw(xyz, uvw)
        endif
        if(ics.eq.1) then
            uvw = pos
            call uvw_to_xyz(uvw, xyz)
        endif
    endif

    call in_plasma(xyz,inp,0,coeffs)
    if(inp) then
        i = coeffs%i
        j = coeffs%j
        k = coeffs%k
        if(inter_grid%nphi .eq. 1) then
            k2 = min(k+1,inter_grid%nphi)
        else
            k2 = k+1
        endif

        fields = coeffs%b111*equil%fields(i,j,k)    + coeffs%b121*equil%fields(i,j+1,k) +   &
                 coeffs%b112*equil%fields(i,j,k2)   + coeffs%b122*equil%fields(i,j+1,k2) +  &
                 coeffs%b211*equil%fields(i+1,j,k)  + coeffs%b221*equil%fields(i+1,j+1,k) + &
                 coeffs%b212*equil%fields(i+1,j,k2) + coeffs%b222*equil%fields(i+1,j+1,k2)

        phi = atan2(uvw(2),uvw(1))
        s = sin(phi) ; c = cos(phi)

        !Convert cylindrical coordinates to uvw
        uvw_bfield(1) = c*fields%br - s*fields%bt
        uvw_bfield(2) = s*fields%br + c*fields%bt
        uvw_bfield(3) = fields%bz
        uvw_efield(1) = c*fields%er - s*fields%et
        uvw_efield(2) = s*fields%er + c*fields%et
        uvw_efield(3) = fields%ez

        if(ocs.eq.0) then
            !Represent fields in beam grid coordinates
            xyz_bfield = matmul(beam_grid%inv_basis,uvw_bfield)
            xyz_efield = matmul(beam_grid%inv_basis,uvw_efield)
            fields%pos = xyz
        endif
        if(ocs.eq.1) then
            xyz_bfield = uvw_bfield
            xyz_efield = uvw_efield
            fields%pos = uvw
        endif

        !Calculate field directions and magnitudes
        fields%b_abs = norm2(xyz_bfield)
        fields%e_abs = norm2(xyz_efield)
        if(fields%b_abs.gt.0.d0) fields%b_norm = xyz_bfield/fields%b_abs
        if(fields%e_abs.gt.0.d0) fields%e_norm = xyz_efield/fields%e_abs

        call calc_perp_vectors(fields%b_norm,fields%a_norm,fields%c_norm)

        fields%uvw = uvw
        fields%in_plasma = .True.
        fields%coords = ocs
        fields%b = coeffs
    endif

end subroutine get_fields

subroutine get_distribution(fbeam, denf, pos, ind, coeffs)
    !+ Gets Guiding Center distribution at position `pos` or [[libfida:beam_grid]] indices `ind`
    real(Float64), dimension(:,:), intent(out)         :: fbeam
        !+ Guiding Center Fast-ion distribution at `pos`/`ind`: F(E,p)
    real(Float64), intent(out)                         :: denf
        !+ Guiding Center Fast-ion density at `pos`/`ind` [fast-ions/cm^3]
    real(Float64), dimension(3), intent(in), optional  :: pos
        !+ Position in beam grid coordinates
    integer(Int32), dimension(3), intent(in), optional :: ind
        !+ [[libfida:beam_grid]] indices
    type(InterpolCoeffs3D), intent(in), optional       :: coeffs
        !+ Precomputed interpolation coefficients

    real(Float64), dimension(3) :: xyz, uvw
    real(Float64) :: R, Z, Phi
    integer :: err

    if(present(coeffs)) then
        call interpol(fbm%r, fbm%z, fbm%phi, fbm%f, R, Z, Phi, fbeam, err, coeffs)
        call interpol(fbm%r, fbm%z, fbm%phi, fbm%denf, R, Z, Phi, denf, err, coeffs)
    else
        if(present(ind)) call get_position(ind,xyz)
        if(present(pos)) xyz = pos

        !! Convert to machine coordinates
        call xyz_to_uvw(xyz,uvw)
        R = sqrt(uvw(1)*uvw(1) + uvw(2)*uvw(2))
        Z = uvw(3)
        Phi = atan2(uvw(2),uvw(1))

        call interpol(fbm%r, fbm%z, fbm%phi, fbm%f, R, Z, Phi, fbeam, err)
        call interpol(fbm%r, fbm%z, fbm%phi, fbm%denf, R, Z, Phi, denf, err)
    endif

end subroutine get_distribution

subroutine get_ep_denf(energy, pitch, denf, pos, ind, coeffs)
    !+ Get guiding center fast-ion density at given energy and pitch
    !+ at position `pos` or [[libfida:beam_grid]] indices `ind`
    real(Float64), intent(in)                          :: energy
        !+ Energy [keV]
    real(Float64), intent(in)                          :: pitch
        !+ Pitch
    real(Float64), intent(out)                         :: denf
        !+ Fast-ion density [fast-ions/(cm^3*dE*dp)]
    real(Float64), dimension(3), intent(in), optional  :: pos
        !+ Position in beam grid coordinates
    integer(Int32), dimension(3), intent(in), optional :: ind
        !+ [[libfida:beam_grid]] indices
    type(InterpolCoeffs3D), intent(in), optional       :: coeffs
        !+ Precomputed interpolation coefficients

    real(Float64), dimension(3) :: xyz, uvw
    real(Float64), dimension(fbm%nenergy,fbm%npitch)  :: fbeam
    integer(Int32), dimension(2) :: epi
    real(Float64) :: R, Phi, Z
    real(Float64) :: dE, dp
    integer :: err

    epi(1) = minloc(abs(fbm%energy - energy),1)
    epi(2) = minloc(abs(fbm%pitch - pitch),1)
    dE = abs(fbm%energy(epi(1)) - energy)
    dp = abs(fbm%pitch(epi(2)) - pitch)

    if((dE.le.fbm%dE).and.(dp.le.fbm%dp)) then
        if(present(coeffs)) then
            call interpol(inter_grid%r, inter_grid%z, inter_grid%phi, fbm%f, R, Z, Phi, fbeam, err, coeffs)
        else
            if(present(ind)) call get_position(ind,xyz)
            if(present(pos)) xyz = pos

            !! Convert to machine coordinates
            call xyz_to_uvw(xyz,uvw)
            R = sqrt(uvw(1)*uvw(1) + uvw(2)*uvw(2))
            Z = uvw(3)
            Phi = atan2(uvw(2),uvw(1))

            call interpol(inter_grid%r, inter_grid%z, inter_grid%phi, fbm%f, R, Z, Phi, fbeam, err)
        endif
        denf = fbeam(epi(1),epi(2))
    else
        denf = 0.0
    endif

end subroutine get_ep_denf

!=============================================================================
!--------------------------Result Storage Routines----------------------------
!=============================================================================
subroutine init_reservoir(res)
    !+ Initialize reservoir
    type(NeutralParticleReservoir), intent(inout) :: res
        !+ Neutral Particle Reservoir
    !integer, intent(in) :: k
        !+ Size of reservoir

    real(Float64) :: randomu(2)

    if (res%n.lt.1) then
        res%n = 1
        res%k = reservoir_size
        call randu(randomu)
        res%W = exp(log(randomu(1))/res%k)
        res%i = res%k + floor(log(randomu(2))/log(1.0-res%W)) + 1
    endif

end subroutine init_reservoir

subroutine update_reservoir(res, vn, wght)
    !+ Store neutrals in reservoir
    type(NeutralParticleReservoir), intent(inout) :: res
        !+ [[libfida:beam_grid]] indices
    real(Float64), dimension(3), intent(in)       :: vn
        !+ Neutral particle velocity [cm/s]
    real(Float64), intent(in)                     :: wght
        !+ Neutral weight

    integer, dimension(1) :: randomi
    real(Float64), dimension(2) :: randomu

    !print*, res%n, res%k, res%i
    !$OMP CRITICAL(update_reservoir_1)
    if (res%n.le.res%k) then
        if (res%n.le.0) then
            call init_reservoir(res)
        endif
        res%R(res%n)%w = wght
        res%R(res%n)%v = vn
    elseif (res%n.eq.res%i) then
        call randind(res%k, randomi)
        res%R(randomi(1)) = NeutralParticle(wght, vn)
        call randu(randomu)
        res%W = res%W*exp(log(randomu(1))/res%k)
        res%i = res%i + floor(log(randomu(2))/log(1.0-res%W)) + 1
    endif
    res%n = res%n + 1
    !$OMP END CRITICAL(update_reservoir_1)

end subroutine update_reservoir

subroutine merge_reservoirs(res1,res2)
    !+ Merge reservoir 2 into reservoir 1
    type(NeutralParticleReservoir), intent(inout) :: res1
        !+ Reservoir 1
    type(NeutralParticleReservoir), intent(inout) :: res2
        !+ Reservoir 2

    integer :: i
    type(NeutralParticleReservoir) :: res3
    real(Float64) :: p, randomu(1)

    if(res1%n.le.0) then
        call init_reservoir(res1) !res2%k)
    endif
    if(res2%n.le.0) return

    !! If neither reservoir is full just fill res1 from res2 normally
    if ((res1%n.le.res1%k).and.(res2%n.le.res2%k)) then
        do i=1,res2%n
            call update_reservoir(res1, res2%R(i)%v, res2%R(i)%w)
        enddo
        res2 = res1
        return
    endif

    !! Switch res1 with res2 if res1 is not full but res2 is
    if (res1%n.le.res1%k) then
        res3 = res1
        res1 = res2
        res2 = res3
    endif

    !! if res2 is not full then fill res1 from res2 normally
    if (res2%n.le.res2%k) then
        do i=1,res2%n
            call update_reservoir(res1, res2%R(i)%v, res2%R(i)%w)
        enddo
    else !! if both are full select particles from res2 according to # particles seen
        p = float(res1%n)/(res1%n + res2%n)
        do i=1,res1%k
            call randu(randomu)
            if(randomu(1).gt.p) then
                res1%R(i) = res2%R(i)
            endif
        enddo
    endif
    !! set res1 to res2
    res2 = res1

end subroutine merge_reservoirs

subroutine init_neutral_population(pop)
    !+ Initialize [[NeutralPopulation]]
    type(NeutralPopulation), intent(inout) :: pop
        !+ Neutral Population to initialize

    integer :: i, j, k, ic, ind(3)

    if(.not.allocated(pop%dens)) then
        allocate(pop%dens(nlevs,beam_grid%nx,beam_grid%ny,beam_grid%nz))
    endif
    if(.not.allocated(pop%res)) then
        allocate(pop%res(beam_grid%nx,beam_grid%ny,beam_grid%nz))
    endif
    pop%dens = 0.d0
    pop%res%n = 0

end subroutine init_neutral_population

subroutine free_neutral_population(pop)
    !+ Initialize [[NeutralPopulation]]
    type(NeutralPopulation), intent(inout) :: pop
        !+ Neutral Population to initialize

    integer :: i, j, k, ic, ind(3)

    if(allocated(pop%dens)) then
        deallocate(pop%dens)
    endif
    if(allocated(pop%res)) then
        deallocate(pop%res)
    endif

end subroutine free_neutral_population

subroutine update_neutrals(pop, ind, vn, dens)
    !+Update [NeutralPopulation]] `pop` at `ind`
    type(NeutralPopulation), intent(inout)   :: pop
        !+ Neutral Population
    integer(Int32), dimension(3), intent(in) :: ind
        !+ [[libfida:beam_grid]] indices
    real(Float64), dimension(3), intent(in)  :: vn
        !+ Neutral particle velocity [cm/s]
    real(Float64), dimension(:), intent(in)  :: dens
        !+ Neutral density [neutrals/cm^3]

    integer :: i, j, k

    i = ind(1) ; j = ind(2) ; k = ind(3)

    !$OMP CRITICAL(update_neutrals_1)
    pop%dens(:,i,j,k) = pop%dens(:,i,j,k) + dens
    !$OMP END CRITICAL(update_neutrals_1)

    call update_reservoir(pop%res(i,j,k), vn, sum(dens))

end subroutine update_neutrals

subroutine store_neutrals(ind, pos, vn, neut_type, dens)
    !Store neutrals in [[libfida:neut]] at indices `ind`
    integer(Int32), dimension(3), intent(in) :: ind
        !+ [[libfida:beam_grid]] indices
    real(Float64), dimension(3), intent(in)  :: pos
        !+ Neutral particle position [cm]
    real(Float64), dimension(3), intent(in)  :: vn
        !+ Neutral particle velocity [cm/s]
    integer, intent(in)                      :: neut_type
        !+ Neutral type
    real(Float64), dimension(nlevs), intent(in)  :: dens
        !+ Neutral density [neutrals/cm^3]

    select case (neut_type)
        case (nbif_type)
            call update_neutrals(neut%full, ind, vn, dens)
        case (nbih_type)
            call update_neutrals(neut%half, ind, vn, dens)
        case (nbit_type)
            call update_neutrals(neut%third, ind, vn, dens)
        case (dcx_type)
            call update_neutrals(neut%dcx, ind, vn, dens)
        case (halo_type)
            call update_neutrals(neut%halo, ind, vn, dens)
        case default
            if(inputs%verbose.ge.0) then
                write(*,'("STORE_NEUTRALS: Unknown neutral type: ",i2)') neut_type
            endif
            stop
    end select

end subroutine store_neutrals

subroutine merge_neutral_populations(pop1, pop2)
    !+ Merge [[NeutralPopulation]]s
    type(NeutralPopulation), intent(inout) :: pop1
        !+ Neutral Population target
    type(NeutralPopulation), intent(inout) :: pop2
        !+ Neutral Population to merge into target

    integer :: i, j, k, ic, ind(3)

    pop1%dens = pop1%dens + pop2%dens
    do ic=1,beam_grid%ngrid
        call ind2sub(beam_grid%dims,ic,ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)
        call merge_reservoirs(pop1%res(i,j,k), pop2%res(i,j,k))
    enddo

end subroutine merge_neutral_populations

#ifdef _MPI

subroutine parallel_merge_reservoirs(res)
    !+ Merge neutral particle reservoirs across MPI processes
    type(NeutralParticleReservoir), intent(inout) :: res
        !+ Neutral Particle Reservoir

    integer :: proc, num_procs, i, j, k
    integer(Int32), dimension(:), allocatable :: nums, inds
    real(Float64), dimension(:), allocatable :: Ws
    real(Float64), dimension(:,:), allocatable :: wghts
    real(Float64), dimension(:,:,:), allocatable :: vels
    type(NeutralParticleReservoir) :: res1, res2

    proc = my_rank() + 1
    num_procs = num_ranks()

    k = res%k

    allocate(nums(num_procs), inds(num_procs), Ws(num_procs))
    nums = 0
    inds = 0
    Ws = 0.d0
    nums(proc) = res%n
    inds(proc) = res%i
    Ws(proc) = res%W

    call parallel_sum(nums)
    call parallel_sum(inds)
    call parallel_sum(Ws)

    allocate(vels(3,reservoir_size,num_procs))
    allocate(wghts(reservoir_size,num_procs))
    vels = 0.d0; wghts = 0.d0

    do i=1,min(k,res%n)
        vels(:,i,proc) = res%R(i)%v
        wghts(i,proc) = res%R(i)%w
    enddo
    call parallel_sum(vels)
    call parallel_sum(wghts)

    call init_reservoir(res1)
    res1%n = nums(1)
    res1%i = inds(1)
    res1%W = Ws(1)
    do i=1,min(k,res1%n)
        res1%R(i) = NeutralParticle(wghts(i,1),vels(:,i,1))
    enddo
    call init_reservoir(res2)
    do i=2,num_procs
        res2%n = nums(i)
        res2%i = inds(i)
        res2%W = Ws(i)
        do j=1,min(k,res2%n)
            res2%R(j) = NeutralParticle(wghts(j,i),vels(:,j,i))
        enddo
        call merge_reservoirs(res1,res2)
    enddo
    deallocate(vels,wghts,nums,inds,Ws)

    res = res1

end subroutine parallel_merge_reservoirs

subroutine parallel_merge_populations(pop)
    !+Update neutral population
    type(NeutralPopulation), intent(inout)   :: pop
        !+ Neutral Population

    integer :: i, j, k, ic, ind(3)

    call parallel_sum(pop%dens)
    do ic=1,beam_grid%ngrid
        call ind2sub(beam_grid%dims,ic,ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)
        call parallel_merge_reservoirs(pop%res(i,j,k))
    enddo

end subroutine parallel_merge_populations

#endif

subroutine store_births(ind, neut_type, dflux)
    !+ Store birth particles/density in [[libfida:birth]]
    integer(Int32), dimension(3), intent(in) :: ind
        !+ [[libfida:beam_grid]] indices
    integer(Int32), intent(in)               :: neut_type
        !+ Neutral type
    real(Float64), intent(in)                :: dflux
        !+ Deposited flux

    !$OMP ATOMIC UPDATE
    birth%dens( neut_type,ind(1),ind(2),ind(3))= &
     birth%dens(neut_type,ind(1),ind(2),ind(3)) + dflux
    !$OMP END ATOMIC
end subroutine store_births

subroutine store_npa(det, ri, rf, vn, flux, orbit_class, passive)
    !+ Store NPA particles in [[libfida:npa]]
    integer, intent(in)                     :: det
        !+ Detector/Channel Number
    real(Float64), dimension(3), intent(in) :: ri
        !+ Birth position in beam coordinates [cm]
    real(Float64), dimension(3), intent(in) :: rf
        !+ Detector position in beam coordinates [cm]
    real(Float64), dimension(3), intent(in) :: vn
        !+ Particle velocity [cm/s]
    real(Float64), intent(in)               :: flux
        !+ Neutral flux [neutrals/s]
    integer, intent(in), optional           :: orbit_class
        !+ Orbit class ID
    logical, intent(in), optional           :: passive
        !+ Indicates whether npa particle is passive

    integer :: iclass, oclass
    type(LocalEMFields) :: fields
    real(Float64), dimension(3) :: uvw_ri, uvw_rf,vn_norm
    real(Float64) :: energy, pitch, dE
    integer(Int32), dimension(1) :: ienergy
    type(NPAParticle), dimension(:), allocatable :: parts
    logical :: pas = .False.

    if(present(orbit_class)) then
        iclass = min(orbit_class,particles%nclass)
        oclass = orbit_class
    else
        iclass = 1
        oclass = 1
    endif

    if(present(passive)) pas = passive

    ! Convert to machine coordinates
    call xyz_to_uvw(ri,uvw_ri)
    call xyz_to_uvw(rf,uvw_rf)

    ! Calculate energy
    energy = beam_mass*v2_to_E_per_amu*dot_product(vn,vn)
    if(pas) then
        dE = pnpa%energy(2)-pnpa%energy(1)
    else
        dE = npa%energy(2)-npa%energy(1)
    endif


    ! Calculate pitch if distribution actually uses pitch
    if(inputs%dist_type.le.2) then
        call get_fields(fields, pos = ri)
        vn_norm = vn/norm2(vn)
        pitch = dot_product(fields%b_norm,vn_norm)
    else
        pitch = 0.d0
    endif

    if(pas) then
        !$OMP CRITICAL(store_npa_1)
        pnpa%npart = pnpa%npart + 1
        if(pnpa%npart.gt.pnpa%nmax) then
            pnpa%nmax = int(pnpa%nmax*2)
            allocate(parts(pnpa%nmax))
            parts(1:(pnpa%npart-1)) = pnpa%part
            deallocate(pnpa%part)
            call move_alloc(parts, pnpa%part)
        endif
        pnpa%part(pnpa%npart)%detector = det
        pnpa%part(pnpa%npart)%class = oclass
        pnpa%part(pnpa%npart)%xi = uvw_ri(1)
        pnpa%part(pnpa%npart)%yi = uvw_ri(2)
        pnpa%part(pnpa%npart)%zi = uvw_ri(3)
        pnpa%part(pnpa%npart)%xf = uvw_rf(1)
        pnpa%part(pnpa%npart)%yf = uvw_rf(2)
        pnpa%part(pnpa%npart)%zf = uvw_rf(3)
        pnpa%part(pnpa%npart)%energy = energy
        pnpa%part(pnpa%npart)%pitch = pitch
        pnpa%part(pnpa%npart)%weight = flux
        ienergy = minloc(abs(pnpa%energy - energy))
        pnpa%flux(ienergy(1),det,iclass) = &
            pnpa%flux(ienergy(1),det,iclass) + flux/dE
        !$OMP END CRITICAL(store_npa_1)
    else
        !$OMP CRITICAL(store_npa_2)
        npa%npart = npa%npart + 1
        if(npa%npart.gt.npa%nmax) then
            npa%nmax = int(npa%nmax*2)
            allocate(parts(npa%nmax))
            parts(1:(npa%npart-1)) = npa%part
            deallocate(npa%part)
            call move_alloc(parts, npa%part)
        endif
        npa%part(npa%npart)%detector = det
        npa%part(npa%npart)%class = oclass
        npa%part(npa%npart)%xi = uvw_ri(1)
        npa%part(npa%npart)%yi = uvw_ri(2)
        npa%part(npa%npart)%zi = uvw_ri(3)
        npa%part(npa%npart)%xf = uvw_rf(1)
        npa%part(npa%npart)%yf = uvw_rf(2)
        npa%part(npa%npart)%zf = uvw_rf(3)
        npa%part(npa%npart)%energy = energy
        npa%part(npa%npart)%pitch = pitch
        npa%part(npa%npart)%weight = flux
        ienergy = minloc(abs(npa%energy - energy))
        npa%flux(ienergy(1),det,iclass) = &
            npa%flux(ienergy(1),det,iclass) + flux/dE
        !$OMP END CRITICAL(store_npa_2)
    endif

end subroutine store_npa

!=============================================================================
!--------------------------Atomic Physics Routines----------------------------
!=============================================================================
subroutine bb_cx_rates(denn, vn, vi, rates)
    !+ Get beam-beam neutralization/cx rates
    real(Float64), dimension(nlevs), intent(in)  :: denn
        !+ Neutral density [cm^-3]
    real(Float64), dimension(3),     intent(in)  :: vn
        !+ Neutral velocity [cm/s]
    real(Float64), dimension(3),     intent(in)  :: vi
        !+ Ion velocity [cm/s]
    real(Float64), dimension(nlevs), intent(out) :: rates
        !+ Reaction rates [1/s]

    real(Float64), dimension(nlevs,nlevs) :: neut  !!rate coeff
    real(Float64) :: eb !! relative Energy
    type(InterpolCoeffs1D) :: c
    real(Float64) :: dlogE, logEmin, logeb
    real(Float64) :: vrel !! relative velocity
    integer :: ebi, neb, err

    !Eeff
    vrel=norm2(vi-vn)
    eb=v2_to_E_per_amu*vrel**2  ! [kev/amu]
    logeb = log10(eb)
    logEmin = tables%H_H_cx_cross%logemin
    dlogE = tables%H_H_cx_cross%dlogE
    neb = tables%H_H_cx_cross%nenergy
    call interpol_coeff(logEmin,dlogE,neb,logeb,c,err)
    ebi = c%i
    if(err.eq.1) then
        if(inputs%verbose.ge.0) then
            write(*,'(a)') "BB_CX_RATES: Eb out of range of H_H_cx table. Using nearest energy value."
            write(*,'("eb = ",ES10.3," [keV]")') eb
        endif
        if(ebi.lt.1) then
            ebi=1
            c%b1=1.0 ; c%b2=0.0
        else
            ebi=neb-1
            c%b1=0.0 ; c%b2=1.0
        endif
    endif

    neut(:,:) = (c%b1*tables%H_H_cx_cross%log_cross(:,:,ebi) + &
                 c%b2*tables%H_H_cx_cross%log_cross(:,:,ebi+1))

    where (neut.lt.tables%H_H_cx_cross%minlog_cross)
        neut = 0.d0
    elsewhere
        neut = exp(neut*log_10)
    end where

    rates=matmul(neut,denn)*vrel

end subroutine bb_cx_rates

subroutine bt_cx_rates(plasma, denn, an, vi, rates)
    !+ Get beam-target neutralization/cx rates
    type(LocalProfiles), intent(in)              :: plasma
        !+ Plasma parameters (for neutral temperature and vrot)
    real(Float64), dimension(nlevs), intent(in)  :: denn
        !+ Neutral density [cm^-3]
    real(Float64), intent(in)                    :: an
        !+ Neutral mass [amu]
    real(Float64), dimension(3),     intent(in)  :: vi
        !+ Ion velocity [cm/s]
    real(Float64), dimension(nlevs), intent(out) :: rates
        !+ Reaction rates [1/s]

    real(Float64) :: logEmin, dlogE, logeb_amu, eb_amu
    real(Float64) :: logTmin, dlogT, logti_amu, vrel
    integer :: neb, nt
    type(InterpolCoeffs2D) :: c
    real(Float64) :: b11, b12, b21, b22
    real(Float64), dimension(nlevs,nlevs) :: H_H_rate
    integer :: ebi, tii, n, err_status

    H_H_rate = 0.d0

    vrel=norm2(vi-plasma%vrot)
    eb_amu=v2_to_E_per_amu*vrel**2  ! [kev/amu]

    logeb_amu = log10(eb_amu)
    logti_amu = log10(plasma%ti/an)

    !!H_H
    err_status = 1
    logEmin = tables%H_H_cx_rate%logemin
    logTmin = tables%H_H_cx_rate%logtmin
    dlogE = tables%H_H_cx_rate%dlogE
    dlogT = tables%H_H_cx_rate%dlogT
    neb = tables%H_H_cx_rate%nenergy
    nt = tables%H_H_cx_rate%ntemp
    call interpol_coeff(logEmin, dlogE, neb, logTmin, dlogT, nt, &
                        logeb_amu, logti_amu, c, err_status)
    ebi = c%i
    tii = c%j
    b11 = c%b11
    b12 = c%b12
    b21 = c%b21
    b22 = c%b22
    if(err_status.eq.1) then
        if(inputs%verbose.ge.0) then
            write(*,'(a)') "BT_CX_RATES: Eb or Ti out of range of H_H_CX table. Setting H_H_CX rates to zero"
            write(*,'("eb/amu = ",ES10.3," [keV/amu]")') eb_amu
            write(*,'("ti/amu = ",ES10.3," [keV/amu]")') plasma%ti/an
        endif
        rates = 0.0
        return
    endif

    H_H_rate = (b11*tables%H_H_cx_rate%log_rate(:,:,ebi,tii)   + &
                b12*tables%H_H_cx_rate%log_rate(:,:,ebi,tii+1) + &
                b21*tables%H_H_cx_rate%log_rate(:,:,ebi+1,tii) + &
                b22*tables%H_H_cx_rate%log_rate(:,:,ebi+1,tii+1))

    where (H_H_rate.lt.tables%H_H_cx_rate%minlog_rate)
        H_H_rate = 0.d0
    elsewhere
        H_H_rate = exp(H_H_rate*log_10) !cm^3/s
    end where

    rates=matmul(H_H_rate,denn) !1/s

end subroutine bt_cx_rates

subroutine get_dd_rate(plasma, eb, rate, branch)
    !+ Gets d(d,p)T rate for a beam with energy `eb` interacting with a target plasma
    type(LocalProfiles), intent(in) :: plasma
        !+ Plasma Paramters
    real(Float64), intent(in)       :: eb
        !+ Beam energy [keV]
    real(Float64), intent(out)      :: rate
        !+ Neutron reaction rate [1/s]
    integer, intent(in), optional   :: branch
        !+ Indicates 1 for cfpd rate and 2 for neutron rate

    integer :: err_status, neb, nt, ebi, tii, is, ib
    real(Float64) :: dlogE, dlogT, logEmin, logTmin
    real(Float64) :: logeb, logti, lograte
    type(InterpolCoeffs2D) :: c
    real(Float64) :: b11, b12, b21, b22

    if(present(branch)) then
        ib = branch
    else
        ib = 2
    endif

    logeb = log10(eb)
    logti = log10(plasma%ti)

    !!D_D
    err_status = 1
    logEmin = tables%D_D%logemin
    logTmin = tables%D_D%logtmin
    dlogE = tables%D_D%dlogE
    dlogT = tables%D_D%dlogT
    neb = tables%D_D%nenergy
    nt = tables%D_D%ntemp
    call interpol_coeff(logEmin, dlogE, neb, logTmin, dlogT, nt, &
                        logeb, logti, c, err_status)
    ebi = c%i
    tii = c%j
    b11 = c%b11
    b12 = c%b12
    b21 = c%b21
    b22 = c%b22
    if(err_status.eq.1) then
        if(inputs%verbose.ge.0) then
            write(*,'(a)') "GET_DD_RATE: Eb or Ti out of range of D_D table. Setting D_D rates to zero"
        endif
        rate = 0.d0
        return
    endif

    lograte = (b11*tables%D_D%log_rate(ebi,tii,ib)   + &
               b12*tables%D_D%log_rate(ebi,tii+1,ib) + &
               b21*tables%D_D%log_rate(ebi+1,tii,ib) + &
               b22*tables%D_D%log_rate(ebi+1,tii+1,ib))

    if (lograte.lt.tables%D_D%minlog_rate) then
        rate = 0.d0
    else
        rate = 0.d0
        do is=1,n_thermal
            if(thermal_mass(is).eq.H2_amu) then
                rate = rate + plasma%deni(is) * exp(lograte*log_10)
            endif
        enddo
    endif

end subroutine get_dd_rate

subroutine get_ddpt_anisotropy(plasma, v1, v3, kappa)
    !+ Gets d(d,p)T anisotropy defecit/enhancement factor for a beam interacting with a target plasma
    type(LocalProfiles), intent(in)         :: plasma
        !+ Plasma Paramters
    real(Float64), dimension(3), intent(in) :: v1
        !+ Beam velocity [cm/s]
    real(Float64), dimension(3), intent(in) :: v3
        !+ Charged Fusion Product velocity [cm/s]
    real(Float64), intent(out)              :: kappa
        !+ Anisotropy factor
    !+ Reference: Eq. (1) and (3) of NIM A236 (1985) 380

    real(Float64), dimension(3,12) :: abc
    real(Float64), dimension(13)   :: bhcor !!! 13?
    real(Float64), dimension(12)   :: e, a, b, c
    real(Float64), dimension(3)    :: vcm, v3cm, vrel
    type(InterpolCoeffs1D) :: c1D

    real(Float64) :: ai, bi, ci, b1, b2, b11, b12, b21, b22, cos_phi, sin_phi
    real(Float64) :: eb, e1com, vnet_square, cos_theta, k, KE, Q, mp, k0, JMeV
    integer :: ei, i, err_status

    !! Calculate effective beam energy
    vrel = v1-plasma%vrot ![cm/s]
    vnet_square=dot_product(vrel, vrel) ![(cm/s)**2]
    eb = v2_to_E_per_amu*fbm%A*vnet_square ![keV]

    !!Calculate anisotropy enhancement/deficit factor
    JMeV = 1.60218d-13 ! Conversion factor from MeV to Joules
    mp = H1_amu*mass_u  ![kg]
    Q = 4.04*JMeV ![J]

    vcm = 0.5*(v1+plasma%vrot) ![cm/s]
    KE = 0.5*mp*vnet_square*1.d-4  ! [J] C-O-M kinetic energy
    k0 = norm2(vcm) * sqrt(2*mp/(3*(Q+KE)))*100.d0 ![(cm/s)**2]
    if ((norm2(vcm)*norm2(v3)).gt.0.d0) then
        cos_phi = dot_product(vcm, v3) / (norm2(vcm)*norm2(v3))
        sin_phi = sin(acos(cos_phi))

        if (abs(k0*sin_phi).le.1) then
            cos_theta = cos_phi*sqrt(1-(k0*sin_phi)**2) - k0*sin_phi**2
        else
            cos_theta = 0.d0
        endif

    else
        cos_theta = 0.d0
    endif

    !Brown-Jarmie coefficients in Table I with prepended isotropic low-energy extrapolated point
    e = [0.0,19.944,29.935,39.927,49.922,59.917,69.914,79.912,89.911,99.909,109.909,116.909]
    a = [0.0,0.0208,0.0886,0.1882,0.3215,0.4636,0.6055,0.7528,0.8976,1.041,1.166,1.243]
    b = [0.0,0.0023,0.0184,0.0659,0.076,0.168,0.251,0.250,0.378,0.367,0.521,0.55]
    c = [0.0,0.0,0.0,0.0,0.048,0.021,0.039,0.162,0.142,0.276,0.203,0.34]
    e(1) = 10.0
    a(1) = 0.00903775/(4*pi)
    abc(1,:) = a ; abc(2,:) = b ; abc(3,:) = c

    !Correction factor to make it consistent with Bosch & Hale
    bhcor=[1.0,1.02614,.98497,1.006,1.015,1.012,1.0115,1.021,1.012,1.012,1.016,1.018,1.012]
    do i=1,12
        abc(:,i) = abc(:,i)*bhcor(i)
    enddo

    e1com=0.5d0*eb
    call interpol_coeff(e, e1com, c1D, err_status)

    ei = c1D%i
    b1 = c1D%b1
    b2 = c1D%b2

    ai = b1*abc(1,ei) + b2*abc(1,ei+1)
    bi = b1*abc(2,ei) + b2*abc(2,ei+1)
    ci = b1*abc(3,ei) + b2*abc(3,ei+1)

    kappa = (ai + bi*cos_theta**2 + ci*cos_theta**4) / (ai+bi/3.d0+ci/5.d0)

end subroutine get_ddpt_anisotropy

subroutine get_pgyro(fields,E3,E1,pitch,plasma,v3_xyz,pgyro,gam0)
    !+ Returns fraction of gyroangles that can produce a reaction with
    !+ given inputs
    type(LocalEMFields), intent(in) :: fields
        !+ Electromagneticfields in beam coordinates
    real(Float64), intent(in)       :: E3
        !+ E3 charged fusion product energy [keV]
    real(Float64), intent(in) :: E1
        !+ E1 fast-ion energy [keV]
    real(Float64), intent(in) :: pitch
        !+ pitch  fast ion pitch relative to the field
    type(LocalProfiles), intent(in)         :: plasma
        !+ Plasma Paramters in beam coordinates
    real(Float64), dimension(3), intent(in) :: v3_xyz
        !+ Charged fusion product velocity in beam coorindates
    real(Float64), intent(out) :: pgyro
        !+ pgyro   DeltaE_3*\partial\gam/\partial E_3/pi
    real(Float64), intent(out) :: gam0
        !+ Gyro angle of fast ion [rad]

    real(Float64), dimension(3) :: a_hat, vrot
    real(Float64) :: JMeV,JkeV,mp,Q,norm_v3,norm_v1,vpar,vperp,vb,va,E3max,E3min
    real(Float64) :: phip,cosphip,cosphib,cosphia,rhs,gammaplus,gammaminus
    real(Float64) :: eps,bracket,ccminus,ccplus,bbminus,bbplus,v3minus,v3plus
    real(Float64) :: E3minus,E3plus,v3plusmag,v3minusmag,v3maxmag,v3minmag
    real(Float64) :: DeltaE3

    pgyro = 0.d0
    gam0 = 0.d0

    ! Preliminaries [SI units]
    JMeV = 1.60218d-13 ! Conversion factor from MeV to Joules
    JkeV = 1.60218d-16 ! Conversion factor from keV to Joules
    mp = H1_amu*mass_u  ![kg]
    Q = 4.04*JMeV ![J]
    norm_v3 = sqrt(E3/(H1_amu*v2_to_E_per_amu)) / 100 ![m/s]
    norm_v1 = sqrt(E1/(beam_mass*v2_to_E_per_amu)) / 100 ![m/s]
    vpar = norm_v1*pitch ![m/s]
    vperp = norm_v1*sqrt(1-pitch**2) ![m/s]
    vrot = plasma%vrot/100.d0 ![m/s]
    a_hat = cross_product(fields%b_norm, v3_xyz) / norm2(v3_xyz) !(b,a,c)

    !! First, check E3 limits are valid for the given inputs
    !! Assumes cos(gamma) = +/- 1 and vrot = 0
    cosphip = dot_product(v3_xyz, fields%b_norm) / norm2(v3_xyz)
    phip = acos(cosphip)

    !! Get plasma rotation components in (b,a,c) coords
    if (all(vrot.eq.0.d0)) then
        vb = 0.d0
        va = 0.d0
    else
        cosphib = dot_product(vrot, fields%b_norm) / norm2(vrot)
        vb = norm2(vrot)*cosphib ![m/s]
        cosphia = dot_product(vrot, a_hat) / norm2(vrot)
        va = norm2(vrot)*cosphia ![m/s]
    endif

    ! LHS coefficient teeny (step #0)
    bracket = vperp*(sin(phip)-2*va/norm_v3)
    if (abs(bracket).lt.1.d-5) then !'Case 0'
        return
    endif

    ! Find E3 limits for these parameters (step #1)
    ccminus = 1.5*Q/mp + 0.5*norm_v1**2 - 2*vpar*vb + 0.5*norm2(vrot) - vperp*va
    ccplus = 1.5*Q/mp + 0.5*norm_v1**2 - 2*vpar*vb + 0.5*norm2(vrot) + vperp*va
    bbminus = -vperp*sin(phip) - (vpar+vb)*cos(phip) - va*sin(phip)
    bbplus = vperp*sin(phip) - (vpar+vb)*cos(phip) - va*sin(phip)
    v3minus = 0.5*(-bbminus + sqrt(bbminus**2+4*ccminus))
    v3plus = 0.5*(-bbplus + sqrt(bbplus**2+4*ccplus))
    E3minus = 0.5*mp*v3minus**2 / JkeV
    E3plus = 0.5*mp*v3plus**2 / JkeV
    E3min = min(E3plus,E3minus) ![keV]
    E3max = max(E3plus,E3minus) ![keV]

    ! Now E3plus and E3minus are charged fusion product energies at edge of bin
    ! 'Case 1'
    E3plus = E3 + 0.5*ctable%dE
    E3minus = E3 - 0.5*ctable%dE

    ! Bin misses gamma curve altogether  (step #2)
    if ((E3plus.le.E3min) .or. (E3minus.ge.E3max)) then !'Case 2'
        return
    endif

    ! Get gammaplus and gammaminus  (step #3)
    v3plusmag = sqrt(2*E3plus*JkeV/mp)
    v3minusmag = sqrt(2*E3minus*JkeV/mp)
    v3maxmag = sqrt(2*E3plus*JkeV/mp)
    v3minmag = sqrt(2*E3minus*JkeV/mp)

    if (E3plus.gt.E3max) then !'Case 3'
        norm_v3 = v3maxmag
        rhs = norm_v3 - 1.5*Q/(mp*norm_v3) - (vpar+vb)*cos(phip) - va*sin(phip) &
              - 0.5*(norm_v1**2 + norm2(vrot))/norm_v3 + 2*vpar*vb/norm_v3
        if (abs(rhs).gt.abs(bracket)) then
            if (rhs*bracket.gt.0.) then
                gammaplus = 0.
            else
                gammaplus = pi
            endif
        else
            gammaplus = acos(rhs/bracket)
        endif
    else
        norm_v3 = v3plusmag
        rhs = norm_v3 - 1.5*Q/(mp*norm_v3) - (vpar+vb)*cos(phip) - va*sin(phip) &
              - 0.5*(norm_v1**2 + norm2(vrot))/norm_v3 + 2*vpar*vb/norm_v3
        if (abs(rhs).gt.abs(bracket)) then
            if (rhs*bracket.gt.0.) then
                gammaplus = 0.
            else
                gammaplus = pi
            endif
        else
            gammaplus = acos(rhs/bracket)
        endif
    endif

    if (E3minus.lt.E3min) then !'Case 3'
        norm_v3 = v3minmag
        rhs = norm_v3 - 1.5*Q/(mp*norm_v3) - (vpar+vb)*cos(phip) - va*sin(phip) &
               - 0.5*(norm_v1**2 + norm2(vrot))/norm_v3 + 2*vpar*vb/norm_v3
        if (abs(rhs).gt.abs(bracket)) then
            if (rhs*bracket.gt.0.) then
                gammaminus = 0.
            else
                gammaminus = pi
            endif
        else
            gammaminus = acos(rhs/bracket)
        endif
    else
        norm_v3 = v3minusmag
        rhs = norm_v3 - 1.5*Q/(mp*norm_v3) - (vpar+vb)*cos(phip) - va*sin(phip) &
              - 0.5*(norm_v1**2 + norm2(vrot))/norm_v3 + 2*vpar*vb/norm_v3
        if (abs(rhs).gt.abs(bracket)) then
            if (rhs*bracket.gt.0.) then
                gammaminus = 0.
            else
                gammaminus = pi
            endif
        else
            gammaminus = acos(rhs/bracket)
        endif
    endif

    gam0 = (gammaplus + gammaminus) / 2.d0 ! return mean gyroangle
    !!! should there be modifications in gamma in if statements below

    if ((E3plus.ge.E3max) .and. (E3minus.le.E3min)) then !'Case 4a'
        pgyro = 1.d0
        return
    else if ((E3plus.lt.E3max) .and. (E3minus.gt.E3min)) then !'Case 4b'
        pgyro = abs(gammaplus-gammaminus)/pi
        return
    else if((E3plus.ge.E3max) .and. (E3minus.gt.E3min)) then
        !'Case 4c'
        if (gammaplus.gt.gammaminus) then
            pgyro = (pi-gammaminus)/pi
            return
        else
            pgyro = gammaminus/pi
            return
        endif
    else if ((E3plus.lt.E3max) .and. (E3minus.le.E3min)) then
        !'Case 4d'
        if (gammaminus.gt.gammaplus) then
            pgyro = (pi-gammaplus)/pi
            return
        else
            pgyro = gammaplus/pi
            return
        endif
    endif


endsubroutine get_pgyro

subroutine neutral_cx_rate(denn, res, v_ion, rates)
    !+ Get probability of a thermal ion charge exchanging with neutral
    !+ population within cell
    real(Float64), dimension(nlevs), intent(in)  :: denn
        !+ Neutral density
    type(NeutralParticleReservoir), intent(in)   :: res
        !+ Neutral Particle Reservoir
    real(Float64), dimension(3), intent(in)      :: v_ion
        !+ Ion velocity [cm/s]
    real(Float64), dimension(nlevs), intent(out) :: rates
        !+ Total charge exchange rate [1/s]

    integer :: i
    real(Float64), dimension(nlevs) :: rates_i
    real(Float64) :: vn(3), totw

    !! Average CX rate over Reservoir
    rates = 0.d0
    if(res%n.eq.0) return

    totw = sum(res%R%w)
    if(totw.eq.0.d0) return

    do i=1, min(res%n,res%k)
        vn = res%R(i)%v
        call bb_cx_rates(denn,vn,v_ion,rates_i)
        rates = rates + res%R(i)%w*rates_i
    enddo
    rates = rates/totw

end subroutine neutral_cx_rate

subroutine get_total_cx_rate(ind, pos, v_ion, types, rate_tot)
    !+ Get probability of a thermal ion charge exchanging with `types` neutrals
    integer(Int32), dimension(3), intent(in)     :: ind
        !+ [[libfida:beam_grid]] indices
    real(Float64), dimension(3), intent(in)      :: pos
        !+ Interaction position in beam grid coordinates
    real(Float64), dimension(3), intent(in)      :: v_ion
        !+ Ion velocity [cm/s]
    integer(Int32), dimension(:), intent(in)     :: types
        !+ Neutral types
    real(Float64), dimension(nlevs), intent(out) :: rate_tot
        !+ Total charge exchange rate [1/s]

    integer :: n, i, j, k, ii
    type(LocalProfiles) :: plasma
    type(NeutralParticleReservoir) :: res
    real(Float64), dimension(nlevs) :: rates, denn
    real(Float64), dimension(3) :: vn

    n = size(types)

    rate_tot = 0.0
    i = ind(1) ; j = ind(2) ; k = ind(3)
    do ii=1,n
        select case (types(ii))
            case (nbif_type)
                denn = neut%full%dens(:,i,j,k)
                res = neut%full%res(i,j,k)
            case (nbih_type)
                denn = neut%half%dens(:,i,j,k)
                res = neut%half%res(i,j,k)
            case (nbit_type)
                denn = neut%third%dens(:,i,j,k)
                res = neut%third%res(i,j,k)
            case (dcx_type)
                denn = neut%dcx%dens(:,i,j,k)
                res = neut%dcx%res(i,j,k)
            case (halo_type)
                denn = neut%halo%dens(:,i,j,k)
                res = neut%halo%res(i,j,k)
        end select
        call neutral_cx_rate(denn, res, v_ion, rates)
        rate_tot = rate_tot + rates
    enddo

end subroutine get_total_cx_rate

subroutine get_rate_matrix(plasma, ab, eb, rmat)
    !+ Gets rate matrix for use in [[libfida:colrad]]
    type(LocalProfiles), intent(in)                    :: plasma
        !+ Plasma parameters
    real(Float64), intent(in)                          :: ab
        !+ "Beam" ion mass [amu]
    real(Float64), intent(in)                          :: eb
        !+ "Beam" ion energy [keV]
    real(Float64), dimension(nlevs,nlevs), intent(out) :: rmat
        !+ Rate matrix

    real(Float64) :: ai
    real(Float64) :: logEmin, dlogE, logeb, logeb_amu
    real(Float64) :: logTmin, dlogT, logti, logti_amu, logte
    integer :: neb, nt, i
    type(InterpolCoeffs2D) :: c
    real(Float64) :: b11, b12, b21, b22, dene, deni(max_species), denimp
    real(Float64), dimension(nlevs,nlevs) :: H_H_pop_i, H_H_pop, H_e_pop, H_Aq_pop
    real(Float64), dimension(nlevs) :: H_H_depop_i, H_H_depop, H_e_depop, H_Aq_depop
    integer :: ebi, tii, tei, n, err_status

    H_H_pop = 0.d0
    H_e_pop = 0.d0
    H_Aq_pop = 0.d0
    H_H_depop = 0.d0
    H_e_depop = 0.d0
    H_Aq_depop = 0.d0

    deni = plasma%deni
    dene = plasma%dene
    denimp = plasma%denimp
    logeb_amu = log10(eb/ab)
    logti = log10(plasma%ti)
    logte = log10(plasma%te)

    !!H_H
    err_status = 1
    logEmin = tables%H_H%logemin
    logTmin = tables%H_H%logtmin
    dlogE = tables%H_H%dlogE
    dlogT = tables%H_H%dlogT
    neb = tables%H_H%nenergy
    nt = tables%H_H%ntemp
    do i=1, n_thermal
        H_H_pop_i = 0.d0
        H_H_depop_i = 0.d0
        logti_amu = log10(plasma%ti/thermal_mass(i))
        call interpol_coeff(logEmin, dlogE, neb, logTmin, dlogT, nt, &
                            logeb_amu, logti_amu, c, err_status)
        ebi = c%i
        tii = c%j
        b11 = c%b11
        b12 = c%b12
        b21 = c%b21
        b22 = c%b22
        if(err_status.eq.1) then
            if(inputs%verbose.ge.0) then
                write(*,'(a)') "GET_RATE_MATRIX: Eb or Ti out of range of H_H table. Setting H_H rates to zero"
                write(*,'("eb/amu = ",ES10.3," [keV/amu]")') eb/ab
                write(*,'("ti/amu = ",ES10.3," [keV/amu]")') plasma%ti/thermal_mass(i)
            endif
        else
            H_H_pop_i = (b11*tables%H_H%log_pop(:,:,ebi,tii)   + &
                         b12*tables%H_H%log_pop(:,:,ebi,tii+1) + &
                         b21*tables%H_H%log_pop(:,:,ebi+1,tii) + &
                         b22*tables%H_H%log_pop(:,:,ebi+1,tii+1))
            where (H_H_pop_i.lt.tables%H_H%minlog_pop)
                H_H_pop_i = 0.d0
            elsewhere
                H_H_pop_i = deni(i) * exp(H_H_pop_i*log_10)
            end where

            H_H_depop_i = (b11*tables%H_H%log_depop(:,ebi,tii)   + &
                           b12*tables%H_H%log_depop(:,ebi,tii+1) + &
                           b21*tables%H_H%log_depop(:,ebi+1,tii) + &
                           b22*tables%H_H%log_depop(:,ebi+1,tii+1))
            where (H_H_depop_i.lt.tables%H_H%minlog_depop)
                H_H_depop_i = 0.d0
            elsewhere
                H_H_depop_i = deni(i) * exp(H_H_depop_i*log_10)
            end where
        endif
        H_H_pop = H_H_pop + H_H_pop_i
        H_H_depop = H_H_depop + H_H_depop_i
    enddo

    !!H_e
    err_status = 1
    logEmin = tables%H_e%logemin
    logTmin = tables%H_e%logtmin
    dlogE = tables%H_e%dlogE
    dlogT = tables%H_e%dlogT
    neb = tables%H_e%nenergy
    nt = tables%H_e%ntemp
    call interpol_coeff(logEmin, dlogE, neb, logTmin, dlogT, nt, &
                        logeb_amu, logte, c, err_status)
    ebi = c%i
    tei = c%j
    b11 = c%b11
    b12 = c%b12
    b21 = c%b21
    b22 = c%b22
    if(err_status.eq.1) then
        if(inputs%verbose.ge.0) then
            write(*,'(a)') "GET_RATE_MATRIX: Eb or Te out of range of H_e table. Setting H_e rates to zero"
            write(*,'("eb/amu = ",ES10.3," [keV/amu]")') eb/ab
            write(*,'("te = ",ES10.3," [keV]")') plasma%te
        endif
    else
        H_e_pop = (b11*tables%H_e%log_pop(:,:,ebi,tei)   + &
                   b12*tables%H_e%log_pop(:,:,ebi,tei+1) + &
                   b21*tables%H_e%log_pop(:,:,ebi+1,tei) + &
                   b22*tables%H_e%log_pop(:,:,ebi+1,tei+1))
        where (H_e_pop.lt.tables%H_e%minlog_pop)
            H_e_pop = 0.d0
        elsewhere
            H_e_pop = dene * exp(H_e_pop*log_10)
        end where

        H_e_depop = (b11*tables%H_e%log_depop(:,ebi,tei)   + &
                     b12*tables%H_e%log_depop(:,ebi,tei+1) + &
                     b21*tables%H_e%log_depop(:,ebi+1,tei) + &
                     b22*tables%H_e%log_depop(:,ebi+1,tei+1))

        where (H_e_depop.lt.tables%H_e%minlog_depop)
            H_e_depop = 0.d0
        elsewhere
            H_e_depop = dene * exp(H_e_depop*log_10)
        end where
    endif

    !!H_Aq
    err_status = 1
    logEmin = tables%H_Aq%logemin
    logTmin = tables%H_Aq%logtmin
    dlogE = tables%H_Aq%dlogE
    dlogT = tables%H_Aq%dlogT
    neb = tables%H_Aq%nenergy
    nt = tables%H_Aq%ntemp
    call interpol_coeff(logEmin, dlogE, neb, logTmin, dlogT, nt, &
                        logeb_amu, logti, c, err_status)
    ebi = c%i
    tii = c%j
    b11 = c%b11
    b12 = c%b12
    b21 = c%b21
    b22 = c%b22
    if(err_status.eq.1) then
        if(inputs%verbose.ge.0) then
            write(*,'(a)') "GET_RATE_MATRIX: Eb or Ti out of range of H_Aq table. Setting H_Aq rates to zero"
            write(*,'("eb/amu = ",ES10.3," [keV/amu]")') eb
            write(*,'("ti = ",ES10.3," [keV]")') plasma%ti
        endif
    else
        H_Aq_pop = (b11*tables%H_Aq%log_pop(:,:,ebi,tii)   + &
                    b12*tables%H_Aq%log_pop(:,:,ebi,tii+1) + &
                    b21*tables%H_Aq%log_pop(:,:,ebi+1,tii) + &
                    b22*tables%H_Aq%log_pop(:,:,ebi+1,tii+1))
        where (H_Aq_pop.lt.tables%H_Aq%minlog_pop)
            H_Aq_pop = 0.d0
        elsewhere
            H_Aq_pop = denimp * exp(H_Aq_pop*log_10)
        end where
        H_Aq_depop = (b11*tables%H_Aq%log_depop(:,ebi,tii)   + &
                      b12*tables%H_Aq%log_depop(:,ebi,tii+1) + &
                      b21*tables%H_Aq%log_depop(:,ebi+1,tii) + &
                      b22*tables%H_Aq%log_depop(:,ebi+1,tii+1))

        where (H_Aq_depop.lt.tables%H_Aq%minlog_depop)
            H_Aq_depop = 0.d0
        elsewhere
            H_Aq_depop = denimp * exp(H_Aq_depop*log_10)
        end where
    endif

    rmat = tables%einstein + H_H_pop + H_e_pop + H_Aq_pop
    do n=1,nlevs
        rmat(n,n) = -sum(tables%einstein(:,n)) - H_H_depop(n) - H_e_depop(n) - H_Aq_depop(n)
    enddo

end subroutine get_rate_matrix

subroutine colrad(plasma,ab,vn,dt,states,dens,photons)
    !+ Evolve density of states in time `dt` via collisional radiative model
    type(LocalProfiles), intent(in)              :: plasma
        !+ Plasma parameters
    real(Float64), intent(in)                    :: ab
        !+ Ion/Neutral mass [amu]
    real(Float64), dimension(:), intent(in)      :: vn
        !+ Neutral velocitiy [cm/s]
    real(Float64), intent(in)                    :: dt
        !+ Time interval [s]
    real(Float64), dimension(:), intent(inout)   :: states
        !+ Density of states
    real(Float64), dimension(nlevs), intent(out) :: dens
        !+ Density of neutrals
    real(Float64), intent(out)                   :: photons
        !+ Emitted photons

    real(Float64), dimension(nlevs,nlevs) :: matrix  !! Matrix
    real(Float64) :: vnet_square    !! net velocity of neutrals squared
    real(Float64) :: eb             !! Energy of the fast neutral

    real(Float64), dimension(nlevs,nlevs) :: eigvec
    real(Float64), dimension(nlevs) :: eigval, coef
    real(Float64), dimension(nlevs) :: exp_eigval_dt
    real(Float64) :: iflux !!Initial total flux
    integer :: n

    photons=0.d0
    dens=0.d0

    iflux=sum(states)

    if(.not.plasma%in_plasma) then
        dens = states*dt
        return
    endif

    vnet_square=dot_product(vn-plasma%vrot,vn-plasma%vrot)  ![cm/s]
    eb = v2_to_E_per_amu*ab*vnet_square ![kev]
    call get_rate_matrix(plasma, ab, eb, matrix)

    call eigen(nlevs,matrix, eigvec, eigval)
    call linsolve(eigvec,states,coef) !coeffs determined from states at t=0
    exp_eigval_dt = exp(eigval*dt)   ! to improve speed (used twice)
    do n=1,nlevs
        if(eigval(n).eq.0.0) eigval(n)=eigval(n)+1 !protect against dividing by zero
    enddo

    states = matmul(eigvec, coef * exp_eigval_dt)  ![neutrals/cm^3/s]!
    dens   = matmul(eigvec,coef*(exp_eigval_dt-1.d0)/eigval)

    where (states.lt.0)
        states = 0.d0
    endwhere

    where (dens.lt.0)
        dens = 0.d0
    endwhere

    photons=dens(initial_state)*tables%einstein(final_state,initial_state) !! - [Ph/(s*cm^3)] - !!

end subroutine colrad

subroutine attenuate(ri, rf, vi, states, dstep_in)
    !+ Attenuate `states` along a trajectory
    real(Float64), dimension(3), intent(in)        :: ri
        !+ Initial position in beam grid coordinates
    real(Float64), dimension(3), intent(in)        :: rf
        !+ Final position in beam grid coordinates
    real(Float64), dimension(3), intent(in)        :: vi
        !+ Initial velocity of neutral
    real(Float64), dimension(nlevs), intent(inout) :: states
        !+ Density of states
    real(Float64), intent(in), optional            :: dstep_in
        !+ Step length [cm]

    type(LocalProfiles) :: plasma
    real(Float64) :: photons, vabs, dt, dstep, dis,max_dis
    real(Float64), dimension(3) :: r0
    real(Float64), dimension(nlevs) :: dens
    logical :: inp
    integer :: ncross

    if(present(dstep_in)) then
        dstep=dstep_in
    else
        dstep = sqrt(inter_grid%da) !cm
    endif

    max_dis = norm2(rf-ri)

    vabs = norm2(vi)
    dt = dstep/vabs

    call get_plasma(plasma,pos=ri)
    r0 = ri
    dis = 0.d0
    ncross = 0
    inp = plasma%in_plasma
    do while (dis.le.max_dis)
        call colrad(plasma,beam_mass,vi,dt,states,dens,photons)
        r0 = r0 + vi*dt
        dis = dis + dstep
        call get_plasma(plasma,pos=r0)
        if(inp.neqv.plasma%in_plasma) then
            ncross = ncross + 1
            inp = plasma%in_plasma
        endif
    enddo

    if(ncross.gt.1) states = 0.0

end subroutine attenuate

subroutine doppler_stark(vecp, vi, fields, lambda0, lambda)
    !+ Calculates doppler shift and stark split wavelengths
    real(Float64), dimension(3), intent(in)        :: vecp
        !+ Vector directing towards optical head
    real(Float64), dimension(3), intent(in)        :: vi
        !+ Particle velocity
    type(LocalEMFields), intent(in)                :: fields
        !+ Electro-magnetic fields
    real(Float64), intent(in)                      :: lambda0
        !+ Reference wavelength [nm]
    real(Float64), dimension(n_stark), intent(out) :: lambda
        !+ Wavelengths [nm]
    real(Float64), dimension(n_stark) :: wavel
    real(Float64), dimension(3) :: vp, vn
    real(Float64), dimension(3) :: bfield, efield
    real(Float64) :: E, lambda_shifted, B, l0, h, q0, q1, m
    !! vector directing towards the optical head
    vp=vecp/norm2(vecp)

    !! Calculate Doppler shift
    vn=vi*0.01d0 ! [m/s]
    lambda_shifted = lambda0*(1.d0 + dot_product(vn,vp)/c0)

    !! Calculate Stark Splitting
    !! Calculate E-field
    bfield = fields%b_norm*fields%b_abs        ! planck constant in SI units
    efield = fields%e_norm*fields%e_abs
    efield(1) = efield(1) +  vn(2)*bfield(3) - vn(3)*bfield(2)
    efield(2) = efield(2) - (vn(1)*bfield(3) - vn(3)*bfield(1))
    efield(3) = efield(3) +  vn(1)*bfield(2) - vn(2)*bfield(1)

    E = norm2(efield)
    B = norm2(bfield)

    m = 9.109384d-31
        ! mass of electron [kg]
    h = 6.62607004d-34
        ! planck constant in SI units
    l0 = lambda0*1d-9
        ! reference wavelength [m]
    ! stark-zeeman corrections to energy of n=2 states are -q0, 0 and q0
    q0 = sqrt((e0*h*B/(4*pi*m))**2 + (3*a_0*e0*E)**2)
    ! stark-zeeman corrections to energy of n=3 states are -q1, -0.5*q1, 0, 0.5*q1, and q1/2
    q1 = sqrt(4*(e0*h*B/(4*pi*m))**2 + 9*(3*a_0*e0*E)**2)
    ! wavelengths calculated from h*c0/lambda =  E_i - E_j for transition from i to j energies
    ! order is small wavelengths to large wavelengths

    if(n_stark.eq.15) then
        wavel(1)  = 2*c0*h*l0/(2*c0*h+(-2*q0*(-1)+q1*(2))*l0)
        wavel(2)  = 2*c0*h*l0/(2*c0*h+(-2*q0*(0)+q1*(2))*l0)
        wavel(3) = 2*c0*h*l0/(2*c0*h+(-2*q0*(-1)+q1*(1))*l0)
        wavel(4) = 2*c0*h*l0/(2*c0*h+(-2*q0*(1)+q1*(2))*l0)
        wavel(5) = 2*c0*h*l0/(2*c0*h+(-2*q0*(0)+q1*(1))*l0)
        wavel(6) = 2*c0*h*l0/(2*c0*h+(-2*q0*(-1)+q1*(0))*l0)
        wavel(7) = 2*c0*h*l0/(2*c0*h+(-2*q0*(1)+q1*(1))*l0)
        wavel(8) = 2*c0*h*l0/(2*c0*h+(-2*q0*(0)+q1*(0))*l0)
        wavel(9) = 2*c0*h*l0/(2*c0*h+(-2*q0*(-1)+q1*(-1))*l0)
        wavel(10) = 2*c0*h*l0/(2*c0*h+(-2*q0*(1)+q1*(0))*l0)
        wavel(11) = 2*c0*h*l0/(2*c0*h+(-2*q0*(0)+q1*(-1))*l0)
        wavel(12) = 2*c0*h*l0/(2*c0*h+(-2*q0*(-1)+q1*(-2))*l0)
        wavel(13) = 2*c0*h*l0/(2*c0*h+(-2*q0*(1)+q1*(-1))*l0)
        wavel(14) = 2*c0*h*l0/(2*c0*h+(-2*q0*(0)+q1*(-2))*l0)
        wavel(15) = 2*c0*h*l0/(2*c0*h+(-2*q0*(1)+q1*(-2))*l0)
    else if(n_stark.eq.3) then
        wavel(1) = c0*h*l0/(c0*h+q0*(-1)*l0)
        wavel(2) = c0*h*l0/(c0*h+q0*(0)*l0)
        wavel(3) = c0*h*l0/(c0*h+q0*(1)*l0)
    endif
    !lambda =  lambda_shifted + E * stark_wavel ![nm] <-old calculation for pure stark effect
    wavel = (wavel-l0)*10**9
    !! Stark Splitting
    lambda = lambda_shifted + wavel
end

subroutine spectrum(vecp, vi, fields, lambda0, sigma_pi, photons, dlength, lambda, intensity, stokes)
    !+ Calculates doppler shift, stark splitting, and intensities
    real(Float64), dimension(3), intent(in)        :: vecp
        !+ Vector directing towards optical head
    real(Float64), dimension(3), intent(in)        :: vi
        !+ Particle velocity
    type(LocalEMFields), intent(in)                :: fields
        !+ Electro-magnetic fields
    real(Float64), intent(in)                      :: lambda0
        !+ Reference wavelength [nm]
    real(Float64), intent(in)                      :: sigma_pi
        !+ Sigma-pi ratio
    real(Float64), intent(in)                      :: photons
        !+ Photon density from [[libfida:colrad]]
    real(Float64), intent(in)                      :: dlength
        !+ LOS intersection length with [[libfida:beam_grid]] cell particle is in
    real(Float64), dimension(n_stark), intent(out) :: lambda
        !+ Wavelengths [nm]
    real(Float64), dimension(n_stark), intent(out) :: intensity
    real(Float64), dimension(n_stark,4), intent(out) :: stokes
        !+ Spectra intensities [Ph/(s cm^2 starkline)]
    real(Float64) :: m, h, normfactor
    real(Float64), dimension(3) :: vp, vn
    real(Float64), dimension(3) :: bfield, efield
    real(Float64) :: E, B, cos_los_Efield, cos_los_Bfield, lambda_shifted, q0, q1, l0, szratio
    integer, dimension(n_stark) :: stark_sign
    real(Float64), dimension(n_stark) :: wavel
    real(Float64), dimension(n_stark) :: circularity
    stark_sign = +1*stark_sigma - 1*stark_pi
    !! vector directing towards the optical head
    vp=vecp/norm2(vecp)

    ! Calculate Doppler shift
    vn=vi*0.01d0 ! [m/s]
    lambda_shifted = lambda0*(1.d0 + dot_product(vn,vp)/c0)

    !! Calculate Stark-Zeeman Splitting
    ! Calculate E-field
    bfield = fields%b_norm*fields%b_abs
    efield = fields%e_norm*fields%e_abs
    efield(1) = efield(1) +  vn(2)*bfield(3) - vn(3)*bfield(2)
    efield(2) = efield(2) - (vn(1)*bfield(3) - vn(3)*bfield(1))
    efield(3) = efield(3) +  vn(1)*bfield(2) - vn(2)*bfield(1)
    E = norm2(efield)
    B = norm2(bfield)
    !Stark-Zeeman Splitting
    h=6.62607004d-34
        ! planck constant in SI units
    m = 9.109384d-31
        ! mass of electron [kg]
    l0 = lambda0*1d-9
        ! reference wavelength [m]
    ! stark-zeeman corrections to energy of n=2 states are -q0, 0 and q0
    q0 = sqrt((e0*h*B/(4*pi*m))**2 + (3*a_0*e0*E)**2)
    ! stark-zeeman corrections to energy of n=3 states are -q1, -0.5*q1, 0, 0.5*q1, and q1/2
    q1 = sqrt(4*(e0*h*B/(4*pi*m))**2 + 9*(3*a_0*e0*E)**2)
    ! szratio is gamma/epsilon factor. can be thought of as roughly zeeman energy shift/ stark energy shift.
    if(E.eq.0.d0)then
    	szratio = 0.d0
    else
    	szratio = (e0*h*B/(4*pi*m))/(3*a_0*e0*E)
    endif
    ! wavelengths calculated from h*c0/lambda =  E_i - E_j for transition from i to j energies
    ! order is small wavelengths to large wavelengths
    if(n_stark.eq.15) then
        wavel(1)  = 2*c0*h*l0/(2*c0*h+(-2*q0*(-1)+q1*(2))*l0)
        wavel(2)  = 2*c0*h*l0/(2*c0*h+(-2*q0*(0)+q1*(2))*l0)
        wavel(3) = 2*c0*h*l0/(2*c0*h+(-2*q0*(-1)+q1*(1))*l0)
        wavel(4) = 2*c0*h*l0/(2*c0*h+(-2*q0*(1)+q1*(2))*l0)
        wavel(5) = 2*c0*h*l0/(2*c0*h+(-2*q0*(0)+q1*(1))*l0)
        wavel(6) = 2*c0*h*l0/(2*c0*h+(-2*q0*(-1)+q1*(0))*l0)
        wavel(7) = 2*c0*h*l0/(2*c0*h+(-2*q0*(1)+q1*(1))*l0)
        wavel(8) = 2*c0*h*l0/(2*c0*h+(-2*q0*(0)+q1*(0))*l0)
        wavel(9) = 2*c0*h*l0/(2*c0*h+(-2*q0*(-1)+q1*(-1))*l0)
        wavel(10) = 2*c0*h*l0/(2*c0*h+(-2*q0*(1)+q1*(0))*l0)
        wavel(11) = 2*c0*h*l0/(2*c0*h+(-2*q0*(0)+q1*(-1))*l0)
        wavel(12) = 2*c0*h*l0/(2*c0*h+(-2*q0*(-1)+q1*(-2))*l0)
        wavel(13) = 2*c0*h*l0/(2*c0*h+(-2*q0*(1)+q1*(-1))*l0)
        wavel(14) = 2*c0*h*l0/(2*c0*h+(-2*q0*(0)+q1*(-2))*l0)
        wavel(15) = 2*c0*h*l0/(2*c0*h+(-2*q0*(1)+q1*(-2))*l0)
        circularity = [1.0/3,4.0/3,5.0/3,-97.0/123,-1.0/2,-17.0/81,-1.0/3,0.0,1.0/3,17.0/81,1.0/2,97.0/123,-5.0/3,-4.0/3,-1.0/3]
        circularity = circularity*szratio
    else if(n_stark.eq.3) then
        wavel(1) = c0*h*l0/(c0*h+q0*(-1)*l0)
        wavel(2) = c0*h*l0/(c0*h+q0*(0)*l0)
        wavel(3) = c0*h*l0/(c0*h+q0*(1)*l0)
        circularity = [1.0,0.0,-1.0]
        circularity = circularity*szratio
    endif
    !lambda =  lambda_shifted + E * stark_wavel ![nm] <-old calculation for pure stark effect
    wavel = (wavel-l0)*10**9    ! converting to [nm]
    lambda = lambda_shifted + wavel
    !Intensities of stark components
    if (E .eq. 0.d0) then
        cos_los_Efield = 0.d0
    else
        cos_los_Efield = dot_product(vp,efield) / E
    endif
    if (B .eq. 0.d0) then
        cos_los_Bfield = 0.d0
    else
        cos_los_Bfield = dot_product(vp,bfield) / B
    endif
    !calculate stokes parameters
    intensity = stark_intens
    stokes(:,1) = intensity*(1.d0+ cos_los_Efield**2)*stark_sigma*sigma_pi +  intensity*(1.d0- cos_los_Efield**2)*stark_pi
    stokes(:,2) = intensity*(1.d0- cos_los_Efield**2)*stark_sigma*sigma_pi - intensity*(1.d0- cos_los_Efield**2)*stark_pi
    stokes(:,3) = 0.d0
    stokes(:,4) = 2*intensity*cos_los_Bfield*circularity
    intensity = stark_intens*(1.d0+ stark_sign* cos_los_Efield**2)
    !! E.g. mirrors may change the pi to sigma intensity ratio
    where (stark_sigma .eq. 1)
        intensity = intensity * sigma_pi
    endwhere
    !! normalize and multiply with photon density from colrad
    normfactor = (1/sum(intensity))*photons*dlength
    intensity = intensity/sum(intensity)*photons*dlength
    stokes = stokes*normfactor
endsubroutine spectrum

subroutine store_photons(pos, vi, lambda0, photons, spectra, stokevec, passive)
    !+ Store photons in `spectra`
    real(Float64), dimension(3), intent(in)      :: pos
        !+ Position of neutral in beam coordinates [machine coordinates for passive case]
    real(Float64), dimension(3), intent(in)      :: vi
        !+ Velocitiy of neutral in beam coordinates [cm/s]
    real(Float64), intent(in)                    :: lambda0
        !+ Reference wavelength [nm]
    real(Float64), intent(in)                    :: photons
        !+ Photons from [[libfida:colrad]] [Ph/(s*cm^3)]
    real(Float64), dimension(:,:,:), intent(inout) :: spectra
    !+ Stark split `spectra`
    real(Float64), dimension(:,:,:,:), intent(inout) :: stokevec
    !+ Stark split `stokes vector`
    logical, intent(in), optional                :: passive
        !+ Indicates whether photon is passive FIDA

    real(Float64), dimension(n_stark) :: lambda, intensity
    real(Float64), dimension(n_stark,4) :: stokes
    real(Float64) :: dlength, sigma_pi
    type(LocalEMFields) :: fields
    integer(Int32), dimension(3) :: ind
    real(Float64), dimension(3) :: pos_xyz, lens_xyz, cyl, vp
    type(LOSinters) :: inter
    integer :: ichan,i,j,bin,nchan, k
    logical :: pas = .False.

    if(present(passive)) pas = passive

    if(pas) then
        cyl(1) = sqrt(pos(1)*pos(1) + pos(2)*pos(2))
        cyl(2) = pos(3)
        cyl(3) = atan2(pos(2), pos(1))
        call get_passive_grid_indices(cyl,ind)
        inter = spec_chords%cyl_inter(ind(1),ind(2),ind(3))
        call uvw_to_xyz(pos, pos_xyz)
    else
        call get_indices(pos,ind)
        inter = spec_chords%inter(ind(1),ind(2),ind(3))
        pos_xyz = pos
    endif

    nchan = inter%nchan
    if(nchan.eq.0) return

    call get_fields(fields,pos=pos_xyz)

    loop_over_channels: do j=1,nchan
        ichan = inter%los_elem(j)%id
        dlength = inter%los_elem(j)%length
        sigma_pi = spec_chords%los(ichan)%sigma_pi
        if(pas) then
            call uvw_to_xyz(spec_chords%los(ichan)%lens_uvw,lens_xyz)
        else
            lens_xyz = spec_chords%los(ichan)%lens
        endif
        vp = pos_xyz - lens_xyz
        call spectrum(vp,vi,fields,lambda0,sigma_pi,photons, &
                      dlength,lambda,intensity, stokes)

        loop_over_stark: do i=1,n_stark
            bin=floor((lambda(i)-inputs%lambdamin)/inputs%dlambda) + 1
            if (bin.lt.1) cycle loop_over_stark
            if (bin.gt.inputs%nlambda) cycle loop_over_stark
            !$OMP CRITICAL(spec_stokes)
            spectra(i,bin,ichan) = spectra(i,bin,ichan) + intensity(i)
            stokevec(i,:,bin,ichan) = stokevec(i,:,bin,ichan) + stokes(i,:)
            !$OMP END CRITICAL(spec_stokes)
        enddo loop_over_stark
    enddo loop_over_channels
end subroutine store_photons

subroutine store_nbi_photons(pos, vi, lambda0, photons, neut_type, store_spatial)
    !+ Store BES photons in [[libfida:spectra]]
    real(Float64), dimension(3), intent(in) :: pos
        !+ Position of neutral in beam grid coordinates
    real(Float64), dimension(3), intent(in) :: vi
        !+ Velocitiy of neutral [cm/s]
    real(Float64), intent(in)               :: lambda0
        !+ Reference wavelength [nm]
    real(Float64), intent(in)               :: photons
        !+ Photons from [[libfida:colrad]] [Ph/(s*cm^3)]
    integer,intent(in)                      :: neut_type
        !+ Neutral type (full,half,third,halo)
    logical, intent(in), optional           :: store_spatial
        !+ Store neutral position and photons in reservoir

    logical :: store_spat = .False.

    if(present(store_spatial)) then
        store_spat = store_spatial
    endif

    select case (neut_type)
           case (nbif_type)
               call store_photons(pos,vi,lambda0,photons,spec%full,spec%fullstokes)
               if(store_spat) then
                   call store_photon_birth(pos, photons, spatres%full)
               endif
           case (nbih_type)
               call store_photons(pos,vi,lambda0,photons,spec%half, spec%halfstokes)
               if(store_spat) then
                   call store_photon_birth(pos, photons, spatres%half)
               endif
           case (nbit_type)
               call store_photons(pos,vi,lambda0,photons,spec%third,spec%thirdstokes)
               if(store_spat) then
                   call store_photon_birth(pos, photons, spatres%third)
               endif
           case default
               if(inputs%verbose.ge.0) then
                   write(*,'("STORE_NBI_PHOTONS: Unknown neutral type: ",i2)') neut_type
               endif
               stop
    end select

end subroutine store_nbi_photons

subroutine store_fida_photons(pos, vi, lambda0, photons, orbit_class, passive)
    !+ Store fida photons in [[libfida:spectra]]
    real(Float64), dimension(3), intent(in) :: pos
        !+ Position of neutral in beam grid coordinates
    real(Float64), dimension(3), intent(in) :: vi
        !+ Velocitiy of neutral [cm/s]
    real(Float64), intent(in)               :: lambda0
        !+ Reference wavelength [nm]
    real(Float64), intent(in)               :: photons
        !+ Photons from [[libfida:colrad]] [Ph/(s*cm^3)]
    integer, intent(in), optional           :: orbit_class
        !+ Orbit class ID
    logical, intent(in), optional           :: passive
        !+ Indicates whether photon is passive FIDA

    integer :: iclass = 1
    logical :: pas = .False.
    type(LOSinters) :: inter

    if(present(orbit_class)) then
        iclass = min(orbit_class,particles%nclass)
    endif

    if(present(passive)) pas = passive

    if(pas) then
        call store_photons(pos, vi, lambda0, photons, spec%pfida(:,:,:,iclass), spec%pfidastokes(:,:,:,:,iclass), passive=.True.)
    else
        call store_photons(pos, vi, lambda0, photons, spec%fida(:,:,:,iclass),spec%fidastokes(:,:,:,:,iclass))
    endif

end subroutine store_fida_photons

subroutine store_photon_birth(pos, photons, res, passive)
    !+ Store neutral birth location of the photon source
    real(Float64), dimension(3), intent(in)       :: pos
        !+ Birth location of the photon source neutral
    real(Float64), intent(in)                     :: photons
        !+ Number of photons
    type(NeutralParticleReservoir), dimension(:), intent(inout) :: res
        !+ reservoir of neutral particles
    logical, intent(in), optional                 :: passive

    integer :: i = 1, ichan
    logical :: pas = .False.
    real(Float64) :: cyl(3), pos_xyz(3)
    integer(Int32) :: ind(3), nchan
    type(LOSinters) :: inter

    if(present(passive)) pas = passive

    if(pas) then
        cyl(1) = sqrt(pos(1)*pos(1) + pos(2)*pos(2))
        cyl(2) = pos(3)
        cyl(3) = atan2(pos(2), pos(1))
        call get_passive_grid_indices(cyl,ind)
        inter = spec_chords%cyl_inter(ind(1),ind(2),ind(3))
    else
        call get_indices(pos,ind)
        inter = spec_chords%inter(ind(1),ind(2),ind(3))
    endif

    nchan = inter%nchan
    if(nchan.eq.0) return

    loop_over_channels: do i=1,nchan
        ichan = inter%los_elem(i)%id
        call update_reservoir(res(ichan), pos, photons)
    enddo loop_over_channels

end subroutine store_photon_birth

subroutine store_neutrons(rate, orbit_class)
    !+ Store neutron rate in [[libfida:neutron]]
    real(Float64), intent(in)     :: rate
        !+ Neutron rate [neutrons/sec]
    integer, intent(in), optional :: orbit_class
        !+ Orbit class ID

    integer :: iclass

    if(present(orbit_class)) then
        iclass = min(orbit_class,particles%nclass)
    else
        iclass = 1
    endif

    !$OMP ATOMIC UPDATE
    neutron%rate(iclass)= neutron%rate(iclass) + rate
    !$OMP END ATOMIC

end subroutine store_neutrons

subroutine store_fw_photons_at_chan(ichan,eind,pind,vp,vi,lambda0,fields,dlength,sigma_pi,denf,photons)
    !+ Store FIDA weight photons in [[libfida:fweight]] for a specific channel
    integer, intent(in)                     :: ichan
        !+ Channel index
    integer, intent(in)                     :: eind
        !+ Energy index
    integer, intent(in)                     :: pind
        !+ Pitch index
    real(Float64), dimension(3), intent(in) :: vp
        !+ Vector pointing toward optical head
    real(Float64), dimension(3), intent(in) :: vi
        !+ Velocity of neutral [cm/s]
    real(Float64), intent(in)               :: lambda0
        !+ Reference Wavelength [nm]
    type(LocalEMFields), intent(in)         :: fields
        !+ Electro-magnetic fields
    real(Float64), intent(in)               :: dlength
        !+ LOS intersection length with [[libfida:beam_grid]] cell particle is in
    real(Float64), intent(in)               :: sigma_pi
        !+ Sigma-pi ratio for channel
    real(Float64), intent(in)               :: denf
        !+ Fast-ion density [cm^-3]
    real(Float64), intent(in)               :: photons
        !+ Photons from [[libfida:colrad]] [Ph/(s*cm^3)]

    real(Float64), dimension(n_stark) :: lambda,intensity
    real(Float64), dimension(n_stark,4) :: stokes
    real(Float64) :: dlambda,intens_fac
    integer :: i,bin


    dlambda=(inputs%lambdamax_wght-inputs%lambdamin_wght)/inputs%nlambda_wght
    intens_fac = (1.d0)/(4.d0*pi*dlambda)
    call spectrum(vp,vi,fields, lambda0, sigma_pi,photons, &
                  dlength,lambda,intensity, stokes = stokes)

    !$OMP CRITICAL(fida_wght)
    loop_over_stark: do i=1,n_stark
        bin=floor((lambda(i) - inputs%lambdamin_wght)/dlambda) + 1
        if (bin.lt.1)                   cycle loop_over_stark
        if (bin.gt.inputs%nlambda_wght) cycle loop_over_stark
        !fida(bin,ichan)= fida(bin,ichan) + &
        !  (denf*intens_fac*1.d4)*intensity(i) !ph/(s*nm*sr*m^2)
        fweight%weight(bin,eind,pind,ichan) = &
          fweight%weight(bin,eind,pind,ichan) + intensity(i)*intens_fac !(ph*cm)/(s*nm*sr*fast-ion*dE*dp)
    enddo loop_over_stark
    if(denf.gt.0.d0) then
        fweight%mean_f(eind,pind,ichan) = fweight%mean_f(eind,pind,ichan) + &
                                          (denf*intens_fac)*sum(intensity)
    endif
    !$OMP END CRITICAL(fida_wght)

end subroutine store_fw_photons_at_chan

subroutine store_fw_photons(eind, pind, pos, vi, lambda0, denf, photons)
    !+ Store FIDA weight photons in [[libfida:fweight]]
    integer, intent(in)                     :: eind
        !+ Energy index
    integer, intent(in)                     :: pind
        !+ Pitch index
    real(Float64), dimension(3), intent(in) :: pos
        !+ Position of neutral
    real(Float64), dimension(3), intent(in) :: vi
        !+ Velocity of neutral [cm/s]
    real(Float64), intent(in)               :: lambda0
        !+ Reference Wavelength [nm]
    real(Float64), intent(in)               :: denf
        !+ Fast-ion density [cm^-3]
    real(Float64), intent(in)               :: photons
        !+ Photons from [[libfida:colrad]] [Ph/(s*cm^3)]

    real(Float64) :: dlength, sigma_pi
    type(LocalEMFields) :: fields
    integer(Int32), dimension(3) :: ind
    real(Float64), dimension(3) :: vp
    type(LOSInters) :: inter
    integer :: ichan,nchan,i

    call get_indices(pos,ind)
    inter = spec_chords%inter(ind(1),ind(2),ind(3))
    nchan = inter%nchan
    if(nchan.eq.0) return

    call get_fields(fields,pos=pos)

    loop_over_channels: do i=1,nchan
        ichan = inter%los_elem(i)%id
        dlength = inter%los_elem(i)%length
        sigma_pi = spec_chords%los(ichan)%sigma_pi
        vp = pos - spec_chords%los(ichan)%lens
        call store_fw_photons_at_chan(ichan, eind, pind, &
             vp, vi, lambda0, fields, dlength, sigma_pi, denf, photons)
    enddo loop_over_channels

end subroutine store_fw_photons

!=============================================================================
!---------------------------Monte Carlo Routines------------------------------
!=============================================================================
subroutine get_nlaunch(nr_markers,papprox, nlaunch)
    !+ Sets the number of MC markers launched from each [[libfida:beam_grid]] cell
    integer(Int64), intent(in)                    :: nr_markers
        !+ Approximate total number of markers to launch
    real(Float64), dimension(:,:,:), target, intent(in)   :: papprox
        !+ [[libfida:beam_grid]] cell weights
    integer(Int32), dimension(:,:,:), intent(out) :: nlaunch
        !+ Number of mc markers to launch for each cell: nlaunch(x,y,z)

    logical, dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz) :: mask
    real(Float64), dimension(beam_grid%ngrid) :: cdf
    integer  :: c, i, j, k, nc, nm, ind(3)
    integer  :: nmin = 5
    integer, dimension(1) :: randomi
    type(rng_type) :: r
    real(Float64), pointer :: papprox_ptr(:)

    !! Fill in minimum number of markers per cell
    nlaunch = 0
    mask = papprox.gt.0.0
    where(mask)
        nlaunch = nmin
    endwhere

    !! If there are any left over distribute according to papprox
    nc = count(mask)
    if(nr_markers.gt.(nmin*nc)) then
        nm = nr_markers - nmin*nc

        !! precalculate cdf to save time
        call c_f_pointer(c_loc(papprox), papprox_ptr, [beam_grid%ngrid])
        call cumsum(papprox_ptr, cdf)

        !! use the same seed for all processes
        call rng_init(r, 932117)
        do c=1, nm
            call randind_cdf(r, cdf, randomi)
            call ind2sub(beam_grid%dims, randomi(1), ind)
            i = ind(1) ; j = ind(2) ; k = ind(3)
            nlaunch(i,j,k) = nlaunch(i,j,k) + 1
        enddo
    endif

end subroutine get_nlaunch

subroutine get_nlaunch_pass_grid(nr_markers,papprox, nlaunch)
    !+ Sets the number of MC markers launched from each [[libfida:pass_grid]] cell
    integer(Int64), intent(in)                    :: nr_markers
        !+ Approximate total number of markers to launch
    real(Float64), dimension(:,:,:), intent(in)   :: papprox
        !+ [[libfida:pass_grid]] cell weights
    integer(Int32), dimension(:,:,:), intent(out) :: nlaunch
        !+ Number of mc markers to launch for each cell: nlaunch(r,z,phi)

    logical, dimension(pass_grid%nr,pass_grid%nz,pass_grid%nphi) :: mask
    real(Float64), dimension(pass_grid%ngrid) :: cdf
    integer, dimension(1) :: randomi
    type(rng_type) :: r
    integer :: c, i, j, k, nc, nm, ind(3)
    integer :: nmin = 5

    !! Fill in minimum number of markers per cell
    nlaunch = 0
    mask = papprox.gt.0.0
    where(mask)
        nlaunch = nmin
    endwhere

    !! If there are any left over distribute according to papprox
    nc = count(mask)
    if(nr_markers.gt.(nmin*nc)) then
        nm = nr_markers - nmin*nc

        !! precalculate cdf to save time
        call cumsum(reshape(papprox,[pass_grid%ngrid]), cdf)
        !! use the same seed for all processes
        call rng_init(r, 932117)
        do c=1, nm
            call randind_cdf(r, cdf, randomi)
            call ind2sub(pass_grid%dims, randomi(1), ind)
            i = ind(1) ; j = ind(2) ; k = ind(3)
            nlaunch(i,j,k) = nlaunch(i,j,k) + 1
        enddo
    endif

end subroutine get_nlaunch_pass_grid

subroutine pitch_to_vec(pitch, gyroangle, fields, vi_norm)
    !+ Calculates velocity vector from pitch, gyroangle and fields
    real(Float64), intent(in)                :: pitch
        !+ Pitch
    real(Float64), intent(in)                :: gyroangle
        !+ Gyroangle [radians]
    type(LocalEMFields), intent(in)          :: fields
        !+ Electromagnetic fields
    real(Float64), dimension(3), intent(out) :: vi_norm
        !+ Normalized velocity vector

    real(Float64) :: sinus

    sinus = sqrt(max(1.d0-pitch**2,0.d0))
    vi_norm = (sinus*cos(gyroangle)*fields%a_norm + &
               pitch*fields%b_norm + &
               sinus*sin(gyroangle)*fields%c_norm)

end subroutine pitch_to_vec

subroutine gyro_step(vi, fields, Ai, r_gyro)
    !+ Calculates gyro-step
    !+
    !+###References
    !+ Belova, E. V., N. N. Gorelenkov, and C. Z. Cheng. "Self-consistent equilibrium model of low aspect-
    !+ ratio toroidal plasma with energetic beam ions." Physics of Plasmas (1994-present) 10.8 (2003):
    !+ 3240-3251. Appendix A: Last equation
    real(Float64), dimension(3), intent(in)  :: vi
        !+ Ion velocity
    type(LocalEMFields), intent(in)          :: fields
        !+ Electro-magnetic fields
    real(Float64), intent(in)                :: Ai
        !+ Atomic mass of ion
    real(Float64), dimension(3), intent(out) :: r_gyro
        !+ Gyro-step
        !+ Gyro-radius vector from particle position to guiding center

    real(Float64), dimension(3) :: vxB, rg_uvw, uvw, cuvrxb, b_rtz, grad_B, rg_rtz
    real(Float64) :: one_over_omega, phi, R, vpar, term1, term2

    if(inputs%flr.ge.1) then
        uvw = fields%uvw
        R = sqrt(uvw(1)**2 + uvw(2)**2)
        phi = atan2(uvw(2),uvw(1))
        one_over_omega=Ai*mass_u/(fields%b_abs*e0)
        vxB = cross_product(vi,fields%b_norm)
        vpar =  dot_product(vi,fields%b_norm)
        r_gyro = vxB*one_over_omega !points towards gyrocenter, in beam coordinates

        if(inputs%flr.ge.2) then
            !! convert the r_gyro vector to machine coordiantes
            if(fields%coords.eq.0) then
                rg_uvw=  matmul(beam_grid%basis, r_gyro)
            endif
            if(fields%coords.eq.1) then
                rg_uvw = r_gyro
            endif

            b_rtz(1) = fields%br/fields%b_abs
            b_rtz(2) = fields%bt/fields%b_abs
            b_rtz(3) = fields%bz/fields%b_abs
            cuvrxb(1) = (1./R*fields%dbz_dphi-fields%dbt_dz)/fields%b_abs
            cuvrxb(2) = (fields%dbr_dz - fields%dbz_dr)/fields%b_abs
            cuvrxb(3) = (1.0/R*fields%bt + fields%dbt_dr - 1.0/R*fields%dbr_dphi)/fields%b_abs
            term1 = vpar*one_over_omega*dot_product(b_rtz,cuvrxb)
            grad_B(1) = (fields%br*fields%dbr_dr + fields%bt * fields%dbt_dr + fields%bz*fields%dbz_dr)/&
                        fields%b_abs
            grad_B(2) = 1.0/R*(fields%br*fields%dbr_dphi + fields%bt * fields%dbt_dphi + fields%bz*fields%dbz_dphi)/&
                        fields%b_abs
            grad_B(3) = (fields%br*fields%dbr_dz + fields%bt * fields%dbt_dz + fields%bz*fields%dbz_dz)/&
                        fields%b_abs
            !convert rg_uvw vector to cylindrical coordiantes
            rg_rtz(1) = rg_uvw(1)*cos(phi) + rg_uvw(2)*sin(phi)
            rg_rtz(2) = -rg_uvw(1)*sin(phi) + rg_uvw(2)*cos(phi)
            rg_rtz(3) = rg_uvw(3)
            term2 = -1.0 / (2.0 * fields%b_abs)*dot_product(rg_rtz,grad_B)
        else
            term1 = 0.0
            term2 = 0.0
        endif

        r_gyro = r_gyro * (1.0 - term1 - term2)
        if ((1.0 - term1 - term2 .le. 0.0) .or. (1.0 - term1 - term2 .ge. 2.0) ) then
            write(*,*) 'GYRO_STEP: Gyro correction results in negative distances or too large shift: ', &
                          1.0-term1-term2
            stop
        endif
    else
        r_gyro = 0.d0
    endif

end subroutine gyro_step

subroutine gyro_correction(fields, energy, pitch, Ai, rp, vp, theta_in)
    !+ Calculates gyro correction for Guiding Center MC distribution calculation
    type(LocalEMFields), intent(in)          :: fields
        !+ Electromagnetic fields at guiding center
    real(Float64), intent(in)                :: energy
        !+ Energy of particle
    real(Float64), intent(in)                :: pitch
        !+ Particle pitch w.r.t the magnetic field
    real(Float64), intent(in)                :: Ai
        !+ Atomic mass of Ion
    real(Float64), dimension(3), intent(out) :: rp
        !+ Particle position
    real(Float64), dimension(3), intent(out) :: vp
        !+ Particle velocity
    real(Float64), intent(in), optional      :: theta_in
        !+ Gyro-angle
    real(Float64), dimension(3) :: vi_norm, r_step
    real(Float64), dimension(1) :: randomu
    real(Float64) :: vabs, theta

    vabs  = sqrt(energy/(v2_to_E_per_amu*Ai))

    if(present(theta_in)) then
        theta = theta_in
    else
        !! Sample gyroangle
        call randu(randomu)
        theta = 2*pi*randomu(1)
    endif

    !! Calculate velocity vector
    call pitch_to_vec(pitch, theta, fields, vi_norm)
    vp = vabs*vi_norm

    !! Move to particle location
    call gyro_step(vp, fields, Ai, r_step)
    rp = fields%pos - r_step

end subroutine gyro_correction

function gyro_radius(fields, energy, pitch, Ai) result (gyro_rad)
    !+ Calculates mean gyro-radius
    type(LocalEMFields), intent(in)          :: fields
        !+ Electromagnetic fields at guiding center
    real(Float64), intent(in)                :: energy
        !+ Energy of particle
    real(Float64), intent(in)                :: pitch
        !+ Particle pitch w.r.t the magnetic field
    real(Float64), intent(in)                :: Ai
        !+ Atomic Mass of ion
    real(Float64) :: gyro_rad
        !+ Mean gyro-radius

    real(Float64), dimension(3) :: vi_norm, r_step
    real(Float64) :: vabs, phi
    integer :: i,n

    vabs  = sqrt(energy/(v2_to_E_per_amu*Ai))

    gyro_rad = 0.d0
    n = 6
    do i=1,n
        phi = i*2*pi/n
        call pitch_to_vec(pitch, phi, fields, vi_norm)
        call gyro_step(vabs*vi_norm, fields, Ai, r_step)
        gyro_rad = gyro_rad + norm2(r_step)/n
    enddo

end function gyro_radius

subroutine mc_beam_grid(ind, ri)
    !+ Sample uniformally in beam grid cell
    integer, dimension(3), intent(in)        :: ind
        !+ Cell index
    real(Float64), dimension(3), intent(out) :: ri
        !+ Position in cell

    real(Float64), dimension(3) :: random3

    call randu(random3)
    ri(1) = beam_grid%xc(ind(1)) + beam_grid%dr(1)*(random3(1) - 0.5)
    ri(2) = beam_grid%yc(ind(2)) + beam_grid%dr(2)*(random3(2) - 0.5)
    ri(3) = beam_grid%zc(ind(3)) + beam_grid%dr(3)*(random3(3) - 0.5)

end subroutine mc_beam_grid

subroutine mc_passive_grid(ind, ri)
    !+ Sample uniformally in passive grid cell
    integer, dimension(3), intent(in)        :: ind
        !+ Cell index
    real(Float64), dimension(3), intent(out) :: ri
        !+ Position in cell

    real(Float64) :: rmin, rmax, zmin, phimin
    real(Float64), dimension(3) :: randomu3, ri_cyl

    call randu(randomu3)
    rmin = pass_grid%r(ind(1))
    rmax = rmin + pass_grid%dr
    zmin = pass_grid%z(ind(2))
    phimin = pass_grid%phi(ind(3))

    ! Sample uniformally in annulus
    ri_cyl(1) = sqrt(randomu3(1)*(rmax**2 - rmin**2) + rmin**2)
    ri_cyl(2) = zmin + randomu3(2)*pass_grid%dz
    ri_cyl(3) = phimin + randomu3(3)*pass_grid%dphi
    call cyl_to_uvw(ri_cyl, ri)

end subroutine mc_passive_grid

subroutine mc_fastion(ind,fields,eb,ptch,denf)
    !+ Samples a Guiding Center Fast-ion distribution function at a given [[libfida:beam_grid]] index
    integer, dimension(3), intent(in) :: ind
        !+ [[libfida:beam_grid]] index
    type(LocalEMFields), intent(out)       :: fields
        !+ Electromagnetic fields at the guiding center
    real(Float64), intent(out)             :: eb
        !+ Energy of the fast ion
    real(Float64), intent(out)             :: ptch
        !+ Pitch of the fast ion
    real(Float64), intent(out)             :: denf
        !+ Fast-ion density at guiding center

    real(Float64), dimension(fbm%nenergy,fbm%npitch) :: fbeam
    real(Float64), dimension(3) :: rg
    real(Float64), dimension(3) :: randomu3
    integer, dimension(2,1) :: ep_ind

    call mc_beam_grid(ind, rg)

    denf=0.d0

    call get_fields(fields,pos=rg)
    if(.not.fields%in_plasma) return

    call get_distribution(fbeam,denf,pos=rg, coeffs=fields%b)
    call randind(fbeam,ep_ind)
    call randu(randomu3)
    eb = fbm%energy(ep_ind(1,1)) + fbm%dE*(randomu3(1)-0.5)
    ptch = fbm%pitch(ep_ind(2,1)) + fbm%dp*(randomu3(2)-0.5)

end subroutine mc_fastion

subroutine mc_fastion_pass_grid(ind,fields,eb,ptch,denf,output_coords)
    !+ Samples a Guiding Center Fast-ion distribution function at a given [[libfida:pass_grid]] index
    integer, dimension(3), intent(in)      :: ind
        !+ [[libfida:pass_grid]] index
    type(LocalEMFields), intent(out)       :: fields
        !+ Electromagnetic fields at the guiding center
    real(Float64), intent(out)             :: eb
        !+ Energy of the fast ion
    real(Float64), intent(out)             :: ptch
        !+ Pitch of the fast ion
    real(Float64), intent(out)             :: denf
        !+ Fast-ion density at guiding center
    integer, intent(in), optional          :: output_coords
        !+ Indicates coordinate system of `fields`. Beam grid (0), machine (1) and cylindrical (2)

    real(Float64), dimension(fbm%nenergy,fbm%npitch) :: fbeam
    real(Float64), dimension(3) :: rg
    real(Float64), dimension(3) :: randomu3
    integer, dimension(2,1) :: ep_ind
    integer :: ocs

    if(present(output_coords)) then
        ocs = output_coords
    else
        ocs = 0
    endif

    denf=0.d0
    call mc_passive_grid(ind, rg)

    call get_fields(fields,pos=rg,input_coords=1,output_coords=ocs)
    if(.not.fields%in_plasma) return

    call get_distribution(fbeam,denf,coeffs=fields%b)
    call randind(fbeam,ep_ind)
    call randu(randomu3)
    eb = fbm%energy(ep_ind(1,1)) + fbm%dE*(randomu3(1)-0.5)
    ptch = fbm%pitch(ep_ind(2,1)) + fbm%dp*(randomu3(2)-0.5)

end subroutine mc_fastion_pass_grid

subroutine mc_halo(plasma, ai, vhalo)
    !+ Generate velocity vector from a thermal Maxwellian distribution
    type(LocalProfiles), intent(in) :: plasma
        !+ Plasma Parameters
    real(Float64), intent(in)       :: ai
        !+ Ion mass [amu]
    real(Float64), dimension(3), intent(out)           :: vhalo
        !+ Velocity [cm/s]

    real(Float64), dimension(3) :: random3

    call randn(random3)

    vhalo = plasma%vrot + sqrt(plasma%ti*0.5/(v2_to_E_per_amu*ai))*random3 !![cm/s]

end subroutine mc_halo

subroutine mc_nbi(vnbi,efrac,rnbi,err)
    !+ Generates a neutral beam particle trajectory
    integer, intent(in)                      :: efrac
        !+ Beam neutral type (1,2,3)
    real(Float64), dimension(3), intent(out) :: vnbi
        !+ Velocity [cm/s]
    real(Float64), dimension(3), intent(out) :: rnbi
        !+ Starting position on [[libfida:beam_grid]]
    logical, intent(out)                     :: err
        !+ Error Code

    real(Float64), dimension(3) :: r_exit
    real(Float64), dimension(3) :: uvw_src    !! Start position on ion source
    real(Float64), dimension(3) :: xyz_src    !! Start position on ion source
    real(Float64), dimension(3) :: uvw_ray    !! NBI velocity in uvw coords
    real(Float64), dimension(3) :: xyz_ray    !! NBI velocity in xyz coords
    real(Float64), dimension(3) :: xyz_ape    !! Aperture plane intersection point
    real(Float64), dimension(2) :: randomu    !! uniform random numbers
    real(Float64), dimension(2) :: randomn    !! normal random numbers
    real(Float64) :: length, sqrt_rho, theta
    integer :: i, j
    logical :: inp, valid_trajectory

    err = .False.
    valid_trajectory = .False.
    rejection_loop: do i=1,1000
        call randu(randomu)
        select case (nbi%shape)
            case (1)
                ! Uniformally sample in rectangle
                xyz_src(1) =  0.d0
                xyz_src(2) =  nbi%widy * 2.d0*(randomu(1)-0.5d0)
                xyz_src(3) =  nbi%widz * 2.d0*(randomu(2)-0.5d0)
            case (2)
                ! Uniformally sample in ellipse
                sqrt_rho = sqrt(randomu(1))
                theta = 2*pi*randomu(2)
                xyz_src(1) = 0.d0
                xyz_src(2) = nbi%widy*sqrt_rho*cos(theta)
                xyz_src(3) = nbi%widz*sqrt_rho*sin(theta)
        end select

        !! Create random velocity vector
        call randn(randomn)
        xyz_ray(1)= 1.d0
        xyz_ray(2)=(-xyz_src(2)/nbi%focy + tan(nbi%divy(efrac)*randomn(1)))
        xyz_ray(3)=(-xyz_src(3)/nbi%focz + tan(nbi%divz(efrac)*randomn(2)))

        aperture_loop: do j=1,nbi%naperture
            xyz_ape = xyz_ray*nbi%adist(j) + xyz_src
            select case (nbi%ashape(j))
                case (1)
                    if ((abs(xyz_ape(2) - nbi%aoffy(j)).gt.nbi%awidy(j)).or.&
                        (abs(xyz_ape(3) - nbi%aoffz(j)).gt.nbi%awidz(j))) then
                        cycle rejection_loop
                    endif
                case (2)
                    if ((((xyz_ape(2) - nbi%aoffy(j))*nbi%awidz(j))**2 + &
                         ((xyz_ape(3) - nbi%aoffz(j))*nbi%awidy(j))**2).gt. &
                         (nbi%awidy(j)*nbi%awidz(j))**2) then
                        cycle rejection_loop
                    endif
            end select
        enddo aperture_loop
        valid_trajectory = .True.

        !! Convert to beam centerline coordinates to beam grid coordinates
        uvw_src = matmul(nbi%basis,xyz_src) + nbi%src
        uvw_ray = matmul(nbi%basis,xyz_ray)
        exit rejection_loop
    enddo rejection_loop

    !Set Default trajectory in case rejection sampling fails
    if(.not.valid_trajectory)then
        if(inputs%verbose.ge.0) then
            write(*,'(a)') "MC_NBI: Failed to find trajectory though aperture(s). Using beam centerline."
        endif
        uvw_src = nbi%src
        uvw_ray = nbi%axis
    endif
    vnbi = uvw_ray/norm2(uvw_ray)

    !! Determine start position on beam grid
    call grid_intersect(uvw_src,vnbi,length,rnbi,r_exit)
    if(length.le.0.0)then
        err = .True.
        nbi_outside = nbi_outside + 1
    endif

    !! Check if start position is in the plasma
    call in_plasma(rnbi,inp)
    if(inp)then
        if(inputs%verbose.ge.0) then
            write(*,'(a)') "MC_NBI: A beam neutral has started inside the plasma."
            write(*,'(a)') "Move the beam grid closer to the source to fix"
        endif
        stop
    endif

    !! Determine velocity of neutrals corrected by efrac
    vnbi = vnbi*nbi%vinj/sqrt(real(efrac))
end subroutine mc_nbi

subroutine mc_nbi_cell(ind, neut_type, vnbi, weight)
    !+ Generates a neutral beam velocity vector
    !+ that passes through cell at `ind` with weight `weight`
    integer, dimension(3), intent(in)        :: ind
        !+ Cell index
    integer, intent(in)                      :: neut_type
        !+ Neutral Type (1=Full,2=Half,3=Third)
    real(Float64), dimension(3), intent(out) :: vnbi
        !+ Normalized Velocity
    real(Float64), intent(out)               :: weight
        !+ Weigth/probability of trajectory

    real(Float64), dimension(3) :: rc         !! Center of cell in uvw coords
    real(Float64), dimension(3) :: uvw_rf     !! End position in xyz coords
    real(Float64), dimension(3) :: xyz_rf     !! End position in xyz coords
    real(Float64), dimension(3) :: uvw_src    !! Start position on ion source
    real(Float64), dimension(3) :: xyz_src    !! Start position on ion source
    real(Float64), dimension(3) :: uvw_ray    !! NBI velocity in uvw coords
    real(Float64), dimension(3) :: xyz_ray    !! NBI velocity in xyz coords
    real(Float64), dimension(3) :: xyz_ape    !! Aperture plane intersection point
    real(Float64), dimension(3) :: randomu    !! uniform random numbers
    real(Float64) :: sqrt_rho, theta, vy, vz, theta_y, theta_z, py, pz
    integer :: i, j
    logical :: valid_trajectory

    rc = [beam_grid%xc(ind(1)), beam_grid%yc(ind(2)), beam_grid%zc(ind(3))]
    valid_trajectory = .False.
    rejection_loop: do i=1,1000
        call randu(randomu)
        select case (nbi%shape)
            case (1)
                ! Uniformally sample in rectangle
                xyz_src(1) =  0.d0
                xyz_src(2) =  nbi%widy * 2.d0*(randomu(1)-0.5d0)
                xyz_src(3) =  nbi%widz * 2.d0*(randomu(2)-0.5d0)
            case (2)
                ! Uniformally sample in ellipse
                sqrt_rho = sqrt(randomu(1))
                theta = 2*pi*randomu(2)
                xyz_src(1) = 0.d0
                xyz_src(2) = nbi%widy*sqrt_rho*cos(theta)
                xyz_src(3) = nbi%widz*sqrt_rho*sin(theta)
        end select

        !! Create random position in the cell
        call randu(randomu)
        uvw_rf = rc + (randomu - 0.5)*beam_grid%dr
        xyz_rf = matmul(nbi%inv_basis, uvw_rf - nbi%src)

        xyz_ray = xyz_rf - xyz_src
        xyz_ray = xyz_ray/norm2(xyz_ray)

        aperture_loop: do j=1,nbi%naperture
            xyz_ape = xyz_ray*nbi%adist(j) + xyz_src
            select case (nbi%ashape(j))
                case (1)
                    if ((abs(xyz_ape(2) - nbi%aoffy(j)).gt.nbi%awidy(j)).or.&
                        (abs(xyz_ape(3) - nbi%aoffz(j)).gt.nbi%awidz(j))) then
                        cycle rejection_loop
                    endif
                case (2)
                    if ((((xyz_ape(2) - nbi%aoffy(j))*nbi%awidz(j))**2 + &
                         ((xyz_ape(3) - nbi%aoffz(j))*nbi%awidy(j))**2).gt. &
                         (nbi%awidy(j)*nbi%awidz(j))**2) then
                        cycle rejection_loop
                    endif
            end select
        enddo aperture_loop
        valid_trajectory = .True.

        !! Convert to beam centerline coordinates to beam grid coordinates
        uvw_src = matmul(nbi%basis,xyz_src) + nbi%src
        uvw_ray = matmul(nbi%basis,xyz_ray)
        vnbi = nbi%vinj*uvw_ray/norm2(uvw_ray)/sqrt(real(neut_type))

        exit rejection_loop
    enddo rejection_loop

    !Set Default trajectory in case rejection sampling fails
    if(.not.valid_trajectory)then
        call randu(randomu)
        uvw_rf = rc + (randomu - 0.5)*beam_grid%dr
        uvw_ray = uvw_rf - nbi%src
        vnbi = nbi%vinj*uvw_ray/norm2(uvw_ray)/sqrt(real(neut_type))
    endif

    !! Find probability of trajectory
    vy = xyz_ray(2)/xyz_ray(1)
    vz = xyz_ray(3)/xyz_ray(1)
    theta_y = atan(vy + xyz_src(2)/nbi%focy)
    theta_z = atan(vz + xyz_src(3)/nbi%focz)

    py = (1.0/(1.0 + (vy + xyz_src(2)/nbi%focy)**2)) * &
         exp(-(theta_y**2)/(2*nbi%divy(neut_type)**2)) / &
         sqrt(2*nbi%divy(neut_type)**2)

    pz = (1.0/(1.0 + (vz + xyz_src(3)/nbi%focy)**2)) * &
         exp(-(theta_z**2)/(2*nbi%divz(neut_type)**2)) / &
         sqrt(2*nbi%divz(neut_type)**2)

    weight = py*pz

end subroutine mc_nbi_cell

!=============================================================================
!------------------------Primary Simulation Routines--------------------------
!=============================================================================
subroutine ndmc
    !+ Calculates neutral beam deposition and spectra
    integer :: neut_type !! full half third energy
    real(Float64)  :: nlaunch   !! nr. of markers
    real(Float64)  :: nneutrals !! # NBI particles
    real(Float64), dimension(3) :: vnbi !! velocities(full..)
    real(Float64), dimension(3) :: rnbi !! initial position

    integer(Int64) :: jj, ii, kk, cnt
    integer :: ntrack, ic
    type(ParticleTrack), dimension(beam_grid%ntrack) :: tracks
    type(LocalProfiles) :: plasma
    type(LocalEMFields) :: fields
    real(Float64), dimension(nlevs) :: states, dens
    real(Float64) :: photons, iflux,flux_tot
    integer(Int32), dimension(3) :: ind
    real(Float64), dimension(3) :: ri,ri_gc,r_gyro
    real(Float64), dimension(1) :: randomu
    integer, dimension(1) :: randi
    logical :: err

    !! Initialize Neutral Population
    call init_neutral_population(neut%full)
    call init_neutral_population(neut%half)
    call init_neutral_population(neut%third)

    if(inputs%verbose.ge.1) then
        write(*,'(T6,"# of markers: ",i10)') inputs%n_nbi
        if(inputs%calc_birth.ge.1) then
            write(*,'(T6,"# of birth markers: 3 x",i10)') int(inputs%n_nbi*inputs%n_birth)
        endif
    endif

    !! # of injected neutrals = NBI power/energy_per_particle
    nneutrals=1.d6*nbi%pinj/ (1.d3*nbi%einj*e0 &
         *( nbi%current_fractions(1)      &
         +  nbi%current_fractions(2)/2.d0 &
         +  nbi%current_fractions(3)/3.d0 ) )

    nlaunch=real(inputs%n_nbi)
    !$OMP PARALLEL DO schedule(guided) &
    !$OMP& private(vnbi,rnbi,tracks,ntrack,plasma,fields,randi,flux_tot, &
    !$OMP& states,dens,iflux,photons,neut_type,jj,ii,kk,ind,err,ri,ri_gc,r_gyro)
    loop_over_markers: do ii=istart,inputs%n_nbi,istep
        energy_fractions: do neut_type=1,3
            !! (type = 1: full energy, =2: half energy, =3: third energy
            if(nbi%current_fractions(neut_type).eq.0.d0) cycle energy_fractions

            call mc_nbi(vnbi,neut_type,rnbi,err)
            if(err) cycle energy_fractions

            call track(rnbi,vnbi,tracks,ntrack)
            if(ntrack.eq.0) cycle energy_fractions

            !! Solve collisional radiative model along track
            flux_tot = 0.d0
            states=0.d0
            states(1)=nneutrals*nbi%current_fractions(neut_type)/beam_grid%dv
            loop_along_track: do jj=1,ntrack
                iflux = sum(states)
                ind = tracks(jj)%ind
                call get_plasma(plasma,pos=tracks(jj)%pos)

                call colrad(plasma,beam_mass,vnbi,tracks(jj)%time,states,dens,photons)
                dens = dens/nlaunch
                call store_neutrals(ind,tracks(jj)%pos,vnbi,neut_type,dens)
                tracks(jj)%flux = (iflux - sum(states))/nlaunch
                flux_tot = flux_tot + tracks(jj)%flux*beam_grid%dv

                if(inputs%calc_birth.ge.1) then
                    call store_births(ind,neut_type,tracks(jj)%flux)
                endif

                if((photons.gt.0.d0).and.(inputs%calc_bes.ge.1)) then
                    call store_nbi_photons(tracks(jj)%pos,vnbi,beam_lambda0, photons/nlaunch, neut_type, inputs%calc_res.eq.1)
                endif
            enddo loop_along_track
            if((inputs%calc_birth.ge.1).and.(flux_tot.gt.0.d0)) then
                !! Sample according to deposited flux along neutral trajectory
                !$OMP CRITICAL(ndmc_birth)
                do kk=1,inputs%n_birth
                    call randind(tracks(1:ntrack)%flux,randi)
                    call randu(randomu)
                    birth%part(birth%cnt)%neut_type = neut_type
                    birth%part(birth%cnt)%energy = nbi%einj/real(neut_type)
                    birth%part(birth%cnt)%weight = flux_tot/inputs%n_birth
                    birth%part(birth%cnt)%ind = tracks(randi(1))%ind
                    birth%part(birth%cnt)%vi = vnbi
                    ri = tracks(randi(1))%pos + vnbi*(tracks(randi(1))%time*(randomu(1)-0.5))
                    birth%part(birth%cnt)%ri = ri

                    call get_fields(fields,pos=ri)
                    birth%part(birth%cnt)%pitch = dot_product(fields%b_norm,vnbi/norm2(vnbi))
                    call gyro_step(vnbi,fields,beam_mass,r_gyro)
                    birth%part(birth%cnt)%ri_gc = ri + r_gyro
                    birth%cnt = birth%cnt + 1
                enddo
                !$OMP END CRITICAL(ndmc_birth)
            endif
        enddo energy_fractions
    enddo loop_over_markers
    !$OMP END PARALLEL DO

#ifdef _MPI
    !! Combine beam neutrals
    call parallel_merge_populations(neut%full)
    call parallel_merge_populations(neut%half)
    call parallel_merge_populations(neut%third)

    call parallel_sum(nbi_outside)
    if(inputs%calc_birth.ge.1) then
        call parallel_sum(birth%dens)
    endif
    !! Combine spectra
    if(inputs%calc_bes.ge.1) then
        call parallel_sum(spec%full)
        call parallel_sum(spec%half)
        call parallel_sum(spec%third)
        if(inputs%calc_res.ge.1) then
            do jj=1,spec_chords%nchan
                call parallel_merge_reservoirs(spatres%full(jj))
                call parallel_merge_reservoirs(spatres%half(jj))
                call parallel_merge_reservoirs(spatres%third(jj))
            enddo
        endif
    endif
#endif

    if(nbi_outside.gt.0)then
        if(inputs%verbose.ge.1) then
            write(*,'(T4,a, f6.2,a)') 'Percent of markers outside the grid: ', &
                                  100.*nbi_outside/(3.*inputs%n_nbi),'%'
        endif
        if(sum(neut%full%dens).eq.0) stop 'Beam does not intersect the grid!'
    endif

end subroutine ndmc

subroutine dcx
    !+ Calculates Direct Charge Exchange (DCX) neutral density and spectra
    integer :: ic,i,j,k,ncell,is
    integer(Int64) :: idcx !! counter
    real(Float64), dimension(3) :: ri    !! start position
    real(Float64), dimension(3) :: vihalo
    integer,dimension(3) :: ind
    integer,dimension(3) :: neut_types = [1,2,3]
    !! Determination of the CX probability
    type(LocalProfiles) :: plasma
    real(Float64), dimension(nlevs) :: denn    !!  neutral dens (n=1-4)
    real(Float64), dimension(nlevs) :: rates    !!  CX rates
    !! Collisiional radiative model along track
    real(Float64), dimension(nlevs) :: states  ! Density of n-states
    integer :: ntrack
    type(ParticleTrack), dimension(beam_grid%ntrack) :: tracks  !! Particle tracks
    integer :: jj       !! counter along track
    real(Float64):: tot_denn, photons  !! photon flux
    integer, dimension(beam_grid%ngrid) :: cell_ind
    real(Float64), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz) :: papprox
    integer(Int32), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz) :: nlaunch
    real(Float64) :: fi_correction, dcx_dens

    !! Initialized Neutral Population
    call init_neutral_population(neut%dcx)

    papprox=0.d0
    tot_denn=0.d0
    do ic=1,beam_grid%ngrid
        call ind2sub(beam_grid%dims,ic,ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)
        call get_plasma(plasma,ind=ind)
        if(.not.plasma%in_plasma) cycle
        tot_denn = sum(neut%full%dens(:,i,j,k)) + &
                   sum(neut%half%dens(:,i,j,k)) + &
                   sum(neut%third%dens(:,i,j,k))
        papprox(i,j,k)= tot_denn*sum(plasma%deni)
    enddo

    ncell = 0
    do ic=1,beam_grid%ngrid
        call ind2sub(beam_grid%dims,ic,ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)
        if(papprox(i,j,k).gt.0.0) then
            ncell = ncell + 1
            cell_ind(ncell) = ic
        endif
    enddo

    call get_nlaunch(inputs%n_dcx,papprox,nlaunch)

    if(inputs%verbose.ge.1) then
       write(*,'(T6,"# of markers: ",i10)') sum(nlaunch)
    endif
    !$OMP PARALLEL DO schedule(dynamic,1) private(i,j,k,ic,is,idcx,ind,vihalo, &
    !$OMP& ri,tracks,ntrack,rates,denn,states,jj,photons,plasma,fi_correction)
    loop_over_cells: do ic = istart, ncell, istep
        call ind2sub(beam_grid%dims,cell_ind(ic),ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)
        loop_over_species: do is=1, n_thermal !This loop has to come first
            !! Loop over the markers
            loop_over_dcx: do idcx=1, nlaunch(i,j,k)
                !! Calculate ri,vhalo and track
                call mc_beam_grid(ind, ri)
                call get_plasma(plasma, pos=ri)
                call mc_halo(plasma, thermal_mass(is), vihalo)
                call track(ri,vihalo,tracks,ntrack)
                if(ntrack.eq.0) cycle loop_over_dcx

                !! Calculate CX probability
                call get_total_cx_rate(tracks(1)%ind, ri, vihalo, neut_types, rates)
                if(sum(rates).le.0.) cycle loop_over_dcx

                !! Solve collisional radiative model along track
                call get_plasma(plasma,pos=tracks(1)%pos)

                !! Weight CX rates by ion source density
                if(beam_mass.eq.thermal_mass(is)) then
                    states = rates*(plasma%deni(is) + plasma%denf)
                    if(sum(states).eq.0) cycle loop_over_dcx
                    fi_correction = max(plasma%deni(is)/(plasma%deni(is)+plasma%denf),0.d0)
                else
                    states = rates*plasma%deni(is)
                    if(sum(states).eq.0) cycle loop_over_dcx
                    fi_correction = 1.d0
                endif

                loop_along_track: do jj=1,ntrack
                    call get_plasma(plasma,pos=tracks(jj)%pos)
                    if(.not.plasma%in_plasma) exit loop_along_track
                    call colrad(plasma,thermal_mass(is),vihalo,tracks(jj)%time,states,denn,photons)
                    call store_neutrals(tracks(jj)%ind,tracks(jj)%pos,vihalo,dcx_type,denn/nlaunch(i,j,k))

                    photons = fi_correction*photons !! Correct for including fast-ions in states

                    if((photons.gt.0.d0).and.(inputs%calc_dcx.ge.1)) then
                        call store_photons(tracks(jj)%pos,vihalo, thermal_lambda0(is), photons/nlaunch(i,j,k),&
                            &spec%dcx(:,:,:,is),spec%dcxstokes(:,:,:,:,is))
                    endif

                    if((photons.gt.0.d0).and.(inputs%calc_res.ge.1)) then
                        call store_photon_birth(tracks(1)%pos, photons/nlaunch(i,j,k), spatres%dcx)
                    endif
                enddo loop_along_track
            enddo loop_over_dcx
        enddo loop_over_species
    enddo loop_over_cells
    !$OMP END PARALLEL DO

#ifdef _MPI
    !! Combine densities
    call parallel_merge_populations(neut%dcx)
    if(inputs%calc_dcx.ge.1) then
        call parallel_sum(spec%dcx)
    endif
    if(inputs%calc_res.ge.1) then
        do jj=1,spec_chords%nchan
            call parallel_merge_reservoirs(spatres%dcx(jj))
        enddo
    endif
#endif

    dcx_dens = sum(neut%dcx%dens)
    if((dcx_dens.eq.0).or.isnan(dcx_dens)) then
        write(*,*) 'DCX: DCX density is zero or nan: ', dcx_dens
        stop
    endif

end subroutine dcx

subroutine halo
    !+ Calculates halo neutral density and spectra
    integer :: ic,i,j,k,ncell
    integer(Int64) :: ihalo !! counter
    real(Float64), dimension(3) :: ri    !! start position
    real(Float64), dimension(3) :: vihalo!! velocity bulk plasma ion
    integer,dimension(3) :: ind, tind    !! actual cell
    !! Determination of the CX probability
    type(LocalProfiles) :: plasma
    real(Float64), dimension(nlevs) :: denn    !!  neutral dens (n=1-4)
    real(Float64), dimension(nlevs) :: rates    !! CX rates
    !! Collisiional radiative model along track
    real(Float64), dimension(nlevs) :: states  ! Density of n-states
    integer :: ntrack
    type(ParticleTrack), dimension(beam_grid%ntrack) :: tracks  !! Particle Tracks
    integer :: ii,jj,kk,it,is
    real(Float64) :: tot_denn, photons  !! photon flux
    integer, dimension(beam_grid%ngrid) :: cell_ind
    real(Float64), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz) :: papprox
    integer(Int32), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz) :: nlaunch
    real(Float64) :: local_iter_dens

    type(NeutralPopulation) :: cur_pop, prev_pop
    !! Halo iteration
    integer(Int64) :: hh, n_halo !! counters
    real(Float64) :: dcx_dens, halo_iter_dens(2), seed_dcx
    integer :: prev_type = 1  ! previous iteration
    integer :: cur_type = 2 ! current iteration
    real(Float64) :: fi_correction

    !! Initialize Neutral Population
    call init_neutral_population(neut%halo)

    dcx_dens = sum(neut%dcx%dens)

    halo_iter_dens(prev_type) = dcx_dens
    if(dcx_dens.eq.0) then
        if(inputs%verbose.ge.0) then
            write(*,'(a)') 'HALO: Density of DCX-neutrals is zero'
        endif
        stop
    endif
    n_halo = inputs%n_halo

    !! Allocate previous neutral populations
    call init_neutral_population(prev_pop)
    prev_pop = neut%dcx

    seed_dcx = 1.0
    iterations: do hh=1,200
        !! Allocate/Reallocate current population
        call init_neutral_population(cur_pop)
        halo_iter_dens(cur_type) = 0.d0

        !! Calculate how many mc markers to launch in each cell
        papprox = 0.d0
        cell_ind = 0
        ncell = 0
        do ic=1,beam_grid%ngrid
            call ind2sub(beam_grid%dims,ic,ind)
            i = ind(1) ; j = ind(2) ; k = ind(3)

            call get_plasma(plasma,ind=ind)
            if(.not.plasma%in_plasma) cycle
            tot_denn = sum(prev_pop%dens(:,i,j,k))
            papprox(i,j,k)= tot_denn*sum(plasma%deni)
            if(papprox(i,j,k).gt.0.0) then
                ncell = ncell + 1
                cell_ind(ncell) = ic
            endif
        enddo
        call get_nlaunch(n_halo, papprox, nlaunch)

        if(inputs%verbose.ge.1) then
            write(*,'(T6,"# of markers: ",i10," --- Seed/DCX: ",f5.3)') sum(nlaunch), seed_dcx
        endif

        !$OMP PARALLEL DO schedule(dynamic,1) private(i,j,k,ic,ihalo,ii,jj,kk,it,is,ind,vihalo, &
        !$OMP& ri,tracks,ntrack,rates,denn,states,photons,plasma,tind,fi_correction)
        loop_over_cells: do ic=istart,ncell,istep
            call ind2sub(beam_grid%dims,cell_ind(ic),ind)
            i = ind(1) ; j = ind(2) ; k = ind(3)
            loop_over_species: do is=1,n_thermal !This loop has to come first
                !! Loop over the markers
                loop_over_halos: do ihalo=1, nlaunch(i,j,k)
                    !! Calculate ri,vhalo and track
                    call mc_beam_grid(ind, ri)
                    call get_plasma(plasma, pos=ri)
                    call mc_halo(plasma, thermal_mass(is), vihalo)
                    call track(ri,vihalo,tracks,ntrack)
                    if(ntrack.eq.0)cycle loop_over_halos

                    !! Get plasma parameters at particle location
                    call get_plasma(plasma, pos=ri)

                    !! Calculate CX probability
                    tind = tracks(1)%ind
                    ii = tind(1); jj = tind(2); kk = tind(3)
                    call neutral_cx_rate(prev_pop%dens(:,ii,jj,kk), prev_pop%res(ii,jj,kk), vihalo, rates)
                    if(sum(rates).le.0.)cycle loop_over_halos

                    !! Get plasma parameters at mean point in cell
                    call get_plasma(plasma, pos=tracks(1)%pos)

                    !! Weight CX rates by ion source density
                    if(beam_mass.eq.thermal_mass(is)) then
                        states = rates*(plasma%deni(is) + plasma%denf)
                        if(sum(states).eq.0) cycle loop_over_halos
                        fi_correction = max(plasma%deni(is)/(plasma%deni(is)+plasma%denf),0.d0)
                    else
                        states = rates*plasma%deni(is)
                        if(sum(states).eq.0) cycle loop_over_halos
                        fi_correction = 1.d0
                    endif

                    loop_along_track: do it=1,ntrack
                        call get_plasma(plasma,pos=tracks(it)%pos)
                        if(.not.plasma%in_plasma) exit loop_along_track
                        call colrad(plasma,thermal_mass(is),vihalo,tracks(it)%time,states,denn,photons)

                        !! Store Neutrals
                        call update_neutrals(cur_pop, tracks(it)%ind, vihalo, denn/nlaunch(i,j,k))

                        photons = fi_correction*photons !! Correct for including fast-ions in states

                        if((photons.gt.0.d0).and.(inputs%calc_halo.ge.1)) then
                            call store_photons(tracks(it)%pos,vihalo,thermal_lambda0(is), &
                                               photons/nlaunch(i,j,k),spec%halo(:,:,:,is),spec%halostokes(:,:,:,:,is))
                        endif

                        if((photons.gt.0.d0).and.(inputs%calc_res.ge.1)) then
                            call store_photon_birth(tracks(1)%pos, photons/nlaunch(i,j,k), spatres%halo)
                        endif
                    enddo loop_along_track
                enddo loop_over_halos
            enddo loop_over_species
        enddo loop_over_cells
        !$OMP END PARALLEL DO

#ifdef _MPI
        !! Combine densities
        call parallel_merge_populations(cur_pop)
#endif
        halo_iter_dens(cur_type) = sum(cur_pop%dens)

        if(halo_iter_dens(cur_type)/halo_iter_dens(prev_type).gt.1.0) then
            write(*,'(a)') "HALO: Halo generation density exceeded seed density. This shouldn't happen."
            exit iterations
        endif

        !! Set current generation to previous generation
        halo_iter_dens(prev_type) = halo_iter_dens(cur_type)
        prev_pop = cur_pop

        !! merge current population with halo population
        call merge_neutral_populations(neut%halo, cur_pop)

        seed_dcx = halo_iter_dens(cur_type)/dcx_dens
        if(isnan(seed_dcx)) then
            write(*,'(a)') "NAN (Halo generation)/DCX. Exiting..."
            stop
        endif
        n_halo=int(inputs%n_halo*seed_dcx,Int64)

        if(seed_dcx.lt.0.01) exit iterations
    enddo iterations

    call free_neutral_population(cur_pop)
    call free_neutral_population(prev_pop)

#ifdef _MPI
    !! Combine Spectra
    if(inputs%calc_halo.ge.1) then
        call parallel_sum(spec%halo)
    endif
    if(inputs%calc_res.ge.1) then
        do jj=1,spec_chords%nchan
            call parallel_merge_reservoirs(spatres%halo(jj))
        enddo
    endif
#endif

end subroutine halo

subroutine nbi_spec
    !+ Calculates approximate neutral beam emission (full, half, third)
    !+ from user supplied neutrals file
    integer :: ic, i, j, k, it
    real(Float64), dimension(3) :: ri, vnbi, random3, rc
    integer,dimension(3) :: ind
    !! Determination of the CX probability
    real(Float64) :: nbif_photons, nbih_photons, nbit_photons
    real(Float64) :: f_wght, h_wght, t_wght
    real(Float64) :: f_tot, h_tot, t_tot
    real(Float64), dimension(n_stark,inputs%nlambda,spec_chords%nchan) :: full, half, third
    real(Float64), dimension(:,:,:,:), allocatable :: fullstokes, halfstokes, thirdstokes
    logical :: inp
    integer :: n = 10000

    !$OMP PARALLEL DO schedule(dynamic,1) private(i,j,k,ic,ind, &
    !$OMP& nbif_photons, nbih_photons, nbit_photons, rc, ri,inp, vnbi,&
    !$OMP& random3,f_tot,h_tot,t_tot,f_wght,h_wght,t_wght,&
    !$OMP& full,half,third, &
    !$OMP& fullstokes, halfstokes, thirdstokes)
    loop_over_cells: do ic = istart, spec_chords%ncell, istep
        call ind2sub(beam_grid%dims,spec_chords%cell(ic),ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)

        nbif_photons = neut%full%dens(3,i,j,k)*tables%einstein(2,3)
        nbih_photons = neut%half%dens(3,i,j,k)*tables%einstein(2,3)
        nbit_photons = neut%third%dens(3,i,j,k)*tables%einstein(2,3)
        if((nbif_photons + nbih_photons + nbit_photons).le.0.0) then
            cycle loop_over_cells
        endif

        rc = [beam_grid%xc(i), beam_grid%yc(j), beam_grid%zc(k)]

        !Find a point in cell that is also in the plasma
        ri = rc
        call in_plasma(ri, inp)
        do while (.not.inp)
            call randu(random3)
            ri = rc + beam_grid%dr*(random3 - 0.5)
            call in_plasma(ri,inp)
        enddo
        allocate(fullstokes(n_stark,4,inputs%nlambda,spec_chords%nchan))
        allocate(halfstokes(n_stark,4,inputs%nlambda,spec_chords%nchan))
        allocate(thirdstokes(n_stark,4,inputs%nlambda,spec_chords%nchan))
        f_tot = 0.0 ; h_tot = 0.0 ; t_tot = 0.0
        full  = 0.0 ; half  = 0.0 ; third = 0.0
        fullstokes  = 0.0 ; halfstokes  = 0.0 ; thirdstokes = 0.0
        do it=1, n
            !! Full Spectra
            call mc_nbi_cell(ind, nbif_type, vnbi, f_wght)
            f_tot = f_tot + f_wght
            call store_photons(ri, vnbi, beam_lambda0, f_wght*nbif_photons, full, fullstokes)

            !! Half Spectra
            call mc_nbi_cell(ind, nbih_type, vnbi, h_wght)
            h_tot = h_tot + h_wght
            call store_photons(ri, vnbi, beam_lambda0, h_wght*nbih_photons, half, halfstokes)

            !! Third Spectra
            call mc_nbi_cell(ind, nbit_type, vnbi, t_wght)
            t_tot = t_tot + t_wght
            call store_photons(ri, vnbi, beam_lambda0, t_wght*nbit_photons, third, thirdstokes)
        enddo
        !$OMP CRITICAL(nbi_spec_1)
        spec%full = spec%full + full/f_tot
        spec%half = spec%half + half/h_tot
        spec%third = spec%third + third/t_tot
        spec%fullstokes = spec%fullstokes + fullstokes/f_tot
        spec%halfstokes = spec%halfstokes + halfstokes/h_tot
        spec%thirdstokes = spec%thirdstokes + thirdstokes/t_tot
        !$OMP END CRITICAL(nbi_spec_1)
        deallocate(fullstokes)
        deallocate(halfstokes)
        deallocate(thirdstokes)
    enddo loop_over_cells
    !$OMP END PARALLEL DO

#ifdef _MPI
    !! Combine Spectra
    call parallel_sum(spec%full)
    call parallel_sum(spec%half)
    call parallel_sum(spec%third)
#endif
end subroutine nbi_spec

subroutine dcx_spec
    !+ Calculates DCX emission from user supplied neutrals file
    integer :: ic, i, j, k, it, is
    real(Float64), dimension(3) :: ri, vhalo, random3, rc
    integer,dimension(3) :: ind
    !! Determination of the CX probability
    type(LocalProfiles) :: plasma
    real(Float64) :: dcx_photons, wght
    logical :: inp
    integer :: n = 10000

    !$OMP PARALLEL DO schedule(dynamic,1) private(i,j,k,ic,is,ind, &
    !$OMP& dcx_photons, wght, rc, ri, inp, vhalo, random3, plasma)
    loop_over_cells: do ic = istart, spec_chords%ncell, istep
        call ind2sub(beam_grid%dims,spec_chords%cell(ic),ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)

        dcx_photons  = neut%dcx%dens(3,i,j,k)*tables%einstein(2,3)
        if(dcx_photons.le.0.0) cycle loop_over_cells

        rc = [beam_grid%xc(i), beam_grid%yc(j), beam_grid%zc(k)]

        !Find a point in cell that is also in the plasma
        ri = rc
        call in_plasma(ri, inp)
        do while (.not.inp)
            call randu(random3)
            ri = rc + beam_grid%dr*(random3 - 0.5)
            call in_plasma(ri,inp)
        enddo

        call get_plasma(plasma, pos=ri)
        do is=1, n_thermal
            wght = plasma%deni(is)/sum(plasma%deni)
            do it=1, n
                !! DCX Spectra
                call mc_halo(plasma, thermal_mass(is), vhalo)
                call store_photons(ri, vhalo, thermal_lambda0(is), wght*dcx_photons/n, &
                    &spec%dcx(:,:,:,is),spec%dcxstokes(:,:,:,:,is))
            enddo
        enddo
    enddo loop_over_cells
    !$OMP END PARALLEL DO

#ifdef _MPI
    !! Combine Spectra
    call parallel_sum(spec%dcx)
#endif

end subroutine dcx_spec

subroutine halo_spec
    !+ Calculates halo emission from user supplied neutrals file
    integer :: ic, i, j, k, it, is
    real(Float64), dimension(3) :: ri, vhalo, random3, rc
    integer,dimension(3) :: ind
    !! Determination of the CX probability
    type(LocalProfiles) :: plasma
    real(Float64) :: halo_photons, wght
    logical :: inp
    integer :: n = 10000

    !$OMP PARALLEL DO schedule(dynamic,1) private(i,j,k,ic,ind,is, &
    !$OMP& halo_photons, wght, rc, ri, inp, vhalo, random3, plasma)
    loop_over_cells: do ic = istart, spec_chords%ncell, istep
        call ind2sub(beam_grid%dims,spec_chords%cell(ic),ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)

        halo_photons = neut%halo%dens(3,i,j,k)*tables%einstein(2,3)
        if(halo_photons.le.0.0) cycle loop_over_cells

        rc = [beam_grid%xc(i), beam_grid%yc(j), beam_grid%zc(k)]

        !Find a point in cell that is also in the plasma
        ri = rc
        call in_plasma(ri, inp)
        do while (.not.inp)
            call randu(random3)
            ri = rc + beam_grid%dr*(random3 - 0.5)
            call in_plasma(ri,inp)
        enddo

        call get_plasma(plasma, pos=ri)
        do is=1, n_thermal
            wght = plasma%deni(is)/sum(plasma%deni)
            do it=1, n
                !! Halo Spectra
                call mc_halo(plasma, thermal_mass(is), vhalo)
                call store_photons(ri, vhalo, thermal_lambda0(is), wght*halo_photons/n, &
                    &spec%halo(:,:,:,is),spec%halostokes(:,:,:,:,is))
            enddo
        enddo
    enddo loop_over_cells
    !$OMP END PARALLEL DO

#ifdef _MPI
    !! Combine Spectra
    call parallel_sum(spec%halo)
#endif

end subroutine halo_spec

subroutine cold_spec
    !+ Calculates cold D-alpha emission
    integer :: ic, i, j, k, it, is
    real(Float64), dimension(3) :: ri, vhalo, ri_cyl
    integer,dimension(3) :: ind
    !! Determination of the CX probability
    type(LocalProfiles) :: plasma
    real(Float64) :: cold_photons
    integer :: n = 10000

    !$OMP PARALLEL DO schedule(dynamic,1) private(i,j,k,ic,is,ind, &
    !$OMP& cold_photons, ri, vhalo, plasma, ri_cyl)
    loop_over_cells: do ic = istart, spec_chords%cyl_ncell, istep
        call ind2sub(pass_grid%dims,spec_chords%cyl_cell(ic),ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)
        ri_cyl(1) = pass_grid%r(i) + 0.5d0*pass_grid%dr
        ri_cyl(2) = pass_grid%z(j) + 0.5d0*pass_grid%dz
        ri_cyl(3) = pass_grid%phi(k) + 0.5d0*pass_grid%dphi
        call cyl_to_uvw(ri_cyl, ri)

        call get_plasma(plasma, pos=ri, input_coords=1)
        do is = 1, n_thermal
            cold_photons = plasma%denn(3,is)*tables%einstein(2,3)
            if(cold_photons.le.0.0) cycle loop_over_cells

            do it=1, n
                !! Cold Spectra
                call mc_halo(plasma, thermal_mass(is), vhalo)
                call store_photons(ri, vhalo, thermal_lambda0(is), cold_photons/n, &
                    &spec%cold(:,:,:,is), spec%coldstokes(:,:,:,:,is), passive=.True.)
            enddo
        enddo
    enddo loop_over_cells
    !$OMP END PARALLEL DO

#ifdef _MPI
    !! Combine Spectra
    call parallel_sum(spec%cold)
#endif

end subroutine cold_spec

subroutine bremsstrahlung
    !+ Calculates bremsstrahlung
    type(LocalProfiles) :: plasma
    integer :: i, ichan, nc, ic
    real(Float64) :: dlength, dlambda, gaunt, max_length
    real(Float64) :: spot_size, theta, sqrt_rho
    real(Float64), dimension(2) :: randomu
    real(Float64), dimension(3) :: vi, xyz, r0
    real(Float64), dimension(3,3) :: basis
    real(Float64), dimension(:), allocatable :: lambda_arr,brems

    allocate(lambda_arr(inputs%nlambda))
    allocate(brems(inputs%nlambda))

    do i=1,inputs%nlambda
        lambda_arr(i)= 10*((i-0.5)*inputs%dlambda+inputs%lambdamin) ! [A]
    enddo
    dlambda = 10*inputs%dlambda ![A]

    dlength = 0.3 !cm
    !! $OMP PARALLEL DO schedule(guided) private(ichan,xyz,vi,basis,spot_size, &
    !! $OMP& max_length, ic, nc,randomu,sqrt_rho,theta,r0,plasma,gaunt,brems)
    loop_over_channels: do ichan=istart,spec_chords%nchan,istep
        xyz = spec_chords%los(ichan)%lens
        vi = spec_chords%los(ichan)%axis
        vi = vi/norm2(vi)
        spot_size = spec_chords%los(ichan)%spot_size
        call line_basis(xyz,vi,basis)

        if(spot_size.le.0.d0) then
            nc = 1
        else
            nc = 100
        endif

        loop_over_los: do ic=1,nc
            call randu(randomu)
            sqrt_rho = sqrt(randomu(1))
            theta = 2*pi*randomu(2)
            r0(1) = 0.d0
            r0(2) = spot_size*sqrt_rho*cos(theta)
            r0(3) = spot_size*sqrt_rho*sin(theta)
            r0 = matmul(basis,r0) + xyz

            ! Find edge of plasma
            call get_plasma(plasma,pos=r0)
            max_length=0.0
            do while (.not.plasma%in_plasma)
                r0 = r0 + vi*dlength ! move dlength
                call get_plasma(plasma,pos=r0)
                max_length = max_length + dlength
                ! If plasma is not within max_length [1000 cm] of LOS lens then no brems is calculated for the channel
                ! 1000 cm is assumed to cover the limit for most devices simulated in FIDASIM as of March 2023
                if(max_length.gt.1000) cycle loop_over_los
            enddo

            ! Calculate bremsstrahlung along los
            do while (plasma%in_plasma)
                if(plasma%te.gt.0.0) then
                    gaunt = 5.542-(3.108-log(plasma%te))*(0.6905-0.1323/plasma%zeff)
                    brems = (7.57d-9)*gaunt*((plasma%dene**2)*plasma%zeff/(lambda_arr &
                            *sqrt(plasma%te*1000.0)))*exp(-h_planck*c0/((1.d-10)*lambda_arr*plasma%te*1.d3)) &
                            *dlambda*(4.d0*pi)*1.d-4

                    spec%brems(:,ichan)= spec%brems(:,ichan) + (brems*dlength*1.d-2)/nc
                endif
                ! Take a step
                r0 = r0 + vi*dlength
                call get_plasma(plasma,pos=r0)
            enddo
        enddo loop_over_los
    enddo loop_over_channels
    !! $OMP END PARALLEL DO

#ifdef _MPI
    !! Combine Brems
    call parallel_sum(spec%brems)
#endif

    deallocate(lambda_arr,brems)

end subroutine bremsstrahlung

subroutine fida_f
    !+ Calculate Active FIDA emission using a Fast-ion distribution function F(E,p,r,z)
    integer :: i,j,k,ic,ncell
    integer(Int64) :: iion
    real(Float64), dimension(3) :: ri      !! start position
    real(Float64), dimension(3) :: vi      !! velocity of fast ions
    real(Float64) :: denf !! fast-ion density
    integer, dimension(3) :: ind      !! new actual cell
    integer, dimension(5) :: neut_types=[1,2,3,4,5]
    logical :: los_intersect
    !! Determination of the CX probability
    type(LocalEMFields) :: fields
    type(LocalProfiles) :: plasma
    real(Float64), dimension(nlevs) :: rates !! CX rates

    !! Collisiional radiative model along track
    integer :: ntrack
    integer :: jj      !! counter along track
    type(ParticleTrack),dimension(beam_grid%ntrack) :: tracks

    real(Float64) :: photons !! photon flux
    real(Float64), dimension(nlevs) :: states  !! Density of n-states
    real(Float64), dimension(nlevs) :: denn

    !! Number of particles to launch
    real(Float64) :: eb, ptch
    integer, dimension(beam_grid%ngrid) :: cell_ind
    real(Float64), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz) :: papprox
    integer(Int32), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz) :: nlaunch

    !! Estimate how many particles to launch in each cell
    papprox=0.d0
    do ic=1,beam_grid%ngrid
        call ind2sub(beam_grid%dims,ic,ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)
        call get_plasma(plasma,ind=ind)
        if(.not.plasma%in_plasma) cycle
        papprox(i,j,k) = (sum(neut%full%dens(:,i,j,k)) + &
                          sum(neut%half%dens(:,i,j,k)) + &
                          sum(neut%third%dens(:,i,j,k)) + &
                          sum(neut%dcx%dens(:,i,j,k)) + &
                          sum(neut%halo%dens(:,i,j,k)))* &
                          plasma%denf
    enddo

    ncell = 0
    do ic=1,beam_grid%ngrid
        call ind2sub(beam_grid%dims,ic,ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)
        if(papprox(i,j,k).gt.0.0) then
            ncell = ncell + 1
            cell_ind(ncell) = ic
        endif
    enddo

    call get_nlaunch(inputs%n_fida, papprox, nlaunch)
    if(inputs%verbose.ge.1) then
        write(*,'(T6,"# of markers: ",i10)') sum(nlaunch)
    endif

    !! Loop over all cells that have neutrals
    !$OMP PARALLEL DO schedule(dynamic,1) private(ic,i,j,k,ind,iion,vi,ri,fields, &
    !$OMP tracks,ntrack,jj,plasma,rates,denn,states,photons,denf,eb,ptch,los_intersect)
    loop_over_cells: do ic = istart, ncell, istep
        call ind2sub(beam_grid%dims,cell_ind(ic),ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)
        loop_over_fast_ions: do iion=1, nlaunch(i, j, k)
            !! Sample fast ion distribution for velocity and position
            call mc_fastion(ind, fields, eb, ptch, denf)
            if(denf.le.0.0) cycle loop_over_fast_ions

            !! Correct for gyro motion and get particle position and velocity
            call gyro_correction(fields, eb, ptch, fbm%A, ri, vi)

            !! Find the particles path through the beam grid
            call track(ri, vi, tracks, ntrack,los_intersect)
            if(.not.los_intersect) cycle loop_over_fast_ions
            if(ntrack.eq.0) cycle loop_over_fast_ions

            !! Calculate CX probability with beam and halo neutrals
            call get_total_cx_rate(tracks(1)%ind, ri, vi, neut_types, rates)
            if(sum(rates).le.0.) cycle loop_over_fast_ions

            !! Weight CX rates by ion source density
            states=rates*denf

            !! Calculate the spectra produced in each cell along the path
            loop_along_track: do jj=1,ntrack
                call get_plasma(plasma,pos=tracks(jj)%pos)

                call colrad(plasma, fbm%A, vi, tracks(jj)%time, states, denn, photons)

                call store_fida_photons(tracks(jj)%pos, vi, beam_lambda0, photons/nlaunch(i,j,k))
                if(inputs%calc_res.ge.1) call store_photon_birth(tracks(1)%pos, photons/nlaunch(i,j,k), spatres%fida)

            enddo loop_along_track
        enddo loop_over_fast_ions
    enddo loop_over_cells
    !$OMP END PARALLEL DO

#ifdef _MPI
    call parallel_sum(spec%fida)
    if(inputs%calc_res.ge.1) then
        do jj=1,spec_chords%nchan
            call parallel_merge_reservoirs(spatres%fida(jj))
        enddo
    endif
#endif

end subroutine fida_f

subroutine pfida_f
    !+ Calculate Passive FIDA emission using a Fast-ion distribution function F(E,p,r,z)
    integer :: i,j,k,ic,ncell,is
    integer(Int64) :: iion
    real(Float64), dimension(3) :: ri      !! start position
    real(Float64), dimension(3) :: vi      !! velocity of fast ions
    real(Float64), dimension(3) :: xyz_vi
    real(Float64) :: denf !! fast-ion density
    integer, dimension(3) :: ind      !! new actual cell
    logical :: los_intersect
    !! Determination of the CX probability
    type(LocalEMFields) :: fields
    type(LocalProfiles) :: plasma
    real(Float64), dimension(nlevs) :: rates, rates_is !! CX rates

    !! Collisiional radiative model along track
    integer :: ntrack
    integer :: jj      !! counter along track
    type(ParticleTrack),dimension(pass_grid%ntrack) :: tracks

    real(Float64) :: photons !! photon flux
    real(Float64), dimension(nlevs) :: states  !! Density of n-states
    real(Float64), dimension(nlevs) :: denn

    !! Number of particles to launch
    real(Float64) :: max_papprox, eb, ptch
    integer, dimension(pass_grid%ngrid) :: cell_ind
    real(Float64), dimension(pass_grid%nr,pass_grid%nz,pass_grid%nphi) :: papprox
    integer(Int32), dimension(pass_grid%nr,pass_grid%nz,pass_grid%nphi) :: nlaunch

    !! Estimate how many particles to launch in each cell
    papprox=0.d0
    do ic=1,pass_grid%ngrid
        call ind2sub(pass_grid%dims,ic,ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)
        call get_plasma(plasma,ind=ind,input_coords=2)
        if(.not.plasma%in_plasma) cycle
        papprox(i,j,k) = sum(plasma%denn)*plasma%denf
    enddo
    max_papprox = maxval(papprox)
    where (papprox.lt.(max_papprox*1.d-6))
        papprox = 0.0
    endwhere

    ncell = 0
    do ic=1,pass_grid%ngrid
        call ind2sub(pass_grid%dims,ic,ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)
        if(papprox(i,j,k).gt.0.0) then
            ncell = ncell + 1
            cell_ind(ncell) = ic
        endif
    enddo

    call get_nlaunch_pass_grid(inputs%n_pfida, papprox, nlaunch)
    if(inputs%verbose.ge.1) then
        write(*,'(T6,"# of markers: ",i10)') sum(nlaunch)
    endif

    !! Loop over all cells that have neutrals
    !$OMP PARALLEL DO schedule(dynamic,1) private(ic,i,j,k,ind,iion,vi,xyz_vi,ri,fields, &
    !$OMP tracks,ntrack,jj,plasma,rates,rates_is,is,denn,states,photons,denf,eb,ptch,los_intersect)
    loop_over_cells: do ic = istart, ncell, istep
        call ind2sub(pass_grid%dims,cell_ind(ic),ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)
        loop_over_fast_ions: do iion=1, nlaunch(i, j, k)
            !! Sample fast ion distribution for velocity and position
            call mc_fastion_pass_grid(ind, fields, eb, ptch, denf, output_coords=1)
            if(denf.le.0.0) cycle loop_over_fast_ions

            !! Correct for gyro motion and get particle position and velocity
            call gyro_correction(fields, eb, ptch, fbm%A, ri, vi)
            xyz_vi = matmul(beam_grid%inv_basis,vi)

            !! Find the particles path through the interpolation grid
            call track_cylindrical(ri, vi, tracks, ntrack,los_intersect)
            if(.not.los_intersect) cycle loop_over_fast_ions
            if(ntrack.eq.0) cycle loop_over_fast_ions

            !! Calculate CX probability with edge neutrals
            call get_plasma(plasma, pos=ri, input_coords=1)
            rates = 0.0
            do is=1,n_thermal
                call bt_cx_rates(plasma, plasma%denn(:,is), thermal_mass(is), xyz_vi, rates_is)
                rates = rates + rates_is
            enddo
            if(sum(rates).le.0.) cycle loop_over_fast_ions

            !! Weight CX rates by ion source density
            states=rates*denf

            !! Calculate the spectra produced in each cell along the path
            loop_along_track: do jj=1,ntrack
                call get_plasma(plasma,pos=tracks(jj)%pos,input_coords=1)

                call colrad(plasma, fbm%A, xyz_vi, tracks(jj)%time, states, denn, photons)

                call store_fida_photons(tracks(jj)%pos, xyz_vi, beam_lambda0, photons/nlaunch(i,j,k), passive=.True.)
                if(inputs%calc_res.ge.1) then 
                    call store_photon_birth(tracks(1)%pos, photons/nlaunch(i,j,k), spatres%pfida, passive=.True.)
                endif
            enddo loop_along_track
        enddo loop_over_fast_ions
    enddo loop_over_cells
    !$OMP END PARALLEL DO

#ifdef _MPI
    call parallel_sum(spec%pfida)
    if(inputs%calc_res.ge.1) then
        do jj=1,spec_chords%nchan
            call parallel_merge_reservoirs(spatres%pfida(jj))
        enddo
    endif
#endif

end subroutine pfida_f

subroutine fida_mc
    !+ Calculate Active FIDA emission using a Monte Carlo Fast-ion distribution
    integer :: iion,igamma,ngamma
    type(FastIon) :: fast_ion
    type(LocalEMFields) :: fields
    type(LocalProfiles) :: plasma
    real(Float64) :: phi
    real(Float64), dimension(3) :: ri      !! start position
    real(Float64), dimension(3) :: vi      !! velocity of fast ions
    !! Determination of the CX probability
    real(Float64), dimension(nlevs) :: denn    !!  neutral dens (n=1-4)
    real(Float64), dimension(nlevs) :: rates    !! CX rates
    !! Collisiional radiative model along track
    real(Float64), dimension(nlevs) :: states  ! Density of n-states
    integer :: ntrack
    type(ParticleTrack), dimension(beam_grid%ntrack) :: tracks
    logical :: los_intersect
    integer :: jj      !! counter along track
    real(Float64) :: photons !! photon flux
    integer, dimension(5) :: neut_types=[1,2,3,4,5]
    real(Float64), dimension(3) :: uvw, uvw_vi
    real(Float64)  :: s, c
    real(Float64), dimension(1) :: randomu

    ngamma = 1
    if(particles%axisym.or.(inputs%dist_type.eq.2)) then
        ngamma = ceiling(dble(inputs%n_fida)/particles%nparticle)
    endif

    if(inputs%verbose.ge.1) then
        write(*,'(T6,"# of markers: ",i10)') int(particles%nparticle*ngamma,Int64)
    endif

    !$OMP PARALLEL DO schedule(dynamic,1) private(iion,igamma,fast_ion,vi,ri,phi,tracks,s,c, &
    !$OMP& randomu,plasma,fields,uvw,uvw_vi,ntrack,jj,rates,denn,los_intersect,states,photons)
    loop_over_fast_ions: do iion=istart,particles%nparticle,istep
        fast_ion = particles%fast_ion(iion)
        if(fast_ion%vabs.eq.0) cycle loop_over_fast_ions
        if(.not.fast_ion%beam_grid_cross_grid) cycle loop_over_fast_ions
        gamma_loop: do igamma=1,ngamma
            if(particles%axisym) then
                !! Pick random toroidal angle
                call randu(randomu)
                phi = fast_ion%beam_grid_phi_enter + fast_ion%delta_phi*randomu(1)
            else
                phi = fast_ion%phi
            endif
            s = sin(phi)
            c = cos(phi)

            !! Calculate position in machine coordinates
            uvw(1) = fast_ion%r*c
            uvw(2) = fast_ion%r*s
            uvw(3) = fast_ion%z

            !! Convert to beam grid coordinates
            call uvw_to_xyz(uvw, ri)

            if(inputs%dist_type.eq.2) then
                !! Get electomagnetic fields
                call get_fields(fields, pos=ri)

                !! Correct for gyro motion and get particle position and velocity
                call gyro_correction(fields, fast_ion%energy, fast_ion%pitch, fast_ion%A, ri, vi)
            else !! Full Orbit
                !! Calculate velocity vector
                uvw_vi(1) = c*fast_ion%vr - s*fast_ion%vt
                uvw_vi(2) = s*fast_ion%vr + c*fast_ion%vt
                uvw_vi(3) = fast_ion%vz
                vi = matmul(beam_grid%inv_basis,uvw_vi)
            endif

            !! Track particle through grid
            call track(ri, vi, tracks, ntrack, los_intersect)
            if(.not.los_intersect) cycle gamma_loop
            if(ntrack.eq.0)cycle gamma_loop

            !! Calculate CX probability
            call get_total_cx_rate(tracks(1)%ind, ri, vi, neut_types, rates)
            if(sum(rates).le.0.)cycle gamma_loop

            !! Weight CX rates by ion source density
            states=rates*fast_ion%weight*(fast_ion%delta_phi/(2*pi))/beam_grid%dv/ngamma

            !! Calculate the spectra produced in each cell along the path
            loop_along_track: do jj=1,ntrack
                call get_plasma(plasma,pos=tracks(jj)%pos)

                call colrad(plasma, fast_ion%A, vi, tracks(jj)%time, states, denn, photons)

                call store_fida_photons(tracks(jj)%pos, vi, beam_lambda0, photons, fast_ion%class)
                if(inputs%calc_res.ge.1) call store_photon_birth(tracks(1)%pos, photons, spatres%fida)
            enddo loop_along_track
        enddo gamma_loop
    enddo loop_over_fast_ions
    !$OMP END PARALLEL DO

#ifdef _MPI
    call parallel_sum(spec%fida)
    if(inputs%calc_res.ge.1) then
        do jj=1,spec_chords%nchan
            call parallel_merge_reservoirs(spatres%fida(jj))
        enddo
    endif
#endif

end subroutine fida_mc

subroutine pfida_mc
    !+ Calculate Passive FIDA emission using a Monte Carlo Fast-ion distribution
    integer :: iion,igamma,ngamma
    type(FastIon) :: fast_ion
    type(LocalEMFields) :: fields
    type(LocalProfiles) :: plasma
    real(Float64) :: phi
    real(Float64), dimension(3) :: ri      !! start position
    real(Float64), dimension(3) :: vi      !! velocity of fast ions
    real(Float64), dimension(3) :: xyz_vi
    !! Determination of the CX probability
    real(Float64), dimension(nlevs) :: denn    !!  neutral dens (n=1-4)
    real(Float64), dimension(nlevs) :: rates, rates_is    !! CX rates
    !! Collisiional radiative model along track
    real(Float64), dimension(nlevs) :: states  ! Density of n-states
    type(ParticleTrack), dimension(pass_grid%ntrack) :: tracks
    integer :: ntrack
    logical :: los_intersect
    integer :: jj, is
    real(Float64) :: photons
    real(Float64)  :: s, c
    real(Float64), dimension(1) :: randomu

    ngamma = 1
    if(particles%axisym.or.(inputs%dist_type.eq.2)) then
        ngamma = ceiling(dble(inputs%n_pfida)/particles%nparticle)
    endif

    if(inputs%verbose.ge.1) then
        write(*,'(T6,"# of markers: ",i10)') int(particles%nparticle*ngamma,Int64)
    endif

    !$OMP PARALLEL DO schedule(dynamic,1) private(iion,igamma,is,fast_ion,vi,ri,phi,tracks,s,c,&
    !$OMP& randomu,plasma,fields,ntrack,jj,rates,rates_is,denn,los_intersect,states,photons,xyz_vi)
    loop_over_fast_ions: do iion=istart,particles%nparticle,istep
        fast_ion = particles%fast_ion(iion)
        if(fast_ion%vabs.eq.0) cycle loop_over_fast_ions
        gamma_loop: do igamma=1,ngamma
            if(particles%axisym) then
                !! Pick random toroidal angle
                call randu(randomu)
                phi = pass_grid%phi(1)+pass_grid%nphi*pass_grid%dphi*randomu(1)
            else
                phi = fast_ion%phi
            endif
            s = sin(phi)
            c = cos(phi)

            !! Calculate position in machine coordinates
            ri(1) = fast_ion%r*c
            ri(2) = fast_ion%r*s
            ri(3) = fast_ion%z

            if(inputs%dist_type.eq.2) then
                !! Get electomagnetic fields
                call get_fields(fields, pos=ri, input_coords=1, output_coords=1)

                !! Correct for gyro motion and get particle position and velocity
                call gyro_correction(fields, fast_ion%energy, fast_ion%pitch, fast_ion%A, ri, vi)
            else !! Full Orbit
                !! Calculate velocity vector
                vi(1) = c*fast_ion%vr - s*fast_ion%vt
                vi(2) = s*fast_ion%vr + c*fast_ion%vt
                vi(3) = fast_ion%vz
            endif
            xyz_vi = matmul(beam_grid%inv_basis,vi)

            !! Track particle through grid
            call track_cylindrical(ri, vi, tracks, ntrack, los_intersect)
            if(.not.los_intersect) cycle gamma_loop
            if(ntrack.eq.0) cycle gamma_loop

            !! Calculate CX probability
            call get_plasma(plasma, pos=ri, input_coords=1)
            rates = 0.d0
            do is=1,n_thermal
                call bt_cx_rates(plasma, plasma%denn(:,is), thermal_mass(is), xyz_vi, rates_is)
                rates = rates + rates_is
            enddo
            if(sum(rates).le.0.) cycle gamma_loop

            !! Weight CX rates by ion source density
            if(particles%axisym) then
                states=rates*fast_ion%weight*(pass_grid%nphi*pass_grid%dphi/(2*pi))  &
                       /(fast_ion%r*pass_grid%dv)/ngamma

            else
                states=rates*fast_ion%weight/(fast_ion%r*pass_grid%dv)/ngamma
            endif

            !! Calculate the spectra produced in each cell along the path
            loop_along_track: do jj=1,ntrack
                call get_plasma(plasma,pos=tracks(jj)%pos,input_coords=1)

                call colrad(plasma, fast_ion%A, xyz_vi, tracks(jj)%time, states, denn, photons)

                call store_fida_photons(tracks(jj)%pos, xyz_vi, beam_lambda0, photons, fast_ion%class,passive=.True.)
                if(inputs%calc_res.ge.1) call store_photon_birth(tracks(1)%pos, photons, spatres%pfida, passive = .True.)
            enddo loop_along_track
        enddo gamma_loop
    enddo loop_over_fast_ions
    !$OMP END PARALLEL DO

#ifdef _MPI
    call parallel_sum(spec%pfida)
    if(inputs%calc_res.ge.1) then
        do jj=1,spec_chords%nchan
            call parallel_merge_reservoirs(spatres%pfida(jj))
        enddo
    endif
#endif

end subroutine pfida_mc

subroutine npa_f
    !+ Calculate Active NPA flux using a fast-ion distribution function F(E,p,r,z)
    integer :: i,j,k,det,ic
    integer(Int64) :: iion
    real(Float64), dimension(3) :: rg,ri,rf,vi
    integer, dimension(3) :: ind,pind
    real(Float64) :: denf
    type(LocalProfiles) :: plasma
    type(LocalEMFields) :: fields
    type(GyroSurface) :: gs
    real(Float64), dimension(2,4) :: gyrange
    integer, dimension(5) :: neut_types=[1,2,3,4,5]
    real(Float64), dimension(nlevs) :: rates
    real(Float64), dimension(nlevs) :: states
    real(Float64) :: flux, theta, dtheta, eb, ptch

    integer :: inpa,ichan,nrange,ir,npart,ncell
    integer, dimension(beam_grid%ngrid) :: cell_ind
    real(Float64), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz) :: papprox
    integer(Int32), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz) :: nlaunch

    papprox=0.d0
    do ic=1,beam_grid%ngrid
        call ind2sub(beam_grid%dims,ic,ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)
        call get_plasma(plasma,ind=ind)
        if(.not.plasma%in_plasma) cycle
        papprox(i,j,k)=(sum(neut%full%dens(:,i,j,k)) + &
                        sum(neut%half%dens(:,i,j,k)) + &
                        sum(neut%third%dens(:,i,j,k)) + &
                        sum(neut%dcx%dens(:,i,j,k)) + &
                        sum(neut%halo%dens(:,i,j,k)))* &
                        plasma%denf
    enddo

    ncell = 0
    do ic=1,beam_grid%ngrid
        call ind2sub(beam_grid%dims,ic,ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)
        if(papprox(i,j,k).gt.0.0) then
            ncell = ncell + 1
            cell_ind(ncell) = ic
        endif
    enddo

    call get_nlaunch(inputs%n_npa, papprox, nlaunch)
    if(inputs%verbose.ge.1) then
        write(*,'(T6,"# of markers: ",i12)') sum(nlaunch)
    endif

    !! Loop over all cells that can contribute to NPA signal
    !$OMP PARALLEL DO schedule(dynamic,1) private(ic,i,j,k,ind,iion,ichan,fields,nrange,gyrange, &
    !$OMP& pind,vi,ri,rf,det,plasma,rates,states,flux,denf,eb,ptch,gs,ir,theta,dtheta)
    loop_over_cells: do ic = istart, ncell, istep
        call ind2sub(beam_grid%dims,cell_ind(ic),ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)
        loop_over_fast_ions: do iion=1, nlaunch(i, j, k)
            !! Sample fast ion distribution for energy and pitch
            call mc_fastion(ind, fields, eb, ptch, denf)
            if(denf.le.0.0) cycle loop_over_fast_ions

            call gyro_surface(fields, eb, ptch, fbm%A, gs)

            detector_loop: do ichan=1,npa_chords%nchan
                call npa_gyro_range(ichan, gs, gyrange, nrange)
                if(nrange.eq.0) cycle detector_loop
                gyro_range_loop: do ir=1,nrange
                    dtheta = gyrange(2,ir)
                    theta = gyrange(1,ir) + 0.5*dtheta
                    call gyro_trajectory(gs, theta, ri, vi)

                    !! Check if particle hits a NPA detector
                    call hit_npa_detector(ri, vi ,det, rf, ichan)
                    if(det.ne.ichan) then
                        if (inputs%verbose.ge.0)then
                            write(*,*) "NPA_F: Missed Detector ",ichan
                        endif
                        cycle gyro_range_loop
                    endif

                    !! Get beam grid indices at ri
                    call get_indices(ri,pind)

                    !! Calculate CX probability with beam and halo neutrals
                    call get_total_cx_rate(pind, ri, vi, neut_types, rates)
                    if(sum(rates).le.0.) cycle gyro_range_loop

                    !! Weight CX rates by ion source density
                    states=rates*denf

                    !! Attenuate states as the particle move through plasma
                    call attenuate(ri,rf,vi,states)

                    !! Store NPA Flux
                    flux = (dtheta/(2*pi))*sum(states)*beam_grid%dv/nlaunch(i,j,k)
                    call store_npa(det,ri,rf,vi,flux)
                enddo gyro_range_loop
            enddo detector_loop
        enddo loop_over_fast_ions
    enddo loop_over_cells
    !$OMP END PARALLEL DO

    npart = npa%npart
#ifdef _MPI
    call parallel_sum(npart)
    call parallel_sum(npa%flux)
#endif

    if(inputs%verbose.ge.1) then
        write(*,'(T4,"Number of Active NPA particles that hit a detector: ",i8)') npart
    endif

end subroutine npa_f

subroutine pnpa_f
    !+ Calculate Passive NPA flux using a fast-ion distribution function F(E,p,r,z)
    integer :: i,j,k,det,ic
    integer(Int64) :: iion
    real(Float64), dimension(3) :: rg,ri,rf,vi,ri_uvw
    integer, dimension(3) :: ind
    real(Float64) :: denf,r
    type(LocalProfiles) :: plasma
    type(LocalEMFields) :: fields
    type(GyroSurface) :: gs
    real(Float64), dimension(2,4) :: gyrange
    real(Float64), dimension(nlevs) :: rates, rates_is
    real(Float64), dimension(nlevs) :: states
    real(Float64) :: flux, theta, dtheta, eb, ptch,max_papprox

    integer :: inpa,ichan,nrange,ir,npart,ncell, is
    integer, dimension(pass_grid%ngrid) :: cell_ind
    real(Float64), dimension(pass_grid%nr,pass_grid%nz,pass_grid%nphi) :: papprox
    integer(Int32), dimension(pass_grid%nr,pass_grid%nz,pass_grid%nphi) :: nlaunch

    papprox=0.d0
    do ic=1,pass_grid%ngrid
        call ind2sub(pass_grid%dims,ic,ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)
        call get_plasma(plasma,ind=ind,input_coords=2)
        if(.not.plasma%in_plasma) cycle
        papprox(i,j,k) = sum(plasma%denn)*plasma%denf
    enddo
    max_papprox = maxval(papprox)
    where (papprox.lt.(max_papprox*1.d-6))
        papprox = 0.0
    endwhere

    ncell = 0
    do ic=1,pass_grid%ngrid
        call ind2sub(pass_grid%dims,ic,ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)
        if(papprox(i,j,k).gt.0.0) then
            ncell = ncell + 1
            cell_ind(ncell) = ic
        endif
    enddo

    call get_nlaunch_pass_grid(inputs%n_pnpa, papprox, nlaunch)
    if(inputs%verbose.ge.1) then
        write(*,'(T6,"# of markers: ",i12)') sum(nlaunch)
    endif

    !! Loop over all cells that can contribute to NPA signal
    !$OMP PARALLEL DO schedule(dynamic,1) private(ic,i,j,k,ind,iion,ichan,fields,nrange,gyrange, &
    !$OMP& vi,ri,rf,det,plasma,rates,is, rates_is,states,flux,denf,eb,ptch,gs,ir,theta,dtheta,r,ri_uvw)
    loop_over_cells: do ic = istart, ncell, istep
        call ind2sub(pass_grid%dims,cell_ind(ic),ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)
        loop_over_fast_ions: do iion=1, nlaunch(i, j, k)
            !! Sample fast ion distribution for energy and pitch
            call mc_fastion_pass_grid(ind, fields, eb, ptch, denf)
            if(denf.le.0.0) cycle loop_over_fast_ions

            call gyro_surface(fields, eb, ptch, fbm%A, gs)

            detector_loop: do ichan=1,npa_chords%nchan
                call npa_gyro_range(ichan, gs, gyrange, nrange)
                if(nrange.eq.0) cycle detector_loop
                gyro_range_loop: do ir=1,nrange
                    dtheta = gyrange(2,ir)
                    theta = gyrange(1,ir) + 0.5*dtheta
                    call gyro_trajectory(gs, theta, ri, vi)

                    !! Check if particle hits a NPA detector
                    call hit_npa_detector(ri, vi ,det, rf, ichan)
                    if(det.ne.ichan) then
                        if (inputs%verbose.ge.0)then
                            write(*,*) "PNPA_F: Missed Detector ",ichan
                        endif
                        cycle gyro_range_loop
                    endif

                    !! Calculate CX probability with beam and halo neutrals
                    call get_plasma(plasma, pos=ri)
                    rates = 0.d0
                    do is=1,n_thermal
                        call bt_cx_rates(plasma, plasma%denn(:,is), thermal_mass(is), vi, rates_is)
                        rates = rates + rates_is
                    enddo
                    if(sum(rates).le.0.) cycle gyro_range_loop

                    !! Weight CX rates by ion source density
                    states=rates*denf

                    !! Attenuate states as the particle move through plasma
                    call attenuate(ri,rf,vi,states)

                    !! Store NPA Flux
                    flux = (dtheta/(2*pi))*sum(states)*pass_grid%r(i)*pass_grid%dv/nlaunch(i,j,k)
                    call store_npa(det,ri,rf,vi,flux,passive=.True.)
                enddo gyro_range_loop
            enddo detector_loop
        enddo loop_over_fast_ions
    enddo loop_over_cells
    !$OMP END PARALLEL DO

    npart = pnpa%npart
#ifdef _MPI
    call parallel_sum(npart)
    call parallel_sum(pnpa%flux)
#endif

    if(inputs%verbose.ge.1) then
        write(*,'(T4,"Number of Passive NPA particles that hit a detector: ",i8)') npart
    endif

end subroutine pnpa_f

subroutine npa_mc
    !+ Calculate Active NPA flux using a Monte Carlo fast-ion distribution
    integer :: iion,igamma,ngamma,npart
    type(FastIon) :: fast_ion
    real(Float64) :: phi,theta,dtheta
    real(Float64), dimension(3) :: ri, rf, rg, vi
    integer :: det,ichan,ir,nrange,it
    type(LocalEMFields) :: fields
    type(GyroSurface) :: gs
    real(Float64), dimension(nlevs) :: rates
    real(Float64), dimension(nlevs) :: states
    real(Float64) :: flux
    integer, dimension(5) :: neut_types=[1,2,3,4,5]
    integer, dimension(3) :: ind
    real(Float64), dimension(3) :: uvw, uvw_vi
    real(Float64), dimension(2,4) :: gyrange
    real(Float64) :: s,c
    real(Float64), dimension(1) :: randomu

    ngamma = 1
    if(particles%axisym.or.(inputs%dist_type.eq.2)) then
        ngamma = ceiling(dble(inputs%n_npa)/particles%nparticle)
    endif

    if(inputs%verbose.ge.1) then
        write(*,'(T6,"# of markers: ",i10)') int(particles%nparticle*ngamma,Int64)
    endif

    !$OMP PARALLEL DO schedule(guided) private(iion,igamma,ind,fast_ion,vi,ri,rf,phi,s,c,ir,it, &
    !$OMP& randomu,rg,fields,uvw,uvw_vi,rates,states,flux,det,ichan,gs,nrange,gyrange,theta,dtheta)
    loop_over_fast_ions: do iion=istart,particles%nparticle,istep
        fast_ion = particles%fast_ion(iion)
        if(fast_ion%vabs.eq.0) cycle loop_over_fast_ions
        if(.not.fast_ion%beam_grid_cross_grid) cycle loop_over_fast_ions
        gamma_loop: do igamma=1,ngamma
            if(particles%axisym) then
                !! Pick random toroidal angle
                call randu(randomu)
                phi = fast_ion%beam_grid_phi_enter + fast_ion%delta_phi*randomu(1)
            else
                phi = fast_ion%phi
            endif
            s = sin(phi)
            c = cos(phi)

            !! Calculate position in machine coordinates
            uvw(1) = fast_ion%r*c
            uvw(2) = fast_ion%r*s
            uvw(3) = fast_ion%z

            if(inputs%dist_type.eq.2) then
                !! Convert to beam grid coordinates
                call uvw_to_xyz(uvw, rg)

                !! Get electomagnetic fields
                call get_fields(fields, pos=rg)

                !! Correct for gyro motion and get position and velocity
                call gyro_surface(fields, fast_ion%energy, fast_ion%pitch, fast_ion%A, gs)

                detector_loop: do ichan=1,npa_chords%nchan
                    call npa_gyro_range(ichan, gs, gyrange, nrange)
                    if(nrange.eq.0) cycle detector_loop

                    gyro_range_loop: do ir=1,nrange
                        dtheta = gyrange(2,ir)
                        theta = gyrange(1,ir) + 0.5*dtheta
                        call gyro_trajectory(gs, theta, ri, vi)

                        !! Check if particle hits a NPA detector
                        call hit_npa_detector(ri, vi ,det, rf, det=ichan)
                        if(det.ne.ichan) then
                            if (inputs%verbose.ge.0)then
                                write(*,*) "NPA_MC: Missed Detector ",ichan
                            endif
                            cycle gyro_range_loop
                        endif

                        !! Get beam grid indices at ri
                        call get_indices(ri,ind)

                        !! Calculate CX probability with beam and halo neutrals
                        call get_total_cx_rate(ind, ri, vi, neut_types, rates)
                        if(sum(rates).le.0.) cycle gyro_range_loop

                        !! Weight CX rates by ion source density
                        states=rates*fast_ion%weight*(fast_ion%delta_phi/(2*pi))/beam_grid%dv/ngamma

                        !! Attenuate states as the particle move through plasma
                        call attenuate(ri,rf,vi,states)

                        !! Store NPA Flux
                        flux = (dtheta/(2*pi))*sum(states)*beam_grid%dv
                        spread_loop: do it=1,25
                            theta = gyrange(1,ir) + (it-0.5)*dtheta/25
                            call gyro_trajectory(gs, theta, ri, vi)
                            call hit_npa_detector(ri, vi ,det, rf, det=ichan)
                            if(det.ne.ichan) then
                                if (inputs%verbose.ge.0)then
                                    write(*,*) "NPA_MC: Missed Detector ",ichan
                                endif
                                cycle spread_loop
                            endif
                            call store_npa(det,ri,rf,vi,flux/25,fast_ion%class)
                        enddo spread_loop
                    enddo gyro_range_loop
                enddo detector_loop
            else !! Full Orbit
                !! Convert to beam grid coordinates
                call uvw_to_xyz(uvw, ri)

                !! Calculate velocity vector
                uvw_vi(1) = c*fast_ion%vr - s*fast_ion%vt
                uvw_vi(2) = s*fast_ion%vr + c*fast_ion%vt
                uvw_vi(3) = fast_ion%vz
                vi = matmul(beam_grid%inv_basis,uvw_vi)

                !! Check if particle hits a NPA detector
                call hit_npa_detector(ri, vi ,det, rf)
                if(det.eq.0) cycle gamma_loop

                !! Get beam grid indices at ri
                call get_indices(ri,ind)

                !! Calculate CX probability with beam and halo neutrals
                call get_total_cx_rate(ind, ri, vi, neut_types, rates)
                if(sum(rates).le.0.) cycle gamma_loop

                !! Weight CX rates by ion source density
                states=rates*fast_ion%weight*(fast_ion%delta_phi/(2*pi))/beam_grid%dv/ngamma

                !! Attenuate states as the particle moves though plasma
                call attenuate(ri,rf,vi,states)

                !! Store NPA Flux
                flux = sum(states)*beam_grid%dv
                call store_npa(det,ri,rf,vi,flux,fast_ion%class)
            endif
        enddo gamma_loop
    enddo loop_over_fast_ions
    !$OMP END PARALLEL DO

    npart = npa%npart
#ifdef _MPI
    call parallel_sum(npart)
    call parallel_sum(npa%flux)
#endif

    if(inputs%verbose.ge.1) then
        write(*,'(T4,"Number of Active NPA particles that hit a detector: ",i8)') npart
    endif

end subroutine npa_mc

subroutine pnpa_mc
    !+ Calculate Passive NPA flux using a Monte Carlo fast-ion distribution
    integer :: iion,igamma,ngamma,npart
    type(FastIon) :: fast_ion
    real(Float64) :: phi,theta,dtheta
    real(Float64), dimension(3) :: ri, rf, rg, vi
    integer :: det,j,ichan,ir,nrange,it,is
    type(LocalEMFields) :: fields
    type(LocalProfiles) :: plasma
    type(GyroSurface) :: gs
    real(Float64), dimension(nlevs) :: rates, rates_is
    real(Float64), dimension(nlevs) :: states
    real(Float64) :: flux
    integer, dimension(3) :: ind
    real(Float64), dimension(3) :: uvw, uvw_vi
    real(Float64), dimension(2,4) :: gyrange
    real(Float64) :: s,c
    real(Float64), dimension(1) :: randomu

    ngamma = 1
    if(particles%axisym.or.(inputs%dist_type.eq.2)) then
        ngamma = ceiling(dble(inputs%n_pnpa)/particles%nparticle)
    endif

    if(inputs%verbose.ge.1) then
        write(*,'(T6,"# of markers: ",i10)') int(particles%nparticle*ngamma,Int64)
    endif

    !$OMP PARALLEL DO schedule(guided) private(iion,igamma,ind,fast_ion,vi,ri,rf,phi,s,c,ir,it,plasma, &
    !$OMP& randomu,rg,fields,uvw,uvw_vi,rates,rates_is,is,states,flux,det,ichan,gs,nrange,gyrange,theta,dtheta)
    loop_over_fast_ions: do iion=istart,particles%nparticle,istep
        fast_ion = particles%fast_ion(iion)
        if(fast_ion%vabs.eq.0) cycle loop_over_fast_ions
        gamma_loop: do igamma=1,ngamma
            if(particles%axisym) then
                !! Pick random toroidal angle
                call randu(randomu)
                phi = pass_grid%phi(1)+pass_grid%nphi*pass_grid%dphi*randomu(1)
            else
                phi = fast_ion%phi
            endif
            s = sin(phi)
            c = cos(phi)

            !! Calculate position in machine coordinates
            uvw(1) = fast_ion%r*c
            uvw(2) = fast_ion%r*s
            uvw(3) = fast_ion%z

            if(inputs%dist_type.eq.2) then
                !! Get electomagnetic fields
                call get_fields(fields, pos=uvw, input_coords=1)

                !! Correct for gyro motion and get position and velocity
                call gyro_surface(fields, fast_ion%energy, fast_ion%pitch, fast_ion%A, gs)

                detector_loop: do ichan=1,npa_chords%nchan
                    call npa_gyro_range(ichan, gs, gyrange, nrange)
                    if(nrange.eq.0) cycle detector_loop

                    gyro_range_loop: do ir=1,nrange
                        dtheta = gyrange(2,ir)
                        theta = gyrange(1,ir) + 0.5*dtheta
                        call gyro_trajectory(gs, theta, ri, vi)

                        !! Check if particle hits a NPA detector
                        call hit_npa_detector(ri, vi ,det, rf, det=ichan)
                        if(det.ne.ichan) then
                            if (inputs%verbose.ge.0)then
                                write(*,*) "PNPA_MC: Missed Detector ",ichan
                            endif
                            cycle gyro_range_loop
                        endif

                        !! Calculate CX probability with beam and halo neutrals
                        call get_plasma(plasma, pos=ri)
                        rates = 0.d0
                        do is=1,n_thermal
                            call bt_cx_rates(plasma, plasma%denn(:,is), thermal_mass(is), vi, rates_is)
                            rates = rates + rates_is
                        enddo
                        if(sum(rates).le.0.) cycle gyro_range_loop

                        !! Weight CX rates by ion source density
                        if(particles%axisym) then
                            states=rates*fast_ion%weight*(pass_grid%nphi*pass_grid%dphi/(2*pi)) &
                                   /(fast_ion%r*pass_grid%dv)/ngamma
                        else
                            states=rates*fast_ion%weight/(fast_ion%r*pass_grid%dv)/ngamma
                        endif

                        !! Attenuate states as the particle move through plasma
                        call attenuate(ri,rf,vi,states)

                        !! Store NPA Flux
                        flux = (dtheta/(2*pi))*sum(states)*(fast_ion%r*pass_grid%dv)
                        spread_loop: do it=1,25
                            theta = gyrange(1,ir) + (it-0.5)*dtheta/25
                            call gyro_trajectory(gs, theta, ri, vi)
                            call hit_npa_detector(ri, vi ,det, rf, det=ichan)
                            if(det.ne.ichan) then
                                if (inputs%verbose.ge.0)then
                                    write(*,*) "PNPA_MC: Missed Detector ",ichan
                                endif
                                cycle spread_loop
                            endif
                            call store_npa(det,ri,rf,vi,flux/25,fast_ion%class, passive=.True.)
                        enddo spread_loop
                    enddo gyro_range_loop
                enddo detector_loop
            else !! Full Orbit
                !! Convert to beam grid coordinates
                call uvw_to_xyz(uvw, ri)

                !! Calculate velocity vector
                uvw_vi(1) = c*fast_ion%vr - s*fast_ion%vt
                uvw_vi(2) = s*fast_ion%vr + c*fast_ion%vt
                uvw_vi(3) = fast_ion%vz
                vi = matmul(beam_grid%inv_basis,uvw_vi)

                !! Check if particle hits a NPA detector
                call hit_npa_detector(ri, vi ,det, rf)
                if(det.eq.0) cycle gamma_loop

                !! Calculate CX probability with edge neutrals
                call get_plasma(plasma, pos=ri)
                rates = 0.d0
                do is=1, n_thermal
                    call bt_cx_rates(plasma, plasma%denn(:,is), thermal_mass(is), vi, rates_is)
                    rates = rates + rates_is
                enddo
                if(sum(rates).le.0.) cycle gamma_loop

                !! Weight CX rates by ion source density
                if(particles%axisym) then
                    states=rates*fast_ion%weight*(pass_grid%nphi*pass_grid%dphi/(2*pi)) &
                           /(fast_ion%r*pass_grid%dv)/ngamma
                else
                    states=rates*fast_ion%weight/(fast_ion%r*pass_grid%dv)/ngamma
                endif

                !! Attenuate states as the particle moves though plasma
                call attenuate(ri,rf,vi,states)

                !! Store NPA Flux
                flux = sum(states)*(fast_ion%r*pass_grid%dv)
                call store_npa(det,ri,rf,vi,flux,fast_ion%class,passive=.True.)
            endif
        enddo gamma_loop
    enddo loop_over_fast_ions
    !$OMP END PARALLEL DO

    npart = pnpa%npart
#ifdef _MPI
    call parallel_sum(npart)
    call parallel_sum(pnpa%flux)
#endif

    if(inputs%verbose.ge.1) then
        write(*,'(T4,"Number of Passive NPA particles that hit a detector: ",i8)') npart
    endif

end subroutine pnpa_mc

subroutine neutron_f
    !+ Calculate neutron emission rate using a fast-ion distribution function F(E,p,r,z)
    integer :: ir, iphi, iz, ie, ip, igamma, ngamma
    type(LocalProfiles) :: plasma
    type(LocalEMFields) :: fields
    real(Float64) :: eb,pitch
    real(Float64) :: erel, rate
    real(Float64), dimension(3) :: ri
    real(Float64), dimension(3) :: vi
    real(Float64), dimension(3) :: uvw, uvw_vi
    real(Float64)  :: vnet_square, factor
    real(Float64)  :: s, c

    if(.not.any(thermal_mass.eq.H2_amu)) then
        write(*,'(T2,a)') 'NEUTRON_F: Thermal Deuterium is not present in plasma'
        return
    endif
    if(any(thermal_mass.eq.H3_amu)) then
        write(*,'(T2,a)') 'NEUTRON_F: D-T neutron production is not implemented'
    endif
    if(beam_mass.ne.H2_amu) then
        write(*,'(T2,a)') 'NEUTRON_F: Fast-ion species is not Deuterium'
        return
    endif

    if(inputs%calc_neutron.ge.2) then
        allocate(neutron%weight(fbm%nenergy,fbm%npitch,fbm%nr,fbm%nz,fbm%nphi))
        neutron%weight = 0.d0
        allocate(neutron%emis(fbm%nr,fbm%nz,fbm%nphi))
        neutron%emis = 0.d0
    endif

    ngamma = 20
    rate = 0
    !$OMP PARALLEL DO schedule(guided) private(fields,vi,ri,pitch,eb,&
    !$OMP& ir,iphi,iz,ie,ip,igamma,plasma,factor,uvw,uvw_vi,vnet_square,rate,erel,s,c)
    z_loop: do iz = istart, fbm%nz, istep
        r_loop: do ir=1, fbm%nr
            phi_loop: do iphi = 1, fbm%nphi
                !! Calculate position
                if(fbm%nphi.eq.1) then
                    s = 0.d0
                    c = 1.d0
                else
                    s = sin(fbm%phi(iphi))
                    c = cos(fbm%phi(iphi))
                endif

                uvw(1) = fbm%r(ir)*c
                uvw(2) = fbm%r(ir)*s
                uvw(3) = fbm%z(iz)

                !! Get fields
                call get_fields(fields,pos=uvw,input_coords=1)
                if(.not.fields%in_plasma) cycle r_loop

                factor = fbm%r(ir)*fbm%dE*fbm%dp*fbm%dr*fbm%dz*fbm%dphi/ngamma
                !! Loop over energy/pitch/gamma
                pitch_loop: do ip = 1, fbm%npitch
                    pitch = fbm%pitch(ip)
                    energy_loop: do ie =1, fbm%nenergy
                        eb = fbm%energy(ie)
                        gyro_loop: do igamma=1, ngamma
                            call gyro_correction(fields,eb,pitch,fbm%A,ri,vi)

                            !! Get plasma parameters at particle position
                            call get_plasma(plasma,pos=ri)
                            if(.not.plasma%in_plasma) cycle gyro_loop

                            !! Calculate effective beam energy
                            vnet_square=dot_product(vi-plasma%vrot,vi-plasma%vrot)  ![cm/s]
                            erel = v2_to_E_per_amu*fbm%A*vnet_square ![kev]

                            !! Get neutron production rate
                            call get_dd_rate(plasma, erel, rate, branch=2)
                            if(inputs%calc_neutron.ge.2) then
                                neutron%weight(ie,ip,ir,iz,iphi) = neutron%weight(ie,ip,ir,iz,iphi) &
                                                                 + rate * factor
                                !$OMP CRITICAL(neutron_emis)
                                neutron%emis(ir,iz,iphi) = neutron%emis(ir,iz,iphi) &
                                                                 + rate * fbm%f(ie,ip,ir,iz,iphi) &
                                                                 * factor
                                !$OMP END CRITICAL(neutron_emis)
                            endif

                            rate = rate*fbm%f(ie,ip,ir,iz,iphi)*factor
                            !! Store neutrons
                            call store_neutrons(rate)
                        enddo gyro_loop
                    enddo energy_loop
                enddo pitch_loop
            enddo phi_loop
        enddo r_loop
    enddo z_loop
    !$OMP END PARALLEL DO

#ifdef _MPI
    call parallel_sum(neutron%rate)
    if(inputs%calc_neutron.ge.2) then
        call parallel_sum(neutron%weight)
        call parallel_sum(neutron%emis)
    endif
#endif

    if(inputs%verbose.ge.1) then
        write(*,'(T4,A,ES14.5," [neutrons/s]")') 'Rate:   ',sum(neutron%rate)
        write(*,'(30X,a)') ''
        write(*,*) 'write neutrons:    ' , time_string(time_start)
    endif

#ifdef _MPI
    if(my_rank().eq.0) call write_neutrons()
#else
    call write_neutrons()
#endif


end subroutine neutron_f

subroutine cfpd_f
    !+ Calculate charged fusion product count rate and weight function using a fast-ion distribution function F(E,p,r,z)
    real(Float64), dimension(3) :: vi, vi_norm, v3_xyz, xyz, r_gyro
    type(LocalProfiles) :: plasma
    type(LocalEMFields) :: fields
    real(Float64) :: pgyro, vnet_square, factor, vabs
    real(Float64) :: eb, pitch, erel, rate, kappa, gyro, fbm_denf
    integer :: ie, ip, ich, ie3, iray, ist, cnt

    if(.not.any(thermal_mass.eq.H2_amu)) then
        write(*,'(T2,a)') 'CFPD_F: Thermal Deuterium is not present in plasma'
        return
    endif
    if(any(thermal_mass.eq.H3_amu)) then
        write(*,'(T2,a)') 'CFPD_F: D-T cfpd production is not implemented'
    endif
    if(beam_mass.ne.H2_amu) then
        write(*,'(T2,a)') 'CFPD_F: Fast-ion species is not Deuterium'
        return
    endif

    allocate(cfpd%flux(ctable%nenergy, ctable%nchan))
    allocate(cfpd%prob(ctable%nenergy, ctable%nchan))
    allocate(cfpd%gam(ctable%nenergy, ctable%nchan))
    allocate(cfpd%weight(ctable%nenergy, ctable%nchan, fbm%nenergy, fbm%npitch))
    cfpd%flux = 0.d0
    cfpd%prob = 0.d0
    cfpd%gam = 0.d0
    cfpd%weight = 0.d0

    rate = 0.d0
    factor = 0.5d0*fbm%dE*fbm%dp*ctable%dl !0.5 for TRANSP-pitch (E,p) space factor
    !$OMP PARALLEL DO schedule(guided) private(vi,vi_norm,v3_xyz,xyz,r_gyro,plasma,fields,pgyro,&
    !$OMP& vnet_square,vabs,eb,pitch,erel,rate,kappa,gyro,fbm_denf,ie,ip,ich,ie3,iray,ist,cnt)
    channel_loop: do ich=1, ctable%nchan
        E3_loop: do ie3=1, ctable%nenergy
            cnt = 0
            ray_loop: do iray=1, ctable%nrays
                step_loop: do ist=1, ctable%nsteps
                    if (ist.gt.ctable%nactual(ie3,iray,ich)) cycle ray_loop

                    !! Calculate position and velocity in beam coordinates
                    call convert_sightline_to_xyz(ie3, ist, iray, ich, xyz, v3_xyz)

                    !! Get fields at sightline position
                    call get_fields(fields, pos=xyz)
                    if(.not.fields%in_plasma) cycle step_loop

                    !! Get plasma parameters at sightline position
                    call get_plasma(plasma, pos=xyz)
                    if(.not.plasma%in_plasma) cycle step_loop

                    !! Loop over energy/pitch/gamma
                    pitch_loop: do ip = 1, fbm%npitch
                        pitch = fbm%pitch(ip)
                        energy_loop: do ie =1, fbm%nenergy
                            eb = fbm%energy(ie)

                            !! Get the probability factor
                            call get_pgyro(fields,ctable%earray(ie3),eb,pitch,plasma,v3_xyz,pgyro,gyro)
                            if (pgyro.le.0.d0) cycle energy_loop
                            cnt = cnt + 1

                            !! Calculate fast-ion velocity
                            call pitch_to_vec(pitch, gyro, fields, vi_norm)
                            vabs = sqrt(eb/(v2_to_E_per_amu*fbm%A))
                            vi = vi_norm*vabs
                            !!Correct for gyro orbit
                            call gyro_step(vi,fields,fbm%A,r_gyro)

                            fbm_denf=0
                            if (inputs%dist_type.eq.1) then
                                !get F at guiding center position
                                call get_ep_denf(eb,pitch,fbm_denf,pos=(xyz+r_gyro))
                            endif
                            if (fbm_denf.ne.fbm_denf) cycle energy_loop

                            !! Calculate effective beam energy
                            vnet_square=dot_product(vi-plasma%vrot,vi-plasma%vrot)  ![cm/s]
                            erel = v2_to_E_per_amu*fbm%A*vnet_square ![kev]

                            !! Get the cfpd production rate and anisotropy term
                            call get_dd_rate(plasma, erel, rate, branch=1)
                            call get_ddpt_anisotropy(plasma, vi, v3_xyz, kappa)

                            !$OMP CRITICAL(cfpd_weight)
                            cfpd%weight(ie3,ich,ie,ip) = cfpd%weight(ie3,ich,ie,ip) &
                                                            + rate * kappa * pgyro &
                                                            * ctable%daomega(ie3,iray,ich) &
                                                            * factor / (fbm%dE*fbm%dp)
                            cfpd%prob(ie3,ich) = cfpd%prob(ie3,ich) + pgyro
                            cfpd%gam(ie3,ich) = cfpd%gam(ie3,ich) + gyro
                            cfpd%flux(ie3,ich) = cfpd%flux(ie3,ich) + rate * kappa * pgyro &
                                                                        * ctable%daomega(ie3,iray,ich) &
                                                                        * fbm_denf * factor
                            !$OMP END CRITICAL(cfpd_weight)

                        enddo energy_loop
                    enddo pitch_loop
                enddo step_loop
            enddo ray_loop
            !$OMP CRITICAL(pweight_cnt)
            cfpd%prob(ie3,ich) = cfpd%prob(ie3,ich) / cnt
            cfpd%gam(ie3,ich) = cfpd%gam(ie3,ich) / cnt
            !$OMP END CRITICAL(pweight_cnt)
        enddo E3_loop
    enddo channel_loop
    !$OMP END PARALLEL DO

#ifdef _MPI
    call parallel_sum(cfpd%flux)
    call parallel_sum(cfpd%weight)
    call parallel_sum(cfpd%prob)
    call parallel_sum(cfpd%gam)
#endif

    if(inputs%verbose.ge.1) then
        write(*,'(30X,a)') ''
        write(*,*) 'write charged fusion products:    ' , time_string(time_start)
    endif

#ifdef _MPI
    if(my_rank().eq.0) call write_cfpd_weights()
#else
    call write_cfpd_weights()
#endif

end subroutine cfpd_f

subroutine neutron_mc
    !+ Calculate neutron flux using a Monte Carlo Fast-ion distribution
    integer :: iion, ngamma, igamma
    type(FastIon) :: fast_ion
    type(LocalProfiles) :: plasma
    type(LocalEMFields) :: fields
    real(Float64) :: eb, rate
    real(Float64), dimension(3) :: ri
    real(Float64), dimension(3) :: vi
    real(Float64), dimension(3) :: uvw, uvw_vi
    real(Float64)  :: vnet_square
    real(Float64)  :: phi, s, c, factor, delta_phi

    if(.not.any(thermal_mass.eq.H2_amu)) then
        write(*,'(T2,a)') 'NEUTRON_MC: Thermal Deuterium is not present in plasma'
        return
    endif
    if(any(thermal_mass.eq.H3_amu)) then
        write(*,'(T2,a)') 'NEUTRON_MC: D-T neutron production is not implemented'
    endif
    if(beam_mass.ne.H2_amu) then
        write(*,'(T2,a)') 'NEUTRON_MC: Fast-ion species is not Deuterium'
        return
    endif


    if(inputs%verbose.ge.1) then
        write(*,'(T6,"# of markers: ",i10)') particles%nparticle
    endif

    !! Correct neutron rate when equilibrium is 3D and MC distribution is 4D
    if(particles%axisym.and.(inter_grid%nphi.gt.1)) then
        delta_phi = inter_grid%phi(inter_grid%nphi)-inter_grid%phi(1)
        delta_phi = delta_phi + delta_phi/(inter_grid%nphi-1)/2 !Add half a cell
        factor = delta_phi/(2*pi) * 2 !Riemann sum below assumes coord's are at midpoint of cell
    else
        factor = 1
    endif

    rate=0.0
    ngamma = 20
    !$OMP PARALLEL DO schedule(guided) private(iion,fast_ion,vi,ri,s,c, &
    !$OMP& plasma,fields,uvw,uvw_vi,vnet_square,rate,eb,igamma,phi)
    loop_over_fast_ions: do iion=istart,particles%nparticle,istep
        fast_ion = particles%fast_ion(iion)
        if(fast_ion%vabs.eq.0.d0) cycle loop_over_fast_ions

        !! Calculate position in machine coordinates
        if(particles%axisym) then
            s = 0.d0
            c = 1.d0
        else
            phi = fast_ion%phi
            s = sin(phi)
            c = cos(phi)
        endif

        uvw(1) = fast_ion%r*c
        uvw(2) = fast_ion%r*s
        uvw(3) = fast_ion%z

        if(inputs%dist_type.eq.2) then
            !! Get electomagnetic fields
            call get_fields(fields, pos=uvw, input_coords=1)
            if(.not.fields%in_plasma) cycle loop_over_fast_ions

            gyro_loop: do igamma=1,ngamma
                !! Correct for Gyro-motion
                call gyro_correction(fields, fast_ion%energy, fast_ion%pitch, fast_ion%A, ri, vi)

                !! Get plasma parameters
                call get_plasma(plasma,pos=ri)
                if(.not.plasma%in_plasma) cycle gyro_loop

                !! Calculate effective beam energy
                vnet_square=dot_product(vi-plasma%vrot,vi-plasma%vrot)  ![cm/s]
                eb = v2_to_E_per_amu*fast_ion%A*vnet_square ![kev]

                !! Get neutron production rate
                call get_dd_rate(plasma, eb, rate, branch=2)
                rate = rate*fast_ion%weight/ngamma*factor

                !! Store neutrons
                call store_neutrons(rate, fast_ion%class)
            enddo gyro_loop
        else
            !! Get plasma parameters
            call get_plasma(plasma,pos=uvw,input_coords=1)
            if(.not.plasma%in_plasma) cycle loop_over_fast_ions

            !! Calculate effective beam energy
            uvw_vi(1) = fast_ion%vr
            uvw_vi(2) = fast_ion%vt
            uvw_vi(3) = fast_ion%vz
            vi = matmul(beam_grid%inv_basis,uvw_vi)
            vnet_square=dot_product(vi-plasma%vrot,vi-plasma%vrot)  ![cm/s]
            eb = v2_to_E_per_amu*fast_ion%A*vnet_square ![kev]

            !! Get neutron production rate
            call get_dd_rate(plasma, eb, rate, branch=2)
            rate = rate*fast_ion%weight*factor

            !! Store neutrons
            call store_neutrons(rate, fast_ion%class)
        endif
    enddo loop_over_fast_ions
    !$OMP END PARALLEL DO

#ifdef _MPI
    call parallel_sum(neutron%rate)
#endif

    if(inputs%verbose.ge.1) then
        write(*,'(T4,A,ES14.5," [neutrons/s]")') 'Rate:   ',sum(neutron%rate)
        write(*,'(30X,a)') ''
        write(*,*) 'write neutrons:    ' , time_string(time_start)
    endif

#ifdef _MPI
    if(my_rank().eq.0) call write_neutrons()
#else
    call write_neutrons()
#endif

end subroutine neutron_mc

subroutine fida_weights_mc
    !+ Calculates FIDA weights
    integer :: i,j,k,ic,ncell
    integer(Int64) :: iion,ip
    real(Float64), dimension(3) :: ri,rg      !! start position
    real(Float64), dimension(3) :: vi     !! velocity of fast ions
    integer, dimension(3) :: ind      !! new actual cell
    integer,dimension(5) :: neut_types=[1,2,3,4,5]
    logical :: los_intersect

    !! Determination of the CX rates
    type(LocalProfiles) :: plasma
    type(LocalEMFields) :: fields
    real(Float64), dimension(nlevs) :: rates !! CX rates

    !! Collisiional radiative model along track
    integer :: ntrack
    integer :: jj      !! counter along track
    type(ParticleTrack),dimension(beam_grid%ntrack) :: tracks

    real(Float64) :: photons !! photon flux
    real(Float64), dimension(nlevs) :: states  !! Density of n-states
    real(Float64), dimension(nlevs) :: denn

    integer :: nwav
    real(Float64) :: etov2, energy, pitch
    real(Float64) :: dE, dP, dEdP
    real(Float64), dimension(:), allocatable :: ebarr, ptcharr
    integer, dimension(1) :: ienergy, ipitch
    real(Float64), dimension(3) :: randomu3

    !! Number of particles to launch
    real(Float64) :: fbm_denf,phase_area
    integer, dimension(beam_grid%ngrid) :: cell_ind
    real(Float64), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz) :: papprox
    integer(Int32), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz) :: nlaunch

    nwav = inputs%nlambda_wght

    !! define arrays
    !! define energy - array
    allocate(ebarr(inputs%ne_wght))
    do i=1,inputs%ne_wght
        ebarr(i)=real(i-0.5)*inputs%emax_wght/real(inputs%ne_wght)
    enddo
    dE = abs(ebarr(2)-ebarr(1))

    !! define pitch - array
    allocate(ptcharr(inputs%np_wght))
    do i=1,inputs%np_wght
        ptcharr(i)=real(i-0.5)*2./real(inputs%np_wght)-1.
    enddo
    dP = abs(ptcharr(2)-ptcharr(1))
    dEdP = dE*dP
    phase_area = dEdP*real(inputs%np_wght)*real(inputs%ne_wght)

    !! allocate storage arrays
    allocate(fweight%weight(nwav,inputs%ne_wght,inputs%np_wght,spec_chords%nchan))
    allocate(fweight%mean_f(inputs%ne_wght,inputs%np_wght,spec_chords%nchan))

    if(inputs%verbose.ge.1) then
        write(*,'(T3,"Number of Channels: ",i5)') spec_chords%nchan
        write(*,'(T3,"Nlambda: ",i4)') nwav
        write(*,'(T3,"Nenergy: ",i3)') inputs%ne_wght
        write(*,'(T3,"Maximum Energy: ",f7.2)') inputs%emax_wght
        write(*,'(T3,"LOS averaged: ",a)') "False"
    endif

    !! zero out arrays
    fweight%weight = 0.d0
    fweight%mean_f = 0.d0

    etov2 = 1.d0/(v2_to_E_per_amu*beam_mass)

    !! Estimate how many particles to launch in each cell
    papprox=0.d0
    do ic=1,beam_grid%ngrid
        call ind2sub(beam_grid%dims,ic,ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)
        call get_plasma(plasma,ind=ind)
        if(.not.plasma%in_plasma) cycle
        papprox(i,j,k)=(sum(neut%full%dens(:,i,j,k)) + &
                        sum(neut%half%dens(:,i,j,k)) + &
                        sum(neut%third%dens(:,i,j,k)) + &
                        sum(neut%dcx%dens(:,i,j,k)) + &
                        sum(neut%halo%dens(:,i,j,k)))
    enddo

    ncell = 0
    do ic=1,beam_grid%ngrid
        call ind2sub(beam_grid%dims,ic,ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)
        if(papprox(i,j,k).gt.0.0) then
            ncell = ncell + 1
            cell_ind(ncell) = ic
        endif
    enddo

    call get_nlaunch(10*inputs%n_fida,papprox, nlaunch)
    if(inputs%verbose.ge.1) then
        write(*,'(T6,"# of markers: ",i10)') sum(nlaunch)
    endif

    !! Loop over all cells that have neutrals
    !$OMP PARALLEL DO schedule(guided) private(ic,i,j,k,ind,iion,vi,ri,rg,ienergy,ipitch, &
    !$OMP tracks,ntrack,jj,plasma,fields,rates,denn,states,photons,energy,pitch, &
    !$OMP los_intersect,randomu3,fbm_denf)
    loop_over_cells: do ic = istart, ncell, istep
        call ind2sub(beam_grid%dims,cell_ind(ic),ind)
        i = ind(1) ; j = ind(2) ; k = ind(3)
        loop_over_fast_ions: do iion=1, nlaunch(i, j, k)
            !! Sample fast ion distribution uniformally
            call randind(inputs%ne_wght, ienergy)
            call randind(inputs%np_wght, ipitch)
            call randu(randomu3)
            energy = ebarr(ienergy(1)) + dE*(randomu3(1)-0.5)
            pitch = ptcharr(ipitch(1)) + dP*(randomu3(2)-0.5)
            if(energy.le.0) cycle loop_over_fast_ions

            call randu(randomu3)
            rg = [beam_grid%xc(i),beam_grid%yc(j),beam_grid%zc(k)] + beam_grid%dr*(randomu3-0.5)

            !! Get velocity
            call get_fields(fields,pos=rg)
            if(.not.fields%in_plasma) cycle loop_over_fast_ions
            call gyro_correction(fields,energy,pitch,fbm%A,ri,vi)

            fbm_denf = 0.0
            if (inputs%dist_type.eq.1) then
                call get_ep_denf(energy,pitch,fbm_denf,coeffs=fields%b)
            endif

            !! Find the particles path through the beam grid
            call track(ri, vi, tracks, ntrack, los_intersect)
            if(.not.los_intersect) cycle loop_over_fast_ions
            if(ntrack.eq.0) cycle loop_over_fast_ions

            !! Calculate CX probability with beam and halo neutrals
            call get_total_cx_rate(tracks(1)%ind, ri, vi, neut_types, rates)
            if(sum(rates).le.0.) cycle loop_over_fast_ions
            states=rates*1.d20

            !! Calculate the spectra produced in each cell along the path
            loop_along_track: do jj=1,ntrack
                call get_plasma(plasma,pos=tracks(jj)%pos)

                call colrad(plasma, beam_mass, vi, tracks(jj)%time, states, denn, photons)

                call store_fw_photons(ienergy(1), ipitch(1), &
                     tracks(jj)%pos, vi, beam_lambda0, fbm_denf, photons/nlaunch(i,j,k))
            enddo loop_along_track
        enddo loop_over_fast_ions
    enddo loop_over_cells
    !$OMP END PARALLEL DO

    fweight%weight = ((1.d-20)*phase_area/dEdP)*fweight%weight
    fweight%mean_f = ((1.d-20)*phase_area/dEdP)*fweight%mean_f

    if(inputs%verbose.ge.1) then
        write(*,*) 'write fida weights:    ' , time_string(time_start)
    endif

#ifdef _MPI
    call parallel_sum(fweight%weight)
    call parallel_sum(fweight%mean_f)
    if(my_rank().eq.0) call write_fida_weights()
#else
    call write_fida_weights()
#endif

end subroutine fida_weights_mc

subroutine fida_weights_los
    !+ Calculates LOS averaged FIDA weights
    type(LocalProfiles) :: plasma, plasma_cell
    type(LocalEMFields) :: fields, fields_cell
    real(Float64) :: denf
    real(Float64) :: wght, wght_tot
    real(Float64) :: photons !! photon flux
    real(Float64) :: length
    type(ParticleTrack), dimension(beam_grid%ntrack) :: tracks
    integer :: nwav
    integer(Int32) :: i, j, k, ienergy, cid, cind
    integer(Int32) :: ipitch, igyro, icell, ichan, is
    real(Float64), dimension(:), allocatable :: ebarr,ptcharr,phiarr
    real(Float64), dimension(:,:), allocatable :: mean_f
    real(Float64), dimension(3) :: vi, vi_norm, vp
    real(Float64), dimension(3) :: vnbi_f, vnbi_h, vnbi_t, vhalo
    real(Float64), dimension(3) :: r_enter, r_exit
    real(Float64) :: vabs, dE, dP
    !! Determination of the CX probability
    real(Float64), dimension(nlevs) :: fdens,hdens,tdens,dcxdens,halodens
    real(Float64), dimension(nlevs) :: rates
    real(Float64), dimension(nlevs) :: states ! Density of n-states
    real(Float64), dimension(nlevs) :: denn  ! Density of n-states
    !! COLRAD
    real(Float64) :: dt, max_dens, dlength, sigma_pi
    type(LOSInters) :: inter
    real(Float64) :: eb, ptch, phi
    !! Solution of differential equation
    integer, dimension(3) :: ind  !!actual cell
    real(Float64), dimension(3) :: ri
    integer(Int32) :: ntrack
    logical :: inp

    real(Float64):: etov2, dEdP

    nwav = inputs%nlambda_wght

    !! Define energy array
    allocate(ebarr(inputs%ne_wght))
    do i=1,inputs%ne_wght
        ebarr(i)=real(i-0.5)*inputs%emax_wght/real(inputs%ne_wght)
    enddo
    dE = abs(ebarr(2)-ebarr(1))

    !! Define pitch array
    allocate(ptcharr(inputs%np_wght))
    do i=1,inputs%np_wght
        ptcharr(i)=real(i-0.5)*2./real(inputs%np_wght)-1.
    enddo
    dP = abs(ptcharr(2)-ptcharr(1))
    dEdP = dE*dP

    !! define gyro - array
    allocate(phiarr(inputs%nphi_wght))
    do i=1,inputs%nphi_wght
        phiarr(i)=real(i-0.5)*2.d0*pi/real(inputs%nphi_wght)
    enddo

    !! allocate storage arrays
    allocate(fweight%mean_f(inputs%ne_wght,inputs%np_wght,spec_chords%nchan))
    allocate(fweight%weight(nwav,inputs%ne_wght,inputs%np_wght,spec_chords%nchan))

    allocate(mean_f(inputs%ne_wght,inputs%np_wght))
    !! zero out arrays
    fweight%weight = 0.d0
    fweight%mean_f = 0.d0
    mean_f = 0.d0

    if(inputs%verbose.ge.1) then
        write(*,'(T3,"Number of Channels: ",i5)') spec_chords%nchan
        write(*,'(T3,"Nlambda: ",i4)') nwav
        write(*,'(T3,"Nenergy: ",i3)') inputs%ne_wght
        write(*,'(T3,"Npitch: ",i3)') inputs%np_wght
        write(*,'(T3,"Ngyro: ", i3)') inputs%nphi_wght
        write(*,'(T3,"Maximum Energy: ",f7.2)') inputs%emax_wght
        write(*,'(T3,"LOS averaged: ",a)') "True"
        write(*,*) ''
    endif

    etov2 = 1.0/(v2_to_E_per_amu*beam_mass)

    chan_loop: do ichan=1,spec_chords%nchan
        fdens = 0.d0 ; hdens = 0.d0 ; tdens = 0.d0
        halodens = 0.d0 ; dcxdens = 0.d0
        plasma = plasma*0.d0
        fields = fields*0.d0
        wght_tot = 0.d0
        mean_f = 0.d0
        do k=1,beam_grid%nz
            do j=1,beam_grid%ny
                x_loop: do i=1,beam_grid%nx
                    inter = spec_chords%inter(i,j,k)
                    cid = 0
                    cind = 0
                    do while (cid.ne.ichan.and.cind.lt.inter%nchan)
                        cind = cind + 1
                        cid = inter%los_elem(cind)%id
                    enddo
                    if(cid.eq.ichan) then
                        ind = [i,j,k]
                        ri = [beam_grid%xc(i), beam_grid%yc(j), beam_grid%zc(k)]
                        call in_plasma(ri,inp)
                        if(.not.inp) cycle x_loop
                        dlength = inter%los_elem(cind)%length
                        fdens = fdens + neut%full%dens(:,i,j,k)*dlength
                        hdens = hdens + neut%half%dens(:,i,j,k)*dlength
                        tdens = tdens + neut%third%dens(:,i,j,k)*dlength
                        dcxdens = dcxdens + neut%dcx%dens(:,i,j,k)*dlength
                        halodens = halodens + neut%halo%dens(:,i,j,k)*dlength
                        wght = (neut%full%dens(3,i,j,k) + neut%half%dens(3,i,j,k) + &
                                neut%third%dens(3,i,j,k) + neut%dcx%dens(3,i,j,k) + &
                                neut%halo%dens(3,i,j,k))*dlength

                        call get_plasma(plasma_cell,pos=ri)
                        call get_fields(fields_cell,pos=ri)
                        plasma = plasma + wght*plasma_cell
                        fields = fields + wght*fields_cell
                        if (inputs%dist_type.eq.1) then
                            do ipitch=1,inputs%np_wght
                                do ienergy=1,inputs%ne_wght
                                    call get_ep_denf(ebarr(ienergy),ptcharr(ipitch),denf,coeffs=fields_cell%b)
                                    mean_f(ienergy,ipitch) = mean_f(ienergy,ipitch) + wght*denf
                                enddo
                            enddo
                        endif
                        wght_tot = wght_tot + wght
                    endif
                enddo x_loop
            enddo
        enddo

        if(wght_tot.le.0) then
            if(inputs%verbose.ge.1) then
                write(*,'(T4,"Skipping channel ",i5,": Neutral density is zero")') ichan
            endif
            cycle chan_loop
        else
            plasma = plasma/wght_tot
            plasma%in_plasma = .True.
            fields = fields/wght_tot
            fields%in_plasma= .True.
            mean_f = mean_f/wght_tot
            if(inputs%verbose.ge.1) then
                write(*,'(T4,"Channel: ",i5)') ichan
                write(*,'(T4,"Radius: ",f7.2)') spec_chords%radius(ichan)
                write(*,'(T4,"Mean Fast-ion Density: ",ES14.5)') sum(mean_f)*dEdP
                write(*,*)''
            endif
        endif

        ri = plasma%pos
        vp = ri - spec_chords%los(ichan)%lens
        vnbi_f = ri - nbi%src
        vnbi_f = vnbi_f/norm2(vnbi_f)*nbi%vinj
        vnbi_h = vnbi_f/sqrt(2.d0)
        vnbi_t = vnbi_f/sqrt(3.d0)
        sigma_pi = spec_chords%los(ichan)%sigma_pi
        dlength = 1.d0

        !$OMP PARALLEL DO schedule(guided) collapse(3) private(eb,vabs,ptch,phi,vi,vi_norm, &
        !$OMP& r_enter,r_exit,length,max_dens,ind,tracks,ntrack,dt,icell,states,rates, &
        !$OMP& vhalo,denn,denf,photons,ienergy,ipitch,igyro,is)
        do ienergy=istart,inputs%ne_wght,istep
            do ipitch=1,inputs%np_wght
                do igyro=1,inputs%nphi_wght
                    eb = ebarr(ienergy)
                    vabs = sqrt(eb*etov2)
                    ptch = ptcharr(ipitch)
                    phi = phiarr(igyro)
                    call pitch_to_vec(ptch,phi,fields,vi_norm)
                    vi = vabs*vi_norm

                    call grid_intersect(ri,vi,length,r_enter,r_exit)
                    call track(r_enter, vi, tracks, ntrack)
                    max_dens = 0.d0
                    do icell=1,ntrack
                        ind = tracks(icell)%ind
                        tracks(icell)%flux = neut%full%dens(3,ind(1),ind(2),ind(3))  + &
                                             neut%half%dens(3,ind(1),ind(2),ind(3))  + &
                                             neut%third%dens(3,ind(1),ind(2),ind(3)) + &
                                             neut%dcx%dens(3,ind(1),ind(2),ind(3))   + &
                                             neut%halo%dens(3,ind(1),ind(2),ind(3))
                        if(tracks(icell)%flux.gt.max_dens) max_dens=tracks(icell)%flux
                    enddo
                    dt = 0.d0
                    do icell=1,ntrack
                        if(tracks(icell)%flux.gt.(0.5*max_dens)) then
                            dt = dt + tracks(icell)%time
                        endif
                    enddo

                    states=0.d0
                    call bb_cx_rates(fdens,vnbi_f,vi,rates)
                    states = states + rates
                    call bb_cx_rates(hdens,vnbi_h,vi,rates)
                    states = states + rates
                    call bb_cx_rates(tdens,vnbi_t,vi,rates)
                    states = states + rates
                    do is=1,n_thermal
                        call bt_cx_rates(plasma, dcxdens + halodens, thermal_mass(is), vi, rates)
                        states = states + rates*plasma%deni(is)/sum(plasma%deni)
                    enddo

                    call colrad(plasma,beam_mass,vi,dt,states,denn,photons)
                    denf = mean_f(ienergy,ipitch)*dEdP
                    photons = photons/real(inputs%nphi_wght)
                    call store_fw_photons_at_chan(ichan, ienergy, ipitch, &
                         vp, vi, beam_lambda0, fields, dlength, sigma_pi, denf, photons)

                enddo
            enddo
        enddo
        !$OMP END PARALLEL DO

    enddo chan_loop

    fweight%mean_f = fweight%mean_f/(dEdP)

    if(inputs%verbose.ge.1) then
        write(*,*) 'write fida weights:    ' , time_string(time_start)
    endif

#ifdef _MPI
    call parallel_sum(fweight%weight)
    call parallel_sum(fweight%mean_f)
    if(my_rank().eq.0) call write_fida_weights()
#else
    call write_fida_weights()
#endif

end subroutine fida_weights_los

subroutine npa_weights
    !+ Calculates NPA weights
    type(LocalEMFields) :: fields
    type(NPAProbability) :: phit
    real(Float64) :: pitch
    real(Float64) :: pcxa
    integer(Int32) :: det
    integer(Int32) :: ii, jj, kk, i, ic   !!indices
    integer,dimension(1) :: ipitch
    real(Float64), dimension(3) :: vi,vi_norm
    real(Float64) :: vabs, fbm_denf, dE, dP
    real(Float64), dimension(nlevs) :: pcx   !! Rate coefficiants for CX
    real(Float64), dimension(nlevs) :: states, states_i  ! Density of n-states
    integer, dimension(5) :: neut_types=[1,2,3,4,5]
    real(Float64), dimension(3) :: pos,dpos,r_gyro
    integer(Int32) :: ichan
    real(Float64), dimension(:), allocatable :: ebarr, ptcharr

    !! define energy - array
    allocate(ebarr(inputs%ne_wght))
    do i=1,inputs%ne_wght
        ebarr(i)=real(i-0.5)*inputs%emax_wght/real(inputs%ne_wght)
    enddo
    dE=abs(ebarr(2)-ebarr(1))

    !! define pitch - array
    allocate(ptcharr(inputs%np_wght))
    do i=1,inputs%np_wght
        ptcharr(i)=real(i-0.5)*2./real(inputs%np_wght)-1.
    enddo
    dP=abs(ptcharr(2)-ptcharr(1))

    if(inputs%verbose.ge.1) then
        write(*,'(T3,"Number of Channels: ",i3)') npa_chords%nchan
        write(*,'(T3,"Nenergy: ",i3)') inputs%ne_wght
        write(*,'(T3,"Npitch: ",i3)') inputs%np_wght
        write(*,'(T3,"Maximum energy: ",f7.2)') inputs%emax_wght
        write(*,*) ''
    endif

    !! define storage arrays
    allocate(nweight%emissivity(beam_grid%nx, &
                                beam_grid%ny, &
                                beam_grid%nz, &
                                npa_chords%nchan))

    allocate(nweight%attenuation(inputs%ne_wght, &
                                 beam_grid%nx, &
                                 beam_grid%ny, &
                                 beam_grid%nz, &
                                 npa_chords%nchan))

    allocate(nweight%cx(inputs%ne_wght, &
                        beam_grid%nx, &
                        beam_grid%ny, &
                        beam_grid%nz, &
                        npa_chords%nchan))

    allocate(nweight%weight(inputs%ne_wght, &
                            inputs%np_wght, &
                            npa_chords%nchan))

    allocate(nweight%flux(inputs%ne_wght, npa_chords%nchan))

    nweight%emissivity = 0.d0
    nweight%attenuation = 0.d0
    nweight%cx = 0.d0
    nweight%weight = 0.d0
    nweight%flux = 0.d0

    loop_over_channels: do ichan=1,npa_chords%nchan
        !$OMP PARALLEL DO schedule(guided) collapse(3) private(ii,jj,kk,fields,phit,&
        !$OMP& ic,det,pos,dpos,r_gyro,pitch,ipitch,vabs,vi,pcx,pcxa,states,states_i,vi_norm,fbm_denf)
        loop_along_z: do kk=1,beam_grid%nz
            loop_along_y: do jj=1,beam_grid%ny
                loop_along_x: do ii=1,beam_grid%nx
                    phit = npa_chords%phit(ii,jj,kk,ichan)
                    if(phit%p.gt.0.d0) then
                        pos = [beam_grid%xc(ii), beam_grid%yc(jj), beam_grid%zc(kk)]
                        call get_fields(fields,pos=pos)
                        if(.not.fields%in_plasma) cycle loop_along_x

                        !!Check if it hits a detector just to make sure
                        dpos = phit%eff_rd
                        vi_norm = phit%dir
                        call hit_npa_detector(pos,vi_norm,det)
                        if (det.ne.ichan) then
                            if(inputs%verbose.ge.0) then
                                write(*,'(a)') 'NPA_WEIGHTS: Missed detector'
                            endif
                            cycle loop_along_x
                        endif

                        !! Determine the angle between the B-field and the Line of Sight
                        pitch = phit%pitch
                        ipitch=minloc(abs(ptcharr - pitch))
                        loop_over_energy: do ic = istart, inputs%ne_wght,istep !! energy loop
                            vabs = sqrt(ebarr(ic)/(v2_to_E_per_amu*beam_mass))
                            vi = vi_norm*vabs
                            !!Correct for gyro orbit
                            call gyro_step(vi,fields,fbm%A,r_gyro)

                            fbm_denf=0
                            if (inputs%dist_type.eq.1) then
                                !get dist at guiding center
                                call get_ep_denf(ebarr(ic),pitch,fbm_denf,pos=(pos+r_gyro))
                            endif
                            if (fbm_denf.ne.fbm_denf) cycle loop_over_energy

                            !! -------------- calculate CX probability -------!!
                            call get_total_cx_rate([ii,jj,kk], pos, vi, neut_types, pcx)
                            if(sum(pcx).le.0) cycle loop_over_energy

                            !!Calculate attenuation
                            states = pcx*1.0d14 !!needs to be large aribitrary number so colrad works
                            states_i=states
                            call attenuate(pos,dpos,vi,states)
                            pcxa=sum(states)/sum(states_i)

                            !$OMP CRITICAL(npa_wght)
                            nweight%attenuation(ic,ii,jj,kk,ichan) = pcxa
                            nweight%cx(ic,ii,jj,kk,ichan) = sum(pcx)
                            nweight%weight(ic,ipitch(1),ichan) = nweight%weight(ic,ipitch(1),ichan) + &
                                      2*sum(pcx)*pcxa*phit%p*beam_grid%dv/dP

                            nweight%flux(ic,ichan) = nweight%flux(ic,ichan) + &
                                      2*beam_grid%dv*fbm_denf*sum(pcx)*pcxa*phit%p
                                      !Factor of 2 above is to convert fbm to ions/(cm^3 dE (domega/4pi))
                            nweight%emissivity(ii,jj,kk,ichan)=nweight%emissivity(ii,jj,kk,ichan)+ &
                                      2*fbm_denf*sum(pcx)*pcxa*phit%p*dE
                            !$OMP END CRITICAL(npa_wght)
                        enddo loop_over_energy
                    endif
                enddo loop_along_x
            enddo loop_along_y
        enddo loop_along_z
        !$OMP END PARALLEL DO
    enddo loop_over_channels

#ifdef _MPI
    call parallel_sum(nweight%weight)
    call parallel_sum(nweight%flux)
    call parallel_sum(nweight%cx)
    call parallel_sum(nweight%attenuation)
    call parallel_sum(nweight%emissivity)
#endif

    do ichan=1,npa_chords%nchan
        if(inputs%verbose.ge.1) then
            write(*,'(T4,"Channel: ",i3)') ichan
            write(*,'(T4,"Radius: ",f10.3)') npa_chords%radius(ichan)
            write(*,'(T4,A,ES14.5)') 'Flux:   ',sum(nweight%flux(:,ichan))*dE
            write(*,'(T4,A,ES14.5)') 'Weight: ',sum(nweight%weight(:,:,ichan))*dE*dP
            write(*,*) ''
        endif
    enddo

    if(inputs%verbose.ge.1) then
        write(*,*) 'write npa weights:    ' , time_string(time_start)
    endif

#ifdef _MPI
    if(my_rank().eq.0) call write_npa_weights()
#else
    call write_npa_weights()
#endif

end subroutine npa_weights

end module libfida

!=============================================================================
!-------------------------------Main Program----------------------------------
!=============================================================================
program fidasim
    !+ FIDASIM {!../VERSION!}
    use libfida
    use hdf5_utils
#ifdef _OMP
    use omp_lib
#endif
#ifdef _MPI
    use mpi_utils
#endif
    implicit none
    character(3)          :: arg = ''
    integer               :: i,narg,nthreads,max_threads,seed

#ifdef _VERSION
    version = _VERSION
#endif

#ifdef _MPI
    call init_mpi()
    if(my_rank().eq.0) call print_banner()
#else
    call print_banner()
#endif

    narg = command_argument_count()
    if(narg.eq.0) then
#ifdef _MPI
        if(my_rank().eq.0) write(*,'(a)') "usage: mpirun -np [num_processes] ./fidasim namelist_file"
        call cleanup_mpi()
#else
        write(*,'(a)') "usage: ./fidasim namelist_file [num_threads]"
#endif
        stop
    else
        call get_command_argument(1,namelist_file)
    endif

    !! Check if compression is possible
    call check_compression_availability()

    !! measure time
    call date_and_time (values=time_start)

    call read_inputs()

#ifdef _OMP
    max_threads = OMP_get_num_procs()
    if(narg.ge.2) then
        call get_command_argument(2,arg)
        read(arg,'(i3)') nthreads
    else
        nthreads = max_threads
    endif
    max_threads = min(nthreads,max_threads)
    if(inputs%verbose.ge.1) then
        write(*,'(a)') "---- OpenMP settings ----"
        write(*,'(T2,"Number of threads: ",i2)') max_threads
        write(*,*) ''
    endif
    call OMP_set_num_threads(max_threads)
#else
    max_threads = 1
#endif

#ifdef _MPI
    istart = my_rank()+1
    istep = num_ranks()
    if(inputs%verbose.ge.1) then
        write(*,'(a)') "---- MPI settings ----"
        write(*,'(T2,"Number of processes: ",i3)') istep
        write(*,*) ''
    endif
#endif

    !! ----------------------------------------------------------
    !! ------ INITIALIZE THE RANDOM NUMBER GENERATOR  -----------
    !! ----------------------------------------------------------
    allocate(rng(max_threads))
#ifdef _OMP
    do i=1,max_threads
        if(inputs%seed.lt.0) then
            call rng_init(rng(i), inputs%seed)
        else
            call rng_init(rng(i), inputs%seed + i)
        endif
    enddo
#else
    call rng_init(rng(1), inputs%seed)
#endif
    if(inputs%verbose.ge.1) then
        write(*,'(a)') "---- Random Number Generator settings ----"
        write(*,'(T2,"RNG Seed: ",i10)') inputs%seed
        write(*,*) ''
    endif

    !! ----------------------------------------------------------
    !! ------- READ GRIDS, PROFILES, LOS, TABLES, & FBM --------
    !! ----------------------------------------------------------
    call read_plasma()
    call read_tables()
    call read_equilibrium()
    call make_beam_grid()
    if(inputs%calc_beam.ge.1) call read_beam()
    call read_distribution()

    call quasineutrality_check()

    allocate(spec_chords%inter(beam_grid%nx,beam_grid%ny,beam_grid%nz))
    if((inputs%calc_spec.ge.1).or.(inputs%calc_fida_wght.ge.1)) then
        call read_chords()
    endif

    if((inputs%calc_npa.ge.1).or.(inputs%calc_npa_wght.ge.1).or.(inputs%calc_pnpa.ge.1)) then
        call read_npa()
    endif

    if(inputs%calc_cfpd.ge.1) then
        call read_cfpd()
    endif

    call make_diagnostic_grids()

    !! ----------------------------------------------------------
    !! --------------- ALLOCATE THE RESULT ARRAYS ---------------
    !! ----------------------------------------------------------
    if(inputs%calc_birth.ge.1) then
        allocate(birth%dens(3, &
                            beam_grid%nx, &
                            beam_grid%ny, &
                            beam_grid%nz))
        allocate(birth%part(int(3*inputs%n_birth*inputs%n_nbi)))
    endif

    !! Spectra
    if(inputs%calc_spec.ge.1) then
        if(inputs%calc_brems.ge.1) then
            allocate(spec%brems(inputs%nlambda,spec_chords%nchan))
            spec%brems = 0.d0
        endif
        if(inputs%calc_bes.ge.1) then
            allocate(spec%full(n_stark,inputs%nlambda,spec_chords%nchan))
            allocate(spec%half(n_stark,inputs%nlambda,spec_chords%nchan))
            allocate(spec%third(n_stark,inputs%nlambda,spec_chords%nchan))
            spec%full = 0.d0
            spec%half = 0.d0
            spec%third = 0.d0
            allocate(spec%fullstokes(n_stark,4,inputs%nlambda,spec_chords%nchan))
            allocate(spec%halfstokes(n_stark,4,inputs%nlambda,spec_chords%nchan))
            allocate(spec%thirdstokes(n_stark,4,inputs%nlambda,spec_chords%nchan))
            spec%fullstokes = 0.d0
            spec%halfstokes = 0.d0
            spec%thirdstokes = 0.d0
        endif
        if(inputs%calc_dcx.ge.1) then
            allocate(spec%dcx(n_stark,inputs%nlambda,spec_chords%nchan,n_thermal))
            spec%dcx = 0.d0
            allocate(spec%dcxstokes(n_stark,4,inputs%nlambda,spec_chords%nchan,n_thermal))
            spec%dcxstokes = 0.d0
        endif
        if(inputs%calc_halo.ge.1) then
            allocate(spec%halo(n_stark,inputs%nlambda,spec_chords%nchan,n_thermal))
            spec%halo = 0.d0
            allocate(spec%halostokes(n_stark,4,inputs%nlambda,spec_chords%nchan,n_thermal))
            spec%halostokes = 0.d0
        endif
        if(inputs%calc_cold.ge.1) then
            allocate(spec%cold(n_stark,inputs%nlambda,spec_chords%nchan,n_thermal))
            spec%cold = 0.d0
            allocate(spec%coldstokes(n_stark,4,inputs%nlambda,spec_chords%nchan,n_thermal))
            spec%coldstokes = 0.d0
        endif
        if(inputs%calc_fida.ge.1) then
            allocate(spec%fida(n_stark,inputs%nlambda,spec_chords%nchan,particles%nclass))
            spec%fida = 0.d0
            allocate(spec%fidastokes(n_stark,4,inputs%nlambda,spec_chords%nchan,particles%nclass))
            spec%fidastokes = 0.d0
        endif
        if(inputs%calc_pfida.ge.1) then
            allocate(spec%pfida(n_stark,inputs%nlambda,spec_chords%nchan,particles%nclass))
            spec%pfida = 0.d0
            allocate(spec%pfidastokes(n_stark,4,inputs%nlambda,spec_chords%nchan,particles%nclass))
            spec%pfidastokes = 0.d0
        endif
    endif

    if(inputs%calc_res.ge.1) then
        if(inputs%calc_bes.ge.1) then
            allocate(spatres%full(spec_chords%nchan))
            allocate(spatres%half(spec_chords%nchan))
            allocate(spatres%third(spec_chords%nchan))
        endif
        if(inputs%calc_dcx.ge.1) allocate(spatres%dcx(spec_chords%nchan))
        if(inputs%calc_halo.ge.1) allocate(spatres%halo(spec_chords%nchan))
        if(inputs%calc_fida.ge.1) allocate(spatres%fida(spec_chords%nchan))
        if(inputs%calc_pfida.ge.1) allocate(spatres%pfida(spec_chords%nchan))
    endif

    if(inputs%calc_npa.ge.1)then
        npa%nchan = npa_chords%nchan
        allocate(npa%part(npa%nmax))
        if(inputs%dist_type.eq.1) then
            npa%nenergy = fbm%nenergy
            allocate(npa%energy(npa%nenergy))
            npa%energy = fbm%energy
        else
            allocate(npa%energy(npa%nenergy))
            do i=1,npa%nenergy
                npa%energy(i)=real(i-0.5)
            enddo
        endif
        allocate(npa%flux(npa%nenergy,npa%nchan,particles%nclass))
        npa%flux = 0.0
    endif

    if(inputs%calc_pnpa.ge.1)then
        pnpa%nchan = npa_chords%nchan
        allocate(pnpa%part(pnpa%nmax))
        if(inputs%dist_type.eq.1) then
            pnpa%nenergy = fbm%nenergy
            allocate(pnpa%energy(pnpa%nenergy))
            pnpa%energy = fbm%energy
        else
            allocate(pnpa%energy(pnpa%nenergy))
            do i=1,pnpa%nenergy
                pnpa%energy(i)=real(i-0.5)
            enddo
        endif
        allocate(pnpa%flux(pnpa%nenergy,pnpa%nchan,particles%nclass))
        pnpa%flux = 0.0
    endif

    if(inputs%calc_neutron.ge.1)then
        allocate(neutron%rate(particles%nclass))
        neutron%rate = 0.d0
    endif

    !! -----------------------------------------------------------------------
    !! --------------- CALCULATE/LOAD the BEAM and HALO DENSITY---------------
    !! -----------------------------------------------------------------------
    if(inputs%load_neutrals.eq.1) then
        call read_neutrals()

        if(inputs%calc_bes.ge.1) then
            if(inputs%verbose.ge.1) then
                write(*,*) 'nbi:     ' , time_string(time_start)
            endif
            call nbi_spec()
            if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
        endif

        if(inputs%calc_dcx.ge.1) then
            if(inputs%verbose.ge.1) then
                write(*,*) 'dcx:     ' , time_string(time_start)
            endif
            call dcx_spec()
            if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
        endif

        if(inputs%calc_halo.ge.1) then
            if(inputs%verbose.ge.1) then
                write(*,*) 'halo:    ' , time_string(time_start)
            endif
            call halo_spec()
            if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
        endif
    else
        if(inputs%calc_beam.ge.1) then
            !! ----------- BEAM NEUTRALS ---------- !!
            if(inputs%calc_nbi_dens.ge.1) then
                if(inputs%verbose.ge.1) then
                    write(*,*) 'nbi:     ' , time_string(time_start)
                endif
                call ndmc()
                if(inputs%verbose.ge.1) write(*,'(30X,a)') ''

                if(inputs%calc_birth.eq.1)then
                    if(inputs%verbose.ge.1) then
                        write(*,*) 'write birth:    ' , time_string(time_start)
                    endif
                    call write_birth_profile()
                    if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
                endif
            endif

            !! ---------- DCX (Direct charge exchange) ---------- !!
            if(inputs%calc_dcx_dens.ge.1) then
                if(inputs%verbose.ge.1) then
                    write(*,*) 'dcx:     ' , time_string(time_start)
                endif
                call dcx()
                if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
            endif

            !! ---------- HALO ---------- !!
            if(inputs%calc_halo_dens.ge.1) then
                if(inputs%verbose.ge.1) then
                    write(*,*) 'halo:    ' , time_string(time_start)
                endif
                call halo()
                if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
            endif

            !! ---------- WRITE NEUTRALS ---------- !!
            if(inputs%verbose.ge.1.and.inputs%write_neutrals.ge.1) then
                write(*,*) 'write neutrals:    ' , time_string(time_start)
            endif
#ifdef _MPI
            if(my_rank().eq.0.and.(inputs%write_neutrals.eq.1)) call write_neutrals()
#else
            if(inputs%write_neutrals.eq.1) call write_neutrals()
#endif
            if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
        endif
    endif

    !! -----------------------------------------------------------------------
    !!------------------------------ COLD D-ALPHA ----------------------------
    !! -----------------------------------------------------------------------
    if(inputs%calc_cold.ge.1) then
        if(inputs%verbose.ge.1) then
            write(*,*) 'cold:    ' ,time_string(time_start)
        endif
        call cold_spec()
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
    endif

    !! -----------------------------------------------------------------------
    !!----------------------------- BREMSSTRAHLUNG ---------------------------
    !! -----------------------------------------------------------------------
    if(inputs%calc_brems.ge.1) then
        if(inputs%verbose.ge.1) then
            write(*,*) 'bremsstrahlung:    ' ,time_string(time_start)
        endif
        call bremsstrahlung()
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
    endif

    !! -----------------------------------------------------------------------
    !! --------------------- CALCULATE the FIDA RADIATION --------------------
    !! -----------------------------------------------------------------------
    if(inputs%calc_fida.ge.1)then
        if(inputs%verbose.ge.1) then
            write(*,*) 'fida:    ' ,time_string(time_start)
        endif
        if(inputs%dist_type.eq.1) then
            call fida_f()
        else
            call fida_mc()
        endif
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
    endif

    if(inputs%calc_pfida.ge.1)then
        if(inputs%verbose.ge.1) then
            write(*,*) 'pfida:   ' ,time_string(time_start)
        endif
        if(inputs%dist_type.eq.1) then
            call pfida_f()
        else
            call pfida_mc()
        endif
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
    endif

    if(inputs%calc_spec.ge.1) then
        if(inputs%verbose.ge.1) then
            write(*,*) 'write spectra:    ' , time_string(time_start)
        endif
#ifdef _MPI
        if(my_rank().eq.0) call write_spectra()
#else
        call write_spectra()
#endif
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
    endif

    !! -----------------------------------------------------------------------
    !! ----------------------- CALCULATE the NPA FLUX ------------------------
    !! -----------------------------------------------------------------------
    if(inputs%calc_npa.ge.1)then
        if(inputs%verbose.ge.1) then
            write(*,*) 'npa:     ' ,time_string(time_start)
        endif
        if(inputs%dist_type.eq.1) then
            call npa_f()
        else
            call npa_mc()
        endif
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
    endif

    if(inputs%calc_pnpa.ge.1)then
        if(inputs%verbose.ge.1) then
            write(*,*) 'pnpa:     ' ,time_string(time_start)
        endif
        if(inputs%dist_type.eq.1) then
            call pnpa_f()
        else
            call pnpa_mc()
        endif
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
    endif

    if((inputs%calc_npa.ge.1).or.(inputs%calc_pnpa.ge.1)) then
        if(inputs%verbose.ge.1) then
            write(*,*) 'write npa:    ' , time_string(time_start)
        endif
        call write_npa()
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
    endif

    !! -------------------------------------------------------------------
    !! ------------------- Calculation of neutron flux -------------------
    !! -------------------------------------------------------------------
    if(inputs%calc_neutron.ge.1) then
        if(inputs%verbose.ge.1) then
            write(*,*) 'neutron rate:    ', time_string(time_start)
        endif
        if(inputs%dist_type.eq.1) then
            call neutron_f()
        else
            call neutron_mc()
        endif
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
    endif

    if(inputs%calc_cfpd.ge.1) then
        if(inputs%verbose.ge.1) then
            write(*,*) 'charged fusion products:    ', time_string(time_start)
        endif
        call cfpd_f()
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
    endif

    !! -------------------------------------------------------------------
    !! ----------- Calculation of weight functions -----------------------
    !! -------------------------------------------------------------------
    if(inputs%calc_fida_wght.ge.1) then
        if(inputs%verbose.ge.1) then
            write(*,*) 'fida weight function:    ', time_string(time_start)
        endif
        if(inputs%calc_fida_wght.eq.1) then
            call fida_weights_los()
        else
            call fida_weights_mc()
        endif
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
    endif

    if(inputs%calc_npa_wght.ge.1) then
        if(inputs%verbose.ge.1) then
            write(*,*) 'npa weight function:    ', time_string(time_start)
        endif
        call npa_weights()
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
    endif

#ifdef _MPI
    call cleanup_mpi()
#endif

    if(inputs%verbose.ge.1) then
        write(*,*) 'END: hour:minute:second ', time_string(time_start)
    endif

end program fidasim
