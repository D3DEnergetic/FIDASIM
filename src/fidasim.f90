!+ This file contains the main routines for FIDASIM {!../VERSION!}
module libfida
!+ Main FIDASIM library
USE H5LT !! High level HDF5 Interface
USE HDF5 !! Base HDF5
USE hdf5_extra !! Additional HDF5 routines
USE eigensystem, ONLY : eigen, matinv
USE parallel_rng

implicit none

character(30) :: version = ''
    !+ FIDASIM version number
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
integer, parameter :: halo_type  = 4
    !+ Identifier for halo neutral interaction
integer, parameter :: fida_type  = 5
    !+ Identifier for fida neutral interaction
integer, parameter :: brems_type = 6
    !+ Identifier for bremsstrahlung interaction. Acts as dummy type
integer, parameter :: ntypes     = 6
    !+ Number of different types of neutrals

integer, parameter :: beam_ion = 1
    !+ Identifier for a beam ion
integer, parameter :: thermal_ion = 2
    !+ Identifier for a thermal ion

!! Physical units
real(Float64), parameter :: e_amu = 5.485799093287202d-4
    !+ Atomic mass of an electron [amu]
real(Float64), parameter :: H_1_amu = 1.00782504d0
    !+ Atomic mass of Hydrogen-1 [amu]
real(Float64), parameter :: H_2_amu = 2.0141017778d0
    !+ Atomic mass of Hydrogen-2 [amu]
real(Float64), parameter :: B5_amu = 10.81d0
    !+ Atomic mass of Boron [amu]
real(Float64), parameter :: C6_amu = 12.011d0
    !+ Atomic mass of Carbon [amu]
real(Float64), parameter :: mass_u    = 1.6605402d-27
    !+ Atomic mass unit [kg]
real(Float64), parameter :: e0        = 1.60217733d-19
    !+ Electron charge [C]
real(Float64), parameter :: pi        = 3.14159265358979323846264d0
    !+ Pi
real(Float64), parameter :: c0        = 2.99792458d+08
    !+ Speed of light [m/s]
real(Float64), parameter :: h_planck  = 4.135667516d-15
    !+ Planck's constant [eV*s]
real(Float64), parameter :: lambda0   = 656.1d0
    !+ D-alpha emission line [nm]
real(Float64), parameter :: v2_to_E_per_amu = mass_u/(2.*e0*1.d3)*1.d-4
    !+ \(cm^2/s^2\) to keV conversion factor

integer, parameter ::n_stark = 15
    !+ Number of Stark lines
real(Float64), parameter, dimension(n_stark) :: stark_wavel = &
     [-2.20200d-07,-1.65200d-07,-1.37700d-07,-1.10200d-07, &
      -8.26400d-08,-5.51000d-08,-2.75600d-08, 0.00000d0,   &
       2.75700d-08, 5.51500d-08, 8.27400d-08, 1.10300d-07, &
       1.38000d-07, 1.65600d-07, 2.20900d-07               ]
    !+ Stark wavelengths [nm*m/V]
real(Float64), parameter, dimension(n_stark) :: stark_intens= &
     [ 1.000d0, 18.00d0, 16.00d0, 1681.d0, 2304.d0, &
       729.0d0, 1936.d0, 5490.d0, 1936.d0, 729.0d0, &
       2304.d0, 1681.d0, 16.00d0, 18.00d0, 1.000d0  ]
    !+ Stark Intensities
integer, parameter, dimension(n_stark) :: stark_pi= &
     [1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1]
    !+ Pi line indicators
integer, parameter, dimension(n_stark) :: stark_sigma=1 - stark_pi
    !+ Sigma line indicators

!!Numerical Settings
integer, parameter :: nlevs=6
    !+ Number of atomic energy levels
real(Float64), parameter :: n_halo_neutrate=20.
    !+ Number of times to average halo H_H_cx
real(Float64) :: colrad_threshold=1.d6
    !+ colrad threshold
real(Float64), dimension(ntypes) :: halo_iter_dens = 0.d0
    !+ Keeps track of how of each generations halo density
integer :: nbi_outside = 0
    !+ Keeps track of how many beam neutrals do not hit the [[libfida:beam_grid]]

type InterpolCoeffs1D
    !+ Linear Interpolation Coefficients and indices
    integer :: i = 0
        !+ Index of position right before `xout`
    real(Float64) :: b1 = 0.d0
        !+ Coefficient for y(i) term
    real(Float64) :: b2 = 0.d0
        !+ Coefficient for y(i+1) term
end type InterpolCoeffs1D

type InterpolCoeffs2D
    !+ 2D Linear Interpolation Coefficients and indices
    integer :: i = 0
        !+ Index of abscissa before `xout`
    integer :: j = 0
        !+ Index of ordinate before `yout`
    real(Float64) :: b11 = 0.d0
        !+ Coefficient for z(i,j) term
    real(Float64) :: b12 = 0.d0
        !+ Coefficient for z(i,j+1) term
    real(Float64) :: b21 = 0.d0
        !+ Coefficient for z(i+1,j) term
    real(Float64) :: b22 = 0.d0
        !+ Coefficient for z(i+1,j+1) term
end type InterpolCoeffs2D

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
    !+ Defines a 2D R-Z grid for interpolating plasma parameters and fields
    integer(Int32) :: nr
        !+ Number of Radii
    integer(Int32) :: nz
        !+ Number of Z values
    real(Float64)  :: dr
        !+ Radial spacing [cm]
    real(Float64)  :: dz
        !+ Vertical spacing [cm]
    real(Float64)  :: da
        !+ Grid element area [\(cm^2\)]
    real(Float64), dimension(:),   allocatable :: r
        !+ Radii values [cm]
    real(Float64), dimension(:),   allocatable :: z
        !+ Z values [cm]
    real(Float64), dimension(:,:), allocatable :: r2d
        !+ 2D R grid [cm]
    real(Float64), dimension(:,:), allocatable :: z2d
        !+ 2D Z grid [cm]
end type InterpolationGrid

type Profiles
    !+ Torodial symmetric plasma parameters at a given R-Z
    real(Float64) :: dene = 0.d0
        !+ Electron density [\(cm^{-3}\)]
    real(Float64) :: denp = 0.d0
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
end type Profiles

type, extends( Profiles ) :: LocalProfiles
    !+ Plasma parameters at given position
    logical :: in_plasma = .False.
        !+ Indicates whether plasma parameters are valid/known
    logical :: machine_coords = .False.
        !+ Indicates whether vectors are in machine coordinates
    real(Float64), dimension(3) :: pos = 0.d0
        !+ Position in beam grid coordinates
    real(Float64), dimension(3) :: uvw = 0.d0
        !+ Position in machine coordinates
    real(Float64), dimension(3) :: vrot = 0.d0
        !+ Plasma rotation in beam grid coordinates
    type(InterpolCoeffs2D) :: c
        !+ Linear Interpolation Coefficients and indicies for interpolation at `pos`
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
end type EMFields

type, extends( EMFields ) :: LocalEMFields
    !+ Electro-magnetic fields at given position
    logical       :: in_plasma = .False.
        !+ Indicates whether fields are valid/known
    logical :: machine_coords = .False.
        !+ Indicates whether vectors are in machine coordinates
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
    type(InterpolCoeffs2D) :: c
        !+ Linear Interpolation Coefficients and indicies for interpolation at `pos`
end type LocalEMFields

type Equilibrium
    !+MHD Equilbrium
    type(EMFields), dimension(:,:), allocatable :: fields
        !+ Electro-magnetic fields at points defined in [[libfida:inter_grid]]
    type(Profiles), dimension(:,:), allocatable :: plasma
        !+ Plasma parameters at points defined in [[libfida:inter_grid]]
    real(Float64), dimension(:,:), allocatable  :: mask
        !+ Indicates whether fields and plasma are well-defined at points defined in [[libfida:inter_grid]]
end type Equilibrium

type FastIonDistribution
    !+ Defines a Guiding Center Fast-ion Distribution Function: F(E,p,R,Z)
    integer(Int32) :: nenergy
        !+ Number of energies
    integer(Int32) :: npitch
        !+ Number of pitches
    integer(Int32) :: nr
        !+ Number of radii
    integer(Int32) :: nz
        !+ Number of z values
    real(Float64)  :: dE
        !+ Energy spacing [keV]
    real(Float64)  :: dp
        !+ Pitch spacing
    real(Float64)  :: dr
        !+ Radial spacing [cm]
    real(Float64)  :: dz
        !+ Z spacing [cm]
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
    real(Float64), dimension(:), allocatable       :: energy
        !+ Energy values [keV]
    real(Float64), dimension(:), allocatable       :: pitch
        !+ Pitch w.r.t. the magnetic field
    real(Float64), dimension(:), allocatable       :: r
        !+ Radius [cm]
    real(Float64), dimension(:), allocatable       :: z
        !+ Z [cm]
    real(Float64), dimension(:,:), allocatable     :: denf
        !+ Fast-ion density defined on the [[libfida:inter_grid]]: denf(R,Z)
    real(Float64), dimension(:,:,:,:), allocatable :: f
        !+ Fast-ion distribution function defined on the [[libfida:inter_grid]]: F(E,p,R,Z)
end type FastIonDistribution

type FastIon
    !+ Defines a fast-ion
    logical        :: cross_grid = .False.
        !+ Indicates whether the fast-ion crosses the [[libfida:beam_grid]]
    real(Float64)  :: r = 0.d0
        !+ Radial position of fast-ion [cm]
    real(Float64)  :: z = 0.d0
        !+ Vertical position of fast-ion [cm]
    real(Float64)  :: phi_enter = 0.d0
        !+ Torodial/phi position where fast-ion enters the [[libfida:beam_grid]] [radians]
    real(Float64)  :: delta_phi = 0.d0
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
    !+ Defines an atomic table for populating and de-populating reaction rates
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
    real(Float64), dimension(2) :: ab = 0.d0
        !+ Atomic mass of beam and thermal ions respectively [amu]
    real(Float64), dimension(:,:,:,:,:), allocatable :: log_pop
        !+ Log-10 reaction rates for populating transistions
    real(Float64), dimension(:,:,:,:), allocatable   :: log_depop
        !+ Log-10 reaction rates for de-populating transistions
end type AtomicRates

type AtomicTables
    !+ Atomic tables for various types of interactions
    type(AtomicCrossSection) :: H_H_cx
        !+ Hydrogen-Hydrogen charge exchange n/m-resolved cross sections
    type(AtomicRates)        :: H_H
        !+ Hydrogen-Hydrogen reaction rates
    type(AtomicRates)        :: H_e
        !+ Hydrogen-Electron reaction rates
    type(AtomicRates)        :: H_Aq
        !+ Hydrogen-Impurity reaction rates
    real(Float64), dimension(nlevs,nlevs) :: einstein
        !+ Einstein coefficients for spontaneous emission
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
end type LineOfSight

type SpectralChords
    !+ Defines an spectral diagnostic system
    integer :: nchan = 0
        !+ Number of channels
    type(LineOfSight), dimension(:), allocatable   :: los
        !+ Line of sight array
    real(Float64), dimension(:), allocatable       :: radius
        !+ Radius of each line of sight
    logical, dimension(:,:,:), allocatable         :: los_inter
        !+ Indicates whether a [[libfida:beam_grid]] cell intersects a LOS
    real(Float32), dimension(:,:,:,:), allocatable :: dlength
        !+ [[libfida:beam_grid]] cell - LOS intersection length: dlength(x,y,z,chan)
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
    real(Float64), dimension(3) :: eff_rd = 0.d0
        !+ Effective position of detector
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
    integer(Int32) :: nloop = 1000
        !+ Increases number of MC markers by factor of `nloop`
    type(NPAParticle), dimension(:), allocatable  :: part
        !+ Array of NPA particles
    real(Float64), dimension(:), allocatable      :: energy
        !+ Energy array [keV]
    real(Float64), dimension(:,:), allocatable    :: flux
        !+ Neutral particle flux [neutrals/(s*dE)]
end type NPAResults

type BirthProfile
    !+ Birth profile structure
    integer :: cnt = 1
        !+ Particle counter
    integer, dimension(:), allocatable             :: neut_type
        !+ Particle birth type (1=Full, 2=Half, 3=Third)
    real(Float64), dimension(:,:), allocatable     :: ri
        !+ Particle birth position [cm]
    real(Float64), dimension(:,:), allocatable     :: vi
        !+ Particle birth velocity [cm/s]
    integer, dimension(:,:), allocatable           :: ind
        !+ Particle [[libfida:beam_grid]] indices
    real(Float64), dimension(:,:,:,:), allocatable :: dens
        !+ Birth density: dens(neutral_type,x,y,z) [fast-ions/(s*cm^3)]
end type BirthProfile

type Spectra
    !+ Spectra storage structure
    real(Float64), dimension(:,:), allocatable   :: brems
        !+ Bremsstruhlung: brems(lambda,chan)
    real(Float64), dimension(:,:,:), allocatable :: bes
        !+ Beam emission: bes(lambda,chan,neutral_type)
    real(Float64), dimension(:,:,:), allocatable :: fida
        !+ FIDA emission: fida(lambda,chan,orbit_type)
end type Spectra

type NeutralDensity
    !+ Neutral density structure
    real(Float64), dimension(:,:,:,:,:), allocatable :: dens
        !+ Neutral density: dens(lev,neutral_type,x,y,z)
end type NeutralDensity

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

    !! Monte Carlo settings
    integer(Int32) :: n_fida
        !+ Number of FIDA mc markers
    integer(Int32) :: n_npa
        !+ Number of NPA mc markers
    integer(Int32) :: n_nbi
        !+ Number of neutral beam mc markers
    integer(Int32) :: n_dcx
        !+ Number of direct charge exchange (DCX) mc markers
    integer(Int32) :: n_halo
        !+ Number of halo mc markers
    integer(Int32) :: n_birth
        !+ Number of birth particles per [[SimulationInputs:n_nbi]]

    !! Simulation switches
    integer(Int32) :: calc_spec
        !+ Calculate spectra: 0 = off, 1=on
    integer(Int32) :: calc_brems
        !+ Calculate bremmstruhlung: 0 = off, 1=on
    integer(Int32) :: calc_bes
        !+ Calculate BES: 0 = off, 1=on
    integer(Int32) :: calc_fida
        !+ Calculate FIDA: 0 = off, 1=on
    integer(Int32) :: load_neutrals
        !+ Load neutrals from file: 0 = off, 1=on
    integer(Int32) :: calc_npa
        !+ Calculate NPA: 0 = off, 1=on, 2=on++
    integer(Int32) :: calc_fida_wght
        !+ Calculate FIDA weight: 0 = off, 1=on, 2=on++
    integer(Int32) :: calc_npa_wght
        !+ Calculate NPA weights: 0 = off, 1=on, 2=on++
    integer(Int32) :: calc_birth
        !+ Calculate birth profile: 0 = off, 1=on
    integer(Int32) :: no_flr
        !+ Turns off Finite Larmor Radius effects: 0=off, 1=on
    integer(Int32) :: dump_dcx
        !+ Output DCX density and spectra: 0 = off, 1=on
    integer(Int32) :: verbose
        !+ Verbosity: 0 = off, 1=on, 2=on++

    !! Neutral Beam Settings
    real(Float64)    :: ab
        !+ Atomic mass of beam neutrals

    !! Plasma parameters
    integer(Int32) :: impurity_charge
        !+ Impurity proton number
    real(Float64)  :: ai
        !+ Atomic mass of thermal ions

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

interface assignment(=)
    !+ Allows for assigning [[Profiles]],[[LocalProfiles]],
    !+ [[EMFields]],[[LocalEMFields]],[[FastIon]], and [[NPAParticle]]
    module procedure pp_assign, lpp_assign, plp_assign, lplp_assign, &
                     ff_assign, lff_assign, flf_assign, lflf_assign, &
                     fast_ion_assign,npa_part_assign
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

interface interpol_coeff
    !+ Calculates linear interpolation coefficients
    module procedure interpol1D_coeff, interpol1D_coeff_arr
    module procedure interpol2D_coeff, interpol2D_coeff_arr
end interface

interface interpol
    !+ Performs linear/bilinear interpolation
    module procedure interpol1D_arr
    module procedure interpol2D_arr, interpol2D_2D_arr
end interface

!! definition of the structures:
type(BeamGrid), save            :: beam_grid
    !+ Variable containing beam grid definition
type(InterpolationGrid), save   :: inter_grid
    !+ Variable containing interpolation grid definition
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
    !+ Variable for storing the calculated NPA results
type(SpectralChords), save      :: spec_chords
    !+ Variable containing the spectral system definition
type(NPAChords), save           :: npa_chords
    !+ Variable containing the NPA system definition
type(SimulationInputs), save    :: inputs
    !+ Variable containing the simulation inputs
type(BirthProfile), save        :: birth
    !+ Variable for storing the calculated birth profile
type(NeutralDensity), save      :: neut
    !+ Variable for storing the calculated beam density
type(Spectra), save             :: spec
    !+ Variable for storing the calculated spectra
type(FIDAWeights), save         :: fweight
    !+ Variable for storing the calculated FIDA weights
type(NPAWeights), save          :: nweight
    !+ Variable for storing the calculated NPA weights

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

#ifdef _OMP
#else
    write(*,'(a)') "########################### ATTENTION ###########################"
    write(*,'(a)') "#              OpenMP threading has been disabled               #"
    write(*,'(a)') "#################################################################"
    write(*,'(a)') ""
#endif

end subroutine print_banner

!============================================================================
!---------------------------Operator Overloading-----------------------------
!============================================================================
subroutine fast_ion_assign(p1, p2)
    !+ Defines how to assign [[FastIon]] types to eachother
    type(FastIon), intent(in)  :: p2
    type(FastIon), intent(out) :: p1

    p1%cross_grid = p2%cross_grid
    p1%r          = p2%r
    p1%z          = p2%z
    p1%phi_enter  = p2%phi_enter
    p1%delta_phi  = p2%delta_phi
    p1%energy     = p2%energy
    p1%pitch      = p2%pitch
    p1%vabs       = p2%vabs
    p1%vr         = p2%vr
    p1%vt         = p2%vt
    p1%vz         = p2%vz
    p1%weight     = p2%weight
    p1%class      = p2%class

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

subroutine pp_assign(p1, p2)
    !+ Defines how to assign [[Profiles]] types to eachother
    type(Profiles), intent(in)    :: p2
    type(Profiles), intent(inout) :: p1

    p1%dene   = p2%dene
    p1%ti     = p2%ti
    p1%te     = p2%te
    p1%denp   = p2%denp
    p1%denf   = p2%denf
    p1%denimp = p2%denimp
    p1%zeff   = p2%zeff
    p1%vr     = p2%vr
    p1%vt     = p2%vt
    p1%vz     = p2%vz

end subroutine pp_assign

subroutine lpp_assign(p1, p2)
    !+ Defines how to assign a [[Profiles]] type to a [[LocalProfiles]] type
    type(Profiles), intent(in)         :: p2
    type(LocalProfiles), intent(inout) :: p1

    p1%dene   = p2%dene
    p1%ti     = p2%ti
    p1%te     = p2%te
    p1%denp   = p2%denp
    p1%denf   = p2%denf
    p1%denimp = p2%denimp
    p1%zeff   = p2%zeff
    p1%vr     = p2%vr
    p1%vt     = p2%vt
    p1%vz     = p2%vz

end subroutine lpp_assign

subroutine plp_assign(p1, p2)
    !+ Defines how to assign a [[LocalProfiles]] type to a [[Profiles]] type
    type(LocalProfiles), intent(in) :: p2
    type(Profiles), intent(inout)   :: p1

    p1%dene   = p2%dene
    p1%ti     = p2%ti
    p1%te     = p2%te
    p1%denp   = p2%denp
    p1%denf   = p2%denf
    p1%denimp = p2%denimp
    p1%zeff   = p2%zeff
    p1%vr     = p2%vr
    p1%vt     = p2%vt
    p1%vz     = p2%vz

end subroutine plp_assign

subroutine lplp_assign(p1, p2)
    !+ Defines how to assign [[LocalProfiles]] types to eachother
    type(LocalProfiles), intent(in)    :: p2
    type(LocalProfiles), intent(inout) :: p1

    p1%pos    = p2%pos
    p1%uvw    = p2%uvw
    p1%dene   = p2%dene
    p1%ti     = p2%ti
    p1%te     = p2%te
    p1%denp   = p2%denp
    p1%denf   = p2%denf
    p1%denimp = p2%denimp
    p1%zeff   = p2%zeff
    p1%vr     = p2%vr
    p1%vt     = p2%vt
    p1%vz     = p2%vz
    p1%vrot   = p2%vrot

end subroutine lplp_assign

subroutine ff_assign(p1, p2)
    !+ Defines how to assign [[EMFields]] types to eachother
    type(EMFields), intent(in)    :: p2
    type(EMFields), intent(inout) :: p1

    p1%br   = p2%br
    p1%bt   = p2%bt
    p1%bz   = p2%bz
    p1%er   = p2%er
    p1%et   = p2%et
    p1%ez   = p2%ez

end subroutine ff_assign

subroutine lff_assign(p1, p2)
    !+ Defines how to assign a [[EMFields]] type to a [[LocalEMFields]] type
    type(EMFields), intent(in)         :: p2
    type(LocalEMFields), intent(inout) :: p1

    p1%br   = p2%br
    p1%bt   = p2%bt
    p1%bz   = p2%bz
    p1%er   = p2%er
    p1%et   = p2%et
    p1%ez   = p2%ez

end subroutine lff_assign

subroutine flf_assign(p1, p2)
    !+ Defines how to assign a [[LocalEMFields]] type to a [[EMFields]] type
    type(LocalEMFields), intent(in) :: p2
    type(EMFields), intent(inout)   :: p1

    p1%br   = p2%br
    p1%bt   = p2%bt
    p1%bz   = p2%bz
    p1%er   = p2%er
    p1%et   = p2%et
    p1%ez   = p2%ez

end subroutine flf_assign

subroutine lflf_assign(p1, p2)
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

end subroutine lflf_assign

function pp_add(p1, p2) result (p3)
    !+ Defines how to add two [[Profiles]] types
    type(Profiles), intent(in) :: p1,p2
    type(Profiles)             :: p3

    p3%dene   = p1%dene   + p2%dene
    p3%ti     = p1%ti     + p2%ti
    p3%te     = p1%te     + p2%te
    p3%denp   = p1%denp   + p2%denp
    p3%denf   = p1%denf   + p2%denf
    p3%denimp = p1%denimp + p2%denimp
    p3%zeff   = p1%zeff   + p2%zeff
    p3%vr     = p1%vr     + p2%vr
    p3%vt     = p1%vt     + p2%vt
    p3%vz     = p1%vz     + p2%vz

end function pp_add

function pp_subtract(p1, p2) result (p3)
    !+ Defines how to subtract two [[Profiles]] types
    type(Profiles), intent(in) :: p1,p2
    type(Profiles)             :: p3

    p3%dene   = p1%dene   - p2%dene
    p3%ti     = p1%ti     - p2%ti
    p3%te     = p1%te     - p2%te
    p3%denp   = p1%denp   - p2%denp
    p3%denf   = p1%denf   - p2%denf
    p3%denimp = p1%denimp - p2%denimp
    p3%zeff   = p1%zeff   - p2%zeff
    p3%vr     = p1%vr     - p2%vr
    p3%vt     = p1%vt     - p2%vt
    p3%vz     = p1%vz     - p2%vz

end function pp_subtract

function lplp_add(p1, p2) result (p3)
    !+ Defines how to add two [[LocalProfiles]] types
    type(LocalProfiles), intent(in) :: p1,p2
    type(LocalProfiles)             :: p3

    p3%pos    = p1%pos    + p2%pos
    p3%uvw    = p1%uvw    + p2%uvw
    p3%dene   = p1%dene   + p2%dene
    p3%ti     = p1%ti     + p2%ti
    p3%te     = p1%te     + p2%te
    p3%denp   = p1%denp   + p2%denp
    p3%denf   = p1%denf   + p2%denf
    p3%denimp = p1%denimp + p2%denimp
    p3%zeff   = p1%zeff   + p2%zeff
    p3%vr     = p1%vr     + p2%vr
    p3%vt     = p1%vt     + p2%vt
    p3%vz     = p1%vz     + p2%vz
    p3%vrot   = p1%vrot   + p2%vrot

end function lplp_add

function lplp_subtract(p1, p2) result (p3)
    !+ Defines how to subtract two [[LocalProfiles]] types
    type(LocalProfiles), intent(in) :: p1,p2
    type(LocalProfiles)             :: p3

    p3%pos    = p1%pos    - p2%pos
    p3%uvw    = p1%uvw    - p2%uvw
    p3%dene   = p1%dene   - p2%dene
    p3%ti     = p1%ti     - p2%ti
    p3%te     = p1%te     - p2%te
    p3%denp   = p1%denp   - p2%denp
    p3%denf   = p1%denf   - p2%denf
    p3%denimp = p1%denimp - p2%denimp
    p3%zeff   = p1%zeff   - p2%zeff
    p3%vr     = p1%vr     - p2%vr
    p3%vt     = p1%vt     - p2%vt
    p3%vz     = p1%vz     - p2%vz
    p3%vrot   = p1%vrot   - p2%vrot

end function lplp_subtract

function ps_multiply(p1, real_scalar) result (p3)
    !+ Defines how to multiply [[Profiles]] types by a scalar
    type(Profiles), intent(in) :: p1
    real(Float64), intent(in)  :: real_scalar
    type(Profiles)             :: p3

    p3%dene   = p1%dene   * real_scalar
    p3%ti     = p1%ti     * real_scalar
    p3%te     = p1%te     * real_scalar
    p3%denp   = p1%denp   * real_scalar
    p3%denf   = p1%denf   * real_scalar
    p3%denimp = p1%denimp * real_scalar
    p3%zeff   = p1%zeff   * real_scalar
    p3%vr     = p1%vr     * real_scalar
    p3%vt     = p1%vt     * real_scalar
    p3%vz     = p1%vz     * real_scalar

end function ps_multiply

function sp_multiply(real_scalar, p1) result (p3)
    !+ Defines how to multiply [[Profiles]] types by a scalar
    type(Profiles), intent(in) :: p1
    real(Float64), intent(in)  :: real_scalar
    type(Profiles)             :: p3

    p3 = p1*real_scalar

end function sp_multiply

function ps_divide(p1, real_scalar) result (p3)
    !+ Defines how to divide [[Profiles]] types by a scalar
    type(Profiles), intent(in) :: p1
    real(Float64), intent(in)  :: real_scalar
    type(Profiles)             :: p3

    p3 = p1*(1.d0/real_scalar)

end function ps_divide

function lps_multiply(p1, real_scalar) result (p3)
    !+ Defines how to multiply [[LocalProfiles]] types by a scalar
    type(LocalProfiles), intent(in) :: p1
    real(Float64), intent(in)       :: real_scalar
    type(LocalProfiles)             :: p3

    p3%pos    = p1%pos    * real_scalar
    p3%uvw    = p1%uvw    * real_scalar
    p3%dene   = p1%dene   * real_scalar
    p3%ti     = p1%ti     * real_scalar
    p3%te     = p1%te     * real_scalar
    p3%denp   = p1%denp   * real_scalar
    p3%denf   = p1%denf   * real_scalar
    p3%denimp = p1%denimp * real_scalar
    p3%zeff   = p1%zeff   * real_scalar
    p3%vr     = p1%vr     * real_scalar
    p3%vt     = p1%vt     * real_scalar
    p3%vz     = p1%vz     * real_scalar
    p3%vrot   = p1%vrot   * real_scalar

end function lps_multiply

function slp_multiply(real_scalar, p1) result (p3)
    !+ Defines how to multiply [[LocalProfiles]] types by a scalar
    type(LocalProfiles), intent(in) :: p1
    real(Float64), intent(in)       :: real_scalar
    type(LocalProfiles)             :: p3

    p3 = p1*real_scalar

end function slp_multiply

function lps_divide(p1, real_scalar) result (p3)
    !+ Defines how to divide [[LocalProfiles]] types by a scalar
    type(LocalProfiles), intent(in) :: p1
    real(Float64), intent(in)       :: real_scalar
    type(LocalProfiles)             :: p3

    p3 = p1*(1.d0/real_scalar)

end function lps_divide

function ff_add(p1, p2) result (p3)
    !+ Defines how to add two [[EMFields]] types
    type(EMFields), intent(in) :: p1,p2
    type(EMFields)             :: p3

    p3%br   = p1%br   + p2%br
    p3%bt   = p1%bt   + p2%bt
    p3%bz   = p1%bz   + p2%bz
    p3%er   = p1%er   + p2%er
    p3%et   = p1%et   + p2%et
    p3%ez   = p1%ez   + p2%ez

end function ff_add

function ff_subtract(p1, p2) result (p3)
    !+ Defines how to subtract two [[EMFields]] types
    type(EMFields), intent(in) :: p1,p2
    type(EMFields)             :: p3

    p3%br   = p1%br   - p2%br
    p3%bt   = p1%bt   - p2%bt
    p3%bz   = p1%bz   - p2%bz
    p3%er   = p1%er   - p2%er
    p3%et   = p1%et   - p2%et
    p3%ez   = p1%ez   - p2%ez

end function ff_subtract

function fs_multiply(p1, real_scalar) result (p3)
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

end function fs_multiply

function sf_multiply(real_scalar, p1) result (p3)
    !+ Defines how to multiply [[EMFields]] types by a scalar
    type(EMFields), intent(in) :: p1
    real(Float64), intent(in)  :: real_scalar
    type(EMFields)             :: p3

    p3 = p1*real_scalar

end function sf_multiply

function fs_divide(p1, real_scalar) result (p3)
    !+ Defines how to divide [[EMFields]] types by a scalar
    type(EMFields), intent(in) :: p1
    real(Float64), intent(in)  :: real_scalar
    type(EMFields)             :: p3

    p3 = p1*(1.d0/real_scalar)

end function fs_divide

function lflf_add(p1, p2) result (p3)
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

end function lflf_add

function lflf_subtract(p1, p2) result (p3)
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

end function lflf_subtract

function lfs_multiply(p1, real_scalar) result (p3)
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

end function lfs_multiply

function slf_multiply(real_scalar, p1) result (p3)
    !+ Defines how to multiply [[LocalEMFields]] types by a scalar
    type(LocalEMFields), intent(in) :: p1
    real(Float64), intent(in)       :: real_scalar
    type(LocalEMFields)             :: p3

    p3 = p1*real_scalar

end function slf_multiply

function lfs_divide(p1, real_scalar) result (p3)
    !+ Defines how to divide [[LocalEMFields]] types by a scalar
    type(LocalEMFields), intent(in) :: p1
    real(Float64), intent(in)       :: real_scalar
    type(LocalEMFields)             :: p3

    p3 = p1*(1.d0/real_scalar)

end function lfs_divide

!============================================================================
!-------------------------------I/O Routines---------------------------------
!============================================================================
subroutine read_inputs
    !+ Reads input namelist file and stores the results into [[libfida:inputs]],
    !+ [[libfida:nbi]], and [[libfida:beam_grid]]
    character(charlim) :: runid,result_dir, tables_file
    character(charlim) :: distribution_file, equilibrium_file
    character(charlim) :: geometry_file, neutrals_file
    integer            :: pathlen
    integer            :: calc_brems,calc_bes,calc_fida,calc_npa
    integer            :: calc_birth,calc_fida_wght,calc_npa_wght
    integer            :: load_neutrals,verbose,dump_dcx,no_flr
    integer(Int32)     :: shot,n_fida,n_npa,n_nbi,n_halo,n_dcx,n_birth
    integer(Int32)     :: nlambda,ne_wght,np_wght,nphi_wght,nlambda_wght
    real(Float64)      :: time,lambdamin,lambdamax,emax_wght
    real(Float64)      :: lambdamin_wght,lambdamax_wght
    real(Float64)      :: ai,ab,pinj,einj,current_fractions(3)
    integer(Int32)     :: impurity_charge
    integer(Int32)     :: nx,ny,nz
    real(Float64)      :: xmin,xmax,ymin,ymax,zmin,zmax
    real(Float64)      :: alpha,beta,gamma,origin(3)
    logical            :: exis, error

    NAMELIST /fidasim_inputs/ result_dir, tables_file, distribution_file, &
        geometry_file, equilibrium_file, neutrals_file, shot, time, runid, &
        calc_brems, calc_bes, calc_fida, calc_npa, calc_birth, no_flr, &
        calc_fida_wght, calc_npa_wght, load_neutrals, dump_dcx, verbose, &
        n_fida,n_npa, n_nbi, n_halo, n_dcx, n_birth, &
        ab, pinj, einj, current_fractions, ai, impurity_charge, &
        nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, &
        origin, alpha, beta, gamma, &
        ne_wght, np_wght, nphi_wght, &
        nlambda, lambdamin,lambdamax,emax_wght, &
        nlambda_wght,lambdamin_wght,lambdamax_wght

    inquire(file=namelist_file,exist=exis)
    if(.not.exis) then
        write(*,'(a,a)') 'READ_INPUTS: Input file does not exist: ', trim(namelist_file)
        stop
    endif

    !!Set Defaults
    no_flr = 0

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

    !!Simulation Switches
    if((calc_brems+calc_bes+calc_fida).gt.0) then
        inputs%calc_spec=1
    else
        inputs%calc_spec=0
    endif
    inputs%calc_brems=calc_brems
    inputs%calc_bes=calc_bes
    inputs%calc_fida=calc_fida
    inputs%calc_npa=calc_npa
    inputs%calc_birth=calc_birth
    inputs%calc_fida_wght=calc_fida_wght
    inputs%calc_npa_wght=calc_npa_wght
    inputs%load_neutrals=load_neutrals
    inputs%dump_dcx=dump_dcx
    inputs%verbose=verbose
    inputs%no_flr = no_flr

    !!Monte Carlo Settings
    inputs%n_fida=max(10,n_fida)
    inputs%n_npa=max(10,n_npa)
    inputs%n_nbi=max(10,n_nbi)
    inputs%n_halo=max(10,n_halo)
    inputs%n_dcx=max(10,n_dcx)
    inputs%n_birth= max(1,nint(n_birth/real(n_nbi)))

    !!Plasma Settings
    inputs%ai=ai
    inputs%impurity_charge=impurity_charge

    !!Neutral Beam Settings
    inputs%ab=ab
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
        write(*,'(a,a)') 'READ_INPUTS: Tables file does not exist: ', &
                         trim(inputs%tables_file)
        error = .True.
    endif

    inquire(file=inputs%geometry_file,exist=exis)
    if(exis) then
        if(inputs%verbose.ge.1) then
            write(*,'(T2,"Geometry file: ",a)') trim(inputs%geometry_file)
        endif
    else
        write(*,'(a,a)') 'READ_INPUTS: Geometry file does not exist: ', &
                         trim(inputs%geometry_file)
        error = .True.
    endif

    inquire(file=inputs%equilibrium_file,exist=exis)
    if(exis) then
        if(inputs%verbose.ge.1) then
            write(*,'(T2,"Equilibrium file: ",a)') trim(inputs%equilibrium_file)
        endif
    else
        write(*,'(a,a)') 'READ_INPUTS: Equilibrium file does not exist: ', &
                         trim(inputs%equilibrium_file)
        error = .True.
    endif

    inquire(file=inputs%distribution_file,exist=exis)
    if(exis) then
        if(inputs%verbose.ge.1) then
            write(*,'(T2,"Distribution file: ",a)') trim(inputs%distribution_file)
        endif
    else
        write(*,'(a,a)') 'READ_INPUTS: Distribution file does not exist: ', &
                         trim(inputs%distribution_file)
        error = .True.
    endif

    pathlen = len_trim(inputs%result_dir)+len_trim(inputs%runid) + 20
    !+20 for suffixes and seperators e.g. /, _npa.h5, ...
    if(pathlen.gt.charlim) then
        write(*,'(a,i3,a,i3)') 'READ_INPUTS: Result directory path + runID use too many characters: ', &
                               pathlen-20,'>', charlim-20
        error = .True.
    endif

    if(inputs%verbose.ge.1) then
        write(*,*) ''
    endif

    if(error) then
        stop
    endif

end subroutine read_inputs

subroutine make_beam_grid
    !+ Makes [[libfida:beam_grid] from user defined inputs
    integer(Int32) :: i
    real(Float64) :: dx, dy, dz

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
    beam_grid%ntrack = beam_grid%nx+beam_grid%ny+beam_grid%nz
    beam_grid%ngrid  = beam_grid%nx*beam_grid%ny*beam_grid%nz

    call tb_zyx(beam_grid%alpha,beam_grid%beta,beam_grid%gamma, &
                beam_grid%basis, beam_grid%inv_basis)

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
        write(*,*) ''
    endif

end subroutine make_beam_grid

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

    nbi%vinj=sqrt(2.d0*nbi%einj*1.d3 *e0/(inputs%ab*mass_u))*1.d2 !! [cm/s]
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
    real(Float64), dimension(:), allocatable :: spot_size, sigma_pi
    real(Float64) :: r0(3), v0(3), r_enter(3), r_exit(3)
    real(Float64) :: xyz_lens(3), xyz_axis(3), length
    real(Float64), dimension(3,3) :: basis
    real(Float64), dimension(2) :: randomu
    real(Float64) :: theta, sqrt_rho
    type(ParticleTrack), dimension(beam_grid%ntrack) :: tracks
    character(len=20) :: system = ''

    integer :: i, j, ic, nc, ncell, ind(3)
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
        write(*,'(a)') 'FIDA/BES geometry is not in the geometry file'
        write(*,'(a)') 'Continuing without spectral diagnostics'
        inputs%calc_spec = 0
        inputs%calc_fida = 0
        inputs%calc_bes = 0
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
    allocate(spot_size(spec_chords%nchan))
    allocate(sigma_pi(spec_chords%nchan))
    allocate(spec_chords%los(spec_chords%nchan))
    allocate(spec_chords%radius(spec_chords%nchan))
    allocate(spec_chords%dlength(spec_chords%nchan, &
                                 beam_grid%nx, &
                                 beam_grid%ny, &
                                 beam_grid%nz) )

    spec_chords%dlength = 0.d0

    dims = [3,spec_chords%nchan]
    call h5ltread_dataset_double_f(gid, "/spec/lens", lenses, dims, error)
    call h5ltread_dataset_double_f(gid, "/spec/axis", axes, dims, error)
    call h5ltread_dataset_double_f(gid, "/spec/spot_size", spot_size, dims(2:2), error)
    call h5ltread_dataset_double_f(gid, "/spec/sigma_pi", sigma_pi, dims(2:2), error)
    call h5ltread_dataset_double_f(gid, "/spec/radius", spec_chords%radius, dims(2:2), error)

    !!Close SPEC group
    call h5gclose_f(gid, error)

    !!Close file id
    call h5fclose_f(fid, error)

    !!Close HDF5 interface
    call h5close_f(error)

    chan_loop: do i=1,spec_chords%nchan
        call uvw_to_xyz(lenses(:,i),xyz_lens)
        xyz_axis = matmul(beam_grid%inv_basis,axes(:,i))
        spec_chords%los(i)%lens = xyz_lens
        spec_chords%los(i)%axis = xyz_axis
        spec_chords%los(i)%sigma_pi = sigma_pi(i)
        spec_chords%los(i)%spot_size = spot_size(i)

        r0 = xyz_lens
        v0 = xyz_axis
        v0 = v0/norm2(v0)
        call line_basis(r0,v0,basis)

        call grid_intersect(r0,v0,length,r_enter,r_exit)
        if(length.le.0.d0) then
            WRITE(*,'("Channel ",i3," missed the beam grid")'),i
            cycle chan_loop
        endif

        if(spot_size(i).le.0.d0) then
            nc = 1
        else
            nc = 100
        endif

        !$OMP PARALLEL DO schedule(guided) private(ic,randomu,sqrt_rho,theta,r0, &
        !$OMP& length, r_enter, r_exit, j, tracks, ncell, ind)
        do ic=1,nc
            ! Uniformally sample within spot size
            call randu(randomu)
            sqrt_rho = sqrt(randomu(1))
            theta = 2*pi*randomu(2)
            r0(1) = 0.d0
            r0(2) = spot_size(i)*sqrt_rho*cos(theta)
            r0(3) = spot_size(i)*sqrt_rho*sin(theta)
            r0 = matmul(basis,r0) + xyz_lens

            call grid_intersect(r0, v0, length, r_enter, r_exit)
            call track(r_enter, v0, tracks, ncell)
            track_loop: do j=1, ncell
                ind = tracks(j)%ind
                !inds can repeat so add rather than assign
                !$OMP CRITICAL(read_chords_1)
                spec_chords%dlength(i,ind(1),ind(2),ind(3)) = &
                spec_chords%dlength(i,ind(1),ind(2),ind(3)) + tracks(j)%time/real(nc) !time == distance
                spec_chords%los_inter(ind(1),ind(2),ind(3)) = .True.
                !$OMP END CRITICAL(read_chords_1)
            enddo track_loop
        enddo
        !$OMP END PARALLEL DO
    enddo chan_loop

    if(inputs%verbose.ge.1) then
        write(*,'(T2,"FIDA/BES System: ",a)') trim(adjustl(system))
        write(*,'(T2,"Number of channels: ",i3)') spec_chords%nchan
        write(*,*) ''
    endif

    deallocate(lenses,axes,spot_size,sigma_pi)

end subroutine read_chords

subroutine read_npa
    !+ Reads the NPA geometry and stores the quantities in [[libfida:npa_chords]]
    integer(HID_T) :: fid, gid
    integer(HSIZE_T), dimension(2) :: dims
    logical :: path_valid

    real(Float64), dimension(:,:), allocatable :: a_tedge,a_redge,a_cent
    real(Float64), dimension(:,:), allocatable :: d_tedge,d_redge,d_cent
    integer, dimension(:), allocatable :: a_shape, d_shape
    character(len=20) :: system = ''

    real(Float64), parameter :: inv_4pi = (4.d0*pi)**(-1)
    real(Float64), dimension(3) :: xyz_a_tedge,xyz_a_redge,xyz_a_cent
    real(Float64), dimension(3) :: xyz_d_tedge,xyz_d_redge,xyz_d_cent
    real(Float64), dimension(3) :: eff_rd, rd, rd_d, r0, r0_d, v0
    real(Float64), dimension(3,3) :: basis, inv_basis
    real(Float64), dimension(50) :: xd, yd
    real(Float64), dimension(:,:,:,:,:), allocatable :: effrd
    real(Float64), dimension(:,:,:,:), allocatable :: phit
    real(Float64) :: total_prob, hh, hw, dprob, dx, dy, r
    integer :: ichan,i,j,k,ix,iy,d_index,nd,cnt
    integer :: error

    !!Initialize HDF5 interface
    call h5open_f(error)

    !!Open HDF5 file
    call h5fopen_f(inputs%geometry_file, H5F_ACC_RDWR_F, fid, error)

    !!Check if NPA group exists
    call h5ltpath_valid_f(fid, "/npa", .True., path_valid, error)
    if(.not.path_valid) then
        write(*,'(a)') 'NPA geometry is not in the geometry file'
        write(*,'(a)') 'Continuing without NPA diagnostics'
        inputs%calc_npa = 0
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
    allocate(a_shape(npa_chords%nchan))
    allocate(d_tedge(3, npa_chords%nchan))
    allocate(d_redge(3, npa_chords%nchan))
    allocate(d_cent(3,  npa_chords%nchan))
    allocate(d_shape(npa_chords%nchan))
    allocate(npa_chords%radius(npa_chords%nchan))
    allocate(npa_chords%det(npa_chords%nchan))
    allocate(npa_chords%phit(beam_grid%nx, &
                             beam_grid%ny, &
                             beam_grid%nz, &
                             npa_chords%nchan) )
    allocate(npa_chords%hit(beam_grid%nx, &
                            beam_grid%ny, &
                            beam_grid%nz) )
    npa_chords%hit = .False.

    allocate(effrd(3,beam_grid%nx,&
                     beam_grid%ny,&
                     beam_grid%nz,&
                     npa_chords%nchan) )
    allocate(phit(beam_grid%nx,&
                  beam_grid%ny,&
                  beam_grid%nz,&
                  npa_chords%nchan) )
    effrd = 0.d0
    phit = 0.d0

    dims = [3,spec_chords%nchan]
    call h5ltread_dataset_double_f(gid,"/npa/radius", npa_chords%radius, dims(2:2), error)
    call h5ltread_dataset_int_f(gid, "/npa/a_shape", a_shape, dims(2:2), error)
    call h5ltread_dataset_double_f(gid, "/npa/a_tedge", a_tedge, dims, error)
    call h5ltread_dataset_double_f(gid, "/npa/a_redge", a_redge, dims, error)
    call h5ltread_dataset_double_f(gid, "/npa/a_cent",  a_cent, dims, error)

    call h5ltread_dataset_int_f(gid, "/npa/d_shape", d_shape, dims(2:2), error)
    call h5ltread_dataset_double_f(gid, "/npa/d_tedge", d_tedge, dims, error)
    call h5ltread_dataset_double_f(gid, "/npa/d_redge", d_redge, dims, error)
    call h5ltread_dataset_double_f(gid, "/npa/d_cent",  d_cent, dims, error)

    !!Close NPA group
    call h5gclose_f(gid, error)

    !!Close file id
    call h5fclose_f(fid, error)

    !!Close HDF5 interface
    call h5close_f(error)

    ! Define detector/aperture shape
    npa_chords%det%detector%shape = d_shape
    npa_chords%det%aperture%shape = a_shape

    if(inputs%verbose.ge.1) then
        write(*,'(T2,a)') "Calculating hit probabilities for NPA channels"
    endif
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
        cnt = 0
        ! For each grid point find the probability of hitting the detector given an isotropic source
        !$OMP PARALLEL DO schedule(guided) collapse(3) private(i,j,k,ix,iy,total_prob,eff_rd,r0,r0_d, &
        !$OMP& rd_d,rd,d_index,v0,dprob,r)
        do k=1,beam_grid%nz
            do j=1,beam_grid%ny
                do i=1,beam_grid%nx
                    cnt = cnt+1
                    total_prob = 0.d0
                    eff_rd = eff_rd*0.d0
                    r0 = [beam_grid%xc(i),beam_grid%yc(j),beam_grid%zc(k)]
                    r0_d = matmul(inv_basis,r0-xyz_d_cent)
                    do ix = 1, nd
                        do iy = 1, nd
                            rd_d = [xd(ix),yd(iy),0.d0]
                            rd = matmul(basis,rd_d) + xyz_d_cent
                            v0 = rd - r0
                            d_index = 0
                            call hit_npa_detector(r0,v0,d_index,det=ichan)
                            if(d_index.ne.0) then
                                r = norm2(rd_d - r0_d)**2
                                dprob = (dx*dy) * inv_4pi * r0_d(3)/(r*sqrt(r))
                                eff_rd = eff_rd + dprob*rd
                                total_prob = total_prob + dprob
                            endif
                        enddo !yd loop
                    enddo !xd loop
                    if(total_prob.gt.0.0) then
                        eff_rd = eff_rd/total_prob
                        phit(i,j,k,ichan) = total_prob
                        call xyz_to_uvw(eff_rd,effrd(:,i,j,k,ichan))
                        npa_chords%phit(i,j,k,ichan)%p = total_prob
                        npa_chords%phit(i,j,k,ichan)%eff_rd = eff_rd
                        npa_chords%hit(i,j,k) = .True.
                    endif
                    if(inputs%verbose.ge.2) then
                        WRITE(*,'(T4,"Channel: ",i3," ",f7.2,"% completed",a,$)') &
                                 ichan, cnt/real(beam_grid%ngrid)*100,char(13)
                    endif
                enddo !x loop
            enddo !y loop
        enddo !z loop
        !$OMP END PARALLEL DO

        total_prob = sum(npa_chords%phit(:,:,:,ichan)%p)
        if(total_prob.le.0.d0) then
            WRITE(*,'("Channel ",i3," missed the beam grid")'),ichan
            cycle chan_loop
        endif

    enddo chan_loop
    if(inputs%verbose.ge.1) write(*,'(50X,a)') ""

    deallocate(phit,effrd)
    deallocate(a_shape,a_cent,a_redge,a_tedge)
    deallocate(d_shape,d_cent,d_redge,d_tedge)

end subroutine read_npa

subroutine read_equilibrium
    !+ Reads in the interpolation grid, plasma parameters, and fields
    !+ and stores the quantities in [[libfida:inter_grid]] and [[libfida:equil]]
    integer(HID_T) :: fid, gid
    integer(HSIZE_T), dimension(2) :: dims

    integer :: impc
    integer :: error

    integer, dimension(:,:), allocatable :: p_mask, f_mask

    !!Initialize HDF5 interface
    call h5open_f(error)

    !!Open HDF5 file
    call h5fopen_f(inputs%equilibrium_file, H5F_ACC_RDONLY_F, fid, error)

    !!Open PLASMA group
    call h5gopen_f(fid, "/plasma", gid, error)

    !!Read in interpolation grid
    call h5ltread_dataset_int_scalar_f(gid, "/plasma/nr", inter_grid%nr, error)
    call h5ltread_dataset_int_scalar_f(gid, "/plasma/nz", inter_grid%nz, error)

    allocate(inter_grid%r(inter_grid%nr),inter_grid%z(inter_grid%nz))
    allocate(inter_grid%r2d(inter_grid%nr,inter_grid%nz))
    allocate(inter_grid%z2d(inter_grid%nr,inter_grid%nz))
    allocate(p_mask(inter_grid%nr,inter_grid%nz))
    allocate(f_mask(inter_grid%nr,inter_grid%nz))

    dims = [inter_grid%nr, inter_grid%nz]
    call h5ltread_dataset_double_f(gid, "/plasma/r", inter_grid%r, dims(1:1), error)
    call h5ltread_dataset_double_f(gid, "/plasma/z", inter_grid%z, dims(2:2), error)
    call h5ltread_dataset_double_f(gid, "/plasma/r2d", inter_grid%r2d, dims, error)
    call h5ltread_dataset_double_f(gid, "/plasma/z2d", inter_grid%z2d, dims, error)

    inter_grid%dr = abs(inter_grid%r(2)-inter_grid%r(1))
    inter_grid%dz = abs(inter_grid%z(2)-inter_grid%z(1))
    inter_grid%da = inter_grid%dr*inter_grid%dz

    if(inputs%verbose.ge.1) then
        write(*,'(a)') '---- Interpolation grid settings ----'
        write(*,'(T2,"Nr: ",i3)') inter_grid%nr
        write(*,'(T2,"Nz: ",i3)') inter_grid%nz
        write(*,'(T2,"dA: ", f4.2," [cm^2]")') inter_grid%da
        write(*,*) ''
    endif

    !!Read in plasma parameters
    allocate(equil%plasma(inter_grid%nr,inter_grid%nz))

    call h5ltread_dataset_double_f(gid, "/plasma/dene", equil%plasma%dene, dims, error)
    call h5ltread_dataset_double_f(gid, "/plasma/te", equil%plasma%te, dims, error)
    call h5ltread_dataset_double_f(gid, "/plasma/ti", equil%plasma%ti, dims, error)
    call h5ltread_dataset_double_f(gid, "/plasma/zeff", equil%plasma%zeff, dims, error)
    call h5ltread_dataset_double_f(gid, "/plasma/vr", equil%plasma%vr, dims, error)
    call h5ltread_dataset_double_f(gid, "/plasma/vt", equil%plasma%vt, dims, error)
    call h5ltread_dataset_double_f(gid, "/plasma/vz", equil%plasma%vz, dims, error)
    call h5ltread_dataset_int_f(gid, "/plasma/mask", p_mask, dims,error)

    impc = inputs%impurity_charge
    where(equil%plasma%zeff.lt.1.0)
        equil%plasma%zeff = 1
    endwhere

    where(equil%plasma%zeff.gt.impc)
        equil%plasma%zeff = impc
    endwhere

    where(equil%plasma%dene.lt.0.0)
        equil%plasma%dene = 0.0
    endwhere

    where(equil%plasma%te.lt.0.0)
        equil%plasma%te = 0.0
    endwhere

    where(equil%plasma%ti.lt.0.0)
        equil%plasma%ti = 0.0
    endwhere

    equil%plasma%denimp = ((equil%plasma%zeff-1.d0)/(impc*(impc-1.d0)))*equil%plasma%dene
    equil%plasma%denp = equil%plasma%dene - impc*equil%plasma%denimp

    !!Close PLASMA group
    call h5gclose_f(gid, error)

    !!Open FIELDS group
    call h5gopen_f(fid, "/fields", gid, error)

    allocate(equil%fields(inter_grid%nr,inter_grid%nz))

    !!Read in electromagnetic fields
    call h5ltread_dataset_double_f(gid, "/fields/br", equil%fields%br, dims, error)
    call h5ltread_dataset_double_f(gid, "/fields/bt", equil%fields%bt, dims, error)
    call h5ltread_dataset_double_f(gid, "/fields/bz", equil%fields%bz, dims, error)
    call h5ltread_dataset_double_f(gid, "/fields/er", equil%fields%er, dims, error)
    call h5ltread_dataset_double_f(gid, "/fields/et", equil%fields%et, dims, error)
    call h5ltread_dataset_double_f(gid, "/fields/ez", equil%fields%ez, dims, error)
    call h5ltread_dataset_int_f(gid, "/fields/mask", f_mask, dims,error)

    !!Close FIELDS group
    call h5gclose_f(gid, error)

    !!Close file
    call h5fclose_f(fid, error)

    !!Close HDF5 interface
    call h5close_f(error)

    allocate(equil%mask(inter_grid%nr,inter_grid%nz))
    equil%mask = 0.d0
    where ((p_mask.eq.1).and.(f_mask.eq.1)) equil%mask = 1.d0

end subroutine read_equilibrium

subroutine read_f(fid, error)
    !+ Reads in the fast-ion distribution function and stores the quantities in [[libfida:fbm]]
    integer(HID_T), intent(inout) :: fid
        !+ HDF5 file ID
    integer, intent(out)          :: error
        !+ Error code

    integer(HSIZE_T), dimension(4) :: dims
    real(Float64) :: dummy(1)

    if(inputs%verbose.ge.1) then
        write(*,'(a)') '---- Fast-ion distribution settings ----'
    endif

    call h5ltread_dataset_int_scalar_f(fid,"/nenergy", fbm%nenergy, error)
    call h5ltread_dataset_int_scalar_f(fid,"/npitch", fbm%npitch, error)
    call h5ltread_dataset_int_scalar_f(fid,"/nr", fbm%nr, error)
    call h5ltread_dataset_int_scalar_f(fid,"/nz", fbm%nz, error)

    if((fbm%nr.ne.inter_grid%nr).or.(fbm%nz.ne.inter_grid%nz)) then
        write(*,'(a)') "READ_F: Distribution file has incompatable grid dimensions"
        stop
    endif

    allocate(fbm%energy(fbm%nenergy), fbm%pitch(fbm%npitch), fbm%r(fbm%nr), fbm%z(fbm%nz))
    allocate(fbm%denf(fbm%nr, fbm%nz))
    allocate(fbm%f(fbm%nenergy, fbm%npitch, fbm%nr, fbm%nz))

    dims = [fbm%nenergy, fbm%npitch, fbm%nr, fbm%nz]
    call h5ltread_dataset_double_f(fid, "/energy", fbm%energy, dims(1:1), error)
    call h5ltread_dataset_double_f(fid, "/pitch", fbm%pitch, dims(2:2), error)
    call h5ltread_dataset_double_f(fid, "/r", fbm%r, dims(3:3), error)
    call h5ltread_dataset_double_f(fid, "/z", fbm%z, dims(4:4), error)
    call h5ltread_dataset_double_f(fid, "/denf",fbm%denf, dims(3:4), error)
    call h5ltread_dataset_double_f(fid, "/f", fbm%f, dims, error)
    equil%plasma%denf = fbm%denf

    fbm%dE = abs(fbm%energy(2) - fbm%energy(1))
    fbm%dp = abs(fbm%pitch(2) - fbm%pitch(1))
    fbm%dr = abs(fbm%r(2) - fbm%r(1))
    fbm%dz = abs(fbm%z(2) - fbm%z(1))

    dummy = minval(fbm%energy)
    fbm%emin = dummy(1)
    dummy = maxval(fbm%energy)
    fbm%emax = dummy(1)
    fbm%e_range = fbm%emax - fbm%emin
    dummy = minval(fbm%pitch)
    fbm%pmin = dummy(1)
    dummy = maxval(fbm%pitch)
    fbm%pmax = dummy(1)
    fbm%p_range = fbm%pmax - fbm%pmin

    if(inputs%verbose.ge.1) then
        write(*,'(T2,"Distribution type: ",a)') "Fast-ion Density Function F(energy,pitch,R,Z)"
        write(*,'(T2,"Nenergy = ",i3)'),fbm%nenergy
        write(*,'(T2,"Npitch  = ",i3)'),fbm%npitch
        write(*,'(T2,"Energy range = [",f5.2,",",f6.2,"]")'),fbm%emin,fbm%emax
        write(*,'(T2,"Pitch  range = [",f5.2,",",f5.2,"]")'),fbm%pmin,fbm%pmax
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
    integer(Int32) :: i,j,ii,ir,iz
    real(Float64) :: phi,phi_enter,phi_exit,delta_phi,r_ratio
    real(Float64), dimension(3) :: uvw,ri,vi,e1_xyz,e2_xyz,C_xyz
    integer(Int32), dimension(1) :: minpos
    real(Float64), dimension(:), allocatable :: weight
    real(Float64), dimension(:), allocatable :: r, z, vr, vt, vz
    real(Float64), dimension(:), allocatable :: energy, pitch
    integer(Int32), dimension(:), allocatable :: orbit_class
    type(LocalEMFields) :: fields
    integer :: cnt,num
    logical :: inp
    integer(Int32) :: npart,nrep
    character(len=32) :: dist_type_name = ''

    if(inputs%verbose.ge.1) then
        write(*,'(a)') '---- Fast-ion distribution settings ----'
    endif

    call h5ltread_dataset_int_scalar_f(fid, "/nparticle", npart, error)
    call h5ltread_dataset_int_scalar_f(fid, "/nclass", particles%nclass, error)

    nrep = ceiling(dble(inputs%n_fida)/npart)
    particles%nparticle = int(nrep*npart)

    !!ALLOCATE SPACE
    allocate(particles%fast_ion(particles%nparticle))
    allocate(r(npart))
    allocate(z(npart))
    allocate(vr(npart))
    allocate(vt(npart))
    allocate(vz(npart))
    allocate(weight(npart))
    allocate(orbit_class(npart))

    dims(1) = npart
    call h5ltread_dataset_double_f(fid, "/r", r, dims, error)
    call h5ltread_dataset_double_f(fid, "/z", z, dims, error)
    call h5ltread_dataset_double_f(fid, "/weight", weight, dims, error)
    call h5ltread_dataset_int_f(fid, "/class", orbit_class, dims, error)

    if(inputs%dist_type.eq.2) then
        dist_type_name = "Guiding Center Monte Carlo"
        allocate(energy(npart))
        allocate(pitch(npart))
        call h5ltread_dataset_double_f(fid, "/energy", energy, dims, error)
        call h5ltread_dataset_double_f(fid, "/pitch", pitch, dims, error)
    else
        dist_type_name = "Full Orbit Monte Carlo"
        call h5ltread_dataset_double_f(fid, "/vr", vr, dims, error)
        call h5ltread_dataset_double_f(fid, "/vt", vt, dims, error)
        call h5ltread_dataset_double_f(fid, "/vz", vz, dims, error)
    endif

    cnt=0
    e1_xyz = matmul(beam_grid%inv_basis,[1.0,0.0,0.0])
    e2_xyz = matmul(beam_grid%inv_basis,[0.0,1.0,0.0])
    !$OMP PARALLEL DO schedule(guided) private(i,ii,j,ir,iz,minpos,fields,uvw,phi,ri,vi, &
    !$OMP& delta_phi,phi_enter,phi_exit,r_ratio,C_xyz)
    particle_loop: do i=1,npart
        if(inputs%verbose.ge.2) then
            WRITE(*,'(f7.2,"% completed",a,$)') cnt/real(npart)*100,char(13)
        endif
        call in_plasma([r(i),0.d0,z(i)],inp,machine_coords=.True.)
        if(.not.inp) cycle particle_loop
        do ii=1,nrep
            j = int((i-1)*nrep + ii)
            if(inputs%dist_type.eq.2) then
                !! Transform to full orbit
                uvw = [r(i), 0.d0, z(i)]
                call get_fields(fields,pos = uvw, machine_coords=.True.)
                call gyro_correction(fields, uvw, energy(i), pitch(i), ri, vi)
                particles%fast_ion(j)%r = sqrt(ri(1)**2 + ri(2)**2)
                particles%fast_ion(j)%z = ri(3)
                !r_ratio = particles%fast_ion(j)%r/r(i)
                r_ratio = 1.0
                phi = atan2(ri(2),ri(1))
                particles%fast_ion(j)%vr =  vi(1)*cos(phi) + vi(2)*sin(phi)
                particles%fast_ion(j)%vt = -vi(1)*sin(phi) + vi(2)*cos(phi)
                particles%fast_ion(j)%vz =  vi(3)
                particles%fast_ion(j)%pitch =  pitch(i)
            else
                particles%fast_ion(j)%r = r(i)
                particles%fast_ion(j)%z = z(i)
                particles%fast_ion(j)%vr = vr(i)
                particles%fast_ion(j)%vt = vt(i)
                particles%fast_ion(j)%vz = vz(i)
                r_ratio = 1.0
            endif
            particles%fast_ion(j)%vabs = sqrt(particles%fast_ion(j)%vr**2 + &
                                              particles%fast_ion(j)%vt**2 + &
                                              particles%fast_ion(j)%vz**2)
            particles%fast_ion(j)%energy = v2_to_E_per_amu*inputs%ab*particles%fast_ion(j)%vabs**2
            particles%fast_ion(j)%class = orbit_class(i)

            phi_enter = 0.0
            phi_exit = 0.0
            call uvw_to_xyz([0.d0, 0.d0, particles%fast_ion(j)%z], C_xyz)
            call circle_grid_intersect(C_xyz,e1_xyz,e2_xyz,particles%fast_ion(j)%r,phi_enter,phi_exit)
            delta_phi = phi_exit-phi_enter
            if(delta_phi.gt.0) then
                particles%fast_ion(j)%cross_grid = .True.
            else
                particles%fast_ion(j)%cross_grid = .False.
            endif
            particles%fast_ion(j)%phi_enter = phi_enter
            particles%fast_ion(j)%delta_phi = delta_phi
            particles%fast_ion(j)%weight = (r_ratio*weight(i)/nrep) * &
                                           (delta_phi/(2*pi))/beam_grid%dv

            minpos = minloc(abs(inter_grid%r - particles%fast_ion(j)%r))
            ir = minpos(1)
            minpos = minloc(abs(inter_grid%z - particles%fast_ion(j)%z))
            iz = minpos(1)
            !$OMP CRITICAL(mc_denf)
            equil%plasma(ir,iz)%denf = equil%plasma(ir,iz)%denf + &
                                       (r_ratio*weight(i)/nrep) / &
                                       (2*pi*particles%fast_ion(j)%r*inter_grid%da)
            !$OMP END CRITICAL(mc_denf)
        enddo
        cnt=cnt+1
    enddo particle_loop
    !$OMP END PARALLEL DO

    num = count(particles%fast_ion%cross_grid)
    if(num.le.0) then
        write(*,'(a)') 'READ_MC: No mc particles in beam grid'
        stop
    endif

    if(inputs%verbose.ge.1) then
        write(*,'(T2,"Distribution type: ",a)') dist_type_name
        write(*,'(T2,"Number of mc particles: ",i9)') npart
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

subroutine read_cross(fid, grp, cross)
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

    dim3 = [n_max, m_max, tables%H_H_cx%nenergy]
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

end subroutine read_cross

subroutine read_rates(fid, grp, b_amu, t_amu, rates)
    !+ Reads in a reaction rates table from file
    !+ and puts it into a [[AtomicRates]] type
    integer(HID_T), intent(in)              :: fid
        !+ HDF5 file ID
    character(len=*), intent(in)            :: grp
        !+ HDF5 group to read from
    real(Float64), dimension(2), intent(in) :: b_amu
        !+ Atomic masses of "beam" species (beam ion and thermal ion)
    real(Float64), intent(in)               :: t_amu
        !+ Atomic mass of "target" species (thermal ion)
    type(AtomicRates), intent(inout)        :: rates
        !+ Atomic reaction rates

    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(4) :: dim4
    integer(HSIZE_T), dimension(5) :: dim5
    logical :: path_valid
    integer :: i, j, n, n_max, m_max, error
    integer :: n_bt_amu, tt_ind, bt_ind, drank
    real(Float64) :: emin,emax,tmin,tmax,rmin
    real(Float64) :: bt_min, tt_min, tt_dum, bt_dum
    real(Float64), dimension(2) :: bt_amu, tt_amu
    real(Float64), dimension(:,:), allocatable :: dummy2
    real(Float64), dimension(:,:,:,:), allocatable :: dummy4
    real(Float64), dimension(:,:,:,:,:), allocatable :: dummy5

    call h5ltread_dataset_int_scalar_f(fid, grp//"/n_bt_amu", n_bt_amu, error)
    allocate(dummy2(2, n_bt_amu))
    dim2 = [2, n_bt_amu]
    call h5ltread_dataset_double_f(fid, grp//"/bt_amu", dummy2, dim2, error)

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

    bt_ind = 1
    tt_ind = 1
    bt_amu = [b_amu(1), t_amu]
    tt_amu = [b_amu(2), t_amu]
    bt_min = norm2(bt_amu - dummy2(:,1))
    tt_min = norm2(tt_amu - dummy2(:,1))
    do i=2, n_bt_amu
        bt_dum = norm2(bt_amu - dummy2(:,i))
        tt_dum = norm2(tt_amu - dummy2(:,i))
        if(bt_dum.lt.bt_min) then
            bt_min = bt_dum
            bt_ind = i
        endif
        if(tt_dum.lt.tt_min) then
            tt_min = tt_dum
            tt_ind = i
        endif
    enddo
    rates%ab(1) = dummy2(1,bt_ind)
    rates%ab(2) = dummy2(1,tt_ind)

    deallocate(dummy2)

    allocate(rates%log_pop(&
                    rates%m_max, &
                    rates%n_max, &
                    rates%nenergy, &
                    rates%ntemp, 2))
    allocate(rates%log_depop(&
                    rates%n_max, &
                    rates%nenergy, &
                    rates%ntemp, 2))
    rates%log_pop = 0.d0
    rates%log_depop = 0.d0

    !!Read CX
    call h5ltpath_valid_f(fid, grp//"/cx", .True., path_valid, error)
    if(path_valid) then
        call h5ltget_dataset_ndims_f(fid, grp//"/cx", drank, error)
        if(drank.eq.4) then
            allocate(dummy4(n_max, &
                           rates%nenergy, &
                           rates%ntemp, n_bt_amu))
            dim4 = [n_max, rates%nenergy, rates%ntemp,n_bt_amu]
            call h5ltread_dataset_double_f(fid, grp//"/cx", dummy4, dim4, error)
            do j=1,rates%ntemp
                do i=1,rates%nenergy
                    do n=1,rates%n_max
                        rates%log_depop(n,i,j,1) = dummy4(n,i,j,bt_ind)
                        rates%log_depop(n,i,j,2) = dummy4(n,i,j,tt_ind)
                    enddo
                enddo
            enddo
            deallocate(dummy4)
        endif
        if(drank.eq.5) then
            allocate(dummy5(n_max, m_max, &
                           rates%nenergy, &
                           rates%ntemp, n_bt_amu))
            dim5 = [n_max, m_max, rates%nenergy, rates%ntemp,n_bt_amu]
            call h5ltread_dataset_double_f(fid, grp//"/cx", dummy5, dim5, error)
            do j=1,rates%ntemp
                do i=1,rates%nenergy
                    do n=1,rates%n_max
                        rates%log_depop(n,i,j,1) = sum(dummy5(n,:,i,j,bt_ind))
                        rates%log_depop(n,i,j,2) = sum(dummy5(n,:,i,j,tt_ind))
                    enddo
                enddo
            enddo
            deallocate(dummy5)
        endif
    endif

    !!Read ionization
    call h5ltpath_valid_f(fid, grp//"/ionization", .True., path_valid, error)
    if(path_valid) then
        allocate(dummy4(n_max, &
                       rates%nenergy, &
                       rates%ntemp, n_bt_amu))
        dim4 = [n_max, rates%nenergy, rates%ntemp,n_bt_amu]
        call h5ltread_dataset_double_f(fid, grp//"/ionization", dummy4, dim4, error)
        do j=1,rates%ntemp
            do i=1,rates%nenergy
                do n=1,rates%n_max
                    rates%log_depop(n,i,j,1) = rates%log_depop(n,i,j,1) + &
                                               dummy4(n,i,j,bt_ind)
                    rates%log_depop(n,i,j,2) = rates%log_depop(n,i,j,2) + &
                                               dummy4(n,i,j,tt_ind)
                enddo
            enddo
        enddo
        deallocate(dummy4)
    endif

    !!Read excitation
    call h5ltpath_valid_f(fid, grp//"/excitation", .True., path_valid, error)
    if(path_valid) then
        allocate(dummy5(n_max, m_max,&
                       rates%nenergy, &
                       rates%ntemp, n_bt_amu))
        dim5 = [n_max, m_max, rates%nenergy, rates%ntemp,n_bt_amu]
        call h5ltread_dataset_double_f(fid, grp//"/excitation", dummy5, dim5, error)
        do j=1,rates%ntemp
            do i=1,rates%nenergy
                rates%log_pop(:,:,i,j,1) = transpose(dummy5(1:nlevs,1:nlevs,i,j,bt_ind))
                rates%log_pop(:,:,i,j,2) = transpose(dummy5(1:nlevs,1:nlevs,i,j,tt_ind))
                do n=1,rates%n_max
                    rates%log_depop(n,i,j,1) = rates%log_depop(n,i,j,1) + &
                                               sum(dummy5(n,:,i,j,bt_ind))
                    rates%log_depop(n,i,j,2) = rates%log_depop(n,i,j,2) + &
                                               sum(dummy5(n,:,i,j,tt_ind))
                enddo
            enddo
        enddo
        deallocate(dummy5)
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

end subroutine read_rates

subroutine read_tables
    !+ Reads in atomic tables from file and stores them in [[libfida:tables]]
    integer(HID_T) :: fid
    integer(HSIZE_T), dimension(2) :: dim2
    integer :: error

    integer :: n_max, m_max
    character(len=4) :: impname
    real(Float64) :: imp_amu
    real(Float64), dimension(2) :: b_amu
    real(Float64), dimension(:,:), allocatable :: dummy2

    if(inputs%verbose.ge.1) then
        write(*,'(a)') "---- Atomic tables settings ----"
    endif

    !!Initialize HDF5 interface
    call h5open_f(error)

    !!Open HDF5 file
    call h5fopen_f(inputs%tables_file, H5F_ACC_RDONLY_F, fid, error)

    !!Read Hydrogen-Hydrogen CX Cross Sections
    call read_cross(fid,"/cross/H_H",tables%H_H_cx)

    !!Read Hydrogen-Hydrogen Rates
    b_amu = [inputs%ab, inputs%ai]
    call read_rates(fid,"/rates/H_H",b_amu, inputs%ai, tables%H_H)
    inputs%ab = tables%H_H%ab(1)
    inputs%ai = tables%H_H%ab(2)

    !!Read Hydrogen-Electron Rates
    call read_rates(fid,"/rates/H_e",b_amu, e_amu, tables%H_e)

    !!Read Hydrogen-Impurity rates
    impname = ''
    select case (inputs%impurity_charge)
        case (5)
            impname = "B5"
            imp_amu = B5_amu
        case (6)
            impname = "C6"
            imp_amu = C6_amu
        case DEFAULT
            impname = "Aq"
            imp_amu = 2.d0*inputs%impurity_charge
    end select

    call read_rates(fid,"/rates/H_"//trim(adjustl(impname)), b_amu, imp_amu, tables%H_Aq)

    !!Read Einstein coefficients
    call h5ltread_dataset_int_scalar_f(fid,"/rates/spontaneous/n_max", n_max, error)
    call h5ltread_dataset_int_scalar_f(fid,"/rates/spontaneous/m_max", m_max, error)
    allocate(dummy2(n_max,m_max))
    dim2 = [n_max, m_max]
    call h5ltread_dataset_double_f(fid,"/rates/spontaneous/einstein",dummy2, dim2, error)
    tables%einstein(:,:) = transpose(dummy2(1:nlevs,1:nlevs))
    deallocate(dummy2)

    !!Close file
    call h5fclose_f(fid, error)

    !!Close HDF5 interface
    call h5close_f(error)

    if(inputs%verbose.ge.1) then
        write(*,'(T2,"Maximum n/m: ",i2)') nlevs
        write(*,'(T2,"Beam/Fast-ion mass: ",f6.3," [amu]")') inputs%ab
        write(*,'(T2,"Thermal/Bulk-ion mass: ",f6.3," [amu]")') inputs%ai
        write(*,'(T2,"Impurity mass: ",f6.3," [amu]")') imp_amu
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
    integer :: error, i, npart

    character(charlim) :: filename
    real(Float64), dimension(:,:), allocatable :: ri
    real(Float64), dimension(:,:), allocatable :: vi
    real(Float64), dimension(3) :: xyz,uvw,v_uvw


    npart = birth%cnt-1
    allocate(ri(3,npart))
    allocate(vi(3,npart))

    do i=1,npart
        ! Convert position to rzphi
        xyz = birth%ri(:,i)
        call xyz_to_uvw(xyz,uvw)
        ri(1,i) = sqrt(uvw(1)*uvw(1) + uvw(2)*uvw(2))
        ri(2,i) = uvw(3)
        ri(3,i) = atan2(uvw(2),uvw(1))

        ! Convert velocity to rzphi
        v_uvw = matmul(beam_grid%basis, birth%vi(:,i))
        vi(1,i) = v_uvw(1)*cos(ri(3,i)) + v_uvw(2)*sin(ri(3,i))
        vi(2,i) = v_uvw(3)
        vi(3,i) = -v_uvw(1)*sin(ri(3,i)) + v_uvw(2)*cos(ri(3,i))
    enddo

    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_birth.h5"

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
    call h5ltmake_compressed_dataset_double_f(fid,"/ri", 2, dim2, ri, error)
    call h5ltmake_compressed_dataset_double_f(fid,"/vi", 2, dim2, vi, error)
    call h5ltmake_compressed_dataset_int_f(fid,"/ind", 2, dim2, birth%ind, error)
    call h5ltmake_compressed_dataset_int_f(fid,"/type", 1, dim2(2:2), birth%neut_type, error)

    !Add attributes
    call h5ltset_attribute_string_f(fid, "/n_birth","description", &
         "Number of birth mc particles deposited", error)
    call h5ltset_attribute_string_f(fid, "/dens", "description", &
         "Birth density: dens(beam_component,x,y,z)", error)
    call h5ltset_attribute_string_f(fid, "/dens", "units", &
         "fast-ions/(s*cm^3)", error)
    call h5ltset_attribute_string_f(fid, "/ri", "description", &
         "Fast-ion birth position in R-Z-Phi: ri([r,z,phi],particle)", error)
    call h5ltset_attribute_string_f(fid, "/ri", "units", "cm, radians", error)
    call h5ltset_attribute_string_f(fid, "/vi", "description", &
         "Fast-ion birth velocity in R-Z-Phi: vi([r,z,phi],particle)", error)
    call h5ltset_attribute_string_f(fid, "/vi", "units", "cm/s", error)
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

    deallocate(ri,vi)
    if(inputs%verbose.ge.1) then
        write(*,'(T4,a,a)') 'birth profile written to: ',trim(filename)
    endif

end subroutine write_birth_profile

subroutine write_dcx
    !+ Writes the direct charge exchange (DCX) neutrals and spectra to a HDF5 file
    integer(HID_T) :: fid
    integer(HSIZE_T), dimension(4) :: dims
    integer(HSIZE_T), dimension(1) :: d
    integer :: error

    character(charlim) :: filename
    character(15) :: spec_str
    integer :: i
    real(Float64), dimension(:),   allocatable :: lambda_arr
    real(Float64), dimension(:,:), allocatable :: dcx_spec

    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_dcx.h5"
    spec_str = ""
    if(inputs%calc_spec.ge.1) then
        spec_str = " spectra and"
        allocate(lambda_arr(inputs%nlambda))
        do i=1,inputs%nlambda
            lambda_arr(i) = (i-0.5)*inputs%dlambda + inputs%lambdamin ! [nm]
        enddo

        allocate(dcx_spec(inputs%nlambda,spec_chords%nchan))
        !! convert [Ph/(s*wavel_bin*cm^2*all_directions)] to [Ph/(s*nm*sr*m^2)]!
        dcx_spec = spec%bes(:,:,halo_type)/(inputs%dlambda)/(4.d0*pi)*1.d4
    endif

    !Open HDF5 interface
    call h5open_f(error)

    !Create file overwriting any existing file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, fid, error)

    !Write variables
    call write_beam_grid(fid, error)

    d(1) =1
    call h5ltmake_dataset_int_f(fid,"/nlevel", 0, d, [nlevs], error)
    dims = [nlevs, beam_grid%nx, beam_grid%ny, beam_grid%nz ]
    call h5ltmake_compressed_dataset_double_f(fid, "/dens", 4, dims, &
         neut%dens(:,halo_type,:,:,:), error)

    if(inputs%calc_spec.ge.1) then
        call h5ltmake_dataset_int_f(fid, "/nchan", 0, d, [spec_chords%nchan], error)
        call h5ltmake_dataset_int_f(fid, "/nlambda", 0, d, [inputs%nlambda], error)
        dims(1) = inputs%nlambda
        dims(2) = spec_chords%nchan
        call h5ltmake_compressed_dataset_double_f(fid, "/spec", 2, dims(1:2), &
             dcx_spec, error)
        call h5ltmake_compressed_dataset_double_f(fid, "/lambda", 1, dims(1:1), &
             lambda_arr, error)
        call h5ltmake_compressed_dataset_double_f(fid, "/radius", 1, dims(2:2), &
             spec_chords%radius, error)
    endif

    !Add attributes
    call h5ltset_attribute_string_f(fid,"/nlevel","description", &
         "Number of atomic energy levels", error)
    call h5ltset_attribute_string_f(fid,"/dens", "description", &
         "Direct Charge Exchange (DCX) neutral density: dcx(level,x,y,z)", error)
    call h5ltset_attribute_string_f(fid,"/dens","units","neutrals*cm^-3", error)

    if(inputs%calc_spec.ge.1) then
        call h5ltset_attribute_string_f(fid,"/nchan", "description", &
             "Number of channels", error)
        call h5ltset_attribute_string_f(fid,"/nlambda", "description", &
             "Number of wavelengths", error)
        call h5ltset_attribute_string_f(fid,"/spec","description", &
             "Direct Charge Exchange (DCX) beam emission: spec(lambda, chan)", error)
        call h5ltset_attribute_string_f(fid,"/spec","units","Ph/(s*nm*sr*m^2)",error)
        call h5ltset_attribute_string_f(fid,"/lambda","description", &
             "Wavelength array", error)
        call h5ltset_attribute_string_f(fid,"/lambda","units","nm", error)
        call h5ltset_attribute_string_f(fid,"/radius", "description", &
             "Line of sight radius at midplane or tangency point", error)
        call h5ltset_attribute_string_f(fid,"/radius","units","cm", error)
    endif

    call h5ltset_attribute_string_f(fid, "/", "version", version, error)
    call h5ltset_attribute_string_f(fid,"/","description", &
         "Direct Charge Exchange (DCX)"//trim(spec_str)//" neutral density calculated by FIDASIM", error)

    !Close file
    call h5fclose_f(fid, error)

    !Close HDF5 interface
    call h5close_f(error)

    if(inputs%calc_spec.ge.1) then
        deallocate(dcx_spec,lambda_arr)
    endif

    if(inputs%verbose.ge.1) then
        write(*,'(T4,a,a)') 'dcx written to: ',trim(filename)
    endif

end subroutine write_dcx

subroutine write_neutrals
    !+ Writes [[libfida:neut]] to a HDF5 file
    integer(HID_T) :: fid
    integer(HSIZE_T), dimension(4) :: dims
    integer(HSIZE_T), dimension(1) :: d
    integer :: error

    character(charlim) :: filename

    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_neutrals.h5"

    !Open HDF5 interface
    call h5open_f(error)

    !Create file overwriting any existing file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, fid, error)

    !Write variables
    call write_beam_grid(fid, error)
    dims = [nlevs, beam_grid%nx, beam_grid%ny, beam_grid%nz]
    d(1) =1
    call h5ltmake_dataset_int_f(fid,"/nlevel", 0, d, [nlevs], error)
    call h5ltmake_compressed_dataset_double_f(fid, "/fdens", 4, dims, &
         neut%dens(:,nbif_type,:,:,:), error)
    call h5ltmake_compressed_dataset_double_f(fid, "/hdens", 4, dims, &
         neut%dens(:,nbih_type,:,:,:), error)
    call h5ltmake_compressed_dataset_double_f(fid, "/tdens", 4, dims, &
         neut%dens(:,nbit_type,:,:,:), error)
    call h5ltmake_compressed_dataset_double_f(fid, "/halodens", 4, dims, &
         neut%dens(:,halo_type,:,:,:), error)

    !Write attributes
    call h5ltset_attribute_string_f(fid,"/nlevel","description", &
         "Number of atomic energy levels", error)
    call h5ltset_attribute_string_f(fid,"/fdens","description", &
         "Neutral density for the full energy component of the beam: fdens(level,x,y,z)", error)
    call h5ltset_attribute_string_f(fid,"/fdens","units","neutrals*cm^-3",error)

    call h5ltset_attribute_string_f(fid,"/hdens","description", &
         "Neutral density for the half energy component of the beam: hdens(level,x,y,z)", error)
    call h5ltset_attribute_string_f(fid,"/hdens","units","neutrals*cm^-3",error)

    call h5ltset_attribute_string_f(fid,"/tdens","description", &
         "Neutral density for the third energy component of the beam: tdens(level,x,y,z)", error)
    call h5ltset_attribute_string_f(fid,"/tdens","units","neutrals*cm^-3",error)

    call h5ltset_attribute_string_f(fid,"/halodens","description", &
         "Neutral density of the beam halo(including dcx): halodens(level,x,y,z)", error)
    call h5ltset_attribute_string_f(fid,"/halodens","units","neutrals*cm^-3",error)

    call h5ltset_attribute_string_f(fid, "/", "version", version, error)
    call h5ltset_attribute_string_f(fid,"/","description", &
         "Beam neutral density calculated by FIDASIM", error)

    !Close file
    call h5fclose_f(fid, error)

    !Close HDF5 interface
    call h5close_f(error)

    if(inputs%verbose.ge.1) then
        write(*,'(T4,a,a)') 'neutral density written to: ',trim(filename)
    endif

end subroutine write_neutrals

subroutine write_npa
    !+ Writes [[libfida:npa]] to a HDF5 file
    integer(HID_T) :: fid, gid
    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(1) :: d
    integer :: error

    integer, dimension(:), allocatable :: dcount
    real(Float64), dimension(:,:), allocatable :: ri, rf
    integer :: i, n
    character(charlim) :: filename = ''

    allocate(dcount(npa_chords%nchan))
    do i=1,npa_chords%nchan
        dcount(i) = count(npa%part%detector.eq.i)
    enddo

    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_npa.h5"

    !Open HDF5 interface
    call h5open_f(error)

    !Create file overwriting any existing file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, fid, error)

    !Write variables
    d(1) = 1
    dim2 = [fbm%nenergy, npa%nchan]
    call h5ltmake_dataset_int_f(fid,"/nenergy", 0, d, [fbm%nenergy], error)
    call h5ltmake_dataset_int_f(fid,"/nchan", 0, d, [npa%nchan], error)
    call h5ltmake_compressed_dataset_double_f(fid,"/energy",1,dim2(1:1),&
         npa%energy, error)
    call h5ltmake_compressed_dataset_double_f(fid,"/radius",1,dim2(2:2),&
         npa_chords%radius, error)
    call h5ltmake_compressed_dataset_double_f(fid,"/flux",2,dim2,npa%flux, error)
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
    call h5ltset_attribute_string_f(fid,"/flux", "description", &
         "Neutral flux: flux(energy,chan)", error)
    call h5ltset_attribute_string_f(fid,"/flux", "units","neutrals/(s*dE)", error)
    call h5ltset_attribute_string_f(fid,"/count","description", &
         "Number of particles that hit the detector: count(chan)", error)

    deallocate(dcount)

    if((npa%npart.ne.0).and.(inputs%calc_npa.ge.2)) then
        n = npa%npart
        allocate(ri(3,n),rf(3,n))
        ri(1,:) = npa%part(1:n)%xi
        ri(2,:) = npa%part(1:n)%yi
        ri(3,:) = npa%part(1:n)%zi
        rf(1,:) = npa%part(1:n)%xf
        rf(2,:) = npa%part(1:n)%yf
        rf(3,:) = npa%part(1:n)%zf

        !Create Group
        call h5gcreate_f(fid,"/particles",gid, error)
        call h5ltmake_dataset_int_f(gid, "nparticle", 0, d, [npa%npart], error)
        d(1) = npa%npart
        dim2 = [3, n]
        call h5ltmake_compressed_dataset_double_f(gid,"ri",2,dim2, ri, error)
        call h5ltmake_compressed_dataset_double_f(gid,"rf",2,dim2, rf, error)
        call h5ltmake_compressed_dataset_double_f(gid,"pitch",1,d, &
             npa%part(1:n)%pitch, error)
        call h5ltmake_compressed_dataset_double_f(gid,"energy",1,d,&
             npa%part(1:n)%energy, error)
        call h5ltmake_compressed_dataset_double_f(gid,"weight",1,d,&
             npa%part(1:n)%weight, error)
        call h5ltmake_compressed_dataset_int_f(gid,"detector",1,d,&
             npa%part(1:n)%detector, error)

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

        call h5ltset_attribute_string_f(fid,"/particles","coordinate_system", &
             "Right-handed cartesian",error)
        call h5ltset_attribute_string_f(fid,"/particles","description", &
             "Monte Carlo particles",error)

        !Close group
        call h5gclose_f(gid, error)
        deallocate(ri,rf)
    endif

    !Close file
    call h5fclose_f(fid, error)

    !Close HDF5 interface
    call h5close_f(error)

    if(inputs%verbose.ge.1) then
        write(*,'(T4,a,a)') 'NPA data written to: ',trim(filename)
    endif

end subroutine write_npa

subroutine write_spectra
    !+ Writes [[libfida:spectra]] to a HDF5 file
    integer(HID_T) :: fid
    integer(HSIZE_T), dimension(3) :: dims
    integer(HSIZE_T), dimension(1) :: d
    integer :: error

    character(charlim) :: filename
    integer :: i
    real(Float64), dimension(:), allocatable :: lambda_arr

    allocate(lambda_arr(inputs%nlambda))
    do i=1,inputs%nlambda
        lambda_arr(i) = (i-0.5)*inputs%dlambda + inputs%lambdamin
    enddo

    !! convert [Ph/(s*wavel_bin*cm^2*all_directions)] to [Ph/(s*nm*sr*m^2)]!
    spec%brems=spec%brems/(inputs%dlambda)/(4.d0*pi)*1.d4
    spec%bes=spec%bes/(inputs%dlambda)/(4.d0*pi)*1.d4
    spec%fida=spec%fida/(inputs%dlambda)/(4.d0*pi)*1.d4

    !! write to file
    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_spectra.h5"

    !Open HDF5 interface
    call h5open_f(error)

    !Create file overwriting any existing file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, fid, error)

    !Write variables
    d(1) = 1
    call h5ltmake_dataset_int_f(fid, "/nchan", 0, d, [spec_chords%nchan], error)
    call h5ltmake_dataset_int_f(fid, "/nlambda", 0, d, [inputs%nlambda], error)
    dims(1) = inputs%nlambda
    dims(2) = spec_chords%nchan
    dims(3) = particles%nclass
    call h5ltmake_compressed_dataset_double_f(fid, "/lambda", 1, dims(1:1), &
         lambda_arr, error)
    call h5ltmake_compressed_dataset_double_f(fid, "/radius", 1, dims(2:2), &
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
        !Write variables
        call h5ltmake_compressed_dataset_double_f(fid, "/brems", 2, &
             dims(1:2), spec%brems, error)
        !Add attributes
        call h5ltset_attribute_string_f(fid,"/brems","description", &
             "Visible Bremsstrahlung: brems(lambda,chan)", error)
        call h5ltset_attribute_string_f(fid,"/brems","units",&
             "Ph/(s*nm*sr*m^2)",error )
    endif

    if(inputs%calc_bes.ge.1) then
        !Write variables
        call h5ltmake_compressed_dataset_double_f(fid, "/full", 2, dims(1:2), &
             spec%bes(:,:,nbif_type), error)
        call h5ltmake_compressed_dataset_double_f(fid, "/half", 2, dims(1:2), &
             spec%bes(:,:,nbih_type), error)
        call h5ltmake_compressed_dataset_double_f(fid, "/third", 2, dims(1:2),&
             spec%bes(:,:,nbit_type), error)
        call h5ltmake_compressed_dataset_double_f(fid, "/halo", 2, dims(1:2), &
             spec%bes(:,:,halo_type), error)
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
        call h5ltset_attribute_string_f(fid,"/halo","description", &
             "Halo component of the beam emmision (includes dcx): halo(lambda,chan)", error)
        call h5ltset_attribute_string_f(fid,"/halo","units","Ph/(s*nm*sr*m^2)",error )
    endif

    if(inputs%calc_fida.ge.1) then
        !Write variables
        if(particles%nclass.le.1) then
            call h5ltmake_compressed_dataset_double_f(fid, "/fida", 2, &
                 dims(1:2), spec%fida(:,:,1), error)
            !Add attributes
            call h5ltset_attribute_string_f(fid,"/fida","description", &
                 "Fast-ion D-alpha (FIDA) emmision: fida(lambda,chan)", error)
        else
            call h5ltmake_dataset_int_f(fid,"/nclass", 0, d, [particles%nclass], error)
            call h5ltmake_compressed_dataset_double_f(fid, "/fida", 3, &
                 dims, spec%fida, error)
            !Add attributes
            call h5ltset_attribute_string_f(fid,"/fida","description", &
                 "Fast-ion D-alpha (FIDA) emmision: fida(lambda,chan,class)", error)
       endif
        call h5ltset_attribute_string_f(fid,"/fida","units","Ph/(s*nm*sr*m^2)",error )
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
    vpa_grid = 100*sqrt((((2.0d3)*e0)/(mass_u*inputs%ab))*e_grid)*p_grid ! [cm/s]
    vpe_grid = 100*sqrt((((2.0d3)*e0)/(mass_u*inputs%ab))*e_grid*(1.0-p_grid**2)) ![cm/s]

    !! define jacobian to convert between E-p to velocity
    allocate(jacobian(inputs%ne_wght,inputs%np_wght))
    jacobian = ((inputs%ab*mass_u)/(e0*1.0d3)) *vpe_grid/sqrt(vpa_grid**2 + vpe_grid**2)

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

subroutine read_neutrals
    !+ Reads neutral density from file and puts it in [[libfida:neut]]
    integer(HID_T) :: fid, gid
    integer(HSIZE_T), dimension(4) :: dims
    integer :: error,nx,ny,nz
    logical :: exis

    if(inputs%verbose.ge.1) then
        write(*,'(a)') '---- loading neutrals ----'
    endif

    inquire(file=inputs%neutrals_file,exist=exis)
    if(exis) then
        write(*,'(T2,"Neutrals file: ",a)') trim(inputs%neutrals_file)
        write(*,*) ''
    else
        write(*,'(a,a)') 'READ_NEUTRALS: Neutrals file does not exist: ',inputs%neutrals_file
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

    if((nx.ne.beam_grid%nx).or. &
       (ny.ne.beam_grid%ny).or. &
       (nz.ne.beam_grid%nz)) then
        write(*,'(a)') 'READ_NEUTRALS: Neutrals file has incompatable grid dimensions'
        stop
    endif

    dims = [nlevs, nx, ny, nz]
    call h5ltread_dataset_double_f(fid,"/fdens", &
         neut%dens(:,nbif_type,:,:,:), dims, error)
    call h5ltread_dataset_double_f(fid,"/hdens", &
         neut%dens(:,nbih_type,:,:,:), dims, error)
    call h5ltread_dataset_double_f(fid,"/tdens", &
         neut%dens(:,nbit_type,:,:,:), dims, error)
    call h5ltread_dataset_double_f(fid,"/halodens", &
         neut%dens(:,halo_type,:,:,:), dims, error)

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

subroutine plane_intercept(l0, l, p0, n, p, t)
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
        !+ Line-plane intercept point
    real(Float64), intent(out)               :: t
        !+ "time" to intercept

    t = dot_product(p0 - l0, n)/dot_product(l, n)

    p = l0 + t*l

end subroutine plane_intercept

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
            write(*,'("IN_BOUNDARY: Unknown boundary shape: ",i2)'),bplane%shape
            stop
    END SELECT

end function in_boundary

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
        call plane_intercept(r0,v0,npa_chords%det(i)%detector%origin, &
             npa_chords%det(i)%detector%basis(:,3),d,t_d)

        !! Find where trajectory crosses aperture plane
        call plane_intercept(r0,v0,npa_chords%det(i)%aperture%origin, &
             npa_chords%det(i)%aperture%basis(:,3),a,t_a)

        !! If both points are in plane boundaries and the
        !! particle is heading toward the detector then its a hit
        if(in_boundary(npa_chords%det(i)%aperture,a) .and. &
           in_boundary(npa_chords%det(i)%detector,d) .and. &
           (t_d.gt.0.0) ) then
            d_index = i
            if(present(rd)) rd = d
            exit detector_loop
        endif
    enddo detector_loop

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

subroutine grid_intersect(r0, v0, length, r_enter, r_exit, center_in, lwh_in)
    !+ Calculates a particles intersection length with the [[libfida:beam_grid]]
    real(Float64), dimension(3), intent(in)           :: r0
        !+ Initial position of particle [cm]
    real(Float64), dimension(3), intent(in)           :: v0
        !+ Velocity of particle [cm/s]
    real(Float64), intent(out)                        :: length
        !+ Intersection length [cm]
    real(Float64), dimension(3), intent(out)          :: r_enter
        !+ Point where particle enters [[libfida:beam_grid]]
    real(Float64), dimension(3), intent(out)          :: r_exit
        !+ Point where particle exits [[libfida:beam_grid]]
    real(Float64), dimension(3), intent(in), optional :: center_in
        !+ Alternative grid center
    real(Float64), dimension(3), intent(in), optional :: lwh_in
        !+ Alternative grid [length,width,height]

    real(Float64), dimension(3,6) :: ipnts
    real(Float64), dimension(3) :: vi
    real(Float64), dimension(3) :: center
    real(Float64), dimension(3) :: lwh
    integer, dimension(6) :: side_inter
    integer, dimension(2) :: ind
    integer :: i, j, nunique, ind1, ind2

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

subroutine circle_grid_intersect(r0, e1, e2, radius, phi_enter, phi_exit)
    !+ Calculates the intersection arclength of a circle with the [[libfida:beam_grid]]
    real(Float64), dimension(3), intent(in) :: r0
        !+ Position of center enter of the circle in beam grid coordinates [cm]
    real(Float64), dimension(3), intent(in) :: e1
        !+ Unit vector pointing towards (R, 0) (r,phi) position of the circle in beam grid coordinates
    real(Float64), dimension(3), intent(in) :: e2
        !+ Unit vector pointing towards (R, pi/2) (r,phi) position of the circle in beam grid coordinates
    real(Float64), intent(in)               :: radius
        !+ Radius of circle [cm]
    real(Float64), intent(out)              :: phi_enter
        !+ Phi value where the circle entered the [[libfida:beam_grid]] [rad]
    real(Float64), intent(out)              :: phi_exit
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

    phi_enter = 0.d0
    phi_exit = 0.d0
    if (count(inter).gt.2) return
    if(any(inter)) then
        phi_enter = minval(phi,inter)
        phi_exit = maxval(phi,inter)
        if(r0_ing.and.any(count(inter,1).ge.2)) then
            if((phi_exit - phi_enter) .lt. pi) then
                tmp = phi_enter
                phi_enter = phi_exit
                phi_exit = tmp + 2*pi
            endif
        else
            if((phi_exit - phi_enter) .gt. pi) then
                tmp = phi_enter
                phi_enter = phi_exit
                phi_exit = tmp + 2*pi
            endif
        endif
        if(approx_eq(phi_exit-phi_enter,pi,tol).and.r0_ing) then
            phi_enter = 0.0
            phi_exit = 2*pi
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
                phi_enter = 0.d0
                phi_exit = 2.d0*pi
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

subroutine get_position(ind, pos)
    !+ Get position `pos` given [[libfida:beam_grid]] indices `ind`
    integer(Int32), dimension(3), intent(in) :: ind
        !+ [[libfida:beam_grid]] indices
    real(Float64), dimension(3), intent(out) :: pos
        !+ Position [cm]

    pos(1) = beam_grid%xc(ind(1))
    pos(2) = beam_grid%yc(ind(2))
    pos(3) = beam_grid%zc(ind(3))

end subroutine get_position

subroutine track(rin, vin, tracks, ncell, los_intersect)
    !+ Computes the path of a neutral through the [[libfida:beam_grid]]
    real(Float64), dimension(3), intent(in)          :: rin
        !+ Initial position of particle
    real(Float64), dimension(3), intent(in)          :: vin
        !+ Initial velocity of particle
    type(ParticleTrack), dimension(:), intent(inout) :: tracks
        !+ Array of [[ParticleTrack]] type
    integer(Int32), intent(out)                      :: ncell
        !+ Number of cells that a particle crosses
    logical, intent(out), optional                   :: los_intersect
        !+ Indicator whether particle intersects a LOS in [[libfida:spec_chords]]

    integer :: cc, i, ii, mind
    integer, dimension(3) :: ind
    logical :: in_plasma1, in_plasma2, in_plasma_tmp, los_inter
    real(Float64) :: dT, dt1, inv_50
    real(Float64), dimension(3) :: dt_arr, dr
    real(Float64), dimension(3) :: vn, inv_vn
    real(Float64), dimension(3) :: ri, ri_tmp, ri_cell
    integer, dimension(3) :: sgn
    integer, dimension(3) :: gdims
    integer, dimension(1) :: minpos

    vn = vin ;  ri = rin ; sgn = 0 ; ncell = 0

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

    dr = beam_grid%dr*sgn
    inv_vn = 1/vn
    inv_50 = 1.0/50.0
    cc=1
    los_inter = .False.
    tracks%time = 0.d0
    tracks%flux = 0.d0
    call in_plasma(ri,in_plasma1)
    track_loop: do i=1,beam_grid%ntrack
        if(cc.gt.beam_grid%ntrack) exit track_loop

        if(spec_chords%los_inter(ind(1),ind(2),ind(3)).and.(.not.los_inter))then
            los_inter = .True.
        endif
        dt_arr = abs(( (ri_cell + 0.5*dr) - ri)*inv_vn)
        minpos = minloc(dt_arr)
        mind = minpos(1)
        dT = dt_arr(mind)
        ri_tmp = ri + dT*vn
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
    enddo track_loop
    ncell = cc-1
    if(present(los_intersect)) then
        los_intersect = los_inter
    endif

end subroutine track

!============================================================================
!---------------------------Interpolation Routines---------------------------
!============================================================================
subroutine interpol1D_coeff(xmin,dx,nx,xout,c,err)
    !+ Linear interpolation coefficients and index for a 1D grid y(x)
    real(Float64), intent(in)           :: xmin
        !+ Minimum abscissa value
    real(Float64), intent(in)           :: dx
        !+ Absissa spacing
    integer, intent(in)                 :: nx
        !+ Number of abscissa
    real(Float64), intent(in)           :: xout
        !+ Abscissa value to interpolate
    type(InterpolCoeffs1D), intent(out) :: c
        !+ Interpolation Coefficients
    integer, intent(out), optional      :: err
        !+ Error code

    real(Float64) :: x1, xp, b1, b2
    integer :: i, err_status

    err_status = 1
    xp = max(xout,xmin)
    i = floor((xp - xmin)/dx)+1

    if ((i.gt.0).and.(i.le.(nx-1))) then
        x1 = xmin + (i-1)*dx

        b2 = (xp - x1)/dx
        b1 = (1.0 - b2)

        c%i = i
        c%b1 = b1
        c%b2 = b2
        err_status = 0
    endif

    if(present(err)) err = err_status

end subroutine interpol1D_coeff

subroutine interpol1D_coeff_arr(x,xout,c,err)
    !+ Linear interpolation coefficients and index for a 1D grid y(x)
    real(Float64), dimension(:), intent(in) :: x
        !+ Abscissa values
    real(Float64), intent(in)               :: xout
        !+ Abscissa value to interpolate
    type(InterpolCoeffs1D), intent(out)     :: c
        !+ Interpolation Coefficients
    integer, intent(out), optional          :: err
        !+ Error code

    real(Float64) :: xmin, dx
    integer :: sx,err_status

    err_status = 1
    sx = size(x)
    xmin = x(1)
    dx = abs(x(2)-x(1))

    call interpol1D_coeff(xmin, dx, sx, xout, c, err_status)

    if(present(err)) err = err_status

end subroutine interpol1D_coeff_arr

subroutine interpol2D_coeff(xmin,dx,nx,ymin,dy,ny,xout,yout,c,err)
    !+ Bilinear interpolation coefficients and indicies for a 2D grid z(x,y)
    real(Float64), intent(in)           :: xmin
        !+ Minimum abscissa
    real(Float64), intent(in)           :: dx
        !+ Abscissa spacing
    integer, intent(in)                 :: nx
        !+ Number of abscissa
    real(Float64), intent(in)           :: ymin
        !+ Minimum ordinate
    real(Float64), intent(in)           :: dy
        !+ Ordinate spacing
    integer, intent(in)                 :: ny
        !+ Number of ordinates points
    real(Float64), intent(in)           :: xout
        !+ Abscissa value to interpolate
    real(Float64), intent(in)           :: yout
        !+ Ordinate value to interpolate
    type(InterpolCoeffs2D), intent(out) :: c
        !+ Interpolation Coefficients
    integer, intent(out), optional      :: err
        !+ Error code

    real(Float64) :: x1, x2, y1, y2, xp, yp
    integer :: i, j, err_status

    err_status = 1
    xp = max(xout,xmin)
    yp = max(yout,ymin)
    i = floor((xp-xmin)/dx)+1
    j = floor((yp-ymin)/dy)+1

    if (((i.gt.0).and.(i.le.(nx-1))).and.((j.gt.0).and.(j.le.(ny-1)))) then
        x1 = xmin + (i-1)*dx
        x2 = x1 + dx
        y1 = ymin + (j-1)*dy
        y2 = y1 + dy

        c%b11 = ((x2 - xp) * (y2 - yp))/(dx*dy)
        c%b21 = ((xp - x1) * (y2 - yp))/(dx*dy)
        c%b12 = ((x2 - xp) * (yp - y1))/(dx*dy)
        c%b22 = ((xp - x1) * (yp - y1))/(dx*dy)
        c%i = i
        c%j = j
        err_status = 0
    endif

    if(present(err)) err = err_status

end subroutine interpol2D_coeff

subroutine interpol2D_coeff_arr(x,y,xout,yout,c,err)
    !!Bilinear interpolation coefficients and indicies for a 2D grid z(x,y)
    real(Float64), dimension(:), intent(in) :: x
        !+ Abscissa values
    real(Float64), dimension(:), intent(in) :: y
        !+ Ordinate values
    real(Float64), intent(in)               :: xout
        !+ Abscissa value to interpolate
    real(Float64), intent(in)               :: yout
        !+ Ordinate value to interpolate
    type(InterpolCoeffs2D), intent(out)     :: c
        !+ Interpolation Coefficients
    integer, intent(out), optional          :: err
        !+ Error code

    real(Float64) :: xmin, ymin, dx, dy
    integer :: sx, sy, err_status

    err_status = 1
    sx = size(x)
    sy = size(y)
    xmin = x(1)
    ymin = y(1)
    dx = abs(x(2)-x(1))
    dy = abs(y(2)-y(1))

    call interpol2D_coeff(xmin, dx, sx, ymin, dy, sy, xout, yout, c, err_status)

    if(present(err)) err = err_status

end subroutine interpol2D_coeff_arr

subroutine interpol1D_arr(x, y, xout, yout, err, coeffs)
    !+ Performs linear interpolation on a uniform 1D grid y(x)
    real(Float64), dimension(:), intent(in)      :: x
        !+ The abscissa values of `y`
    real(Float64), dimension(:), intent(in)      :: y
        !+ Values at abscissa values `x`: y(x)
    real(Float64), intent(in)                    :: xout
        !+ Abscissa value to interpolate
    real(Float64), intent(out)                   :: yout
        !+ Interpolant: y(xout)
    integer, intent(out), optional               :: err
        !+ Error code
    type(InterpolCoeffs1D), intent(in), optional :: coeffs
        !+ Precomputed Linear Interpolation Coefficients

    type(InterpolCoeffs1D) :: c
    integer :: i, err_status

    err_status = 1
    if(present(coeffs)) then
        c = coeffs
        err_status = 0
    else
        call interpol_coeff(x,xout,c,err_status)
    endif

    if(err_status.eq.0) then
        i = c%i
        yout = c%b1*y(i) + c%b2*y(i+1)
    else
        yout = 0.d0
    endif

    if(present(err)) err = err_status

end subroutine interpol1D_arr

subroutine interpol2D_arr(x, y, z, xout, yout, zout, err, coeffs)
    !+ Performs bilinear interpolation on a 2D grid z(x,y)
    real(Float64), dimension(:), intent(in)   :: x
        !+ The abscissa values of `z`
    real(Float64), dimension(:), intent(in)   :: y
        !+ The ordinate values of `z`
    real(Float64), dimension(:,:), intent(in) :: z
        !+ Values at the abscissa/ordinates: z(x,y)
    real(Float64), intent(in)                 :: xout
        !+ The abscissa value to interpolate
    real(Float64), intent(in)                 :: yout
        !+ The ordinate value to interpolate
    real(Float64), intent(out)                :: zout
        !+ Interpolant: z(xout,yout)
    integer, intent(out), optional            :: err
        !+ Error code
    type(InterpolCoeffs2D), intent(in), optional :: coeffs
        !+ Precomputed Linear Interpolation Coefficients

    type(InterpolCoeffs2D) :: c
    integer :: i, j, err_status

    err_status = 1
    if(present(coeffs)) then
        c = coeffs
        err_status = 0
    else
        call interpol_coeff(x,y,xout,yout,c,err_status)
    endif

    if(err_status.eq.0) then
        i = c%i
        j = c%j
        zout = c%b11*z(i,j) + c%b12*z(i,j+1) + c%b21*z(i+1,j) + c%b22*z(i+1,j+1)
    else
        zout = 0.d0
    endif

    if(present(err)) err = err_status

end subroutine interpol2D_arr

subroutine interpol2D_2D_arr(x, y, z, xout, yout, zout, err, coeffs)
    !+ Performs bilinear interpolation on a 2D grid of 2D arrays z(:,:,x,y)
    real(Float64), dimension(:), intent(in)       :: x
        !+ The abscissa values of `z`
    real(Float64), dimension(:), intent(in)       :: y
        !+ The ordinate values of `z`
    real(Float64), dimension(:,:,:,:), intent(in) :: z
        !+ Values at the abscissa/ordinates: z(:,:,x,y)
    real(Float64), intent(in)                     :: xout
        !+ The abscissa value to interpolate
    real(Float64), intent(in)                     :: yout
        !+ The ordinate value to interpolate
    real(Float64), dimension(:,:), intent(out)    :: zout
        !+ Interpolant: z(:,:,xout,yout)
    integer, intent(out), optional                :: err
        !+ Error code
    type(InterpolCoeffs2D), intent(in), optional :: coeffs
        !+ Precomputed Linear Interpolation Coefficients

    type(InterpolCoeffs2D) :: c
    integer :: i, j, err_status

    err_status = 1
    if(present(coeffs)) then
        c = coeffs
        err_status = 0
    else
        call interpol_coeff(x,y,xout,yout,c,err_status)
    endif

    if(err_status.eq.0) then
        i = c%i
        j = c%j
        zout = c%b11*z(:,:,i,j) + c%b12*z(:,:,i,j+1) + c%b21*z(:,:,i+1,j) + c%b22*z(:,:,i+1,j+1)
    else
        zout = 0.0
    endif

    if(present(err)) err = err_status

end subroutine interpol2D_2D_arr

!=============================================================================
!-------------------------Profiles and Fields Routines------------------------
!=============================================================================
subroutine in_plasma(xyz, inp, machine_coords, coeffs, uvw_out)
    !+ Indicator subroutine to determine if a position is in a region where
    !+ the plasma parameter and fields are valid/known
    real(Float64), dimension(3), intent(in) :: xyz
        !+ Position in beam coordinates
    logical, intent(out)                    :: inp
        !+ Indicates whether plasma parameters and fields are valid/known
    logical, intent(in), optional           :: machine_coords
        !+ Indicates that xyz is in machine coordinates
    type(InterpolCoeffs2D), intent(out), optional      :: coeffs
        !+ Linear Interpolation coefficients used in calculation
    real(Float64), dimension(3), intent(out), optional :: uvw_out
        !+ Position in machine coordinates

    real(Float64), dimension(3) :: uvw
    type(InterpolCoeffs2D) :: c
    real(Float64) :: R, W, mask
    logical :: mc
    integer :: i, j, err

    err = 1
    mc = .False.
    if(present(machine_coords)) mc = machine_coords

    if(mc) then
        uvw = xyz
    else
        !! Convert to machine coordinates
        call xyz_to_uvw(xyz,uvw)
    endif

    R = sqrt(uvw(1)*uvw(1) + uvw(2)*uvw(2))
    W = uvw(3)

    !! Interpolate mask value
    call interpol_coeff(inter_grid%r, inter_grid%z, R, W, c, err)

    inp = .False.
    if(err.eq.0) then
        i = c%i
        j = c%j
        mask = c%b11*equil%mask(i,j)   + c%b12*equil%mask(i,j+1) + &
               c%b21*equil%mask(i+1,j) + c%b22*equil%mask(i+1,j+1)

        if((mask.ge.0.5).and.(err.eq.0)) then
            inp = .True.
        endif
    endif

    if(present(coeffs)) coeffs = c
    if(present(uvw_out)) uvw_out = uvw

end subroutine in_plasma

subroutine get_plasma(plasma, pos, ind)
    !+ Gets plasma parameters at position `pos` or [[libfida:beam_grid]] indices `ind`
    type(LocalProfiles), intent(out)                   :: plasma
        !+ Plasma parameters at `pos`/`ind`
    real(Float64), dimension(3), intent(in), optional  :: pos
        !+ Position in beam grid coordinates
    integer(Int32), dimension(3), intent(in), optional :: ind
        !+ [[libfida:beam_grid]] indices

    logical :: inp
    type(InterpolCoeffs2D) :: coeffs
    real(Float64), dimension(3) :: xyz, uvw, vrot_uvw
    real(Float64) :: phi, s, c
    integer :: i, j

    plasma%in_plasma = .False.

    if(present(ind)) call get_position(ind,xyz)
    if(present(pos)) xyz = pos

    call in_plasma(xyz,inp,.False.,coeffs,uvw)
    if(inp) then
        phi = atan2(uvw(2),uvw(1))
        i = coeffs%i
        j = coeffs%j

        plasma = coeffs%b11*equil%plasma(i,j)   + coeffs%b12*equil%plasma(i,j+1) + &
                 coeffs%b21*equil%plasma(i+1,j) + coeffs%b22*equil%plasma(i+1,j+1)

        s = sin(phi) ; c = cos(phi)
        vrot_uvw(1) = plasma%vr*c - plasma%vt*s
        vrot_uvw(2) = plasma%vr*s + plasma%vt*c
        vrot_uvw(3) = plasma%vz
        plasma%vrot = matmul(beam_grid%inv_basis,vrot_uvw)
        plasma%pos = xyz
        plasma%uvw = uvw
        plasma%in_plasma = .True.
        plasma%c = coeffs
    endif

end subroutine get_plasma

subroutine calc_perp_vectors(b, a, c)
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

subroutine get_fields(fields, pos, ind, machine_coords)
    !+ Gets electro-magnetic fields at position `pos` or [[libfida:beam_grid]] indices `ind`
    type(LocalEMFields),intent(out)                    :: fields
        !+ Electro-magnetic fields at `pos`/`ind`
    real(Float64), dimension(3), intent(in), optional  :: pos
        !+ Position in beam grid coordinates
    integer(Int32), dimension(3), intent(in), optional :: ind
        !+ [[libfida:beam_grid]] indices
    logical, intent(in), optional :: machine_coords
        !+ Indicates that pos is machine coordinates

    logical :: inp, mc
    real(Float64), dimension(3) :: xyz, uvw
    real(Float64), dimension(3) :: uvw_bfield, uvw_efield
    real(Float64), dimension(3) :: xyz_bfield, xyz_efield
    real(Float64) :: phi, s, c
    type(InterpolCoeffs2D) :: coeffs
    integer :: i, j

    fields%in_plasma = .False.

    if(present(ind)) call get_position(ind,xyz)
    if(present(pos)) xyz = pos

    mc = .False.
    if(present(machine_coords)) mc = machine_coords

    call in_plasma(xyz,inp,mc,coeffs,uvw)
    if(inp) then
        phi = atan2(uvw(2),uvw(1))
        i = coeffs%i
        j = coeffs%j

        fields = coeffs%b11*equil%fields(i,j) + coeffs%b12*equil%fields(i,j+1) + &
                 coeffs%b21*equil%fields(i+1,j) + coeffs%b22*equil%fields(i+1,j+1)

        phi = atan2(uvw(2),uvw(1))
        s = sin(phi) ; c = cos(phi)

        !Convert cylindrical coordinates to uvw
        uvw_bfield(1) = c*fields%br - s*fields%bt
        uvw_bfield(2) = s*fields%br + c*fields%bt
        uvw_bfield(3) = fields%bz
        uvw_efield(1) = c*fields%er - s*fields%et
        uvw_efield(2) = s*fields%er + c*fields%et
        uvw_efield(3) = fields%ez

        if(mc) then
            xyz_bfield = uvw_bfield
            xyz_efield = uvw_efield
        else
            !Represent fields in beam grid coordinates
            xyz_bfield = matmul(beam_grid%inv_basis,uvw_bfield)
            xyz_efield = matmul(beam_grid%inv_basis,uvw_efield)
        endif

        !Calculate field directions and magnitudes
        fields%b_abs = norm2(xyz_bfield)
        fields%e_abs = norm2(xyz_efield)
        if(fields%b_abs.gt.0.d0) fields%b_norm = xyz_bfield/fields%b_abs
        if(fields%e_abs.gt.0.d0) fields%e_norm = xyz_efield/fields%e_abs

        call calc_perp_vectors(fields%b_norm,fields%a_norm,fields%c_norm)

        fields%pos = xyz
        fields%uvw = uvw
        fields%in_plasma = .True.
        fields%machine_coords = mc
        fields%c = coeffs
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
    type(InterpolCoeffs2D), intent(in), optional       :: coeffs
        !+ Precomputed Linear Interpolation Coefficients

    real(Float64), dimension(3) :: xyz, uvw
    real(Float64) :: R, Z
    integer :: err

    if(present(coeffs)) then
        call interpol(fbm%r, fbm%z, fbm%f, R, Z, fbeam, err, coeffs)
        call interpol(fbm%r, fbm%z, fbm%denf, R, Z, denf, err, coeffs)
    else
        if(present(ind)) call get_position(ind,xyz)
        if(present(pos)) xyz = pos

        !! Convert to machine coordinates
        call xyz_to_uvw(xyz,uvw)
        R = sqrt(uvw(1)*uvw(1) + uvw(2)*uvw(2))
        Z = uvw(3)

        call interpol(fbm%r, fbm%z, fbm%f, R, Z, fbeam, err)
        call interpol(fbm%r, fbm%z, fbm%denf, R, Z, denf, err)
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
    type(InterpolCoeffs2D), intent(in), optional       :: coeffs
        !+ Precomputed Linear Interpolation Coefficients

    real(Float64), dimension(3) :: xyz, uvw
    real(Float64), dimension(fbm%nenergy,fbm%npitch)  :: fbeam
    integer(Int32), dimension(2) :: epi
    integer(Int32), dimension(1) :: dummy
    real(Float64) :: R, Z
    real(Float64) :: dE, dp
    integer :: err

    dummy = minloc(abs(fbm%energy - energy))
    epi(1) = dummy(1)
    dummy = minloc(abs(fbm%pitch - pitch))
    epi(2) = dummy(1)
    dE = abs(fbm%energy(epi(1)) - energy)
    dp = abs(fbm%pitch(epi(2)) - pitch)

    if((dE.le.fbm%dE).and.(dp.le.fbm%dp)) then
        if(present(coeffs)) then
            call interpol(inter_grid%r, inter_grid%z, fbm%f, R, Z, fbeam, err, coeffs)
        else
            if(present(ind)) call get_position(ind,xyz)
            if(present(pos)) xyz = pos

            !! Convert to machine coordinates
            call xyz_to_uvw(xyz,uvw)
            R = sqrt(uvw(1)*uvw(1) + uvw(2)*uvw(2))
            Z = uvw(3)

            call interpol(inter_grid%r, inter_grid%z, fbm%f, R, Z, fbeam, err)
        endif
        denf = fbeam(epi(1),epi(2))
    else
        denf = 0.0
    endif

end subroutine get_ep_denf

!=============================================================================
!--------------------------Result Storage Routines----------------------------
!=============================================================================
subroutine store_neutrals(ind, neut_type, dens, store_iter)
    !Store neutrals in [[libfida:neut]] at indices `ind`
    integer(Int32), dimension(3), intent(in) :: ind
        !+ [[libfida:beam_grid]] indices
    integer, intent(in)                      :: neut_type
        !+ Neutral type
    real(Float64), dimension(:), intent(in)  :: dens
        !+ Neutral density [neutrals/cm^3]
    logical, intent(in), optional            :: store_iter
        !+ Store DCX/Halo iteration density in [[libfida:halo_iter_dens]]
    logical :: iter

    if(present(store_iter)) then
        iter = store_iter
    else
        iter = .False.
    endif

    !$OMP CRITICAL(store_neutrals_1)
    if(iter) halo_iter_dens(neut_type) = halo_iter_dens(neut_type) + sum(dens)

    neut%dens(:,neut_type,ind(1),ind(2),ind(3)) = &
      neut%dens(:,neut_type,ind(1),ind(2),ind(3))+dens ![neutrals/cm^3]
    !$OMP END CRITICAL(store_neutrals_1)
end subroutine store_neutrals

subroutine store_births(ind, neut_type, dflux)
    !+ Store birth particles/density in [[libfida:birth]]
    integer(Int32), dimension(3), intent(in) :: ind
        !+ [[libfida:beam_grid]] indices
    integer(Int32), intent(in)               :: neut_type
        !+ Neutral type
    real(Float64), intent(in)                :: dflux
        !+ Deposited flux

    !$OMP CRITICAL(store_births_1)
    birth%dens( neut_type,ind(1),ind(2),ind(3))= &
     birth%dens(neut_type,ind(1),ind(2),ind(3)) + dflux
    !$OMP END CRITICAL(store_births_1)
end subroutine store_births

subroutine store_npa(det, ri, rf, vn, flux)
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

    type(LocalEMFields) :: fields
    real(Float64), dimension(3) :: uvw_ri, uvw_rf,vn_norm
    real(Float64) :: energy, pitch
    integer(Int32), dimension(1) :: ienergy
    type(NPAParticle), dimension(:), allocatable :: parts

    ! Convert to machine coordinates
    call xyz_to_uvw(ri,uvw_ri)
    call xyz_to_uvw(rf,uvw_rf)

    ! Calculate energy
    energy = inputs%ab*v2_to_E_per_amu*dot_product(vn,vn)

    ! Calculate pitch if distribution actually uses pitch
    if(inputs%dist_type.le.2) then
        call get_fields(fields, pos = ri)
        vn_norm = vn/norm2(vn)
        pitch = dot_product(fields%b_norm,vn_norm)
    else
        pitch = 0.d0
    endif

    !$OMP CRITICAL(store_npa_1)
    npa%npart = npa%npart + 1
    if(npa%npart.gt.npa%nmax) then
        allocate(parts(npa%npart-1))
        parts = npa%part
        deallocate(npa%part)
        npa%nmax = int(npa%nmax*2)
        allocate(npa%part(npa%nmax))
        npa%part(1:(npa%npart-1)) = parts
        deallocate(parts)
    endif
    npa%part(npa%npart)%detector = det
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
    npa%flux(ienergy(1),det) = npa%flux(ienergy(1),det) + flux/fbm%dE
    !$OMP END CRITICAL(store_npa_1)

end subroutine store_npa

!=============================================================================
!--------------------------Atomic Physics Routines----------------------------
!=============================================================================
subroutine neut_rates(denn, vi, vn, rates)
    !+ Get neutralization/cx rates
    real(Float64), dimension(nlevs), intent(in)  :: denn
        !+ Neutral density [cm^-3]
    real(Float64), dimension(3),     intent(in)  :: vi
        !+ Ion velocity [cm/s]
    real(Float64), dimension(3),     intent(in)  :: vn
        !+ Neutral velocity [cm/s]
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
    logEmin = tables%H_H_cx%logemin
    dlogE = tables%H_H_cx%dlogE
    neb = tables%H_H_cx%nenergy
    call interpol_coeff(logEmin,dlogE,neb,logeb,c,err)
    ebi = c%i
    if(err.eq.1) then
        write(*,'(a)') "NEUT_RATES: Eb out of range of H_H_cx table. Using nearest energy value."
        write(*,'("eb = ",f5.3," [keV]")') eb
        if(ebi.lt.1) then
            ebi=1
            c%b1=1.0 ; c%b2=0.0
        else
            ebi=neb
            c%b1=0.0 ; c%b2=1.0
        endif
    endif

    neut(:,:) = (c%b1*tables%H_H_cx%log_cross(:,:,ebi) + &
                 c%b2*tables%H_H_cx%log_cross(:,:,ebi+1))

    where (neut.lt.tables%H_H_cx%minlog_cross)
        neut = 0.d0
    elsewhere
        neut = 10.d0**neut
    end where

    rates=matmul(neut,denn)*vrel

end subroutine neut_rates

subroutine get_beam_cx_prob(ind, pos, v_ion, types, prob)
    !+ Get probability of a thermal ion charge exchanging with `types` neutrals
    integer(Int32), dimension(3), intent(in)     :: ind
        !+ [[libfida:beam_grid]] indices
    real(Float64), dimension(3), intent(in)      :: pos
        !+ Interaction position in beam grid coordinates
    real(Float64), dimension(3), intent(in)      :: v_ion
        !+ Ion velocity [cm/s]
    integer(Int32), dimension(:), intent(in)     :: types
        !+ Neutral types
    real(Float64), dimension(nlevs), intent(out) :: prob
        !+ Charge exchange rate/probability [1/s]

    integer :: ntypes, i, ii
    real(Float64), dimension(nlevs) :: rates, denn
    real(Float64), dimension(3) :: vhalo, vnbi ,vn

    vnbi = pos - nbi%src
    vnbi = vnbi/norm2(vnbi)*nbi%vinj

    ntypes = size(types)
    prob = 0

    do i=1,ntypes
        if((types(i).le.3).and.(types(i).ne.0)) then
            ! CX with full type'th energy NBI neutrals
            denn = neut%dens(:,types(i),ind(1),ind(2),ind(3))
            vn = vnbi/sqrt(real(types(i)))
            call neut_rates(denn,v_ion,vn,rates)
            prob = prob+rates
        else
            denn = neut%dens(:,types(i),ind(1),ind(2),ind(3))
            do ii=1,int(n_halo_neutrate)
                call mc_halo(ind, vhalo)
                call neut_rates(denn,v_ion,vhalo,rates)
                prob=prob+rates/n_halo_neutrate
            enddo
        endif
    enddo

end subroutine get_beam_cx_prob

subroutine get_rate_matrix(plasma, i_type, eb, rmat)
    !+ Gets rate matrix for use in [[libfida:colrad]]
    type(LocalProfiles), intent(in)                    :: plasma
        !+ Plasma parameters
    integer, intent(in)                                :: i_type
        !+ Ion type
    real(Float64), intent(in)                          :: eb
        !+ Ion energy [keV]
    real(Float64), dimension(nlevs,nlevs), intent(out) :: rmat
        !+ Rate matrix

    real(Float64) :: logEmin, dlogE, logeb
    real(Float64) :: logTmin, dlogT, logti, logte
    integer :: neb, nt
    type(InterpolCoeffs2D) :: c
    real(Float64) :: b11, b12, b21, b22, dene, denp, denimp
    real(Float64), dimension(nlevs,nlevs) :: H_H_pop, H_e_pop, H_Aq_pop
    real(Float64), dimension(nlevs) :: H_H_depop, H_e_depop, H_Aq_depop
    integer :: ebi, tii, tei, n, err_status

    H_H_pop = 0.d0
    H_e_pop = 0.d0
    H_Aq_pop = 0.d0
    H_H_depop = 0.d0
    H_e_depop = 0.d0
    H_Aq_depop = 0.d0
    denp = plasma%denp
    dene = plasma%dene
    denimp = plasma%denimp
    logeb = log10(eb)
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
    call interpol_coeff(logEmin, dlogE, neb, logTmin, dlogT, nt, &
                        logeb, logti, c, err_status)
    ebi = c%i
    tii = c%j
    b11 = c%b11
    b12 = c%b12
    b21 = c%b21
    b22 = c%b22
    if(err_status.eq.1) then
        write(*,'(a)') "GET_RATE_MATRIX: Eb or Ti out of range of H_H table. Setting H_H rates to zero"
        write(*,'("eb = ",f7.3," [keV]")') eb
        write(*,'("ti = ",f6.3," [keV]")') plasma%ti
        denp = 0.d0
    endif

    H_H_pop = (b11*tables%H_H%log_pop(:,:,ebi,tii,i_type)   + &
               b12*tables%H_H%log_pop(:,:,ebi,tii+1,i_type) + &
               b21*tables%H_H%log_pop(:,:,ebi+1,tii,i_type) + &
               b22*tables%H_H%log_pop(:,:,ebi+1,tii+1,i_type))
    where (H_H_pop.lt.tables%H_H%minlog_pop)
        H_H_pop = 0.d0
    elsewhere
        H_H_pop = denp * 10.d0**H_H_pop
    end where

    H_H_depop = (b11*tables%H_H%log_depop(:,ebi,tii,i_type)   + &
                 b12*tables%H_H%log_depop(:,ebi,tii+1,i_type) + &
                 b21*tables%H_H%log_depop(:,ebi+1,tii,i_type) + &
                 b22*tables%H_H%log_depop(:,ebi+1,tii+1,i_type))
    where (H_H_depop.lt.tables%H_H%minlog_depop)
        H_H_depop = 0.d0
    elsewhere
        H_H_depop = denp * 10.d0**H_H_depop
    end where

    !!H_e
    err_status = 1
    logEmin = tables%H_e%logemin
    logTmin = tables%H_e%logtmin
    dlogE = tables%H_e%dlogE
    dlogT = tables%H_e%dlogT
    neb = tables%H_e%nenergy
    nt = tables%H_e%ntemp
    call interpol_coeff(logEmin, dlogE, neb, logTmin, dlogT, nt, &
                        logeb, logte, c, err_status)
    ebi = c%i
    tei = c%j
    b11 = c%b11
    b12 = c%b12
    b21 = c%b21
    b22 = c%b22
    if(err_status.eq.1) then
        write(*,'(a)') "GET_RATE_MATRIX: Eb or Te out of range of H_e table. Setting H_e rates to zero"
        write(*,'("eb = ",f7.3," [keV]")') eb
        write(*,'("te = ",f6.3," [keV]")') plasma%te
        dene = 0.d0
    endif

    H_e_pop = (b11*tables%H_e%log_pop(:,:,ebi,tei,i_type)   + &
               b12*tables%H_e%log_pop(:,:,ebi,tei+1,i_type) + &
               b21*tables%H_e%log_pop(:,:,ebi+1,tei,i_type) + &
               b22*tables%H_e%log_pop(:,:,ebi+1,tei+1,i_type))
    where (H_e_pop.lt.tables%H_e%minlog_pop)
        H_e_pop = 0.d0
    elsewhere
        H_e_pop = dene * 10.d0**H_e_pop
    end where

    H_e_depop = (b11*tables%H_e%log_depop(:,ebi,tei,i_type)   + &
                 b12*tables%H_e%log_depop(:,ebi,tei+1,i_type) + &
                 b21*tables%H_e%log_depop(:,ebi+1,tei,i_type) + &
                 b22*tables%H_e%log_depop(:,ebi+1,tei+1,i_type))

    where (H_e_depop.lt.tables%H_e%minlog_depop)
        H_e_depop = 0.d0
    elsewhere
        H_e_depop = dene * 10.d0**H_e_depop
    end where

    !!H_Aq
    err_status = 1
    logEmin = tables%H_Aq%logemin
    logTmin = tables%H_Aq%logtmin
    dlogE = tables%H_Aq%dlogE
    dlogT = tables%H_Aq%dlogT
    neb = tables%H_Aq%nenergy
    nt = tables%H_Aq%ntemp
    call interpol_coeff(logEmin, dlogE, neb, logTmin, dlogT, nt, &
                        logeb, logti, c, err_status)
    ebi = c%i
    tii = c%j
    b11 = c%b11
    b12 = c%b12
    b21 = c%b21
    b22 = c%b22
    if(err_status.eq.1) then
        write(*,'(a)') "GET_RATE_MATRIX: Eb or Ti out of range of H_Aq table. Setting H_Aq rates to zero"
        write(*,'("eb = ",f7.3," [keV]")') eb
        write(*,'("ti = ",f6.3," [keV]")') plasma%ti
        denimp = 0.d0
    endif

    H_Aq_pop = (b11*tables%H_Aq%log_pop(:,:,ebi,tii,i_type)   + &
                b12*tables%H_Aq%log_pop(:,:,ebi,tii+1,i_type) + &
                b21*tables%H_Aq%log_pop(:,:,ebi+1,tii,i_type) + &
                b22*tables%H_Aq%log_pop(:,:,ebi+1,tii+1,i_type))
    where (H_Aq_pop.lt.tables%H_Aq%minlog_pop)
        H_Aq_pop = 0.d0
    elsewhere
        H_Aq_pop = denimp * 10.d0**H_Aq_pop
    end where
    H_Aq_depop = (b11*tables%H_Aq%log_depop(:,ebi,tii,i_type)   + &
                  b12*tables%H_Aq%log_depop(:,ebi,tii+1,i_type) + &
                  b21*tables%H_Aq%log_depop(:,ebi+1,tii,i_type) + &
                  b22*tables%H_Aq%log_depop(:,ebi+1,tii+1,i_type))

    where (H_Aq_depop.lt.tables%H_Aq%minlog_depop)
        H_Aq_depop = 0.d0
    elsewhere
        H_Aq_depop = denimp * 10.d0**H_Aq_depop
    end where

    rmat = tables%einstein + H_H_pop + H_e_pop + H_Aq_pop
    do n=1,nlevs
        rmat(n,n) = -sum(tables%einstein(:,n)) - H_H_depop(n) - H_e_depop(n) - H_Aq_depop(n)
    enddo

end subroutine get_rate_matrix

subroutine colrad(plasma,i_type,vn,dt,states,dens,photons)
    !+ Evolve density of states in time `dt` via collisional radiative model
    type(LocalProfiles), intent(in)              :: plasma
        !+ Plasma parameters
    integer, intent(in)                          :: i_type
        !+ Ion/Neutral type (beam,thermal)
    real(Float64), dimension(:), intent(in)      :: vn
        !+ Neutral velocitiy [cm/s]
    real(Float64), intent(in)                    :: dt
        !+ Time interval [s]
    real(Float64), dimension(:), intent(inout)   :: states
        !+ Density of states
    real(Float64), dimension(nlevs), intent(out) :: dens
        !+ Density of neutrals
    real(Float64), intent(out)                   :: photons
        !+ Emitted photons(3->2)

    real(Float64), dimension(nlevs,nlevs) :: matrix  !! Matrix
    real(Float64) :: b_amu
    real(Float64) :: vnet_square    !! net velocity of neutrals squared
    real(Float64) :: eb             !! Energy of the fast neutral

    real(Float64), dimension(nlevs,nlevs) :: eigvec, eigvec_inv
    real(Float64), dimension(nlevs) :: eigval, coef
    real(Float64), dimension(nlevs) :: exp_eigval_dt
    real(Float64) :: iflux !!Initial total flux
    integer :: n

    photons=0.d0
    dens=0.d0

    iflux=sum(states)
    if(iflux.lt.colrad_threshold .and. inputs%calc_npa.eq.0)then
        return
    endif

    if(.not.plasma%in_plasma) then
        dens = states*dt
        return
    endif

    if(i_type.eq.beam_ion) then
        b_amu = inputs%ab
    else
        b_amu = inputs%ai
    endif
    vnet_square=dot_product(vn-plasma%vrot,vn-plasma%vrot)  ![cm/s]
    eb = v2_to_E_per_amu*b_amu*vnet_square ![kev]
    call get_rate_matrix(plasma, i_type, eb, matrix)

    call eigen(nlevs,matrix, eigvec, eigval)
    call matinv(eigvec, eigvec_inv)
    coef = matmul(eigvec_inv, states)!coeffs determined from states at t=0
    exp_eigval_dt = exp(eigval*dt)   ! to improve speed (used twice)
    do n=1,nlevs
        if(eigval(n).eq.0.0) eigval(n)=eigval(n)+1 !protect against dividing by zero
    enddo

    states = matmul(eigvec, coef * exp_eigval_dt)  ![neutrals/cm^3/s]!
    dens   = matmul(eigvec,coef*(exp_eigval_dt-1.d0)/eigval)

    if ((minval(states).lt.0).or.(minval(dens).lt.0)) then
        do n=1,nlevs
            if(states(n).lt.0) states(n)=0.d0
            if(dens(n).lt.0) dens(n)=0.d0
        enddo
    endif

    photons=dens(3)*tables%einstein(2,3) !! - [Ph/(s*cm^3)] - !!

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
    do while (dis.le.max_dis)
        call colrad(plasma,beam_ion,vi,dt,states,dens,photons)
        r0 = r0 + vi*dt
        dis = dis + dstep
        call get_plasma(plasma,pos=r0)
    enddo

end subroutine attenuate

subroutine spectrum(vecp, vi, fields, sigma_pi, photons, dlength, lambda, intensity)
    !+ Calculates doppler shift and stark splitting
    real(Float64), dimension(3), intent(in)        :: vecp
        !+ Vector directing towards optical head
    real(Float64), dimension(3), intent(in)        :: vi
        !+ Particle velocity
    type(LocalEMFields), intent(in)                :: fields
        !+ Electro-magnetic fields
    real(Float64), intent(in)                      :: sigma_pi
        !+ Sigma-pi ratio
    real(Float64), intent(in)                      :: photons
        !+ Photon density from [[libfida:colrad]]
    real(Float64), intent(in)                      :: dlength
        !+ LOS intersection length with [[libfida:beam_grid]] cell particle is in
    real(Float64), dimension(n_stark), intent(out) :: lambda
        !+ Wavelengths [nm]
    real(Float64), dimension(n_stark), intent(out) :: intensity
        !+ Spectra intensities [Ph/(s cm^2 starkline)]

    real(Float64), dimension(3) :: vp, vn
    real(Float64), dimension(3) :: bfield, efield
    real(Float64) :: E, cos_los_Efield, lambda_shifted
    integer, parameter, dimension(n_stark) :: stark_sign= +1*stark_sigma -1*stark_pi

    !! vector directing towards the optical head
    vp=vecp/norm2(vecp)

    ! Calculate Doppler shift
    vn=vi*0.01d0 ! [m/s]
    lambda_shifted = lambda0*(1.d0 + dot_product(vn,vp)/c0)

    !! Calculate Stark Splitting
    ! Calculate E-field
    bfield = fields%b_norm*fields%b_abs
    efield = fields%e_norm*fields%e_abs
    efield(1) = efield(1) +  vn(2)*bfield(3) - vn(3)*bfield(2)
    efield(2) = efield(2) - (vn(1)*bfield(3) - vn(3)*bfield(1))
    efield(3) = efield(3) +  vn(1)*bfield(2) - vn(2)*bfield(1)
    E = norm2(efield)

    !Stark Splitting
    lambda =  lambda_shifted + E * stark_wavel ![nm]

    !Intensities of stark components
    if (E .eq. 0.d0) then
        cos_los_Efield = 0.d0
    else
        cos_los_Efield = dot_product(vp,efield) / E
    endif

    intensity = stark_intens*(1.d0+ stark_sign* cos_los_Efield**2)
    !! E.g. mirrors may change the pi to sigma intensity ratio
    where (stark_sigma .eq. 1)
        intensity = intensity * sigma_pi
    endwhere

    !! normalize and multiply with photon density from colrad
    intensity = intensity/sum(intensity)*photons*dlength

endsubroutine spectrum

subroutine store_bes_photons(pos, vi, photons, neut_type)
    !+ Store BES photons in [[libfida:spectra]]
    real(Float64), dimension(3), intent(in) :: pos
        !+ Position of neutral in beam grid coordinates
    real(Float64), dimension(3), intent(in) :: vi
        !+ Velocitiy of neutral [cm/s]
    real(Float64), intent(in)               :: photons
        !+ Photons from [[libfida:colrad]] [Ph/(s*cm^3)]
    integer,intent(in)                      :: neut_type
        !+ Neutral type (full,half,third,halo)

    real(Float64), dimension(n_stark) :: lambda, intensity
    real(Float64) :: dlength, sigma_pi
    type(LocalEMFields) :: fields
    integer(Int32), dimension(3) :: ind
    real(Float64), dimension(3) :: vp
    integer :: ichan,i,bin

    call get_indices(pos,ind)
    call get_fields(fields,pos=pos)

    loop_over_channels: do ichan=1,spec_chords%nchan
        dlength = spec_chords%dlength(ichan,ind(1),ind(2),ind(3))
        if(dlength.le.0.0) cycle loop_over_channels
        sigma_pi = spec_chords%los(ichan)%sigma_pi
        vp = pos - spec_chords%los(ichan)%lens
        call spectrum(vp,vi,fields,sigma_pi,photons, &
                      dlength,lambda,intensity)

        loop_over_stark: do i=1,n_stark
            bin=floor((lambda(i)-inputs%lambdamin)/inputs%dlambda) + 1
            if (bin.lt.1) cycle loop_over_stark
            if (bin.gt.inputs%nlambda) cycle loop_over_stark
            !$OMP CRITICAL(bes_spectrum)
            spec%bes(bin,ichan,neut_type)= &
              spec%bes(bin,ichan,neut_type) + intensity(i)
            !$OMP END CRITICAL(bes_spectrum)
        enddo loop_over_stark
    enddo loop_over_channels

end subroutine store_bes_photons

subroutine store_fida_photons(pos, vi, photons, orbit_class)
    !+ Store fida photons in [[libfida:spectra]]
    real(Float64), dimension(3), intent(in) :: pos
        !+ Position of neutral in beam grid coordinates
    real(Float64), dimension(3), intent(in) :: vi
        !+ Velocitiy of neutral [cm/s]
    real(Float64), intent(in)               :: photons
        !+ Photons from [[libfida:colrad]] [Ph/(s*cm^3)]
    integer, intent(in), optional           :: orbit_class
        !+ Orbit class ID

    real(Float64), dimension(n_stark) :: lambda, intensity
    real(Float64) :: dlength, sigma_pi
    type(LocalEMFields) :: fields
    integer(Int32), dimension(3) :: ind
    real(Float64), dimension(3) :: vp
    integer :: ichan, i, bin, iclass

    if(present(orbit_class)) then
        iclass = orbit_class
    else
        iclass = 1
    endif

    call get_indices(pos,ind)
    call get_fields(fields,pos=pos)

    loop_over_channels: do ichan=1,spec_chords%nchan
        dlength = spec_chords%dlength(ichan,ind(1),ind(2),ind(3))
        if(dlength.le.0.0) cycle loop_over_channels
        sigma_pi = spec_chords%los(ichan)%sigma_pi
        vp = pos - spec_chords%los(ichan)%lens
        call spectrum(vp,vi,fields,sigma_pi,photons, &
                      dlength,lambda,intensity)

        loop_over_stark: do i=1,n_stark
            bin=floor((lambda(i)-inputs%lambdamin)/inputs%dlambda) + 1
            if (bin.lt.1) cycle loop_over_stark
            if (bin.gt.inputs%nlambda) cycle loop_over_stark
            !$OMP CRITICAL(fida_spectrum)
            spec%fida(bin,ichan,iclass)= &
              spec%fida(bin,ichan,iclass) + intensity(i)
            !$OMP END CRITICAL(fida_spectrum)
        enddo loop_over_stark
    enddo loop_over_channels

end subroutine store_fida_photons

subroutine store_fw_photons_at_chan(ichan,eind,pind,vp,vi,fields,dlength,sigma_pi,denf,photons)
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
    real(Float64) :: dlambda,intens_fac
    integer :: i,bin


    dlambda=(inputs%lambdamax_wght-inputs%lambdamin_wght)/inputs%nlambda_wght
    intens_fac = (1.d0)/(4.d0*pi*dlambda)
    call spectrum(vp,vi,fields,sigma_pi,photons, &
                  dlength,lambda,intensity)

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

subroutine store_fw_photons(eind, pind, pos, vi, denf, photons)
    !+ Store FIDA weight photons in [[libfida:fweight]]
    integer, intent(in)                     :: eind
        !+ Energy index
    integer, intent(in)                     :: pind
        !+ Pitch index
    real(Float64), dimension(3), intent(in) :: pos
        !+ Position of neutral
    real(Float64), dimension(3), intent(in) :: vi
        !+ Velocity of neutral [cm/s]
    real(Float64), intent(in)               :: denf
        !+ Fast-ion density [cm^-3]
    real(Float64), intent(in)               :: photons
        !+ Photons from [[libfida:colrad]] [Ph/(s*cm^3)]

    real(Float64) :: dlength, sigma_pi
    type(LocalEMFields) :: fields
    integer(Int32), dimension(3) :: ind
    real(Float64), dimension(3) :: vp
    integer :: ichan

    call get_indices(pos,ind)
    call get_fields(fields,pos=pos)

    loop_over_channels: do ichan=1,spec_chords%nchan
        dlength = spec_chords%dlength(ichan,ind(1),ind(2),ind(3))
        if(dlength.le.0.0) cycle loop_over_channels
        sigma_pi = spec_chords%los(ichan)%sigma_pi
        vp = pos - spec_chords%los(ichan)%lens
        call store_fw_photons_at_chan(ichan, eind, pind, &
             vp, vi, fields, dlength, sigma_pi, denf, photons)
    enddo loop_over_channels

end subroutine store_fw_photons

!=============================================================================
!---------------------------Monte Carlo Routines------------------------------
!=============================================================================
subroutine get_nlaunch(nr_markers,papprox,papprox_tot,nlaunch)
    !+ Sets the number of MC markers launched from each [[libfida:beam_grid]] cell
    integer(Int32), intent(in)                   :: nr_markers
        !+ Approximate total number of markers to launch
    real(Float64), dimension(:,:,:), intent(in)  :: papprox
        !+ [[libfida:beam_grid]] cell weights
    real(Float64), intent(in)                    :: papprox_tot
        !+ Total cell weights
    real(Float64), dimension(:,:,:), intent(out) :: nlaunch
        !+ Number of mc markers to launch for each cell: nlaunch(x,y,z)

    integer  :: i, j, k, cc
    real(Float64), dimension(:), allocatable :: randomu

    do i=1,1000
       nlaunch(:,:,:)=papprox(:,:,:)/papprox_tot*nr_markers*(1.+i*0.01)
       if(sum(nlaunch).gt.nr_markers) then
          exit
       endif
    enddo

    allocate(randomu(count(nlaunch.gt.0)))
    call randu(randomu)
    cc=1
    do k = 1, beam_grid%nz
        do j = 1, beam_grid%ny
            do i = 1, beam_grid%nx
                if(nlaunch(i,j,k).gt.0.)then
                    if(mod(nlaunch(i,j,k),1.).gt.randomu(cc))then
                        nlaunch(i,j,k)=nlaunch(i,j,k)+1.
                    endif
                    cc=cc+1
                endif
            enddo
        enddo
    enddo

    do k = 1, beam_grid%nz
        do j = 1, beam_grid%ny
            do i = 1, beam_grid%nx
                nlaunch(i,j,k)=floor(nlaunch(i,j,k))
            enddo
        enddo
    enddo
    deallocate(randomu)

end subroutine get_nlaunch

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

    sinus = sqrt(1.d0-pitch**2)
    vi_norm = (sinus*cos(gyroangle)*fields%a_norm + &
               pitch*fields%b_norm + &
               sinus*sin(gyroangle)*fields%c_norm)

end subroutine pitch_to_vec

subroutine gyro_step(vi, fields, r_gyro)
    !+ Calculates gyro-step
    real(Float64), dimension(3), intent(in)  :: vi
        !+ Ion velocity
    type(LocalEMFields), intent(in)          :: fields
        !+ Electro-magnetic fields
    real(Float64), dimension(3), intent(out) :: r_gyro
        !+ Gyro-step

    real(Float64), dimension(3) :: vxB,rg_uvw,uvw
    real(Float64) :: one_over_omega, phi, R, rg_r

    if(inputs%no_flr.eq.0) then
        one_over_omega=inputs%ab*mass_u/(fields%b_abs*e0)
        vxB = cross_product(vi,fields%b_norm)
        r_gyro = vxB*one_over_omega !points towards gyrocenter

        !! Second order correction approximation derived from
        !! Belova, E. V., N. N. Gorelenkov, and C. Z. Cheng.
        !! "Self-consistent equilibrium model of low aspect-ratio
        !! toroidal plasma with energetic beam ions."
        !! Physics of Plasmas (1994-present) 10.8 (2003): 3240-3251.
        !! Appendix A: Last equation
        uvw = fields%uvw
        R = sqrt(uvw(1)**2 + uvw(2)**2)
        phi = atan2(uvw(2),uvw(1))
        if(fields%machine_coords) then
            rg_uvw = r_gyro
        else
            rg_uvw = matmul(beam_grid%basis,r_gyro)
        endif
        rg_r = rg_uvw(1)*cos(phi) + rg_uvw(2)*sin(phi)
        r_gyro = r_gyro*(1 - rg_r/(2*R))
    else
        r_gyro = 0.d0
    endif

end subroutine gyro_step

subroutine gyro_correction(fields, rg, energy, pitch, rp, vp)
    !+ Calculates gyro correction for Guiding Center MC distribution calculation
    type(LocalEMFields), intent(in)          :: fields
        !+ Electromagnetic fields at guiding center
    real(Float64), dimension(3), intent(in)  :: rg
        !+ Gyro-center position
    real(Float64), intent(in)                :: energy
        !+ Energy of particle
    real(Float64), intent(in)                :: pitch
        !+ Particle pitch w.r.t the magnetic field
    real(Float64), dimension(3), intent(out) :: rp
        !+ Particle position
    real(Float64), dimension(3), intent(out) :: vp
        !+ Particle velocity

    real(Float64), dimension(3) :: vi_norm, r_step
    real(Float64), dimension(1) :: randomu
    real(Float64) :: vabs, phi

    vabs  = sqrt(energy/(v2_to_E_per_amu*inputs%ab))

    !! Sample gyroangle according to radius to counter-act geometric effect
    call randu(randomu)
    phi = 2*pi*randomu(1)

    !! Calculate velocity vector
    call pitch_to_vec(pitch, phi, fields, vi_norm)
    vp = vabs*vi_norm

    !! Move to particle location
    call gyro_step(vp, fields, r_step)
    rp = rg - r_step

end subroutine gyro_correction

subroutine mc_fastion(ind,rp,vi,denf)
    !+ Samples a Guiding Center Fast-ion distribution function at a given [[libfida:beam_grid]] index
    integer, dimension(3), intent(in)        :: ind
        !+ [[libfida:beam_grid]] index
    real(Float64), dimension(3), intent(out) :: rp
        !+ Fast-ion particle position [cm]
    real(Float64), dimension(3), intent(out) :: vi
        !+ Fast-ion particle velocity [cm/s]
    real(Float64), intent(out)               :: denf
        !+ Fast-ion density at guiding center

    type(LocalEMFields) :: fields
    real(Float64), dimension(fbm%nenergy,fbm%npitch) :: fbeam
    real(Float64), dimension(3) :: rg
    real(Float64) :: eb ,ptch
    real(Float64), dimension(3) :: randomu3
    integer, dimension(2,1) :: ep_ind

    call randu(randomu3)
    rg(1) = beam_grid%xc(ind(1)) + beam_grid%dr(1)*(randomu3(1) - 0.5)
    rg(2) = beam_grid%yc(ind(2)) + beam_grid%dr(2)*(randomu3(2) - 0.5)
    rg(3) = beam_grid%zc(ind(3)) + beam_grid%dr(3)*(randomu3(3) - 0.5)
    vi=0.d0
    denf=0.d0

    call get_fields(fields,pos=rg)
    if(.not.fields%in_plasma) return

    call get_distribution(fbeam,denf,pos=rg, coeffs=fields%c)
    call randind(fbeam,ep_ind)
    call randu(randomu3)
    eb = fbm%energy(ep_ind(1,1)) + fbm%dE*(randomu3(1)-0.5)
    ptch = fbm%pitch(ep_ind(2,1)) + fbm%dp*(randomu3(2)-0.5)
    call gyro_correction(fields,rg,eb,ptch,rp,vi)

end subroutine mc_fastion

subroutine mc_halo(ind,vhalo,ri,plasma_in)
    !+ Sample thermal Maxwellian distribution at [[libfida:beam_grid]] indices `ind`
    integer, dimension(3), intent(in)                  :: ind
        !+ [[libfida:beam_grid]] indices
    real(Float64), dimension(3), intent(out)           :: vhalo
        !+ Velocity [cm/s]
    real(Float64), dimension(3), intent(out), optional :: ri
        !+ Position in [[libfida:beam_grid]] cell
    type(LocalProfiles), intent(in), optional          :: plasma_in
        !+ Plasma parameters

    type(LocalProfiles) :: plasma
    real(Float64), dimension(3) :: random3

    if(.not.present(plasma_in)) then
        if(present(ri)) then
            call randu(random3)
            ri(1) = beam_grid%xc(ind(1)) + beam_grid%dr(1)*(random3(1) - 0.5)
            ri(2) = beam_grid%yc(ind(2)) + beam_grid%dr(2)*(random3(2) - 0.5)
            ri(3) = beam_grid%zc(ind(3)) + beam_grid%dr(3)*(random3(3) - 0.5)
            call get_plasma(plasma,pos=ri)
        else
            call get_plasma(plasma,ind=ind)
        endif
    else
        plasma=plasma_in
    endif

    call randn(random3)

    vhalo = plasma%vrot + sqrt(plasma%ti*0.5/(v2_to_E_per_amu*inputs%ai))*random3 !![cm/s]

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
        write(*,'(a)') "MC_NBI: Failed to find trajectory though aperture(s). Using beam centerline."
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
        write(*,'(a)') "MC_NBI: A beam neutral has started inside the plasma."
        write(*,'(a)') "Move the beam grid closer to the source to fix"
        stop
    endif

    !! Determine velocity of neutrals corrected by efrac
    vnbi = vnbi*nbi%vinj/sqrt(real(efrac))
end subroutine mc_nbi

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

    integer :: jj, ii, kk
    integer :: ncell
    type(ParticleTrack), dimension(beam_grid%ntrack) :: tracks
    integer, dimension(3) :: nl_birth
    type(LocalProfiles) :: plasma
    real(Float64), dimension(nlevs) :: states, dens
    real(Float64) :: photons, iflux
    integer(Int32), dimension(3) :: ind
    real(Float64), dimension(1) :: randomu
    integer, dimension(1) :: randi
    logical :: err

    if(inputs%verbose.ge.1) then
        write(*,'(T6,"# of markers: ",i9)') inputs%n_nbi
    endif

    !! # of injected neutrals = NBI power/energy_per_particle
    nneutrals=1.d6*nbi%pinj/ (1.d3*nbi%einj*e0 &
         *( nbi%current_fractions(1)      &
         +  nbi%current_fractions(2)/2.d0 &
         +  nbi%current_fractions(3)/3.d0 ) )

    nlaunch=real(inputs%n_nbi)

    !$OMP PARALLEL DO schedule(guided) &
    !$OMP& private(vnbi,rnbi,tracks,ncell,plasma,nl_birth,randi, &
    !$OMP& states,dens,iflux,photons,neut_type,jj,ii,kk,ind,err)
    loop_over_markers: do ii=1,inputs%n_nbi
        if(inputs%calc_birth.ge.1) then
            !! Determine the type of birth particle (1, 2, or 3)
            nl_birth = 0
            do kk=1,inputs%n_birth
                call randind(nbi%current_fractions,randi)
                nl_birth(randi(1)) = nl_birth(randi(1)) + 1
            enddo
        endif
        energy_fractions: do neut_type=1,3
            !! (type = 1: full energy, =2: half energy, =3: third energy
            call mc_nbi(vnbi,neut_type,rnbi,err)
            if(err) cycle energy_fractions

            call track(rnbi,vnbi,tracks,ncell)
            if(ncell.eq.0) cycle energy_fractions

            !! Solve collisional radiative model along track
            states=0.d0
            states(1)=nneutrals*nbi%current_fractions(neut_type)/beam_grid%dv
            loop_along_track: do jj=1,ncell
                iflux = sum(states)
                ind = tracks(jj)%ind
                call get_plasma(plasma,pos=tracks(jj)%pos)

                call colrad(plasma,beam_ion,vnbi,tracks(jj)%time,states,dens,photons)
                call store_neutrals(ind,neut_type,dens/nlaunch)
                tracks(jj)%flux = (iflux - sum(states))*beam_grid%dv/nlaunch

                if(inputs%calc_birth.ge.1) then
                    call store_births(ind,neut_type,tracks(jj)%flux)
                endif

                if((photons.gt.0.d0).and.(inputs%calc_bes.ge.1)) then
                    call store_bes_photons(tracks(jj)%pos,vnbi,photons/nlaunch,neut_type)
                endif
            enddo loop_along_track
            if(inputs%calc_birth.ge.1) then
                !! Sample according to deposited flux along neutral trajectory
                !$OMP CRITICAL(ndmc_birth)
                do kk=1,nl_birth(neut_type)
                    call randind(tracks(1:ncell)%flux,randi)
                    call randu(randomu)
                    birth%neut_type(birth%cnt) = neut_type
                    birth%ind(:,birth%cnt) = tracks(randi(1))%ind
                    birth%vi(:,birth%cnt) = vnbi
                    birth%ri(:,birth%cnt) = tracks(randi(1))%pos + &
                                            vnbi*(tracks(randi(1))%time*(randomu(1)-0.5))
                    birth%cnt = birth%cnt+1
                enddo
                !$OMP END CRITICAL(ndmc_birth)
            endif
        enddo energy_fractions
    enddo loop_over_markers
    !$OMP END PARALLEL DO

    if(nbi_outside.gt.0)then
         write(*,'(T4,a, f6.2)') 'Percent of markers outside the grid: ', &
                              100.*nbi_outside/(3.*inputs%n_nbi)
         if(sum(neut%dens).eq.0) stop 'Beam does not intersect the grid!'
    endif

end subroutine ndmc

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
    !! $OMP& ic, nc,randomu,sqrt_rho,theta,r0,plasma,gaunt,brems)
    loop_over_channels: do ichan=1,spec_chords%nchan
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

        do ic=1,nc
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
                if(max_length.gt.300) cycle loop_over_channels
            enddo

            ! Calculate bremsstrahlung along los
            do while (plasma%in_plasma)
                if(plasma%te.gt.0.0) then
                    gaunt = 5.542-(3.108-log(plasma%te))*(0.6905-0.1323/plasma%zeff)
                    brems = 7.57d-9*gaunt*plasma%dene**2*plasma%zeff/(lambda_arr &
                            *sqrt(plasma%te*1000.0))*exp(-h_planck*c0/(lambda_arr*plasma%te*1000.0)) &
                            *dlambda*(4.d0*pi)*1.d-4

                    spec%brems(:,ichan)= spec%brems(:,ichan) + (brems*dlength*1.d-2)/nc
                endif
                ! Take a step
                r0 = r0 + vi*dlength
                call get_plasma(plasma,pos=r0)
            enddo
        enddo

        if (inputs%verbose.eq.2)then
            WRITE(*,'(f7.2,"% completed",a,$)') 100*ichan/real(spec_chords%nchan),char(13)
        endif
    enddo loop_over_channels
    !! $OMP END PARALLEL DO

    deallocate(lambda_arr,brems)

end subroutine bremsstrahlung

subroutine dcx
    !+ Calculates Direct Charge Exchange (DCX) neutral density and spectra
    integer :: i,j,k !indices of cells
    integer :: idcx !! counter
    real(Float64), dimension(3) :: ri    !! start position
    real(Float64), dimension(3) :: vihalo
    integer,dimension(3) :: ind    !! actual cell
    integer,dimension(3) :: neut_types = [1,2,3]
    !! Determination of the CX probability
    type(LocalProfiles) :: plasma
    real(Float64), dimension(nlevs) :: denn    !!  neutral dens (n=1-4)
    real(Float64), dimension(nlevs) :: prob    !!  Prob. for CX
    !! Collisiional radiative model along track
    real(Float64), dimension(nlevs) :: states  ! Density of n-states
    integer :: ncell
    type(ParticleTrack), dimension(beam_grid%ntrack) :: tracks  !! Particle tracks
    integer :: jj       !! counter along track
    real(Float64):: tot_denn, photons  !! photon flux
    real(Float64), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz) :: papprox, nlaunch !! approx. density
    real(Float64) :: papprox_tot, ccnt, inv_ng

    halo_iter_dens(halo_type) = 0.d0
    papprox=0.d0
    papprox_tot=0.d0
    tot_denn=0.d0
    do k=1,beam_grid%nz
        do j=1,beam_grid%ny
            do i=1,beam_grid%nx
                ind = [i,j,k]
                call get_plasma(plasma,ind=ind)
                tot_denn = sum(neut%dens(:,nbif_type,i,j,k)) + &
                           sum(neut%dens(:,nbih_type,i,j,k)) + &
                           sum(neut%dens(:,nbit_type,i,j,k))
                papprox(i,j,k)= tot_denn*(plasma%denp-plasma%denf)
                if(plasma%in_plasma) papprox_tot=papprox_tot+papprox(i,j,k)
            enddo
        enddo
    enddo

    call get_nlaunch(inputs%n_dcx,papprox,papprox_tot,nlaunch)

    if(inputs%verbose.ge.1) then
       write(*,'(T6,"# of markers: ",i9)') int(sum(nlaunch))
    endif

    ccnt=0.d0
    inv_ng = 100.0/real(beam_grid%ngrid)
    loop_along_z: do k = 1, beam_grid%nz
        loop_along_y: do j = 1, beam_grid%ny
            loop_along_x: do i = 1, beam_grid%nx
                !! Loop over the markers
                !$OMP PARALLEL DO schedule(guided) private(idcx,ind,vihalo, &
                !$OMP& ri,tracks,ncell,prob,denn,states,jj,photons,plasma)
                loop_over_dcx: do idcx=1,int(nlaunch(i,j,k))
                    !! Calculate ri,vhalo and track
                    ind = [i, j, k]
                    call mc_halo(ind,vihalo,ri)
                    call track(ri,vihalo,tracks,ncell)
                    if(ncell.eq.0) cycle loop_over_dcx

                    !! Calculate CX probability
                    call get_beam_cx_prob(tracks(1)%ind,ri,vihalo,neut_types,prob)
                    if(sum(prob).le.0.) cycle loop_over_dcx

                    !! Solve collisional radiative model along track
                    call get_plasma(plasma,pos=tracks(1)%pos)

                    states = prob*(plasma%denp - plasma%denf)

                    loop_along_track: do jj=1,ncell
                        call get_plasma(plasma,pos=tracks(jj)%pos)

                        call colrad(plasma,thermal_ion,vihalo,tracks(jj)%time,states,denn,photons)
                        call store_neutrals(tracks(jj)%ind,halo_type,denn/nlaunch(i,j,k),plasma%in_plasma)

                        if((photons.gt.0.d0).and.(inputs%calc_bes.ge.1)) then
                          call store_bes_photons(tracks(jj)%pos,vihalo,photons/nlaunch(i,j,k),halo_type)
                        endif

                    enddo loop_along_track
                enddo loop_over_dcx
                !$OMP END PARALLEL DO
                ccnt=ccnt+1
                if (inputs%verbose.eq.2)then
                    WRITE(*,'(f7.2,"% completed",a,$)') ccnt*inv_ng,char(13)
                endif
            enddo loop_along_x
        enddo loop_along_y
    enddo loop_along_z

end subroutine dcx

subroutine halo
    !+ Calculates halo neutral density and spectra
    integer :: i,j,k !indices of cells
    integer :: ihalo !! counter
    real(Float64), dimension(3) :: ri    !! start position
    real(Float64), dimension(3) :: vihalo!! velocity bulk plasma ion
    integer,dimension(3) :: ind    !! actual cell
    !! Determination of the CX probability
    type(LocalProfiles) :: plasma
    real(Float64), dimension(nlevs) :: denn    !!  neutral dens (n=1-4)
    real(Float64), dimension(nlevs) :: prob    !!  Prob. for CX
    !! Collisiional radiative model along track
    real(Float64), dimension(nlevs) :: states  ! Density of n-states
    integer :: ncell
    type(ParticleTrack), dimension(beam_grid%ntrack) :: tracks  !! Particle Tracks
    integer :: jj       !! counter along track
    real(Float64) :: tot_denn, photons  !! photon flux
    real(Float64), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz) :: papprox, nlaunch !! approx. density
    real(Float64) :: papprox_tot, ccnt, inv_ng
    !! Halo iteration
    integer :: hh !! counters
    real(Float64) :: dcx_dens, halo_iteration_dens
    integer :: s1type  ! halo iteration
    integer :: s2type  ! halo iteration

    s1type = fida_type
    s2type = brems_type

    dcx_dens = halo_iter_dens(halo_type)
    if(dcx_dens.eq.0) then
        write(*,'(a)') 'HALO: Density of DCX-neutrals is zero'
        stop
    endif
    inv_ng = 100.0/real(beam_grid%ngrid)
    neut%dens(:,s1type,:,:,:) = neut%dens(:,halo_type,:,:,:)
    iterations: do hh=1,200
        papprox=0.d0
        papprox_tot=0.d0
        tot_denn=0.d0
        halo_iter_dens(s2type) = 0.d0
        do k=1,beam_grid%nz
            do j=1,beam_grid%ny
                do i=1,beam_grid%nx
                    ind = [i,j,k]
                    call get_plasma(plasma,ind=ind)
                    tot_denn = sum(neut%dens(:,s1type,i,j,k))
                    papprox(i,j,k)= tot_denn*(plasma%denp-plasma%denf)

                    if(plasma%in_plasma) papprox_tot=papprox_tot+papprox(i,j,k)
                enddo
            enddo
        enddo

        call get_nlaunch(inputs%n_halo,papprox,papprox_tot,nlaunch)

        if(inputs%verbose.ge.1) then
            write(*,'(T6,"# of markers: ",i9)') int(sum(nlaunch))
        endif
        ccnt=0.d0
        !$OMP PARALLEL DO schedule(guided) collapse(3) private(i,j,k,ihalo,ind,vihalo, &
        !$OMP& ri,tracks,ncell,prob,denn,states,jj,photons,plasma)
        loop_along_z: do k = 1, beam_grid%nz
            loop_along_y: do j = 1, beam_grid%ny
                loop_along_x: do i = 1, beam_grid%nx
                    !! Loop over the markers
                    loop_over_halos: do ihalo=1,int(nlaunch(i,j,k))
                        !! Calculate ri,vhalo and track
                        ind = [i, j, k]
                        call mc_halo(ind,vihalo,ri)
                        call track(ri,vihalo,tracks,ncell)
                        if(ncell.eq.0)cycle loop_over_halos

                        !! Calculate CX probability
                        call get_beam_cx_prob(tracks(1)%ind,ri,vihalo,[s1type],prob)
                        if(sum(prob).le.0.)cycle loop_over_halos

                        !! Solve collisional radiative model along track
                        call get_plasma(plasma,pos=tracks(1)%pos)

                        states = prob*plasma%denp

                        loop_along_track: do jj=1,ncell
                            call get_plasma(plasma,pos=tracks(jj)%pos)

                            call colrad(plasma,thermal_ion,vihalo,tracks(jj)%time,states,denn,photons)
                            call store_neutrals(tracks(jj)%ind,s2type, &
                                 denn/nlaunch(i,j,k),plasma%in_plasma)

                            if((photons.gt.0.d0).and.(inputs%calc_bes.ge.1)) then
                              call store_bes_photons(tracks(jj)%pos,vihalo,photons/nlaunch(i,j,k),halo_type)
                            endif

                        enddo loop_along_track
                    enddo loop_over_halos
                    ccnt=ccnt+1
                    if (inputs%verbose.eq.2)then
                        WRITE(*,'(f7.2,"% completed",a,$)') ccnt*inv_ng,char(13)
                    endif
                enddo loop_along_x
            enddo loop_along_y
        enddo loop_along_z
        !$OMP END PARALLEL DO

        halo_iteration_dens = halo_iter_dens(s2type)
        neut%dens(:,halo_type,:,:,:)= neut%dens(:,halo_type,:,:,:) &
                                           + neut%dens(:,s2type,:,:,:)
        neut%dens(:,s1type,:,:,:)= neut%dens(:,s2type,:,:,:)
        neut%dens(:,s2type,:,:,:)= 0.

        if(halo_iteration_dens/dcx_dens.gt.1)then
            write(*,'(a)') "HALO: Halo generation density exceeded DCX density. This shouldn't happen."
            exit iterations
        endif

        inputs%n_halo=int(inputs%n_dcx*halo_iteration_dens/dcx_dens)

        if(inputs%n_halo.lt.inputs%n_dcx*0.01)exit iterations
    enddo iterations
    !! set the neutral density in s1type(fida_type) and s2type (brems) to 0!
    neut%dens(:,s1type,:,:,:) = 0.d0
    neut%dens(:,s2type,:,:,:) = 0.d0

end subroutine halo

subroutine fida_f
    !+ Calculate FIDA emission using a Fast-ion distribution function F(E,p,r,z)
    integer :: i,j,k  !! indices  x,y,z of cells
    integer(kind=8) :: iion,ip
    real(Float64), dimension(3) :: ri      !! start position
    real(Float64), dimension(3) :: vi      !! velocity of fast ions
    real(Float64) :: denf !! fast-ion density
    integer, dimension(3) :: ind      !! new actual cell
    integer, dimension(4) :: neut_types=[1,2,3,4]
    logical :: los_intersect
    !! Determination of the CX probability
    type(LocalProfiles) :: plasma
    real(Float64), dimension(nlevs) :: prob !! Prob. for CX

    !! Collisiional radiative model along track
    integer :: ncell
    integer :: jj      !! counter along track
    type(ParticleTrack),dimension(beam_grid%ntrack) :: tracks

    real(Float64) :: photons !! photon flux
    real(Float64), dimension(nlevs) :: states  !! Density of n-states
    real(Float64), dimension(nlevs) :: denn

    !! Number of particles to launch
    integer(kind=8) :: pcnt
    real(Float64) :: papprox_tot, inv_maxcnt, cnt
    integer, dimension(3,beam_grid%ngrid) :: pcell
    real(Float64), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz) :: papprox, nlaunch !! approx. density

    !! Estimate how many particles to launch in each cell
    papprox=0.d0
    papprox_tot=0.d0
    pcnt=1
    do k=1,beam_grid%nz
        do j=1,beam_grid%ny
            do i=1,beam_grid%nx
                ind =[i,j,k]
                call get_plasma(plasma,ind=ind)
                papprox(i,j,k)=(sum(neut%dens(:,nbif_type,i,j,k)) + &
                                sum(neut%dens(:,nbih_type,i,j,k)) + &
                                sum(neut%dens(:,nbit_type,i,j,k)) + &
                                sum(neut%dens(:,halo_type,i,j,k)))* &
                                plasma%denf
                if(papprox(i,j,k).gt.0) then
                    pcell(:,pcnt)= ind
                    pcnt=pcnt+1
                endif
                if(plasma%in_plasma) papprox_tot=papprox_tot+papprox(i,j,k)
            enddo
        enddo
    enddo
    pcnt=pcnt-1
    inv_maxcnt=100.0/real(pcnt)
    call get_nlaunch(inputs%n_fida,papprox,papprox_tot,nlaunch)
    if(inputs%verbose.ge.1) then
        write(*,'(T6,"# of markers: ",i9)') int(sum(nlaunch))
    endif

    !! Loop over all cells that have neutrals
    cnt=0.d0
    loop_over_cells: do ip = 1, int(pcnt)
        i = pcell(1,ip)
        j = pcell(2,ip)
        k = pcell(3,ip)
        ind = [i, j, k]
        !$OMP PARALLEL DO schedule(guided) private(ip,iion,vi,ri, &
        !$OMP tracks,ncell,jj,plasma,prob,denn,states,photons,denf)
        loop_over_fast_ions: do iion=1,int8(nlaunch(i, j, k))
            !! Sample fast ion distribution for velocity and position
            call mc_fastion(ind, ri, vi, denf)
            if(denf.eq.0) cycle loop_over_fast_ions

            !! Find the particles path through the beam grid
            call track(ri, vi, tracks, ncell,los_intersect)
            if(.not.los_intersect) cycle loop_over_fast_ions
            if(ncell.eq.0) cycle loop_over_fast_ions

            !! Calculate CX probability with beam and halo neutrals
            call get_beam_cx_prob(tracks(1)%ind, ri, vi, neut_types, prob)
            if(sum(prob).le.0.) cycle loop_over_fast_ions

            !! Calculate initial states of particle
            states=prob*denf

            !! Calculate the spectra produced in each cell along the path
            loop_along_track: do jj=1,ncell
                call get_plasma(plasma,pos=tracks(jj)%pos)

                call colrad(plasma,beam_ion, vi, tracks(jj)%time, states, denn, photons)

                call store_fida_photons(tracks(jj)%pos, vi, photons/nlaunch(i,j,k))
            enddo loop_along_track
        enddo loop_over_fast_ions
        !$OMP END PARALLEL DO
        cnt=cnt+1

        if (inputs%verbose.eq.2)then
            WRITE(*,'(f7.2,"% completed",a,$)') cnt*inv_maxcnt,char(13)
        endif
    enddo loop_over_cells

end subroutine fida_f

subroutine fida_mc
    !+ Calculate FIDA emission using a Monte Carlo Fast-ion distribution
    integer :: iion
    type(FastIon) :: fast_ion
    type(LocalProfiles) :: plasma
    real(Float64) :: phi
    real(Float64), dimension(3) :: ri      !! start position
    real(Float64), dimension(3) :: vi      !! velocity of fast ions
    !! Determination of the CX probability
    real(Float64), dimension(nlevs) :: denn    !!  neutral dens (n=1-4)
    real(Float64), dimension(nlevs) :: prob    !! Prob. for CX
    !! Collisiional radiative model along track
    real(Float64), dimension(nlevs) :: states  ! Density of n-states
    integer :: ncell
    type(ParticleTrack), dimension(beam_grid%ntrack) :: tracks
    logical :: los_intersect
    integer :: jj      !! counter along track
    real(Float64) :: photons !! photon flux
    integer, dimension(4) :: neut_types=[1,2,3,4]
    real(Float64), dimension(3) :: uvw, uvw_vi
    real(Float64)  :: s, c
    real(Float64)  :: maxcnt, inv_maxcnt, cnt
    real(Float64), dimension(2) :: randomu

    maxcnt=particles%nparticle
    inv_maxcnt = 100.d0/maxcnt
    if(inputs%verbose.ge.1) then
        write(*,'(T6,"# of markers: ",i9)') particles%nparticle
    endif

    cnt=0.0
    !$OMP PARALLEL DO schedule(guided) private(iion,fast_ion,vi,ri,phi,tracks,s,c, &
    !$OMP& plasma,randomu,uvw,uvw_vi,ncell,jj,prob,denn,los_intersect,states,photons)
    loop_over_fast_ions: do iion=1,particles%nparticle
        fast_ion = particles%fast_ion(iion)
        cnt=cnt+1
        if(fast_ion%vabs.eq.0) cycle loop_over_fast_ions
        if(fast_ion%cross_grid) then
            !! Pick random torodial angle
            call randu(randomu)
            phi = fast_ion%phi_enter + fast_ion%delta_phi*randomu(1)
            s = sin(phi) ; c = cos(phi)

            !! Calculate position
            uvw(1) = fast_ion%r*c
            uvw(2) = fast_ion%r*s
            uvw(3) = fast_ion%z
            call uvw_to_xyz(uvw, ri)

            !! Calculate velocity vector
            uvw_vi(1) = c*fast_ion%vr - s*fast_ion%vt
            uvw_vi(2) = s*fast_ion%vr + c*fast_ion%vt
            uvw_vi(3) = fast_ion%vz
            vi = matmul(beam_grid%inv_basis,uvw_vi)

            !! Track particle through grid
            call track(ri, vi, tracks, ncell, los_intersect)
            if(.not.los_intersect) cycle loop_over_fast_ions
            if(ncell.eq.0)cycle loop_over_fast_ions

            !! Calculate CX probability
            call get_beam_cx_prob(tracks(1)%ind,ri,vi,neut_types,prob)
            if(sum(prob).le.0.)cycle loop_over_fast_ions

            !! Calculate the spectra produced in each cell along the path
            states=prob*fast_ion%weight
            loop_along_track: do jj=1,ncell
                call get_plasma(plasma,pos=tracks(jj)%pos)

                call colrad(plasma,beam_ion, vi, tracks(jj)%time, states, denn, photons)

                call store_fida_photons(tracks(jj)%pos, vi, photons, fast_ion%class)
            enddo loop_along_track
        endif
        if (inputs%verbose.ge.2)then
          WRITE(*,'(f7.2,"% completed",a,$)') cnt*inv_maxcnt,char(13)
        endif
    enddo loop_over_fast_ions
    !$OMP END PARALLEL DO

end subroutine fida_mc

subroutine npa_f
    !+ Calculate NPA flux using a fast-ion distribution function F(E,p,r,z)
    integer :: i,j,k  !! indices  x,y,z  of cells
    integer :: iion, det, ip
    real(Float64), dimension(3) :: ri      !! start position
    real(Float64), dimension(3) :: rf      !! end position
    real(Float64), dimension(3) :: vi      !! velocity of fast ions
    real(Float64) :: denf                  !! fast-ion density
    integer, dimension(3) :: ind      !! new actual cell
    integer, dimension(3,beam_grid%ngrid) :: pcell

    !! Determination of the CX probability
    type(LocalProfiles) :: plasma
    integer, dimension(4) :: neut_types=[1,2,3,4]
    real(Float64), dimension(nlevs) :: prob    !! Prob. for CX

    !! Collisiional radiative model along track
    real(Float64), dimension(nlevs) :: states  !! Density of n-states
    real(Float64) :: flux !! flux

    integer :: inpa,pcnt
    real(Float64) :: papprox_tot, maxcnt, cnt, inv_maxcnt
    real(Float64), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz) :: papprox, nlaunch !! approx. density

    papprox=0.d0
    papprox_tot=0.d0

    pcnt=1
    do k=1,beam_grid%nz
        do j=1,beam_grid%ny
            do i=1,beam_grid%nx
                ind =[i,j,k]
                call get_plasma(plasma,ind=ind)
                papprox(i,j,k)=(sum(neut%dens(:,nbif_type,i,j,k)) + &
                                sum(neut%dens(:,nbih_type,i,j,k)) + &
                                sum(neut%dens(:,nbit_type,i,j,k)) + &
                                sum(neut%dens(:,halo_type,i,j,k)))* &
                                plasma%denf

                if(papprox(i,j,k).gt.0.and.(npa_chords%hit(i,j,k))) then
                    !the only doing viewable cells is techically wrong
                    !since a guiding center not in a viewable cell
                    !can gyrostep into one but this is ridiculously faster
                    !and should be fine most of the time
                    pcell(:,pcnt)= ind
                    pcnt = pcnt + 1
                endif

                if(plasma%in_plasma) papprox_tot=papprox_tot+papprox(i,j,k)
            enddo
        enddo
    enddo
    pcnt = pcnt - 1
    maxcnt=real(pcnt)
    inv_maxcnt = 100.0/maxcnt

    call get_nlaunch(inputs%n_npa,papprox,papprox_tot,nlaunch)
    if(inputs%verbose.ge.1) then
        write(*,'(T6,"# of markers: ",i12)') int(sum(nlaunch))
    endif

    !! Loop over all cells that can contribute to NPA signal
    cnt=0.d0
    loop_over_cells: do ip = 1, int(pcnt)
        i = pcell(1,ip)
        j = pcell(2,ip)
        k = pcell(3,ip)
        ind = [i, j, k]
        !$OMP PARALLEL DO schedule(guided) private(iion, &
        !$OMP& vi,ri,rf,det,plasma,prob,states,flux,denf)
        loop_over_fast_ions: do iion=1,int(nlaunch(i, j, k)*npa%nloop)
            !! Sample fast ion distribution for velocity and position
            call mc_fastion(ind, ri, vi, denf)
            if(sum(vi).eq.0)cycle loop_over_fast_ions

            !! Check if particle hits a NPA detector
            call hit_npa_detector(ri, vi ,det, rf)
            if(det.eq.0) cycle loop_over_fast_ions

            !! Calculate CX probability with beam and halo neutrals
            call get_beam_cx_prob(ind,ri,vi,neut_types,prob)
            if(sum(prob).le.0.)cycle loop_over_fast_ions

            !! Attenuate states as the particle move through plasma
            states=prob*denf
            call attenuate(ri,rf,vi,states)

            !! Store NPA Flux
            flux = sum(states)*beam_grid%dv/(nlaunch(i,j,k)*npa%nloop)
            call store_npa(det,ri,rf,vi,flux)
        enddo loop_over_fast_ions
        !$OMP END PARALLEL DO
        cnt=cnt+1
        if (inputs%verbose.eq.2)then
            WRITE(*,'(f7.2,"% completed",a,$)') cnt*inv_maxcnt,char(13)
        endif
    enddo loop_over_cells
    write(*,'(T4,"Number of NPA particles that hit a detector: ",i8)') npa%npart

end subroutine npa_f

subroutine npa_mc
    !+ Calculate NPA flux using a Monte Carlo fast-ion distribution
    integer :: iion
    type(FastIon) :: fast_ion
    real(Float64) :: phi
    real(Float64), dimension(3) :: ri, rf      !! positions
    real(Float64), dimension(3) :: vi          !! velocity of fast ions
    integer :: det !! detector
    real(Float64), dimension(nlevs) :: prob    !! Prob. for CX
    real(Float64), dimension(nlevs) :: states  ! Density of n-states
    real(Float64) :: flux
    integer, dimension(4) :: neut_types=[1,2,3,4]
    integer, dimension(3) :: ind
    real(Float64), dimension(3) :: uvw, uvw_vi
    real(Float64) :: s,c
    real(Float64) :: maxcnt, inv_maxcnt, cnt
    real(Float64), dimension(2) :: randomu

    maxcnt=particles%nparticle
    inv_maxcnt = 100.d0/maxcnt

    cnt=0.0
    !$OMP PARALLEL DO schedule(guided) private(iion,ind,fast_ion,vi,ri,rf,phi,s,c, &
    !$OMP& randomu,uvw,uvw_vi,prob,states,flux,det)
    loop_over_fast_ions: do iion=1,particles%nparticle
        fast_ion = particles%fast_ion(iion)
        if(fast_ion%vabs.eq.0)cycle loop_over_fast_ions
        if(fast_ion%cross_grid) then
            !! Pick random torodial angle
            call randu(randomu)
            phi = fast_ion%phi_enter + fast_ion%delta_phi*randomu(1)
            s = sin(phi) ; c = cos(phi)

            !! Calculate position
            uvw(1) = fast_ion%r*c
            uvw(2) = fast_ion%r*s
            uvw(3) = fast_ion%z
            call uvw_to_xyz(uvw,ri)

            !! Calculate velocity
            uvw_vi(1) = c*fast_ion%vr - s*fast_ion%vt
            uvw_vi(2) = s*fast_ion%vr + c*fast_ion%vt
            uvw_vi(3) = fast_ion%vz
            vi = matmul(beam_grid%inv_basis,uvw_vi)

            !! Check if particle hits a NPA detector
            call hit_npa_detector(ri, vi ,det, rf)
            if(det.eq.0) cycle loop_over_fast_ions

            !! Get beam grid indices
            call get_indices(ri,ind)

            !! Calculate CX probability with beam and halo neutrals
            call get_beam_cx_prob(ind,ri,vi,neut_types,prob)
            if(sum(prob).le.0.)cycle loop_over_fast_ions

            !! Attenuate states as the particle moves though plasma
            states=prob*fast_ion%weight
            call attenuate(ri,rf,vi,states)

            !! Store NPA Flux
            flux = sum(states)*beam_grid%dv
            call store_npa(det,ri,rf,vi,flux)
        endif
        cnt=cnt+1
        if (inputs%verbose.ge.2)then
          WRITE(*,'(f7.2,"% completed",a,$)') cnt*inv_maxcnt,char(13)
        endif
    enddo loop_over_fast_ions
    !$OMP END PARALLEL DO

end subroutine npa_mc

subroutine fida_weights_mc
    !+ Calculates FIDA weights
    integer :: i,j,k !! indices  x,y,z of cells
    integer(kind=8) :: iion,ip
    real(Float64), dimension(3) :: ri,rg      !! start position
    real(Float64), dimension(3) :: vi     !! velocity of fast ions
    integer,dimension(3) :: ind      !! new actual cell
    integer,dimension(4) :: neut_types=[1,2,3,4]
    logical :: los_intersect

    !! Determination of the CX probability
    type(LocalProfiles) :: plasma
    type(LocalEMFields) :: fields
    real(Float64), dimension(nlevs) :: prob !! Prob. for CX

    !! Collisiional radiative model along track
    integer :: ncell
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
    integer(kind=8) :: pcnt
    real(Float64) :: papprox_tot,inv_maxcnt,cnt,fbm_denf,phase_area
    integer,dimension(3,beam_grid%ngrid) :: pcell
    real(Float64), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz) :: papprox,nlaunch !! approx. density

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
        write(*,'(T2,"Number of Channels: ",i3)') spec_chords%nchan
        write(*,'(T2,"Nlambda: ",i4)') nwav
        write(*,'(T2,"Nenergy: ",i3)') inputs%ne_wght
        write(*,'(T2,"Maximum Energy: ",f7.2)') inputs%emax_wght
        write(*,'(T2,"LOS averaged: ",a)') "False"
    endif

    !! zero out arrays
    fweight%weight = 0.d0
    fweight%mean_f = 0.d0

    etov2 = 1.d0/(v2_to_E_per_amu*inputs%ab)

    !! Estimate how many particles to launch in each cell
    papprox=0.d0
    papprox_tot=0.d0
    pcnt=1
    do k=1,beam_grid%nz
        do j=1,beam_grid%ny
            do i=1,beam_grid%nx
                ind =[i,j,k]
                call get_plasma(plasma,ind=ind)
                papprox(i,j,k)=(sum(neut%dens(:,nbif_type,i,j,k)) + &
                                sum(neut%dens(:,nbih_type,i,j,k)) + &
                                sum(neut%dens(:,nbit_type,i,j,k)) + &
                                sum(neut%dens(:,halo_type,i,j,k)))
                if(papprox(i,j,k).gt.0) then
                    pcell(:,pcnt)= ind
                    pcnt=pcnt+1
                endif
                if(plasma%in_plasma) papprox_tot=papprox_tot+papprox(i,j,k)
            enddo
        enddo
    enddo
    pcnt=pcnt-1
    inv_maxcnt=100.0/real(pcnt)
    call get_nlaunch(10*inputs%n_fida,papprox,papprox_tot,nlaunch)
    if(inputs%verbose.ge.1) then
        write(*,'(T6,"# of markers: ",i9)') int(sum(nlaunch))
    endif

    !! Loop over all cells that have neutrals
    cnt=0.d0
    loop_over_cells: do ip = 1, int(pcnt)
        i = pcell(1,ip)
        j = pcell(2,ip)
        k = pcell(3,ip)
        ind = [i, j, k]
        !$OMP PARALLEL DO schedule(guided) private(iion,vi,ri,rg,ienergy,ipitch, &
        !$OMP tracks,ncell,jj,plasma,fields,prob,denn,states,photons,energy,pitch, &
        !$OMP los_intersect,randomu3,fbm_denf)
        loop_over_fast_ions: do iion=1,int8(nlaunch(i, j, k))
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
            call gyro_correction(fields,rg,energy,pitch,ri,vi)

            fbm_denf = 0.0
            if (inputs%dist_type.eq.1) then
                call get_ep_denf(energy,pitch,fbm_denf,coeffs=fields%c)
            endif

            !! Find the particles path through the beam grid
            call track(ri, vi, tracks, ncell, los_intersect)
            if(.not.los_intersect) cycle loop_over_fast_ions
            if(ncell.eq.0) cycle loop_over_fast_ions

            !! Calculate CX probability with beam and halo neutrals
            call get_beam_cx_prob(tracks(1)%ind, ri, vi, neut_types, prob)
            if(sum(prob).le.0.) cycle loop_over_fast_ions
            states=prob*1.d20

            !! Calculate the spectra produced in each cell along the path
            loop_along_track: do jj=1,ncell
                call get_plasma(plasma,pos=tracks(jj)%pos)

                call colrad(plasma,beam_ion, vi, tracks(jj)%time, states, denn, photons)

                call store_fw_photons(ienergy(1), ipitch(1), &
                     tracks(jj)%pos, vi, fbm_denf, photons/nlaunch(i,j,k))
            enddo loop_along_track
        enddo loop_over_fast_ions
        !$OMP END PARALLEL DO
        cnt=cnt+1

        if(inputs%verbose.eq.2) then
            WRITE(*,'(f7.2,"% completed",a,$)') cnt*inv_maxcnt,char(13)
        endif
    enddo loop_over_cells

    fweight%weight = ((1.d-20)*phase_area/dEdP)*fweight%weight
    fweight%mean_f = ((1.d-20)*phase_area/dEdP)*fweight%mean_f
    call write_fida_weights()

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
    integer(Int32) :: i, j, k, ienergy
    integer(Int32) :: ipitch, igyro, icell, ichan
    real(Float64), dimension(:), allocatable :: ebarr,ptcharr,phiarr
    real(Float64), dimension(:,:), allocatable :: mean_f
    real(Float64), dimension(3) :: vi, vi_norm, vp
    real(Float64), dimension(3) :: vnbi_f, vnbi_h, vnbi_t, vhalo
    real(Float64), dimension(3) :: r_enter, r_exit
    real(Float64) :: vabs, dE, dP
    !! Determination of the CX probability
    real(Float64), dimension(nlevs) :: fdens,hdens,tdens,halodens
    real(Float64), dimension(nlevs) :: rates
    real(Float64), dimension(nlevs) :: states ! Density of n-states
    real(Float64), dimension(nlevs) :: denn  ! Density of n-states
    !! COLRAD
    real(Float64) :: dt, max_dens, dlength, sigma_pi
    real(Float64) :: eb, ptch, phi
    !! Solution of differential equation
    integer, dimension(3) :: ind  !!actual cell
    real(Float64), dimension(3) :: ri
    integer(Int32) :: ncell

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
        write(*,'(T2,"Number of Channels: ",i3)') spec_chords%nchan
        write(*,'(T2,"Nlambda: ",i4)') nwav
        write(*,'(T2,"Nenergy: ",i3)') inputs%ne_wght
        write(*,'(T2,"Npitch: ",i3)') inputs%np_wght
        write(*,'(T2,"Ngyro: ", i3)') inputs%nphi_wght
        write(*,'(T2,"Maximum Energy: ",f7.2)') inputs%emax_wght
        write(*,'(T2,"LOS averaged: ",a)') "True"
        write(*,*) ''
    endif

    etov2 = 1.0/(v2_to_E_per_amu*inputs%ab)

    chan_loop: do ichan=1,spec_chords%nchan
        fdens = 0.d0 ; hdens = 0.d0 ; tdens = 0.d0 ; halodens = 0.d0
        plasma = plasma*0.d0
        fields = fields*0.d0
        wght_tot = 0.d0
        mean_f = 0.d0
        do k=1,beam_grid%nz
            do j=1,beam_grid%ny
                do i=1,beam_grid%nx
                    if(spec_chords%dlength(ichan,i,j,k).gt.0.0) then
                        ind = [i,j,k]
                        dlength = spec_chords%dlength(ichan,i,j,k)
                        fdens = fdens + neut%dens(:,nbif_type,i,j,k)*dlength
                        hdens = hdens + neut%dens(:,nbih_type,i,j,k)*dlength
                        tdens = tdens + neut%dens(:,nbit_type,i,j,k)*dlength
                        halodens = halodens + neut%dens(:,halo_type,i,j,k)*dlength
                        wght = sum(neut%dens(3,1:4,i,j,k))*dlength

                        call get_plasma(plasma_cell,ind=ind)
                        call get_fields(fields_cell,ind=ind)
                        plasma = plasma + wght*plasma_cell
                        fields = fields + wght*fields_cell
                        if (inputs%dist_type.eq.1) then
                            do ipitch=1,inputs%np_wght
                                do ienergy=1,inputs%ne_wght
                                    call get_ep_denf(ebarr(ienergy),ptcharr(ipitch),denf,coeffs=fields_cell%c)
                                    mean_f(ienergy,ipitch) = mean_f(ienergy,ipitch) + wght*denf
                                enddo
                            enddo
                        endif
                        wght_tot = wght_tot + wght
                    endif
                enddo
            enddo
        enddo

        if(wght_tot.le.0) then
            if(inputs%verbose.ge.1) then
                write(*,'(T4,"Skipping channel ",i3,": Neutral density is zero")') ichan
            endif
            cycle chan_loop
        else
            plasma = plasma/wght_tot
            plasma%in_plasma = .True.
            fields = fields/wght_tot
            fields%in_plasma= .True.
            mean_f = mean_f/wght_tot
            if(inputs%verbose.ge.1) then
                write(*,'(T4,"Channel: ",i3)') ichan
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
        !$OMP& r_enter,r_exit,length,max_dens,ind,tracks,ncell,dt,icell,states,rates, &
        !$OMP& vhalo,denn,denf,photons,ienergy,ipitch,igyro)
        do ienergy=1,inputs%ne_wght
            do ipitch=1,inputs%np_wght
                do igyro=1,inputs%nphi_wght
                    eb = ebarr(ienergy)
                    vabs = sqrt(eb*etov2)
                    ptch = ptcharr(ipitch)
                    phi = phiarr(igyro)
                    call pitch_to_vec(ptch,phi,fields,vi_norm)
                    vi = vabs*vi_norm

                    call grid_intersect(ri,vi,length,r_enter,r_exit)
                    call track(r_enter, vi, tracks, ncell)
                    max_dens = 0.d0
                    do icell=1,ncell
                        ind = tracks(icell)%ind
                        tracks(icell)%flux = sum(neut%dens(3,1:4,ind(1),ind(2),ind(3)))
                        if(tracks(icell)%flux.gt.max_dens) max_dens=tracks(icell)%flux
                    enddo
                    dt = 0.d0
                    do icell=1,ncell
                        if(tracks(icell)%flux.gt.(0.5*max_dens)) then
                            dt = dt + tracks(icell)%time
                        endif
                    enddo

                    states=0.d0
                    call neut_rates(fdens,vi,vnbi_f,rates)
                    states = states + rates
                    call neut_rates(hdens,vi,vnbi_h,rates)
                    states = states + rates
                    call neut_rates(tdens,vi,vnbi_t,rates)
                    states = states + rates
                    do i=1,int(n_halo_neutrate)
                        call mc_halo(ind,vhalo,plasma_in=plasma)
                        call neut_rates(halodens,vi,vhalo,rates)
                        states = states + rates/real(n_halo_neutrate)
                    enddo

                    call colrad(plasma,beam_ion,vi,dt,states,denn,photons)
                    denf = mean_f(ienergy,ipitch)*dEdP
                    photons = photons/real(inputs%nphi_wght)
                    call store_fw_photons_at_chan(ichan, ienergy, ipitch, &
                         vp, vi, fields, dlength, sigma_pi, denf, photons)

                enddo
            enddo
        enddo
        !$OMP END PARALLEL DO
    enddo chan_loop

    fweight%mean_f = fweight%mean_f/(dEdP)
    call write_fida_weights()

end subroutine fida_weights_los

subroutine npa_weights
    !+ Calculates NPA weights
    type(LocalEMFields) :: fields
    real(Float64) :: pitch
    real(Float64) :: pcxa
    integer(Int32) :: det
    integer(Int32) :: ii, jj, kk, i, ic   !!indices
    integer,dimension(1) :: ipitch
    real(Float64), dimension(3) :: vi,vi_norm
    real(Float64) :: vabs, fbm_denf, dE, dP, ccnt
    real(Float64), dimension(nlevs) :: pcx   !! Rate coefficiants for CX
    real(Float64), dimension(nlevs) :: states, states_i  ! Density of n-states
    integer, dimension(4) :: neut_types=[1,2,3,4]
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
        write(*,'(T2,"Number of Channels: ",i3)') npa_chords%nchan
        write(*,'(T2,"Nenergy: ",i3)') inputs%ne_wght
        write(*,'(T2,"Npitch: ",i3)') inputs%np_wght
        write(*,'(T2,"Maximum energy: ",f7.2)') inputs%emax_wght
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
        if(inputs%verbose.ge.1) then
            write(*,'(T4,"Channel: ",i3)') ichan
            write(*,'(T4,"Radius: ",f10.3)') npa_chords%radius(ichan)
        endif

        ccnt=0.d0
        !$OMP PARALLEL DO schedule(guided) collapse(3) private(ii,jj,kk,fields, &
        !$OMP& ic,det,pos,dpos,r_gyro,pitch,ipitch,vabs,vi,pcx,pcxa,states,states_i,vi_norm,fbm_denf)
        loop_along_z: do kk=1,beam_grid%nz
            loop_along_y: do jj=1,beam_grid%ny
                loop_along_x: do ii=1,beam_grid%nx
                    if(npa_chords%phit(ii,jj,kk,ichan)%p.gt.0.d0) then
                        pos = [beam_grid%xc(ii), beam_grid%yc(jj), beam_grid%zc(kk)]
                        call get_fields(fields,pos=pos)
                        if(.not.fields%in_plasma) cycle loop_along_x

                        !!Determine velocity vector
                        dpos = npa_chords%phit(ii,jj,kk,ichan)%eff_rd
                        vi_norm = (dpos - pos)/norm2(dpos - pos)

                        !!Check if it hits a detector just to make sure
                        call hit_npa_detector(pos,vi_norm,det)
                        if (det.ne.ichan) then
                            write(*,'(a)') 'NPA_WEIGHTS: Missed detector'
                            cycle loop_along_x
                        endif

                        !! Determine the angle between the B-field and the Line of Sight
                        pitch = dot_product(fields%b_norm,vi_norm)
                        ipitch=minloc(abs(ptcharr - pitch))
                        loop_over_energy: do ic = 1, inputs%ne_wght !! energy loop
                            vabs = sqrt(ebarr(ic)/(v2_to_E_per_amu*inputs%ab))
                            vi = vi_norm*vabs
                            !!Correct for gyro orbit
                            call gyro_step(vi,fields,r_gyro)

                            fbm_denf=0
                            if (inputs%dist_type.eq.1) then
                                !get dist at guiding center
                                call get_ep_denf(ebarr(ic),pitch,fbm_denf,pos=(pos+r_gyro))
                            endif
                            if (fbm_denf.ne.fbm_denf) cycle loop_over_energy

                            !! -------------- calculate CX probability -------!!
                            call get_beam_cx_prob([ii,jj,kk],pos,vi,neut_types,pcx)
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
                                      2*sum(pcx)*pcxa*npa_chords%phit(ii,jj,kk,ichan)%p*beam_grid%dv/dP

                            nweight%flux(ic,ichan) = nweight%flux(ic,ichan) + &
                                      2*beam_grid%dv*fbm_denf*sum(pcx)*pcxa*npa_chords%phit(ii,jj,kk,ichan)%p
                                      !Factor of 2 above is to convert fbm to ions/(cm^3 dE (domega/4pi))
                            nweight%emissivity(ii,jj,kk,ichan)=nweight%emissivity(ii,jj,kk,ichan)+ &
                                      2*fbm_denf*sum(pcx)*pcxa*npa_chords%phit(ii,jj,kk,ichan)%p*dE
                            !$OMP END CRITICAL(npa_wght)
                        enddo loop_over_energy
                    endif
                    ccnt=ccnt+1
                    if (inputs%verbose.eq.2)then
                        WRITE(*,'(f7.2,"% completed",a,$)') ccnt/real(beam_grid%ngrid)*100,char(13)
                    endif
                enddo loop_along_x
            enddo loop_along_y
        enddo loop_along_z
        !$OMP END PARALLEL DO

       if(inputs%verbose.ge.1) then
           write(*,'(T4,A,ES14.5)'),'Flux:   ',sum(nweight%flux(:,ichan))*dE
           write(*,'(T4,A,ES14.5)'),'Weight: ',sum(nweight%weight(:,:,ichan))*dE*dP
           write(*,*) ''
       endif
    enddo loop_over_channels

    call write_npa_weights()

end subroutine npa_weights

end module libfida

!=============================================================================
!-------------------------------Main Program----------------------------------
!=============================================================================
program fidasim
    !+ FIDASIM {!../VERSION!}
    use libfida
    use hdf5_extra
#ifdef _OMP
    use omp_lib
#endif
    implicit none
    character(3)          :: arg = ''
    integer, dimension(8) :: time_arr,time_start,time_end !Time array
    integer               :: i,narg,nthreads,max_threads
    integer               :: hour,minu,sec

#ifdef _VERSION
    version = _VERSION
#endif

    call print_banner()

    narg = command_argument_count()
    if(narg.eq.0) then
        write(*,'(a)') "usage: ./fidasim namelist_file [num_threads]"
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

    !! ----------------------------------------------------------
    !! ------ INITIALIZE THE RANDOM NUMBER GENERATOR  -----------
    !! ----------------------------------------------------------
    allocate(rng(max_threads))
    do i=1,max_threads
        call rng_init(rng(i),932117 + i)
    enddo

    !! ----------------------------------------------------------
    !! ------- READ GRIDS, PROFILES, LOS, TABLES, & FBM --------
    !! ----------------------------------------------------------
    call make_beam_grid()
    call read_equilibrium()
    call read_beam()
    call read_tables()
    call read_distribution()

    allocate(spec_chords%los_inter(beam_grid%nx,beam_grid%ny,beam_grid%nz))
    spec_chords%los_inter = .False.
    if((inputs%calc_spec.ge.1).or.(inputs%calc_fida_wght.ge.1)) then
        call read_chords()
    endif

    if((inputs%calc_npa.ge.1).or.(inputs%calc_npa_wght.ge.1)) then
        call read_npa()
    endif

    !! ----------------------------------------------------------
    !! --------------- ALLOCATE THE RESULT ARRAYS ---------------
    !! ----------------------------------------------------------
    !! neutral density array!
    allocate(neut%dens(nlevs,ntypes,beam_grid%nx,beam_grid%ny,beam_grid%nz))
    neut%dens = 0.d0

    !! birth profile
    if(inputs%calc_birth.ge.1) then
        allocate(birth%dens(3, &
                            beam_grid%nx, &
                            beam_grid%ny, &
                            beam_grid%nz))
        allocate(birth%neut_type(int(inputs%n_birth*inputs%n_nbi)))
        allocate(birth%ind(3,int(inputs%n_birth*inputs%n_nbi)))
        allocate(birth%ri(3,int(inputs%n_birth*inputs%n_nbi)))
        allocate(birth%vi(3,int(inputs%n_birth*inputs%n_nbi)))
        birth%neut_type = 0
        birth%dens = 0.d0
        birth%ind = 0
        birth%ri = 0.d0
        birth%vi = 0.d0
    endif

    if(inputs%calc_spec.ge.1) then
        allocate(spec%brems(inputs%nlambda,spec_chords%nchan))
        allocate(spec%bes(inputs%nlambda,spec_chords%nchan,4))
        allocate(spec%fida(inputs%nlambda,spec_chords%nchan,particles%nclass))
        spec%brems = 0.d0
        spec%bes = 0.d0
        spec%fida = 0.d0
    endif

    if(inputs%calc_npa.ge.1)then
        npa%nchan = npa_chords%nchan
        allocate(npa%part(npa%nmax))
        allocate(npa%energy(fbm%nenergy))
        allocate(npa%flux(fbm%nenergy,npa%nchan))
        npa%energy = fbm%energy
        npa%flux = 0.0
    endif

    !! -----------------------------------------------------------------------
    !! --------------- CALCULATE/LOAD the BEAM and HALO DENSITY---------------
    !! -----------------------------------------------------------------------
    if(inputs%load_neutrals.eq.1) then
        call read_neutrals()
    else
        !! ----------- BEAM NEUTRALS ---------- !!
        call date_and_time (values=time_arr)
        if(inputs%verbose.ge.1) then
            write(*,'(A,I2,":",I2.2,":",I2.2)') 'ndmc:    ' , &
                  time_arr(5),time_arr(6),time_arr(7)
        endif
        call ndmc
        if(inputs%calc_birth.eq.1)then
            call write_birth_profile()
        endif
        write(*,'(30X,a)') ''

        !! ---------- DCX (Direct charge exchange) ---------- !!
        call date_and_time (values=time_arr)
        if(inputs%verbose.ge.1) then
            write(*,'(A,I2,":",I2.2,":",I2.2)') 'dcx:     ' , &
                  time_arr(5),time_arr(6),time_arr(7)
        endif
        call dcx()
        if(inputs%dump_dcx.eq.1) call write_dcx()
        write(*,'(30X,a)') ''

        !! ---------- HALO ---------- !!
        call date_and_time (values=time_arr)
        if(inputs%verbose.ge.1) then
            write(*,'(A,I2,":",I2.2,":",I2.2)') 'halo:    ' , &
                  time_arr(5),time_arr(6),time_arr(7)
        endif
        call halo()
        !! ---------- WRITE NEUTRALS ---------- !!
        call write_neutrals()
        write(*,'(30X,a)') ''
    endif

    !! -----------------------------------------------------------------------
    !!----------------------------- BREMSSTRAHLUNG ---------------------------
    !! -----------------------------------------------------------------------
    if(inputs%calc_brems.ge.1) then
        call date_and_time (values=time_arr)
        if(inputs%verbose.ge.1) then
            write(*,'(A,I2,":",I2.2,":",I2.2)') 'bremsstrahlung:    ' , &
                  time_arr(5),time_arr(6),time_arr(7)
        endif
        call bremsstrahlung()
        write(*,'(30X,a)') ''
    endif

    !! -----------------------------------------------------------------------
    !! --------------------- CALCULATE the FIDA RADIATION --------------------
    !! -----------------------------------------------------------------------
    if(inputs%calc_fida.ge.1)then
        call date_and_time (values=time_arr)
        if(inputs%verbose.ge.1) then
            write(*,'(A,I2,":",I2.2,":",I2.2)') 'fida:    ' , &
                  time_arr(5),time_arr(6),time_arr(7)
        endif
        if(inputs%dist_type.eq.1) then
            call fida_f()
        else
            call fida_mc()
        endif
        write(*,'(30X,a)') ''
    endif

    if(inputs%calc_spec.ge.1) then
        call write_spectra()
        write(*,'(30X,a)') ''
    endif

    !! -----------------------------------------------------------------------
    !! ----------------------- CALCULATE the NPA FLUX ------------------------
    !! -----------------------------------------------------------------------
    if(inputs%calc_npa.ge.1)then
        call date_and_time (values=time_arr)
        if(inputs%verbose.ge.1) then
            write(*,'(A,I2,":",I2.2,":",I2.2)') 'npa:    ' , &
            time_arr(5),time_arr(6),time_arr(7)
        endif
        if(inputs%dist_type.eq.1) then
            call npa_f()
        else
            call npa_mc()
        endif
        write(*,'(30X,a)') ''
    endif

    if(inputs%calc_npa.ge.1) then
        call write_npa()
        write(*,'(30X,a)') ''
    endif

    !! -------------------------------------------------------------------
    !! ----------- Calculation of weight functions -----------------------
    !! -------------------------------------------------------------------
    if(inputs%calc_fida_wght.ge.1) then
        colrad_threshold=0. !! to speed up simulation!
        call date_and_time (values=time_arr)
        if(inputs%verbose.ge.1) then
            write(*,'(A,I2,":",I2.2,":",I2.2)') 'fida weight function:    ',  &
                  time_arr(5),time_arr(6),time_arr(7)
        endif
        if(inputs%calc_fida_wght.eq.1) then
            call fida_weights_los()
        else
            call fida_weights_mc()
        endif
        write(*,'(30X,a)') ''
    endif

    if(inputs%calc_npa_wght.ge.1) then
        call date_and_time (values=time_arr)
        if(inputs%verbose.ge.1) then
            write(*,'(A,I2,":",I2.2,":",I2.2)') 'npa weight function:    ',  &
                  time_arr(5),time_arr(6),time_arr(7)
        endif
        call npa_weights()
        write(*,'(30X,a)') ''
    endif

    call date_and_time (values=time_arr)
    if(inputs%verbose.ge.1) then
        write(*,'(A,I2,":",I2.2,":",I2.2)') 'END: hour, minute, second: ',&
              time_arr(5),time_arr(6),time_arr(7)
    endif

    call date_and_time (values=time_end)
    hour = time_end(5) - time_start(5)
    minu = time_end(6) - time_start(6)
    sec  = time_end(7) - time_start(7)
    if (minu.lt.0.) then
        minu = minu +60
        hour = hour -1
    endif
    if (sec.lt.0.) then
        sec  = sec +60
        minu = minu -1
    endif

    if(inputs%verbose.ge.1) then
        write(*,'(A,18X,I2,":",I2.2,":",I2.2)') 'duration:',hour,minu,sec
    endif

end program fidasim
