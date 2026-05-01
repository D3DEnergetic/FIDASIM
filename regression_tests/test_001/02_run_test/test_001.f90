module test_001_utils
  use iso_fortran_env, only: Int32 => int32, Int64 => int64, Float32 => real32, Float64 => real64
  implicit none
  ! private
  ! public :: read_test_input

  type :: InputData
    ! Arrays (variable length)
    real(Float64), allocatable :: dene(:)  ! [cm^-3]
    real(Float64), allocatable :: Lx(:)    ! [cm]

    ! Scalars
    integer(Int32) :: nden  = -1
    integer(Int32) :: nx    = -1
    integer(Int32) :: nt    = -1
    integer(Int32) :: neb   = -1

    real(Float64)  :: tmin  = -1.0d0 ! [eV]
    real(Float64)  :: tmax  = -1.0d0 ! [eV]
    real(Float64)  :: emin  = -1.0d0 ! [eV]
    real(Float64)  :: emax  = -1.0d0 ! [eV]
    real(Float64)  :: beam_mass = -1.0d0 ! [amu]
    real(Float64)  :: ion_mass = -1.0d0 ! [amu]

  end type InputData

  ! Defintion of structures:
  type(InputData) :: test_inputs

  contains

    subroutine read_test_input(filename, out)
      implicit none
      character(len=*), intent(in)    :: filename
      type(InputData),  intent(inout) :: out

      integer, parameter :: nmax = 64

      ! --- variables that NAMELIST can read into (must be plain variables) ---
      integer(Int32) :: nden
      real(Float64)  :: dene(nmax), Lx(nmax)
      integer(Int32) :: nx, nt, neb
      real(Float64)  :: tmin, tmax, emin, emax, beam_mass, ion_mass
      ! ----------------------------------------------------------------------

      integer :: iu, ios
      character(len=256) :: iomsg

      ! Names MUST match the file (your file uses: nden, dene, Lx, nx, ...)
      namelist /run_input/ nden, dene, Lx, nx, &
                           tmin, tmax, nt, &
                           emin, emax, neb, beam_mass, ion_mass

      ! Defaults (important)
      nden  = -1
      dene  = 0.0d0
      Lx    = 0.0d0
      nx    = -1
      tmin  = -1.0d0
      tmax  = -1.0d0
      nt    = -1
      emin  = -1.0d0
      emax  = -1.0d0
      neb   = -1
      beam_mass = -1.0d0
      ion_mass = -1.0d0

      open(newunit=iu, file=trim(filename), status='old', action='read', iostat=ios, iomsg=iomsg)
      if (ios /= 0) error stop 'OPEN failed: ' // trim(iomsg)

      read(iu, nml=run_input, iostat=ios, iomsg=iomsg)
      if (ios /= 0) error stop 'NAMELIST read failed: ' // trim(iomsg)
      close(iu)

      ! Validation
      if (nden <= 0) error stop 'nden must be > 0'
      if (nden > nmax) error stop 'nden exceeds buffer size (increase nmax)'
      if (nx <= 0) error stop 'nx must be > 0'
      if (nt <= 0) error stop 'nt must be > 0'
      if (neb <= 0) error stop 'neb must be > 0'
      if (beam_mass <= 0.0d0) error stop 'beam_mass must be > 0'
      if (ion_mass <= 0.0d0) error stop 'ion_mass must be > 0'
      if (tmin <= 0.0d0 .or. tmax <= tmin) error stop 'invalid temperature range'
      if (emin <= 0.0d0 .or. emax <= emin) error stop 'invalid energy range'

      ! Store scalars
      out%nden  = nden
      out%nx    = nx
      out%nt = nt
      out%neb   = neb
      out%tmin  = tmin
      out%tmax  = tmax
      out%emin  = emin
      out%emax  = emax
      out%beam_mass = beam_mass
      out%ion_mass = ion_mass

      ! Replace arrays safely
      if (allocated(out%dene)) deallocate(out%dene)
      if (allocated(out%Lx))   deallocate(out%Lx)

      allocate(out%dene(nden), out%Lx(nden))
      out%dene = dene(1:nden)
      out%Lx   = Lx(1:nden)

      if (any(out%dene <= 0.0d0)) error stop 'All dene values must be > 0'
      if (any(out%Lx   <= 0.0d0)) error stop 'All Lx values must be > 0'

    end subroutine read_test_input

    subroutine get_atomic_tables_path(path)
      implicit none

      character(len=*), intent(out) :: path
      character(len=4096) :: buf
      integer(Int32) :: stat, n

      call get_environment_variable("FIDASIM_DIR", value=buf, length=n, status=stat)

      if (stat /= 0 .or. n <= 0) then
        error stop "FIDASIM_DIR environment variable is not set."
      end if

      if (len(path) < n + len("/tables/atomic_tables.h5")) then
        error stop "Provided path buffer is too small"
      end if

      path = trim(buf(1:n)) // "/tables/atomic_tables.h5"
    end subroutine get_atomic_tables_path

    subroutine linspace(a, b, n, x)
      implicit none

      real(Float64), intent(in)  :: a, b
      integer,      intent(in)  :: n
      real(Float64), allocatable, intent(out) :: x(:)

      integer :: i

      if (n <= 0) error stop "linspace: n must be > 0"

      allocate(x(n))

      if (n == 1) then
        x(1) = a
      else
        do i = 1, n
          x(i) = a + (b - a) * real(i-1, Float64) / real(n-1, Float64)
        end do
      end if

    end subroutine linspace

end module test_001_utils

! =============================================================================
! =============================================================================
! MAIN PROGRAM
! =============================================================================
! =============================================================================

program test_001
  use hdf5_utils
  use test_001_utils
  use libfida
  implicit none

  real(Float64), allocatable :: dene(:), eb(:), Te(:), Ti(:)
  real(Float64), allocatable :: denn(:,:,:,:) ! (nx,nden,neb,nt)
  real(Float64), allocatable :: denn_n(:,:,:,:,:) ! (nlevs,nx,nden,neb,nt)
  integer(Int32) :: nden, nx, nt, neb
  character(len=512) :: atomic_tables_file, h5filename
  logical :: ex
  real(Float64), allocatable :: x(:,:)
  integer :: xx, nn, ee, tt
  real(FLoat64), allocatable :: dx(:)
  real(Float64), dimension(nlevs) :: states, density
  real(Float64) :: vabs, dt, photons
  real(Float64 ), dimension(3) :: vn
  type(LocalProfiles) :: plasma
  integer :: error
  integer(HID_T) :: fid
  integer(HSIZE_T), dimension(5) :: dim5
  integer(HSIZE_T), dimension(4) :: dim4
  integer(HSIZE_T), dimension(2) :: dim2
  integer(HSIZE_T), dimension(1) :: dim1

  ! Read test input namelist;
  call read_test_input('../01_reference/input.nml', test_inputs)

  ! Get path to default FIDASIM table:
  call get_atomic_tables_path(atomic_tables_file)
  inquire(file=trim(atomic_tables_file), exist=ex)
  if (.not. ex) error stop "atomic_tables.h5 not found: "//trim(atomic_tables_file)
  inputs%tables_file = atomic_tables_file

  ! Read in atomic tables:
  inputs%non_thermal_beam_stopping = 0 ! from libfida
  inputs%full_f = 0 ! from libfida
  inputs%calc_neutron = 0 ! from libfida
  inputs%calc_cfpd = 0 ! from libfida
  impurity_charge = 6 ! from libfida, assuming also denimp = 0
  call read_tables()
  nx    = test_inputs%nx
  nden  = test_inputs%nden
  nt = test_inputs%nt
  neb   = test_inputs%neb
  beam_mass = test_inputs%beam_mass ! from libfida
  thermal_mass(1) = test_inputs%ion_mass ! from libfida
  initial_state = 3 ! from libfida
  final_state = 2 ! from libfida

  ! Make grids:
  if (allocated(x)) deallocate(x)
  allocate(x(nx,nden))
  allocate(dx(nden))
  do nn = 1, nden
    dx(nn) = test_inputs%Lx(nn) / real(nx - 1, Float64)
    do xx = 1, nx
      x(xx,nn) = real(xx-1, Float64) * dx(nn) ! [cm]
    end do
  end do

  allocate(dene(nden))
  dene = test_inputs%dene

  allocate(eb(neb))
  call linspace(test_inputs%emin,test_inputs%emax,neb,eb) ! [eV]

  allocate(Te(nt),Ti(nt))
  call linspace(test_inputs%tmin,test_inputs%tmax,nt,Ti) ! [eV]
  Ti = Ti*1e-3 ! [keV]
  Te = Ti ! [keV]

  write(*,*) "test_inputs.dene: ", test_inputs%dene
  write(*,*) "inputs.tables_file: ", inputs%tables_file
  write(*,*) "dx: ", dx
  write(*,*) "eb: ", eb
  write(*,*) "Ti: ", Ti

  ! Allocate and init main output variables:
  allocate(denn(nx,nden,neb,nt))
  allocate(denn_n(nlevs,nx,nden,neb,nt))
  denn = 0.0d0
  denn_n = 0.0d0

  ! Initialize variables:
  plasma%in_plasma = .true.
  vn = 0.0d0

  ! Compute attenuation:
  do nn = 1, nden
    do ee = 1, neb
      do tt = 1, nt
        ! Reset beam states:
        states = 0.0d0
        density = 0.0d0
        photons = 0.0d0
        states(1) = 1

        ! Define plasma density:
        plasma%dene = dene(nn)
        plasma%deni(1) = dene(nn)

        ! Neutral beam velocity:
        vabs = sqrt(2.d0*eb(ee)*e0/(beam_mass*mass_u))*1.d2 ! [cm/s]
        vn(1) = vabs

        ! Time in cell:
        dt = dx(nn)/vabs

        ! Define plasma temperature: [keV]
        plasma%te = Te(tt)
        plasma%ti = Ti(tt)

        ! Attenuate along path:
        do xx = 1, nx
          call colrad(plasma,beam_mass,vn,dt,states,density,photons)
          denn_n(:,xx,nn,ee,tt) = density
          denn(xx,nn,ee,tt) = sum(density)
        end do
      end do
    end do
  end do

  ! Save data to file:
  ! TODO:
  ! - Add dimension integers

  h5filename = "test_001.h5"
  call h5open_f(error)
  call h5fcreate_f(h5filename, H5F_ACC_TRUNC_F, fid, error)

  dim4 = shape(denn)
  call h5ltmake_compressed_dataset_double_f(fid,"/denn",4,dim4,denn,error)
  call h5ltset_attribute_string_f(fid,"/denn","description", &
    "total neutral density profile along 'x' for each (dene,eb,T). The initial density is given as 1 cm^-3",error)
  call h5ltset_attribute_string_f(fid,"/denn","dimensions", "[nx, nden, neb, nt]",error)
  call h5ltset_attribute_string_f(fid,"/denn","units", "cm^-3",error)

  dim2 = [nx, nden]
  call h5ltmake_compressed_dataset_double_f(fid,"/x",2,dim2,x,error)
  call h5ltset_attribute_string_f(fid,"/x","description", &
    "dimension along the path of ray for each density value",error)
  call h5ltset_attribute_string_f(fid,"/x","dimensions", "[nx, nden]",error)
  call h5ltset_attribute_string_f(fid,"/x","units","cm",error)

  dim5 = shape(denn_n)
  call h5ltmake_compressed_dataset_double_f(fid,"/denn_n",5,dim5,denn_n, error)
  call h5ltset_attribute_string_f(fid,"/denn_n","description", &
    "level-resolved neutral density profle along 'x' for different conditions",error)
  call h5ltset_attribute_string_f(fid,"/denn_n","dimensions", "[nlev,nx, nden, neb, nt]",error)
  call h5ltset_attribute_string_f(fid,"/denn_n","units","cm^-3",error)

  dim1 = [nden]
  call h5ltmake_compressed_dataset_double_f(fid,"/dene",1,dim1,dene,error)
  call h5ltset_attribute_string_f(fid,"/dene","description","electron density grid",error)
  call h5ltset_attribute_string_f(fid,"/dene","dimensions", "[nden]",error)
  call h5ltset_attribute_string_f(fid,"/dene","units","cm^-3",error)

  dim1 = [neb]
  call h5ltmake_compressed_dataset_double_f(fid,"/eb",1,dim1,eb,error)
  call h5ltset_attribute_string_f(fid,"/eb","description","beam energy grid",error)
  call h5ltset_attribute_string_f(fid,"/eb","dimensions", "[neb]",error)
  call h5ltset_attribute_string_f(fid,"/eb","units","eV/amu",error)

  dim1 = [nt]
  call h5ltmake_compressed_dataset_double_f(fid,"/temperature",1,dim1,Ti*1e3, error)
  call h5ltset_attribute_string_f(fid,"/temperature","description","plasma temperature grid",error)
  call h5ltset_attribute_string_f(fid,"/temperature","dimensions", "[nt]",error)
  call h5ltset_attribute_string_f(fid,"/temperature","units","eV",error)

  dim1 = [1]
  call h5ltmake_compressed_dataset_double_f(fid,"/beam_mass",1,dim1,[beam_mass], error)
  call h5ltset_attribute_string_f(fid,"/beam_mass","description","Beam atomic mass",error)
  call h5ltset_attribute_string_f(fid,"/beam_mass","dimensions", "[1]",error)
  call h5ltset_attribute_string_f(fid,"/beam_mass","units","amu",error)

  call h5fclose_f(fid, error)
  call h5close_f(error)

end program test_001
