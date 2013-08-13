!FIDASIM2.0 
!Original version from W. W. Heidbrink 2007
!rewritten in F90 by Benedikt Geiger 2012 (IPP-Garching, ASDEX Upgrade)
!The main routine (fidasim) is at the end of the file!
module application
  implicit none
  !!                      Indizes of the different components:
  character(120)        :: result_dir
  character(120)        :: root_dir
  integer , parameter   :: fida_type = 1 ! fida spectra/density
  integer , parameter   :: nbif_type = 2 ! full energy NBI spectra/density
  integer , parameter   :: nbih_type = 3 ! half energy NBI spectra/density
  integer , parameter   :: nbit_type = 4 ! third energy NBI spectra/density
  integer , parameter   :: halo_type = 5 ! halo spectra/density
  integer , parameter   :: s1type    = 6 ! halo iteration
  integer , parameter   :: s2type    = 7 ! halo iteration
  integer , parameter   :: brems_type= 8 ! halo iteration  
  !!                      Definition for the kind of the variables: 
  integer , parameter   :: long      = kind(1)
  integer , parameter   :: float     = kind(1.e0)
  integer , parameter   :: double    = kind(1.d0) 
  real(double),parameter :: ONE=1.d0,TWO=2.d0,ZERO=0.d0
  real(double),parameter :: XMACH_EPS=2.22d-16
  integer , parameter   :: MAXIT=50
  !!                      Physical units:
  real(double),parameter:: mass_u    = 1.6605402d-27  ! [kg]
  real(double),parameter:: e0        = 1.60217733d-19 ! [C]
  real(double),parameter:: pi        = 3.14159265358979323846264d0
  real(double),parameter:: c0        = 2.99792458d+08 !! [m/s]
  real(double),parameter:: h_planck  = 4.135667516d-15 !![eV/s]
  real(double),parameter:: lambda0   = 6561.d0        !!D-alpha [A]
  real(double),parameter:: v_to_E    = mass_u/(2.*e0*1.e3)*1.e-4 !!conversion cm^2/s^2 to keV
  !! ---- Stark splitting, wavelength and intenisty from Schroedinger ---- !!
  real(double),parameter,dimension(9)::stark_wavel =(/-257.56d0,-193.17d0 &
       ,-128.78d0,-64.39d0,0.d0,64.39d0,128.78d0,193.17d0,257.56d0/) &
       *1.d-16*lambda0**2 
  real(double),parameter,dimension(9)::stark_intens=(/ 1681.d0, 2304.d0 &
       , 729.d0, 1936.d0, 5490.d0, 1936.d0, 729.d0, 2304.d0, 1681.d0/)
  integer(long), dimension(1) :: minpos  !! dummy array to determine minloc
  !!Numerical Settings
  integer(long),parameter:: nlevt=6       !!nr of tabulated states 
  integer(long),parameter:: nlevs=6        !!nr of simulated quantum states  
  real(double),parameter :: nr_halo_neutrate=20. !! to average halo neut-rate 
  real(double)           :: colrad_threshold !! to speed up simulation!
  type cell_type_plasma
     !! kinetic profiles
     real(double) :: te   !! Electron temperature
     real(double) :: ti   !! Ion temperature
     real(double) :: dene !! electron density
     real(double) :: denp !! D-ions density
     real(double) :: deni !! Carbon density
     real(double) :: denf !! Fast-ion-density average in energy and pitch
     real(double) :: zeff !! zeff
     real(double), dimension(3) :: E    !! Electric field
     real(double), dimension(3) :: vrot !! Plasma rotation
     real(double), dimension(3) :: B    !! Magnetic field
  end type cell_type_plasma
  type cell_type
     type(cell_type_plasma)                  :: plasma  !! Kinetic profiles
     real(double)                            :: rho     !! Normalized flux coord
     real(float),dimension(:)  ,allocatable  :: los_wght!! Weights of LOS
     real(float),dimension(7,nlevs)          :: neut_dens! Density of neutrals
  end type cell_type
  type nbi_type 
     integer(long)  :: number            !! number of the NBI
     real(double)   :: dv                !! half width in y direction
     real(double)   :: dw                !! half width in z direction
     real(double)   :: focy              !! focal lenght in y direction
     real(double)   :: focz              !! focal lenght in z direction
     real(double), dimension(3)  :: divay!! divergence in y direction
     real(double), dimension(3)  :: divaz!! divergence in z direction
     real(double), dimension(3)   :: species_mix
     real(double), dimension(3)   :: xyz_pos !! position of source
     real(double)                 :: einj    !! NBI voltage  [kV]
     real(double)                 :: pinj    !! NBI power    [MW]
     real(double)                 :: vinj    !! NBI velocity [cm/s]
     real(double), dimension(3,3) :: Arot    !! Rotation matrizes of NBI
     real(double), dimension(3,3) :: Brot  
     real(double), dimension(3,3) :: Crot 
  end type nbi_type
  type grid_type
     real(double), dimension(3) :: dr    !! dx, dy, dz
     real(double)               :: drmin !! min(dx,dy,dz)
     real(float)                :: dv    !! volume of cells
     integer(long)              :: Nx    !! Nr. of cells in x direction
     integer(long)              :: Ny    !! Nr. of cells in y direction
     integer(long)              :: Nz    !! Nr. of cells in z direction
     integer(long)              :: ntrack!! Maximum Nr. of cells for tracking
     real(double), dimension(:), allocatable        :: xx,xxc
     real(double), dimension(:), allocatable        :: yy,yyc 
     real(double), dimension(:), allocatable        :: zz,zzc
  end type grid_type
  type atomic_type
     !! energy array
     real(double)                                  :: d_eb_qp
     real(double)                                  :: d_eb_qi
     real(double)                                  :: d_eb_qe
     real(double)                                  :: d_eb_neut
     integer(long)                                 :: nr_eb_qp
     integer(long)                                 :: nr_eb_qi
     integer(long)                                 :: nr_eb_qe
     integer(long)                                 :: nr_eb_neut
     !! temperature array
     real(double)                                  :: d_ti_qp
     real(double)                                  :: d_ti_qi
     real(double)                                  :: d_te_qe
     real(double)                                  :: d_ti_neut
     integer(long)                                 :: nr_ti_qp
     integer(long)                                 :: nr_ti_qi
     integer(long)                                 :: nr_te_qe
     integer(long)                                 :: nr_ti_neut
     !! TABLES
     real(double), dimension(:,:,:,:) , allocatable :: qp 
     real(double), dimension(:,:,:,:) , allocatable :: qi
     real(double), dimension(:,:,:,:) , allocatable :: qe
     real(double), dimension(:,:,:)   , allocatable :: neut
     real(double),dimension(:,:)     , allocatable :: einstein
  end type atomic_type
  type distri_type
     real(double), dimension(:)        ,allocatable :: energy  !! Energy array
     real(double), dimension(:)        ,allocatable :: pitch   !! Pitch array
     real(float),  dimension(:,:,:,:,:),allocatable :: fbm     !! Distribution
     real(double)  :: emin
     real(double)  :: emax
     real(double)  :: pitchmin
     real(double)  :: pitchmax
     integer(long) :: nenergy
     integer(long) :: npitch
  end type distri_type
  type spec_type
     real(double), dimension(:,:,:),allocatable :: spectra
     real(double),dimension(:,:),  allocatable :: xyzlos
     real(double),dimension(:,:),  allocatable :: xyzhead
     real(double),dimension(:),    allocatable :: headsize
     integer(long) :: nchan
     integer(long) :: nlambda
     real(double)  :: dlambda
     real(double)  :: lambdamin
     real(double)  :: lambdamax
  end type spec_type

  type npa_type
     real(double), dimension(:,:) ,allocatable  :: v    !! velocity array
     real(double), dimension(:,:) ,allocatable  :: ipos !! initial position arra
     real(double), dimension(:,:) ,allocatable  :: fpos !! final position array
     real(double), dimension(:)   ,allocatable  :: wght !! weight
     real(double), dimension(:)   ,allocatable  :: size !! active area of detector
     integer(long)                 :: counter
     logical                       :: at_detector  
  end type npa_type
  type inputs_type
     integer(long) :: shot_number
     real(double)  :: time
     character(15) :: runid
     character(4)  :: diag 
     !! Monte Carlo Settings
     integer(long) :: nr_fida
     integer(long) :: nr_ndmc
     integer(long) :: nr_dcx
     integer(long) :: nr_halo 
     integer(long) :: nr_npa
     !! general settings
     integer(long) :: nospec
     integer(long) :: nofida
     integer(long) :: load_neutrals
     integer(long) :: npa 
     integer(long) :: guidingcenter  !! 0 for full-orbit F
     integer(long) :: f90brems       !! 0 to use IDL v.b.
     integer(long) :: calc_wght
     !! Plasma parameters
     integer(long) :: impurity_charge
     real(double)  :: sigma_pi_ratio
     real(double)  :: btipsign 
     real(double)  :: ai   !! atomic mass of plasma ions
     real(double)  :: ab   !! atomic mass of beam neutrals
     !! Settings for weight function calculation
     integer(long) :: nr_wght
     integer(long) :: ichan_wght
     real(double)  :: emax_wght
     real(double)  :: dwav_wght
     real(double)  :: wavel_start_wght
     real(double)  :: wavel_end_wght
  end type inputs_type
  !! definition of the structures:
  type(cell_type), dimension(:,:,:), allocatable  :: cell   
  type(nbi_type)    :: nbi      
  type(grid_type)   :: grid      
  type(atomic_type) :: atomic    
  type(distri_type) :: distri
  type(npa_type)    :: npa
  type(spec_type)   :: spec
  type(inputs_type) :: inputs
  !! routines:
  public :: ndmc      
  public :: halo
  public :: fida
  !! structures:
  public :: nbi         
  public :: grid         
  public :: atomic      
  public :: distri     
  public :: inputs  
  public :: npa  
  !! routines to read inputs
  public :: read_inputs
  public :: read_los
  public :: read_plasma
  public :: read_atomic
  public :: read_fbm  
  public :: read_neutrals
  public :: write_neutrals 
  public :: write_fida_spectra
  public :: write_nbi_halo_spectra
contains  
  !****************************************************************************
  subroutine read_inputs
    character(120)   :: filename
    integer(long) :: i,j,k
    print*,'---- loading inputs -----' 
    filename=trim(adjustl(result_dir))//"/inputs.dat"
    open(66,form='formatted',file=filename)
    read(66,*) !# FIDASIM input file created...
    read(66,"(A120)") root_dir
    read(66,*) inputs%shot_number
    read(66,*) inputs%time
    read(66,*) inputs%runid
    read(66,*) inputs%diag
    read(66,*) !# general settings:
    read(66,*) inputs%nospec
    read(66,*) inputs%nofida
    read(66,*) inputs%npa
    read(66,*) inputs%load_neutrals
    read(66,*) inputs%guidingcenter
    read(66,*) inputs%f90brems
    read(66,*) inputs%calc_wght
    read(66,*) !# weight function settings
    read(66,*) inputs%nr_wght
    read(66,*) inputs%ichan_wght
    read(66,*) inputs%emax_wght
    read(66,*) inputs%dwav_wght
    read(66,*) inputs%wavel_start_wght
    read(66,*) inputs%wavel_end_wght
    read(66,*) !# Monte Carlo settings:
    read(66,*) inputs%nr_fida
    read(66,*) inputs%nr_ndmc   
    read(66,*) inputs%nr_halo   
    inputs%nr_dcx = inputs%nr_halo
    read(66,*) inputs%impurity_charge 
    read(66,*) !# Location of transp cdf file:   
    read(66,*) !cdf-file    
    read(66,*) !# discharge parameters:         
    read(66,*) inputs%btipsign
    read(66,*) inputs%ab       
    read(66,*) inputs%ai    
    read(66,*) !# wavelength grid:    
    read(66,*) spec%nlambda 
    read(66,*) spec%lambdamin
    read(66,*) spec%lambdamax
    spec%dlambda=(spec%lambdamax-spec%lambdamin)/spec%nlambda
    read(66,*) !# simulation grid: 
    read(66,*) grid%Nx 
    read(66,*) grid%Ny      
    read(66,*) grid%Nz
    allocate(grid%xx(grid%Nx)  &
         ,   grid%yy(grid%Ny)  &
         ,   grid%zz(grid%Nz))
    do i=1,grid%Nx 
       read(66,*) grid%xx(i)
    enddo
    do i=1,grid%Ny 
       read(66,*) grid%yy(i)
    enddo
    do i=1,grid%Nz 
       read(66,*) grid%zz(i)
    enddo
    read(66,*) !# Neutral beam injection
    read(66,*) nbi%dv
    read(66,*) nbi%dw  
    read(66,*) nbi%number
    read(66,*) nbi%divay(1) 
    read(66,*) nbi%divay(2) 
    read(66,*) nbi%divay(3) 
    read(66,*) nbi%divaz(1)  
    read(66,*) nbi%divaz(2)  
    read(66,*) nbi%divaz(3)  
    read(66,*) nbi%focy  
    read(66,*) nbi%focz   
    read(66,*) nbi%einj 
    read(66,*) nbi%pinj 
    read(66,*) !# Species-mix (Particles):
    read(66,*) nbi%species_mix(1)
    read(66,*) nbi%species_mix(2)
    read(66,*) nbi%species_mix(3)
    read(66,*) !#position of NBI source in xyz coords:
    read(66,*) nbi%xyz_pos(1)
    read(66,*) nbi%xyz_pos(2)
    read(66,*) nbi%xyz_pos(3)
    read(66,*) !# 3 rotation matrizes 3x3
    do j=1,3 
       do k=1,3
          read(66,*) nbi%Arot(j,k)
          read(66,*) nbi%Brot(j,k)
          read(66,*) nbi%Crot(j,k)
       enddo
    enddo
    close(66) 
    nbi%vinj=sqrt(2.d0*nbi%einj*1.d3 &
            *e0/(inputs%ab*mass_u))*1.d2 !! [cm/s]
    allocate(cell(grid%Nx,grid%Ny,grid%Nz)) 
    allocate(grid%xxc(grid%Nx) &
         ,   grid%yyc(grid%Ny) &
         ,   grid%zzc(grid%Nz))
    grid%dr(1)=grid%xx(2)-grid%xx(1)
    grid%dr(2)=grid%yy(2)-grid%yy(1) 
    grid%dr(3)=grid%zz(2)-grid%zz(1)
    grid%xxc(:)=grid%xx(:)+0.5d0*grid%dr(1)
    grid%yyc(:)=grid%yy(:)+0.5d0*grid%dr(2)
    grid%zzc(:)=grid%zz(:)+0.5d0*grid%dr(3)
    minpos=minloc(grid%dr)
    grid%drmin=grid%dr(minpos(1))
    grid%dv=real(grid%dr(1)*grid%dr(2)*grid%dr(3),float)
    grid%ntrack=grid%Nx+grid%Ny+grid%Nz 
    print*,                         'Shot  :',inputs%shot_number
    print*,                         'Time:',int(inputs%time*1.e3),'ms'
    print*, 'NBI #',nbi%number+1
    print*,'Injected NBI power   :', nbi%pinj
    print*,'Injected NBI voltage :', nbi%einj
    if(inputs%npa .eq. 1)then
       print*, 'this is a NPA simultation!'
       inputs%nr_fida=inputs%nr_fida*1000
       inputs%nr_npa=1000000
       allocate(npa%v(inputs%nr_npa,3)) 
       allocate(npa%ipos(inputs%nr_npa,3)) 
       allocate(npa%fpos(inputs%nr_npa,3)) 
       allocate(npa%wght(inputs%nr_npa))
       npa%counter=0
  !     npa%size=spec%headsize(1)    
    endif
  end subroutine read_inputs
  !****************************************************************************
 subroutine read_los
    character(120)  :: filename
    integer(long) :: i, j, k,ichan
    !filename="RESULTS/"//trim(adjustl(inputs%runid))//"/los.bin" 
    filename=trim(adjustl(result_dir))//"/los.bin"
    print*,'---- loading detector information ----'
    open(66,form='unformatted',file=filename,access='stream')
    read(66)spec%nchan
    allocate(spec%xyzhead(spec%nchan,3))
    allocate(spec%xyzlos(spec%nchan,3))
    allocate(spec%headsize(spec%nchan))
    do ichan = 1, spec%nchan
       read(66)spec%xyzhead(ichan,1)
       read(66)spec%xyzhead(ichan,2)
       read(66)spec%xyzhead(ichan,3)
       read(66)spec%headsize(ichan)
       read(66) spec%xyzlos(ichan,1)
       read(66) spec%xyzlos(ichan,2)
       read(66) spec%xyzlos(ichan,3)
    enddo
    read(66)inputs%sigma_pi_ratio
    do i = 1,grid%Nx
       do j = 1, grid%Ny
          do k = 1, grid%Nz           
             allocate(cell(i,j,k)%los_wght(spec%nchan))
             do ichan = 1, spec%nchan
                read(66) cell(i,j,k)%los_wght(ichan)
             enddo
          enddo
       enddo
    enddo
    close(66)

    allocate(spec%spectra(spec%nlambda,spec%nchan,8))
    if (inputs%npa.eq.1) then
        allocate(npa%size(spec%nchan))
        npa%size=spec%headsize
    endif
  end subroutine read_los

  !***************************************************************************!
  subroutine read_plasma
    character(120)          :: filename
    integer(long) :: Nx,Ny,Nz,i, j, k 

real(double) :: mte,mti,mdene,mdeni,mvrot,mb,me,mrho,mdenf,mzeff
mte=0.d0
mti=0.d0
mdene=0.d0
mdeni=0.d0
mvrot=0.d0
mb=0.d0
me=0.d0
mrho=0.d0
mdenf=0.d0
mzeff=0.d0

    filename=trim(adjustl(result_dir))//"/plasma.bin"
    print*,'---- loading plasma data from ', filename
    open(66,form='unformatted',file=filename,access='stream')
    read(66)Nx
    read(66)Ny
    read(66)Nz
    do i = 1,Nx
       do j = 1,Ny
          do k = 1,Nz
             read(66) cell(i,j,k)%plasma%te  , cell(i,j,k)%plasma%ti     &
                  ,cell(i,j,k)%plasma%dene   , cell(i,j,k)%plasma%denp   &
                  ,cell(i,j,k)%plasma%deni   , cell(i,j,k)%plasma%vrot(1)&
                  ,cell(i,j,k)%plasma%vrot(2), cell(i,j,k)%plasma%vrot(3)&
                  ,cell(i,j,k)%plasma%B(1)   , cell(i,j,k)%plasma%B(2)   &
                  ,cell(i,j,k)%plasma%B(3)   , cell(i,j,k)%plasma%E(1)   &
                  ,cell(i,j,k)%plasma%E(2)   , cell(i,j,k)%plasma%E(3)   &
                  ,cell(i,j,k)%rho           , cell(i,j,k)%plasma%denf   &
                  ,cell(i,j,k)%plasma%zeff
          enddo
       enddo
    enddo
    close(66)
    do i = 1,Nx
       do j = 1,Ny
          do k = 1,Nz
if(cell(i,j,k)%plasma%te .gt. mte) then 
  mte=cell(i,j,k)%plasma%te
endif
if(cell(i,j,k)%plasma%ti .gt. mti) then 
  mti=cell(i,j,k)%plasma%ti
endif
if(cell(i,j,k)%plasma%dene .gt. mdene) then 
  mdene=cell(i,j,k)%plasma%dene
endif
if(cell(i,j,k)%plasma%deni .gt. mdeni) then 
  mdeni=cell(i,j,k)%plasma%deni
endif
if(cell(i,j,k)%plasma%vrot(1) .gt. mvrot) then 
  mvrot=cell(i,j,k)%plasma%vrot(1)
endif
if(cell(i,j,k)%plasma%b(1) .gt. mb) then 
  mb=cell(i,j,k)%plasma%b(1)
endif
if(cell(i,j,k)%plasma%e(1) .gt. me) then 
  me=cell(i,j,k)%plasma%e(1)
endif
if(cell(i,j,k)%rho .gt. mrho) then 
  mrho=cell(i,j,k)%rho
endif
if(cell(i,j,k)%plasma%denf .gt. mdenf) then 
  mdenf=cell(i,j,k)%plasma%denf
endif
if(cell(i,j,k)%plasma%zeff .gt. mzeff) then 
  mzeff=cell(i,j,k)%plasma%zeff
endif
          enddo
       enddo
    enddo

print*,'Te:    ',mte
print*,'Ti:    ',mti
print*,'dene:  ',mdene
print*,'deni:  ',mdeni
print*,'vrot:  ',mvrot
print*,'B:     ',mb
print*,'E:     ',me
print*,'rho:   ',mrho
print*,'denf:  ',mdenf
print*,'zeff:  ',mzeff

  end subroutine read_plasma

  !****************************************************************************
  ! Read vb written by IDL routine  WWH 6/2013
  subroutine read_bremsstrahlung
    character(120)          :: filename
    integer(long) :: i
    real(double) :: vb
    real(double), dimension(:)  , allocatable :: brems
    filename=trim(adjustl(result_dir))//"/bremsstrahlung.bin"
    print*,'---- loading bremsstrahlung data from ', filename
    open(66,form='unformatted',file=filename,access='stream')
    allocate(brems(spec%nlambda))
    do i = 1,spec%nchan
      read(66) vb
      spec%spectra(:,i,brems_type)=vb
    enddo
    close(66)
    deallocate(brems)
  end subroutine read_bremsstrahlung

  !****************************************************************************
  subroutine read_atomic
    character(120)  :: filename
    integer         :: n,m !! initial/final state
    integer(long)   :: nlev

   !-------------------ELECTRON EXCITATION/IONIZATION TABLE--------
    filename=trim(adjustl(root_dir))//"TABLES/qetable.bin"
    open(66,form='unformatted',file=filename,access='stream')
    read(66) atomic%nr_te_qe
    read(66) atomic%d_te_qe
    read(66) atomic%nr_eb_qe
    read(66) atomic%d_eb_qe
    read(66) nlev
    if(nlev.ne.nlevs) stop 'stop at "read atomics qetable"'
    allocate(atomic%qe(nlevs+1,nlevs,atomic%nr_eb_qe,atomic%nr_te_qe))
    read(66) atomic%qe(:,:,:,:)
    close(66)

    !-------------------Deuterium EXCITATION/IONIZATION/CX TABLE------
    filename=trim(adjustl(root_dir))//"TABLES/qptable.bin"
    open(66,form='unformatted',file=filename,access='stream')
    read(66) atomic%nr_ti_qp
    read(66) atomic%d_ti_qp
    read(66) atomic%nr_eb_qp
    read(66) atomic%d_eb_qp
    read(66) nlev
    if(nlev.ne.nlevs)then
       print*, atomic%nr_ti_qp,atomic%d_ti_qp, atomic%nr_eb_qp,atomic%d_eb_qp
       print*, nlev,nlevs
       stop 'stop at "read qptable"'
    endif
    allocate(atomic%qp(nlevs+1,nlevs,atomic%nr_eb_qp,atomic%nr_ti_qp))
    read(66) atomic%qp(:,:,:,:)
    close(66)

    !------------------ m-resolved CHARGE EXCHANGE cross-sections  ---
    ! H(+) + H(n) --> H(m) + H(+)
    ! energy in keV/amu
    filename=trim(adjustl(root_dir))//"TABLES/neuttable.bin"
    open(66,form='unformatted',file=filename,access='stream')
    read(66) atomic%nr_eb_neut
    read(66) atomic%d_eb_neut
    read(66) nlev
    if(nlev.ne.nlevs) stop 'stop at read tables (neut_rates)'
    allocate(atomic%neut(nlevs,nlevs,atomic%nr_eb_neut))
    read(66) atomic%neut(:,:,:)
    close(66)

    !-------------------Impurity EXCITATION/IONIZATION/CX TABLE--------
   if(inputs%impurity_charge.lt.5.or.inputs%impurity_charge.gt.7) &
         stop 'wrong impurity charge!'
    if(inputs%impurity_charge.eq.5) &
         filename=trim(adjustl(root_dir))//"TABLES/qbtable.bin"
    if(inputs%impurity_charge.eq.6) &
         filename=trim(adjustl(root_dir))//"TABLES/qctable.bin" 
    if(inputs%impurity_charge.eq.7) &
         filename=trim(adjustl(root_dir))//"TABLES/qntable.bin"
    open(66,form='unformatted',file=filename,access='stream')
    read(66) atomic%nr_ti_qi
    read(66) atomic%d_ti_qi
    read(66) atomic%nr_eb_qi
    read(66) atomic%d_eb_qi
    read(66) nlev
    if(nlev.ne.nlevs) stop 'stop at "read atomics qptable"'
    allocate(atomic%qi(nlevs+1,nlevs,atomic%nr_eb_qi,atomic%nr_ti_qi))
    read(66) atomic%qi(:,:,:,:)
    close(66)

    !-------------------EINSTEIN COEFFICIENTS ----------------------
    filename='TABLES/einstein.dat'
    filename=trim(adjustl(root_dir))//"TABLES/einstein.dat"
    open(66,file=filename)
    read(66,*)! 
    read(66,*)!
    read(66,*) nlev
    if(nlev.ne.nlevs) stop 'stop at "read atomics"'
    allocate(atomic%einstein(nlevs,nlevs))
    do n=1,nlevs !! lower level
       do m=1,nlevs  !! upper level
          read(66,*)atomic%einstein(m,n) 
       enddo
    enddo
    close(66)

  end subroutine read_atomic

  !****************************************************************************
  !-------------------FASTION DISTRIBUTION FUNCTION ----------------------
  subroutine read_fbm
    character(120)  :: filename
    integer(long) :: i,j,k
    !! sub_grid
    integer(long) :: l,m,n,nxsub,nysub,nzsub
    real(double)    :: dxsub,dysub,dzsub,rdummy
    real(float)    :: xsub, ysub, rsub, zsub
    !! for transp input grid
    integer(long) :: transp_nzones
    real(float), dimension(:), allocatable     :: transp_r,transp_z
    real(float), dimension(:,:,:), allocatable :: transp_fbm   
    filename=trim(adjustl(result_dir))//"/transp_fbm.bin"
    print*,'---- loading fast ion distribution function ----'
    open(66,form='unformatted',file=filename,access='stream')
    read(66) transp_nzones
    allocate(transp_r(transp_nzones))
    do i=1,transp_nzones
       read(66) transp_r(i)
    enddo
    allocate(transp_z(transp_nzones))
    do i=1,transp_nzones
       read(66) transp_z(i)
    enddo
    read(66) distri%nenergy
    read(66) distri%emin
    read(66) distri%emax
    allocate(distri%energy(distri%nenergy))
    do i=1, distri%nenergy
       read(66)distri%energy(i)
    enddo
    read(66) distri%npitch
    read(66) distri%pitchmin
    read(66) distri%pitchmax
    allocate(distri%pitch(distri%npitch))
    do i=1, distri%npitch
       read(66) rdummy
        distri%pitch(i) = rdummy*inputs%btipsign
    enddo
    allocate(transp_fbm(distri%nenergy,distri%npitch,transp_nzones))
    do i=1, distri%nenergy
       do j=1,distri%npitch
          do k=1,transp_nzones
             read(66) transp_fbm(i,j,k)
          enddo
       enddo
    enddo
    close(66)
    !! map TRANSP velocity space on grid
    !! Use spatial resolution of sub grids with ~1cm
    nxsub=int(grid%dr(1)/1.d0+0.5d0)
    dxsub=grid%dr(1)/nxsub
    nysub=int(grid%dr(2)/1.d0+0.5d0)
    dysub=grid%dr(2)/nysub
    nzsub=int(grid%dr(3)/1.d0+0.5d0)
    dzsub=grid%dr(3)/nzsub
    allocate(distri%fbm(grid%Nx,grid%Ny,grid%Nz,distri%nenergy,distri%npitch)) 
    distri%fbm=0.
    do i=1, grid%nx
       do j=1,grid%ny
          do k=1,grid%nz
             if (cell(i,j,k)%plasma%denf.gt.0. .and. &
                  sum(cell(i,j,k)%neut_dens(:,:)).gt.0.)then
                if(inputs%npa.eq.1)then
                   if(.not.(any(cell(i,j,k)%los_wght.gt.0.)))cycle
                endif
                do l=1,nxsub
                   do m=1,nysub
                      do n=1,nzsub
                         xsub=real(grid%xx(i)+(l-0.5d0)*dxsub,float)
                         ysub=real(grid%yy(j)+(m-0.5d0)*dysub,float)
                         rsub=sqrt(xsub**2+ysub**2)
                         zsub=real(grid%zz(k)+(n-0.5d0)*dzsub,float)
                         minpos=minloc((transp_r-rsub)**2+(transp_z-zsub)**2)
                         distri%fbm(i,j,k,:,:)= distri%fbm(i,j,k,:,:)  &
                              + transp_fbm(:,:,minpos(1))
                      enddo
                   enddo
                enddo
             endif
          enddo
       enddo
    enddo
    distri%fbm= distri%fbm/(nxsub*nysub*nzsub)
    deallocate(transp_r) 
    deallocate(transp_z)  
    deallocate(transp_fbm)  
  end subroutine read_fbm

  subroutine write_neutrals
    integer(long) :: i,j,k,n 
    character(120)  :: filename
    filename=trim(adjustl(result_dir))//"/neutrals.bin"    
    open (66, form='unformatted',file =filename,access='stream')
    write(66)real(inputs%shot_number,float )
    write(66)real(inputs%time)
    write(66)real(grid%Nx,float)
    write(66)real(grid%Ny,float)
    write(66)real(grid%Nz,float) 
    write(66)real(nlevs  ,float) 
    do i = 1, grid%Nx
       do j = 1, grid%Ny
          do k = 1, grid%Nz 
             do n= 1, nlevs
                !!   write(66)real(cell(i,j,k)%neut_dens(fida_type,n),float)
                write(66)real(cell(i,j,k)%neut_dens(nbif_type,n),float)
                write(66)real(cell(i,j,k)%neut_dens(nbih_type,n),float)
                write(66)real(cell(i,j,k)%neut_dens(nbit_type,n),float)
                write(66)real(cell(i,j,k)%neut_dens(halo_type,n),float)
             enddo
          enddo
       enddo
    enddo
    close (66)
    print*, 'neutral density written to:      ',filename
  end subroutine write_neutrals

  subroutine write_npa
    integer(long) :: i
    character(120)  :: filename
    npa%wght(:)=npa%wght(:)/(pi*npa%size(1)**2)
    filename=trim(adjustl(result_dir))//"/npa.bin"      
    open (66, form='unformatted',file =filename,access='stream')
    write(66)real(inputs%shot_number,float )
    write(66)real(inputs%time,float)
    write(66)real(npa%counter, float)
    do i =1, npa%counter
       write(66)real(npa%ipos(i,1),float)
       write(66)real(npa%ipos(i,2),float)
       write(66)real(npa%ipos(i,3),float)
       write(66)real(npa%fpos(i,1),float)
       write(66)real(npa%fpos(i,2),float)
       write(66)real(npa%fpos(i,3),float)
       write(66)real(npa%v(i,1),float)
       write(66)real(npa%v(i,2),float)
       write(66)real(npa%v(i,3),float)
       write(66)real(npa%wght(i),float)
    enddo
    close (66)
    print*, 'NPA data written to: ',filename
  end subroutine write_npa
  
  subroutine write_nbi_halo_spectra
    integer(long)  :: i,ichan
    character(120)  :: filename
    real(float), dimension(:)  , allocatable :: lambda_arr
    !! ------------------------ calculate wavelength array ------------------ !!
    allocate(lambda_arr(spec%nlambda))
    do i=1,spec%nlambda 
       lambda_arr(i)=real((i-0.5)*spec%dlambda*0.1d0 &
            +spec%lambdamin*0.1d0,float)
    enddo 
    !! convert [Ph/(s*wavel_bin*cm^2*all_directions)] to [Ph/(s*nm*sr*m^2)]!
    spec%spectra(:,:,2:5)=spec%spectra(:,:,2:5)/(0.1d0*spec%dlambda)/(4.d0*pi)*1.d4
    !! write to file
    filename=trim(adjustl(result_dir))//"/nbi_halo_spectra.dat"
    open (66, file =filename,ACTION = 'WRITE')
    write(66,*)real( inputs%shot_number ,float)
    write(66,*)real( inputs%time        ,float)
    write(66,*)inputs%diag
    write(66,*)real( spec%nchan        ,float)
    write(66,*)real( spec%nlambda     ,float)
    write(66,*)'--- wavelength array: -----'
    do i=1,spec%nlambda
       write(66,*)real(lambda_arr(i)    ,float)  
    enddo
    write(66,*)'-------- spectra:----------' 
    do i=1,spec%nlambda   
       do ichan=1,spec%nchan 
          write(66,*)real(spec%spectra(i,ichan,nbif_type),float)
          write(66,*)real(spec%spectra(i,ichan,nbih_type),float)
          write(66,*)real(spec%spectra(i,ichan,nbit_type),float)
          write(66,*)real(spec%spectra(i,ichan,halo_type),float)
          write(66,*)real(spec%spectra(i,ichan,brems_type),float)
       enddo
    enddo
    close (66)
    !! result arrays 
    deallocate(lambda_arr)
    print*, 'NBI and HALO spectra written to: ', filename
  end subroutine write_nbi_halo_spectra

 subroutine write_fida_spectra
    integer(long)  :: i,ichan
    character(120)  :: filename
    real(float), dimension(:)  , allocatable :: lambda_arr
    !! ------------------------ calculate wavelength array ------------------ !!
    allocate(lambda_arr(spec%nlambda))
    do i=1,spec%nlambda 
       lambda_arr(i)=real((i-0.5)*spec%dlambda*0.1d0 &
            +spec%lambdamin*0.1d0,float)
    enddo 
    !! convert [Ph/(s*wavel_bin*cm^2*all_directions)] to [Ph/(s*nm*sr*m^2)]!
    spec%spectra(:,:,fida_type)=spec%spectra(:,:,fida_type)/(0.1d0*spec%dlambda)/(4.d0*pi)*1.d4
    !! write to file
    filename=trim(adjustl(result_dir))//"/fida_spectra.dat"  
    open (66, file =filename,ACTION = 'WRITE')
    write(66,*)real( inputs%shot_number ,float)
    write(66,*)real( inputs%time        ,float)
    write(66,*)inputs%diag
    write(66,*)real( spec%nchan        ,float)
    write(66,*)real( spec%nlambda     ,float)
    write(66,*)'--- wavelength array: -----'
    do i=1,spec%nlambda
       write(66,*)real(lambda_arr(i)    ,float)  
    enddo
    write(66,*)'-------- spectra:----------' 
    do i=1,spec%nlambda
       do ichan=1,spec%nchan 
          write(66,*)real(spec%spectra(i,ichan,fida_type),float)
       enddo
    enddo
    close (66)
    !! result arrays 
    deallocate(lambda_arr)
    print*, 'FIDA spectra written to ',filename
  end subroutine write_fida_spectra

  subroutine read_neutrals
    integer(long) :: i,j,k,n 
    character(120)  :: filename
    real(float) :: fdummi
    print*,'---- load neutrals RESULTS/neutrals.bin ----' 
    filename=trim(adjustl(result_dir))//"/neutrals.bin" 
    open (66, form='unformatted',file =filename,access='stream')
    read(66)fdummi
    read(66)fdummi
    read(66)fdummi
    read(66)fdummi
    read(66)fdummi
    read(66)fdummi
    do i = 1, grid%Nx
       do j = 1, grid%Ny
          do k = 1, grid%Nz 
             do n= 1, nlevs
                read(66) fdummi
                cell(i,j,k)%neut_dens(nbif_type,n) = fdummi
                read(66) fdummi
                cell(i,j,k)%neut_dens(nbih_type,n) = fdummi
                read(66) fdummi
                cell(i,j,k)%neut_dens(nbit_type,n) = fdummi
                read(66) fdummi
                cell(i,j,k)%neut_dens(halo_type,n) = fdummi
              enddo
           enddo
        enddo
     enddo
     close (66)
   end subroutine read_neutrals



  
  !*****************************************************************************
  !------------random number generator-----------------------------------------
  !*****************************************************************************
  function ran(idum)
    !!uniform random number generator from NUMERICAL RECEPIES
    integer(long), intent(INOUT) :: idum
    real(float)                  :: ran
    integer(long), parameter     :: IA=16807,IM=2147483647,IQ=127773,IR=2836
    real(float),   save          :: am
    integer(long), save          :: ix=-1,iy=-1,k
    if (idum <= 0 .or. iy < 0) then !Initialize.
       am=nearest(1.0,-1.0)/IM
       iy=ior(ieor(888889999,abs(idum)),1)
       ix=ieor(777755555,abs(idum))
       idum=abs(idum)+1 !Set idum positive.
    endif
    ix=ieor(ix,ishft(ix,13)) !Marsaglia shift sequence with period 232 − 1.
    ix=ieor(ix,ishft(ix,-17))
    ix=ieor(ix,ishft(ix,5))
    k=iy/IQ !Park-Miller sequence by Schrage’s method,
    iy=IA*(iy-k*IQ)-IR*k !period 231 − 2.
    if (iy < 0) iy=iy+IM
    ran=am*ior(iand(IM,ieor(ix,iy)),1) !Combine the two generators
  end function ran 
  subroutine randn(randomn)
    !!Box Mueller Method to calculate normal distribution
    real(double), dimension(:), intent(inout):: randomn
    integer(long)                            :: idum
    integer(long)                            :: nran
    integer(long)                            :: i
    real(double)                              :: x1,x2,w  
    randomn=0.d0
    idum=1
    nran=size(randomn)
    i=1
    !$OMP CRITICAL(randomA)
    do while (i <= nran)
       w=1.d0
       do while (w >= 1.d0) 
          x1=2.d0*ran(idum)-1.d0
          x2=2.d0*ran(idum)-1.d0
          w=x1*x1+x2*x2
       enddo
       w=sqrt((-2.d0*log(w))/w)
       randomn(i)=x1*w
       i=i+1
       if (i > nran) exit
       randomn(i)=x2*w
       i=i+1
    enddo
  !$OMP END CRITICAL(randomA)
  end subroutine randn 
  subroutine randu(randomu)
    real(double), dimension(:), intent(inout):: randomu
    integer(long)                            :: idum  
    integer(long)                            :: nran
    integer(long)                            :: i
    randomu=0.d0
    idum=1
    nran=size(randomu)
    !$OMP CRITICAL(randomB)
    do i=1,nran
       randomu(i)=ran(idum)
    enddo
    !$OMP END CRITICAL(randomB)
  end subroutine randu  
  subroutine randseed(seed)
    integer(long), intent(in) :: seed
    integer(long)             :: idum 
    real(float)               :: dummy
    idum=seed*(-1)
    !$OMP CRITICAL(randomC)
    dummy=ran(idum)
    !$OMP END CRITICAL(randomC)
  end subroutine randseed 
  
  subroutine neut_rates(denn,vi,vn,rates)
    !GET neutralization rate from tables
    real(double),  dimension(nlevs),intent(in):: denn !!density of neutrals cm-3
    real(double), dimension(3),    intent(in) :: vi,vn!!of neutrals/ions (cm/s)
    real(double), dimension(nlevs),intent(out):: rates!! rates
    real(double), dimension(nlevs,nlevs)      :: neut !!rate coeff
    real(double)                :: eb            !! relative Energy
    real(double)                :: vrel          !! relative velocity
    integer(long)               :: ebi      !! indizes
    !Eeff 
    vrel=sqrt(dot_product(vi-vn,vi-vn))
    eb=v_to_E*vrel**2  ! [kev/amu]
    ebi= int(eb/atomic%d_eb_neut+0.5)+1   
    if(ebi.gt.atomic%nr_eb_neut) stop 'EB out of range of neuttable!'
    neut=dble(atomic%neut(1:nlevs,1:nlevs,ebi))
    rates=matmul(neut(:,:),denn(:))*vrel
  end subroutine neut_rates
                 

  !****************************************************************************
  !----------------mc_fastion------------------------------------------------
  !****************************************************************************
  subroutine mc_fastion(ac,vi)
    !!mcbeam computes monte carlo velocity of fast ions from distri%fbm
    integer(long),dimension(3), intent(in) :: ac !ind of actual cell
    real(double), dimension(3), intent(out):: vi !velocity [cm/s]
    real(double), dimension(3)             :: bfield    ! Magnetic field vector
    real(double), dimension(3)             :: b    ! nomralized bfield
    real(double), dimension(3)             :: a,c  ! vectors perpendicular to b
    real(double)                           :: eb ,ptch
    integer(long)                          :: ienergy, ipitch ,counter
    real(double)                           :: vabs, phi, sinus
    real(double), dimension(3)             :: randomu3
    real(double), dimension(1)             :: randomu1
    logical                                :: reject
    bfield=cell(ac(1),ac(2),ac(3))%plasma%B(:)              
    b= bfield/sqrt(dot_product(bfield,bfield)) 
    !! --- calculate vectors a,c that are perpendicular to b -- !!
    if (abs(b(3)).eq.1) then
       a=(/1.d0,0.d0,0.d0/)
       c=(/0.d0,1.d0,0.d0/)
    else 
       if (b(3).eq.0.) then
          a=(/0.d0,0.d0,1.d0/)
          c=(/b(2),-b(1), 0.d0/)/sqrt(b(1)**2+b(2)**2)
       else
          a=(/b(2),-b(1),0.d0/)/sqrt(b(1)**2+b(2)**2)
          c=(/ a(2) , -a(1) , (a(1)*b(2)-a(2)*b(1))/b(3) /)
          c=c/sqrt(dot_product(c,c))
       endif
    endif
    !! -- use a rejection method to determine vi from distrbution function --!!
    reject=.true.
    counter=0
    vi=0.
    do while (reject)
       call randu(randomu3)
       eb   = distri%emin +(distri%emax - distri%emin) *randomu3(1)
       ptch= -1.d0 + 2.d0 *randomu3(2) 
       minpos=minloc(abs(eb   - distri%energy))  
       ienergy= minpos(1)
       minpos=minloc(abs(ptch - distri%pitch ))  
       ipitch = minpos(1)   
       if(dble(distri%fbm(ac(1),ac(2),ac(3),ienergy,ipitch)).gt.randomu3(3))then
          call randu(randomu1)
          vabs          = sqrt(eb/(v_to_E*inputs%ab))
          phi           = 2.d0*pi*randomu1(1)
          sinus         = sqrt(1.d0-ptch**2)
          vi(:) = vabs * (sinus*cos(phi)*a + ptch*b + sinus*sin(phi)*c) 
          reject=.false.
       endif
       counter=counter+1
       if(counter.gt.10000)return
    enddo
  end subroutine mc_fastion
 
  !****************************************************************************
  !----------------mc_halo------------------------------------------------
  !****************************************************************************
  subroutine mc_halo(ac,vhalo)
    integer(long), dimension(3) , intent(in)    :: ac    !! index of actual cell
    real(double),    dimension(3) , intent(out) :: vhalo !! velocity [cm/s]
    real(double),    dimension(3)               :: randomn
    call randn(randomn)  
    vhalo(:)=   cell(ac(1),ac(2),ac(3))%plasma%vrot(:) &
         + sqrt(cell(ac(1),ac(2),ac(3))%plasma%ti      &
         * 0.5/(v_to_E*inputs%ai)) &
         * randomn(:) !![cm/s]   
  end subroutine mc_halo

  !****************************************************************************
  !----------------mc_nbi------------------------------------------------------
  !****************************************************************************
  subroutine rotate(uvw_vec,updown,xyz_vec)
    real(double), dimension(3), intent(in) :: uvw_vec !! vector in uvw coords
    integer(long)             , intent(in) :: updown  !! source has two plates
    real(double), dimension(3), intent(out):: xyz_vec !! vector in xyz coords
    real(double), dimension(3)             :: uvz_vec   
    !! rotate uvw vector in vertical dirction 
    if(updown.lt.0) uvz_vec(:)=matmul(nbi%Arot(:,:),uvw_vec(:))
    if(updown.ge.0) uvz_vec(:)=matmul(nbi%Brot(:,:),uvw_vec(:))
    !! rotate uvz_vec by phi_box onto xyz_coordinates
    xyz_vec=matmul(nbi%Crot(:,:),uvz_vec(:))
  end subroutine rotate
  subroutine mc_nbi(vnbi,efrac,rnbi)
    !!-- mc_nbi computes monte carlo velocity and initial start position of
    !!-- NBI neutrals on the FIDASIM grid
    integer(long)             , intent(in)    :: efrac !! energy fraction
    real(double), dimension(3), intent(out)   :: vnbi  !! velocity [cm/s]
    real(double), dimension(3), intent(out), optional :: rnbi  !! postition
    integer(long)                :: jj, updown
    real(double), dimension(3)   :: uvw_pos    !! Start position on ion source
    real(double), dimension(3)   :: xyz_pos    !! Start position on ion source
    real(double), dimension(3)   :: uvw_ray    !! NBI veloicity in uvw coords
    real(double), dimension(2)   :: randomu    !! uniform random numbers
    real(double), dimension(2)   :: randomn    !! normal random numbers
    !! ------------ Random start postion on ion source grid -------------- !!
    call randu(randomu)
    uvw_pos(1) =  0.d0
    uvw_pos(2) =  nbi%dv * 2.d0*(randomu(1)-0.5d0)
    uvw_pos(3) =  nbi%dw * 2.d0*(randomu(2)-0.5d0)
    if(uvw_pos(3).gt.0)then 
       updown=1
    else
       updown=-1
    endif
    call randn(randomn)
    uvw_ray(1)=-1.d0
    uvw_ray(2)=uvw_ray(1)*(uvw_pos(2)/nbi%focy &
         +tan(nbi%divay(efrac)*randomn(1)))
    uvw_ray(3)=uvw_ray(1)*(uvw_pos(3)/nbi%focz &
         +tan(nbi%divaz(efrac)*randomn(2)))
    call rotate(uvw_ray,updown,vnbi(:))
    vnbi(:)=vnbi(:)/sqrt(dot_product(vnbi(:),vnbi(:)))
    call rotate(uvw_pos,updown,xyz_pos)
    xyz_pos(:)=xyz_pos(:)+nbi%xyz_pos(:)
    !! ----------- Determine start postition on FIDASIM grid --------- !!
    if(present(rnbi)) then
       nbi_track: do jj=1,2000
          xyz_pos(1) = xyz_pos(1) + grid%dr(1) * vnbi(1)
          xyz_pos(2) = xyz_pos(2) + grid%dr(1) * vnbi(2)
          xyz_pos(3) = xyz_pos(3) + grid%dr(1) * vnbi(3)
          if ( xyz_pos(1).gt.grid%xx(1) .and. &
               xyz_pos(1).lt.grid%xx(grid%nx)+grid%dr(1).and. & 
               xyz_pos(2).gt.grid%yy(1) .and. &
               xyz_pos(2).lt.grid%yy(grid%ny)+grid%dr(2).and. &
               xyz_pos(3).gt.grid%zz(1) .and. &
               xyz_pos(3).lt.grid%zz(grid%nz)+grid%dr(3)) then
             exit nbi_track
          endif
       enddo nbi_track
       if (jj.ge.2000) then
          print*, 'NBI marker outside of grid!'
          rnbi(:)=(/-1,0,0/)
          return
       endif
       rnbi(:)=xyz_pos(:)
    endif
    !! ---- Determine velocity of neutrals corrected by efrac ---- !!
    vnbi(:) = vnbi(:)*nbi%vinj/sqrt(dble(efrac))
  end subroutine mc_nbi 


  !****************************************************************************
  !----------------mc_start-------------------------------------------------
  !****************************************************************************
  subroutine mc_start(ac,vi,ri)
    !! determine random start position within a cell and correct for gyro
    !! orbit
    integer(long)  , dimension(3)  , intent(in)    :: ac !ind of actual cell
    real(double)   , dimension(3)  , intent(in)    :: vi !velocity [cm/s]
    real(double)   , dimension(3)  , intent(out)   :: ri !starting position
    real(double)   , dimension(3)    :: B   ! Magnetic field vector
    real(double)   , dimension(3)    :: vxB       ! crossproduct
    real(double)   , dimension(3)    :: r_gyro! gyro-radius
    real(double)                     :: one_over_omega! For gyro-radius
    real(double)   , dimension(3)    :: randomu    
    call randu(randomu)  
    ri(1)=grid%xx(ac(1))+ grid%dr(1)*randomu(1)
    ri(2)=grid%yy(ac(2))+ grid%dr(2)*randomu(2) 
    ri(3)=grid%zz(ac(3))+ grid%dr(3)*randomu(3)  
    if (inputs%guidingcenter.eq.1) then  ! WWH 
    B=cell(ac(1),ac(2),ac(3))%plasma%B(:) 
    one_over_omega=inputs%ab*mass_u/(dot_product(B,B)*e0)*1.d-2    
    vxB(1)= (vi(2) * B(3) - vi(3) * B(2))
    vxB(2)= (vi(3) * B(1) - vi(1) * B(3))
    vxB(3)= (vi(1) * B(2) - vi(2) * B(1))
    r_gyro(:)=vxB(:)*one_over_omega
    ri(:)=ri(:)-r_gyro(:) !! '-'because v x B is towards the gyrocenter    
    end if   ! WWH
  end subroutine mc_start


  !****************************************************************************
  !----------------------------- colrad  ------------------------------
  !****************************************************************************
  ! first subroutines for eigenvalue decomposition 
  subroutine RSWAP(a,b)
    real(double) :: a,b, t
    t=a; a=b; b=t
  end subroutine RSWAP
  subroutine balance(n,     &  !size of matrix         
       mat,   &  !input matrix
       scal,  &  !Scaling data
       low,   &  !first relevant row index
       high )   !last relevant row index                 
    integer(long), intent(in) :: n
    real(double)   :: mat(0:n,0:n),scal(0:n)
    integer(long), intent(out) :: high, low
    integer(long), parameter :: basis = 2
    real(double)  :: b2, r, c, f, g, s
    integer(long) :: m, k, i, j, iter
    !*====================================================================*
    !*  balance balances the matrix so that the rows with zero entries    *
    !*  off the diagonal are isolated and the remaining columns and rows  *
    !*  are resized to have one norm close to 1.                          *
    !*   Input parameters:                                                *
    !*      n        integer;  ( n > 0 )                                  *
    !*               Dimension of mat                                     *
    !*      mat      n x n input matrix                                   *
    !*                                                                    *
    !*   Output parameters:                                               *
    !*      mat      n x n scaled matrix                                  *
    !*      low      integer;                                             *
    !*      high     integer;                                             *
    !*               the rows 0 to low-1 and those from high to n-1       *
    !*               contain isolated eigenvalues (only nonzero entry on  *
    !*               the diagonal)                                        *
    !*      scal     vector of size n                                     *
    !*               the vector scal contains the isolated eigenvalues in *
    !*               the positions 0 to low-1 and high to n-1, its other  *
    !*               components contain the scaling factors for           *
    !*               transforming mat.                                    *
    !*====================================================================*
    scal=0.d0
    b2 = basis * basis
    m = 0
    k = n - 1
    iter=1
    do while(iter==1)
       iter = 0
       do j = k, 0, -1
          r = ZERO
          do i = 0, k
             if (i.ne.j)  r = r + DABS(mat(j,i))
          enddo
          if (r == ZERO) then
             scal(k) = j
             if (j.ne.k) then
                do i = 0, k 
                   call RSWAP(mat(i,j), mat(i,k))
                enddo
                do i = m, n-1 
                   call RSWAP(mat(j,i), mat(k,i))
                enddo
             endif
             k=k-1
             iter = 1
          endif
       enddo !j loop
    enddo !while iter=1
    iter=1
    do while (iter==1)
       iter = 0
       do j = m, k
          c = ZERO
          do i = m, k
             if (i.ne.j)  c = c + DABS(mat(i,j))
          enddo
          if (c == ZERO) then
             scal(m) = j
             if (j.ne.m) then
                do i = 0, k 
                   call RSWAP(mat(i,j), mat(i,m))
                enddo
                do i = m, n-1 
                   call RSWAP(mat(j,i), mat(m,i))
                enddo
             endif
             m = m + 1
             iter = 1
          endif
       enddo !j loop
    enddo !while iter=1
    low = m
    high = k
    do i = m, k 
       scal(i) = ONE
    enddo
    iter=1
    do while (iter==1)
       iter = 0
       do i = m, k
          c=ZERO; r=ZERO
          do j = m, k
             if (j.ne.i) then
                c = c + DABS(mat(j,i))
                r = r + DABS(mat(i,j))
             endif
          enddo
          g = r / basis
          f = ONE
          s = c + r
          do while (c < g)
             f = f * basis
             c = c * b2
          enddo
          g = r * basis
          do while (c >= g)
             f = f / basis
             c = c / b2
          enddo
          if ((c + r) / f < 0.95 * s) then
             g = ONE / f
             scal(i) = scal(i) * f
             iter = 1
             do j = m, n-1 
                mat(i,j) = mat(i,j) * g
             enddo
             do j = 0, k  
                mat(j,i) = mat(j,i) * f
             enddo
          endif
       enddo !i loop
    enddo !while iter=1
    return
  end subroutine balance
  subroutine balback(n,     &  !Dimension of matrix .........
       low,   &  !first nonzero row ...........
       high,  &  !last nonzero row ............
       scal,  &  !Scaling data ................
       eivec )   !Eigenvectors ................
    integer(long),intent(in)   :: high, low
    integer(long),intent(in)   :: n
    real(double), intent(in)   ::  scal(0:n)
    real(double), intent(inout):: eivec(0:n,0:n)
    real(double) :: s
    integer(long) :: i,j,k
    !*====================================================================*
    !*  balback reverses the balancing of balance for the eigenvactors.   *
    !*   Input parameters:                                                *
    !*   ----------------                                                 *
    !*      n        integer;  ( n > 0 )                                  *
    !*               Dimension of mat                                     *
    !*      low      integer;                                             *
    !*      high     integer;   see balance                               *
    !*      eivec    n x n matrix of eigenvectors, as computed in  qr2    *
    !*      scal     vector of size n;                                    *
    !*               Scaling data from  balance                           *
    !*   Output parameter:                                                *
    !*   ----------------                                                 *
    !*      eivec    n x n matrix;                                        *
    !*               Non-normalized eigenvectors of the original matrix   *
    !*====================================================================*
    do i = low, high
       s = scal(i)
       do j = 0, n-1  
            eivec(i,j) = eivec(i,j) * s
       enddo
    enddo
    do i = low-1, 0, -1
       k = Int(scal(i))
       if (k.ne.i) then
          do j = 0, n-1
             call RSWAP(eivec(i,j), eivec(k,j))
          enddo
       endif
    enddo
    do i = high + 1, n-1
       k = Int(scal(i))
       if (k.ne.i) then
          do j = 0, n-1
             call RSWAP(eivec(i,j), eivec(k,j))
          enddo
       endif
    enddo
    return
  end subroutine balback
  subroutine elmhes(n,    &  !Dimension of matrix
       low,  &  !first nonzero row ...........
       high, &  !last nonzero row ............
       mat,  &  !input/output matrix .........
       perm )  !Permutation vector ..........
    integer(long),intent(in)   :: n
    integer(long),intent(in)   :: high, low
    real(double), intent(inout):: mat(0:n,0:n)
    integer(long),intent(out)  :: perm(0:n)
    integer(long) :: i, j, m
    real(double) ::  x, y
    !*====================================================================*
    !*  elmhes transforms the matrix mat to upper Hessenberg form.        *
    !*   Input parameters:                                                *
    !*      n        integer;  ( n > 0 )                                  *
    !*               Dimension of mat                                     *
    !*      low      integer;                                             *
    !*      high     integer; see  balance                                *
    !*      mat      n x n matrix                                         *
    !*   Output parameter:                                                *
    !*      mat      n x n matrix;                                        *
    !*               upper Hessenberg matrix; additional information on   *
    !*               the transformation is stored in the lower triangle   *
    !*      perm     integer vector of size n;                            *
    !*               Permutation vector for elmtrans                      *
    !*====================================================================*
    do m = low + 1, high-1
       i = m
       x = ZERO
       do j = m, high
          if (DABS(mat(j,m-1)) > DABS (x)) then
             x = mat(j,m-1)
             i = j
          endif
       enddo
       perm(m) = i
       if (i.ne.m) then
          do j = m - 1, n-1 
             call RSWAP(mat(i,j), mat(m,j))
          enddo
          do j = 0, high 
             call RSWAP(mat(j,i), mat(j,m))
          enddo
       endif
       if (x.ne.ZERO) then
          do i = m + 1, high
             y = mat(i,m-1)
             if (y.ne.ZERO) then
                y = y / x
                mat(i,m-1) = y
                do j = m, n-1 
                   mat(i,j) = mat(i,j) - y * mat(m,j)
                enddo
                do j = 0, high 
                   mat(j,m) = mat(j,m) + y * mat(j,i)
                enddo
             endif
          enddo !i loop
       endif !x <> ZERO
    enddo !m loop
  end subroutine elmhes
  Subroutine elmtrans(n,    &  !Dimension of matrix .........
       low,  &  !first nonzero row ...........
       high, &  !last nonzero row ............
       mat,  &  !input matrix ................
       perm, &  !row permutations ............
       h  )     !Hessenberg matrix ...........
    integer(long),intent(in)   :: n
    integer(long),intent(in)   :: high, low
    real(double), intent(in)   :: mat(0:n,0:n)
    integer(long),intent(in)   :: perm(0:n)
    real(double),intent(out)   :: h(0:n,0:n)
    integer(long) :: i, j, k
    !*====================================================================*
    !*  Elmtrans copies the Hessenberg matrix stored in mat to h.         *
    !*   Input parameters:                                                *
    !*      n        integer;  ( n > 0 )                                  *
    !*               Dimension of  mat and eivec                          *
    !*      low      integer;                                             *
    !*      high     integer; see  balance                                *
    !*      mat      n x n input matrix                                   *
    !*      perm     Integer vector of size n;                            *
    !*               Permutation data from  elmhes                        *
    !*   Output parameter:                                                *
    !*      h        n x n matrix;                                        *
    !*               Hessenberg matrix                                    *
    !*====================================================================*
    do i = 0, n-1
       do k = 0, n-1 
          h(i,k) = ZERO
       enddo
       h(i,i) = ONE
    enddo
    do i = high - 1, low+1, -1
       j = perm(i)
       do k = i + 1, high 
          h(k,i) = mat(k,i-1)
       enddo
       if (i.ne.j) then
          do k = i, high
             h(i,k) = h(j,k)
             h(j,k) = ZERO
          enddo
          h(j,i) = ONE
       endif
    enddo
  end subroutine elmtrans
  subroutine Comdiv(ar,     &       !Real part of numerator ..........
       ai,     &       !Imaginary part of numerator .....
       br,     &       !Real part of denominator ........
       bi,     &       !Imaginary part of denominator ...
       cr,     &       !Real part of quotient ...........
       ci,     &       !Imaginary part of quotient ......
       rc )            !return code .....................
    real(double) ::  ar,ai,br,bi,cr,ci
    integer(long) :: rc
    real(double) :: tmp
    !*====================================================================*
    !*  Complex division  c = a / b                                       *
    !*   Input parameters:                                                *
    !*   ================                                                 *
    !*      ar,ai    real, imaginary parts of numerator                   *
    !*      br,bi    real, imaginary parts of denominator                 *
    !*   Output parameters:                                               *
    !*   ==================                                               *
    !*      cr,ci     real , imaginary parts of the quotient              *
    !*====================================================================*
    if (br == ZERO.AND.bi == ZERO) then
       rc = 1
       return
    endif
    if (dabs(br) > dabs(bi)) then
       tmp = bi / br
       br  = tmp * bi + br
       cr  = (ar + tmp * ai) / br
       ci  = (ai - tmp * ar) / br
    else
       tmp = br / bi
       bi  = tmp * br + bi
       cr  = (tmp * ar + ai) / bi
       ci  = (tmp * ai - ar) / bi
    endif
    rc = 0
  end subroutine Comdiv !Comdiv
  function comabs(ar,ai)          !Real part ,Imaginary part ................. 
    real(double) :: ar,ai
    real(double) :: comabs
    !*====================================================================*
    !*   Input parameters:                                                *
    !*      ar,ai     Real, imaginary parts of  a                         *
    !*   Return value :                                                   *
    !*      Absolute value of a (real)                                    *
    !*====================================================================*
    if (ar == ZERO.and.ai == ZERO) then
       Comabs = ZERO
       return
    endif
    ar = DABS(ar)
    ai = DABS(ai)
    if (ai > ar) then                                  !Switch  ai and ar
       call RSWAP(ai, ar)
    endif
    if (ai == ZERO) then
       Comabs = ar
    else
       Comabs = ar * DSQRT(ONE + ai / ar * ai / ar)
    endif
  end function comabs
  subroutine  hqrvec(n,     & !Dimension of matrix .......
       low,   & !first nonzero row .........
       high,  & !last nonzero row ..........
       h,     & !upper Hessenberg matrix ...
       wr,    & !Real parts of evalues .....
       wi,    & !Imaginary parts of evalues 
       eivec, & !Eigenvectors ..............
       rc  )   !return code ...............
    integer(long),intent(in)   :: n
    integer(long),intent(in)   :: high, low
    real(double), intent(in)   :: wr(0:n),wi(0:n)
    real(double), intent(out)  :: eivec(0:n,0:n)
    real(double)  :: h(0:n,0:n)
    integer(long) :: rc
    integer(long) :: i, j, m, k, na, l
    integer(long) :: code, en
    real(double)  :: p, q, r, s, t, w, x, y, z, ra, sa, vr, vi, norm, temp
    !*====================================================================*
    !*  hqrvec computes the eigenvectors for the eigenvalues found in hqr2*
    !*   Input parameters:                                                *
    !*   ================                                                 *
    !*      n        int n;  ( n > 0 )                                    *
    !*               Dimension of  mat and eivec, number of eigenvalues.  *
    !*      low      int low;                                             *
    !*      high     int high; see  balance                               *
    !*      h        n x n upper Hessenberg matrix                        *
    !*      wr       vector of size n;                                    *
    !*               Real parts of the n eigenvalues.                     *
    !*      wi       vector of size n;                                    *
    !*               Imaginary parts of the n eigenvalues.                *
    !*   Output parameter:                                                *
    !*   ================                                                 *
    !*      eivec    n x n matrix, whose columns are the eigenvectors     *
    !*====================================================================*
    r=ZERO; s=ZERO; z=ZERO; norm=ZERO
    do i = 0, n-1                               !find norm of h
       do j = i, n-1
          norm = norm + DABS(h(i,j))
       enddo
    enddo
    if (norm == ZERO) then
       rc = 1                                    !zero matrix
       return
    endif
    do en = n-1, 0, -1                          !transform back
       p = wr(en)
       q = wi(en)
       na = en - 1
       if (q == ZERO) then
          m = en
          h(en,en) = ONE
          do i = na, 0, -1
             w = h(i,i) - p
             r = h(i,en)
             do j = m, na 
                r = r + h(i,j) * h(j,en)
             enddo
             if (wi(i) < ZERO) then
                z = w
                s = r
             else
                m = i
                if (wi(i) == ZERO) then
                   if (w.ne.ZERO) then 
                      temp = w 
                   else 
                      temp=XMACH_EPS * norm
                   endif
                   h(i,en) = -r/temp            
                else
                   !Solve the linear system:
                   !| w   x |  | h[i][en]   |   | -r |
                   !|       |  |            | = |    |
                   !| y   z |  | h[i+1][en] |   | -s |
                   x = h(i,i+1)
                   y = h(i+1,i)
                   q = (wr(i) - p)**2 + wi(i)**2
                   h(i,en) = (x * s - z * r) / q
                   t = h(i,en)
                   if (DABS(x) > DABS(z)) then 
                      temp = (-r -w * t) / x 
                   else 
                      temp = (-s -y * t) / z
                   endif
                   h(i+1,en) = temp
                endif
             endif !wi[i] < 0
          enddo !i loop
       else if (q < ZERO) then
          m = na
          if (DABS(h(en,na)) > DABS(h(na,en))) then
             h(na,na) = - (h(en,en) - p) / h(en,na)
             h(na,en) = - q / h(en,na)
          else
             call Comdiv(-h(na,en),0.d0, h(na,na)-p, q, h(na,na), h(na,en),code)
          endif
          h(en,na) = ONE
          h(en,en) = ZERO
          do i = na - 1, 0, -1
             w = h(i,i) - p
             ra = h(i,en)
             sa = ZERO
             do j = m, na
                ra = ra + h(i,j) * h(j,na)
                sa = sa + h(i,j) * h(j,en)
             enddo
             if (wi(i) < ZERO) then
                z = w
                r = ra
                s = sa
             else
                m = i
                if (wi(i) == ZERO) then
                   call Comdiv(-ra, -sa, w, q, h(i,na), h(i,en),code)
                else
            !  solve complex linear system:
                   !| w+i*q     x | | h[i][na] + i*h[i][en]  |   | -ra+i*sa |
                   !|             | |                        | = |          |
            !|   y    z+i*q| | h[i+1][na]+i*h[i+1][en]|   | -r+i*s   |
                   x = h(i,i+1)
                   y = h(i+1,i)
                   vr = (wr(i) - p)**2 + wi(i)**2 - q*q
                   vi = TWO * q * (wr(i) - p)
                   if (vr == ZERO.AND.vi == ZERO) then
                      vr = XMACH_EPS * norm * (DABS(w) + DABS(q)  &
                           + DABS(x) + DABS(y) + DABS(z))
                   endif
                   
                   call Comdiv (x*r-z*ra+q*sa,x*s-z*sa-q*ra &
                        ,vr,vi,h(i,na),h(i,en),code)
                   if (DABS(x) > DABS(z) + DABS(q)) then
                      h(i+1,na) = (-ra - w * h(i,na) + q * h(i,en)) / x
                      h(i+1,en) = (-sa - w * h(i,en) - q * h(i,na)) / x
                   else
                      call Comdiv (-r - y * h(i,na), -s - y * h(i,en) &
                           , z, q, h(i+1,na), h(i+1,en),code)
                   endif
                endif !wi[i] = 0
             endif !wi[i] < 0
          enddo !i loop
       endif !else if q < 0
    enddo !en loop
    do i = 0, n-1                        !Eigenvectors for the evalues for
       if (i < low.or.i > high) then      !rows < low  and rows  > high
          do k = i + 1, n-1
             eivec(i,k) = h(i,k)
          enddo
       endif
    enddo
    j = n-1
    do while (j>=low)
       if(j<=high)then
          m =j 
       else 
          j = high
       endif
       if (wi(j) < ZERO) then
          l=j-1
          do i = low, high
             y=ZERO; z=ZERO
             do k = low, m
                y = y + eivec(i,k) * h(k,l)
                z = z + eivec(i,k) * h(k,j)
             enddo
             eivec(i,l) = y
             eivec(i,j) = z
          enddo
       else
          if (wi(j) == ZERO) then
             do i = low, high
                z = ZERO
                do k = low, m
                   z = z + eivec(i,k) * h(k,j)
                enddo
                eivec(i,j) = z
             enddo
          endif
       endif
       j = j - 1
    enddo !j loop
    rc = 0
  end subroutine hqrvec
  subroutine hqr2(n,    &  !Dimension of matrix .........
       low,  &  !first nonzero row ...........
       high, &  !last nonzero row ............
       h,    &  !Hessenberg matrix ...........
       wr,   &  !Real parts of eigenvalues ...
       wi,   &  !Imaginary parts of evalues ..
       eivec,&  !Matrix of eigenvectors ......
       cnt,  &  !Iteration counter ...........
       rc   )    !return code .................              
    integer(long),intent(in)    :: n
    integer(long),intent(in)    :: high, low
    real(double) ,intent(out)   :: h(0:n,0:n)
    real(double), intent(out)   :: wr(0:n),wi(0:n)
    real(double), intent(out)   :: eivec(0:n,0:n)
    integer(long),intent(out)   :: rc
    integer(long),intent(out)   :: cnt(0:n)
    integer(long) :: en
    integer(long) :: i, j, na, iter, l, ll, m, k
    real(double)  :: p, q, r, s, t, w, x, y, z
    !**********************************************************************
    !* hqr2 computes the eigenvalues and (if vec = True) the eigenvectors *
    !* of an  n * n upper Hessenberg matrix.                              *
    !*   Input parameters:                                                *
    !*   ----------------                                                 *
    !*      n        integer;  ( n > 0 )                                  *
    !*               Dimension of  h and eivec,                           *
    !*               length of the real parts vector  wr and of the       *
    !*               imaginary parts vector  wi of the eigenvalues.       *
    !*      low      integer;                                             *
    !*      high     integer;  see balance                                *
    !*      h        n x n matrix;                                        *
    !*               upper Hessenberg matrix as output of Elmhes          *
    !*               (destroyed in the process).                          *
    !*   Output parameters:                                               *
    !*   -----------------                                                *
    !*      eivec    n x n matrix;  (only if vec = 1)                     *
    !*               Matrix, which for vec = 1 contains the               *
    !*               eigenvectors as follows:                             *
    !*               For real eigebvalues the corresponding column        *
    !*               contains the corresponding eigenvactor, while for    *
    !*               complex eigenvalues the corresponding column contains*
    !*               the real part of the eigenvactor with its imaginary  *
    !*               part is stored in the subsequent column of eivec.    *
    !*               The eigenvactor for the complex conjugate eigenvactor*
    !*               is given by the complex conjugate eigenvactor.       *
    !*      wr       vector of size n;                                    *
    !*               Real part of the n eigenvalues.                      *
    !*      wi       vector of size n;                                    *
    !*               Imaginary parts of the eigenvalues                   *
    !*      cnt      Integer vector of size n;                            *
    !*               vector of iterations used for each eigenvalue.       *
    !*               For a complex conjugate eigenvalue pair the second   *
    !*               entry is negative.                                   *
    !**********************************************************************
    p=ZERO; q=ZERO; r=ZERO 
    do i = 0, n-1
       if (i < low.or.i > high) then
          wr(i) = h(i,i)
          wi(i) = ZERO
          cnt(i) = 0
       endif
    enddo
    en = high
    t = ZERO
    do while (en >= low)
       iter = 0
       na = en - 1
       do while(1<2)
          ll=999                          
          do l = en, low+1, -1                      !search for small
             !subdiagonal element
             if(DABS(h(l,l-1))<=XMACH_EPS*(DABS(h(l-1,l-1))+DABS(h(l,l))))then
                ll=l;      !save current index
                goto 10    !exit l loop
             endif
          enddo
10        if(ll.ne.999)then 
             l=ll 
          else 
             l=0          !restore l
          endif
          x = h(en,en)
          if (l == en) then                         !found one evalue
             wr(en) = x + t
             h(en,en) = x + t
             wi(en) = ZERO
             cnt(en) = iter
             en = en - 1
             goto 15      !exit from loop while(True)
          endif
          y = h(na,na)
          w = h(en,na) * h(na,en)
          if (l == na) then                         !found two evalues
             p = (y - x) * 0.5d0
             q = p * p + w
             z = DSQRT(DABS(q))
             x = x + t
             h(en,en) = x + t
             h(na,na) = y + t
             cnt(en) = -iter
             cnt(na) = iter
             if (q >= ZERO) then                     !real eigenvalues
                if (p<ZERO) then 
                   z=p-z 
                else 
                   z=p+z
                endif
                wr(na) = x + z
                wr(en) = x - w / z
                s = w - w / z
                wi(na) = ZERO
                wi(en) = ZERO
                x = h(en,na)
                r = DSQRT (x * x + z * z)
                p = x / r
                q = z / r
                do j = na, n-1
                   z = h(na,j)
                   h(na,j) = q * z + p * h(en,j)
                   h(en,j) = q * h(en,j) - p * z
                enddo
                do i = 0, en
                   z = h(i,na)
                   h(i,na) = q * z + p * h(i,en)
                   h(i,en) = q * h(i,en) - p * z
                enddo
                do i = low, high
                   z = eivec(i,na)
                   eivec(i,na) = q * z + p * eivec(i,en)
                   eivec(i,en) = q * eivec(i,en) - p * z
                enddo
             else                                  !pair of complex
                wr(na) = x + p
                wr(en) = x + p
                wi(na) =   z
                wi(en) = - z
             endif !if q>=ZERO
             en = en - 2
             goto 15                               !exit while(1<2)
          endif !if l = na
          if (iter >= MAXIT) then
             cnt(en) = MAXIT + 1
             rc = en
             write(*,*) ' stop at iter >= MAXIT.'
             return
          endif
          if (iter.ne.0.and.MOD(iter,10) == 0) then
             t = t + x
             do i = low, en 
                h(i,i) = h(i,i) - x
             enddo
             s = DABS(h(en,na)) + DABS(h(na,en-2))
             x = 0.75d0 * s; y = x
             w = -0.4375d0 * s * s
          endif
          iter = iter + 1
          do m = en - 2, l, -1
             z = h(m,m)
             r = x - z
             s = y - z
             p = ( r * s - w ) / h(m+1,m) + h(m,m+1)
             q = h(m + 1,m + 1) - z - r - s
             r = h(m + 2,m + 1)
             s = DABS(p) + DABS(q) + DABS (r)
             p = p / s
             q = q / s
             r = r / s
             if (m == l)  goto 12
             if (DABS(h(m,m-1)) * (DABS(q) + DABS(r)) <= XMACH_EPS * DABS(p) &
                  * (DABS(h(m-1,m-1)) + DABS(z) + DABS(h(m+1,m+1)))) then
                goto 12                !exit m loop
             endif
          enddo
12        do i = m + 2, en 
             h(i,i-2) = ZERO
          enddo
          do i = m + 3, en 
             h(i,i-3) = ZERO
          enddo
          do k = m, na
             if(k.ne.m)then!double QR step, for rows l to en and columns m to en
                p = h(k,k-1)
                q = h(k+1,k-1)
                if (k.ne.na) then 
                   r = h(k+2,k-1) 
                else 
                   r = ZERO
                endif
                x = DABS(p) + DABS(q) + DABS(r)
                if (x == ZERO) goto 30                  !next k
                p = p / x
                q = q / x
                r = r / x
             endif
             s = DSQRT(p * p + q * q + r * r)
             if (p < ZERO) s = -s
             if (k.ne.m) then
                h(k,k-1) = -s * x
             else if (l.ne.m) then
                h(k,k-1) = -h(k,k-1)
             endif
             p = p + s
             x = p / s
             y = q / s
             z = r / s
             q = q / p
             r = r / p
             do j = k, n-1                          !modify rows
                p = h(k,j) + q * h(k+1,j)
                if (k.ne.na) then
                   p = p + r * h(k+2,j)
                   h(k+2,j) = h(k+2,j) - p * z
                endif
                h(k+1,j) = h(k+1,j) - p * y
                h(k,j)   = h(k,j) - p * x
             enddo
             if (k+3 < en) then 
                j=k+3 
             else 
                j=en
             endif
             do i = 0, j                            !modify columns
                p = x * h(i,k) + y * h(i,k+1)
                if (k.ne.na) then
                   p = p + z * h(i,k+2)
                   h(i,k+2) = h(i,k+2) - p * r
                endif
                h(i,k+1) = h(i,k+1) - p * q
                h(i,k)   = h(i,k) - p
             enddo
             do i = low, high
                p = x * eivec(i,k) + y * eivec(i,k+1)
                if (k.ne.na) then
                   p = p + z * eivec(i,k+2)
                   eivec(i,k+2) = eivec(i,k+2) - p * r
                endif
                eivec(i,k+1) = eivec(i,k+1) - p * q
                eivec(i,k)   = eivec(i,k) - p
             enddo
30           continue
          enddo !k loop
       enddo !while(1<2)
15  continue
    enddo !while en >= low                         All evalues found
    !transform evectors back
    call hqrvec (n, low, high, h, wr, wi, eivec,rc)
  end subroutine hqr2
  
  subroutine eigen (matrix, eigvec, eigval)    
    real(double) ,intent(in),dimension(nlevs,nlevs)  :: matrix
    real(double) ,intent(out),dimension(nlevs,nlevs) :: eigvec
    real(double) ,intent(out),dimension(nlevs)       :: eigval   
    real(double)     :: mat(0:nlevs,0:nlevs)  
    real(double)     :: eivec(0:nlevs,0:nlevs)
    real(double)     :: valre(0:nlevs) !real parts of eigenvalues
    real(double)     :: valim(0:nlevs) !imaginary parts of eigenvalues
    integer(long)    :: rc             !return code
    integer(long)    :: cnt(0:nlevs)   !Iteration counter
    integer(long)    :: high, low
    real(double)     :: d(0:nlevs), scale(0:nlevs)
    integer(long)    :: perm(0:nlevs)
    integer(long)    :: i,j,k, check ! counter
    real(double)     :: w,v,norm
    integer(long)    :: n ! nlevels
    n=nlevs
    check=0
    !**********************************************************************
    !* The subroutine eigen  determines all eigenvalues and (if desired)  *
    !* all eigenvectors of a real square  n * n  matrix via the QR method *
    !* in the version of Martin, Parlett, Peters, Reinsch and Wilkinson.  * 
    !*   Litterature:                                                     *
    !*   -----------                                                      *
    !*      1) Peters, Wilkinson: Eigenvectors of real and complex        *
    !*         matrices by LR and QR triangularisations,                  *
    !*         Num. Math. 16, p.184-204, (1970); [PETE70]; contribution   *
    !*         II/15, p. 372 - 395 in [WILK71].                           *
    !*      2) Martin, Wilkinson: Similarity reductions of a general      *
    !*         matrix to Hessenberg form, Num. Math. 12, p. 349-368,(1968)*
    !*         [MART 68]; contribution II,13, p. 339 - 358 in [WILK71].   *
    !*      3) Parlett, Reinsch: Balancing a matrix for calculations of   *
    !*         eigenvalues and eigenvectors, Num. Math. 13, p. 293-304,   *
    !*         (1969); [PARL69]; contribution II/11, p.315 - 326 in       *
    !*         [WILK71].                                                  *
    !*   Input parameters:                                                *
    !*   ----------------                                                 *
    !*      n        integer; ( n > 0 )                                   *
    !*               size of matrix, number of eigenvalues                *
    !*      mat      n x n matrix;                                        *
    !*               input matrix                                         *
    !*   Output parameters:                                               *
    !*   -----------------                                                *
    !*      eivec    n x n matrix;     (only if vec = 1)                  *
    !*               matrix, if  vec = 1  that holds the eigenvectors     *
    !*               thus :                                               *
    !*               If the jth eigenvalue of the matrix is real then the *
    !*               jth column is the corresponding real eigenvector;    *
    !*               if the jth eigenvalue is complex then the jth column *
    !*               of eivec contains the real part of the eigenvector   *
    !*               while its imaginary part is in column j+1.           *
    !*               (the j+1st eigenvector is the complex conjugate      *
    !*               vector.)                                             *
    !*      valre    vector of size n;                                    *
    !*               Real parts of the eigenvalues.                       *
    !*      valim    vector of size n;                                    *
    !*               Imaginary parts of the eigenvalues                   *
    !*      cnt      Integer vector of size n;                            *
    !*               vector containing the number of iterations for each  *
    !*               eigenvalue. (for a complex conjugate pair the second *
    !*               entry is negative).                                  *
    !**********************************************************************
    cnt=0 ; d=0.d0
    mat(0:n-1,0:n-1)=matrix(1:n,1:n)
    !balance mat for nearly
    call balance(n, mat, scale, low, high)      !equal row and column 
    !reduce mat to upper
    call elmhes(n, low, high, mat, perm)        !reduce mat to upper
    !Hessenberg form
    call elmtrans(n, low, high, mat, perm, eivec)
    !QR algorithm for eigenvalues and eigenvectors
    call hqr2(n, low, high, mat, valre, valim, eivec, cnt,rc)  
    !reverse balancing to determine eigenvectors
    call balback(n, low, high, scale, eivec) 
    if (rc.ne.0) stop 'problem in eigen!'
    eigval(1:n)=valre(0:n-1)
    eigvec(1:n,1:n)=eivec(0:n-1,0:n-1)
    !! check
    if(check.eq.1) then
       mat(0:n-1,0:n-1)=matrix(1:n,1:n)
       norm = ZERO; k=0
       do while (k <= n-1)
          if (valim(k) == ZERO) then
             do i = 0, n-1
                w = ZERO
                do j = 0, n-1
                   w = w + mat(i,j) * eivec(j,k)
                enddo
                w = w - valre(k) * eivec(i,k)
                norm = norm + DABS(w)
             enddo
          else
             do i = 0, n-1
                w = ZERO
                do j = 0, n-1
                   w = w + mat(i,j) * eivec(j,k)
                enddo
                w = w - valre(k) * eivec(i,k) - valim(k) * eivec(i,k+1)
                v = ZERO
                do j = 0, n-1
                   v = v + mat(i,j) * eivec(j,k+1)
                enddo
                v = v - valre(k) * eivec(i,k+1) + valim(k) * eivec(i,k)
                norm = norm + 2.d0 * DSQRT(v*v + w*w)
             enddo
             k=k+1
          endif
          k=k+1 
        enddo !while k<=n-1
        print*, 'norm:', norm
        stop
     endif
  end subroutine eigen

  function outerprod(a,b)
    real(double), dimension(:), intent(IN) :: a,b
    real(double), dimension(size(a),size(b)) :: outerprod
    outerprod = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  end function outerprod
  subroutine swap(a,b)
    real(double), dimension(:), intent(INOUT) :: a,b
    real(double), dimension(size(a)) :: dum
    dum=a
    a=b
    b=dum
  end subroutine swap
  subroutine ludcmp(a,indx,d)
    real(double), dimension(:,:),intent(INOUT):: a
    integer(long),dimension(:),  intent(OUT)  :: indx
    real(double),                intent(OUT)  :: d
    real(double), dimension(size(a,1))        :: vv
    integer(long),dimension(1)                :: imaxloc
    integer(long) :: j,n,imax
    n=size(indx)
    d=1.0
    vv=maxval(abs(a),dim=2)
    if(any(vv.eq.0.))stop 'singular matrix in ludcmp'
    vv=1.d0/vv
    do j=1,n
       imaxloc=maxloc(vv(j:n)*abs(a(j:n,j)))
       imax=(j-1)+imaxloc(1)
       if (j /= imax) then
          call swap(a(imax,:),a(j,:))
          d=-d
          vv(imax)=vv(j)
       endif
       indx(j)=imax
       if (a(j,j) == 0.0) a(j,j)=1.0d-20
       a(j+1:n,j)=a(j+1:n,j)/a(j,j)
       a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
    enddo
  end subroutine ludcmp
  subroutine lubksb(a,indx,b)
    real(double), dimension(:,:),intent(IN)   :: a
    integer(long),dimension(:),  intent(IN)   :: indx
    real(double), dimension(:),  intent(INOUT):: b
    integer(long) :: i,n,ii,ll
    real(double)  :: summ
    n=size(indx)
    ii=0
    do i=1,n
       ll=indx(i)
       summ=b(ll)
       b(ll)=b(i)
       if (ii /= 0) then
          summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
       else if (summ /= 0.0) then
          ii=i
       endif
       b(i)=summ
    enddo
    do i=n,1,-1
       b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
    enddo
  end subroutine lubksb
  subroutine matinv(a, b)
    !! - Matrix inversion with LU-decomposition
    !====================================================
    real(double), dimension(:,:), intent(IN)             :: a
    real(double), dimension(:,:), intent(OUT)            :: b
    real(double), dimension(size(a,dim=1),size(a,dim=2)) :: ah, y
    integer(long)                                        :: i, N
    integer(long), dimension(size(a,dim=1))              :: indx
    real(double)                                         :: d
    N = size(a,dim=1)
    if (N /= size(a,dim=2)) stop 'SUB matinv: ludcmp matrix must be square!'
    ah = a
    y  = 0.
    do i = 1, N
       y(i,i) = 1.d0
    enddo
    call ludcmp(ah,indx,d)
    do i = 1, N
       call lubksb(ah,indx,y(:,i))
    enddo
    b = y
  end subroutine matinv

  subroutine table_interp(coef_matrix,deb,dti,ebi,tii,eb,ti,interp_out)
    !! bilinear interpolation of the effective rate coefficient tables !!
    real(double),dimension(:,:,:,:), intent(in) :: coef_matrix ! (m,n,eb,ti)
    real(double),                intent(in)     :: deb, dti ! eb and ti grid
    integer,                    intent(in)      :: ebi, tii    ! lower indices
    real(double),                intent(in)     :: eb, ti      ! desired values
    real(double),dimension(nlevs+1,nlevs),intent(out) :: interp_out  ! output 
    real(double), dimension(nlevs+1,nlevs)      :: c00, c10, c01, c11, c0, c1
    real(double)    :: eb0, ti0
    eb0=(ebi-1)*deb
    ti0=(tii-1)*dti
    c00 = coef_matrix(:,:,ebi  ,tii  )
    c10 = coef_matrix(:,:,ebi+1,tii  )
    c01 = coef_matrix(:,:,ebi  ,tii+1)
    c11 = coef_matrix(:,:,ebi+1,tii+1)
    !linear interpolation between C00 and C10 to find C0, C01 and C11 to find C1
    c0  = c00 + (c10 - c00) * (eb - eb0) / deb
    c1  = c01 + (c11 - c01) * (eb - eb0) / deb
    interp_out=( c0  + (c1  - c0) * (ti - ti0) / dti)
  end subroutine table_interp

  subroutine colrad(ac,vn,dt,states,photons,neut_type,nlaunch)
    !colrad solves the collisional-radiative balance equations
    real(double) , dimension(:),  intent(in)   :: vn  !!velocitiy (cm/s)
    real(double)               ,  intent(in)   :: dt  !!time interval in cell
    integer(long)              ,  intent(in)   :: neut_type!!type of neutral
    real(double)               ,  intent(in)   :: nlaunch !! nr of markers
    integer(long), dimension(:),  intent(in)   :: ac   !!actual cell
    real(double) , dimension(:),  intent(inout):: states  !!density of states
    real(double),                 intent(out)  :: photons !!emitted photons 
    !! ---- to determine rate coefficients ---- !   
    real(double), dimension(nlevt+1,nlevt)     :: qp !! Proton rate coefficants
    real(double), dimension(nlevt+1,nlevt)     :: qi !! Impurity rate coefficant
    real(double), dimension(nlevt+1,nlevt)     :: qe !! Electron rate coefficant
    real(double), dimension(nlevs,nlevs)       :: matrix  !! Matrix
    real(double), dimension(3)  :: vrot           !! Rotation velocity of plasma
    real(double)                :: vnet_square    !! netto velocity of neutrals 
    real(double)                :: ti,te          !! Ion/electron temperature
    real(double)                :: denp,dene,deni !! P/impurity/electron density
    real(double)                :: eb             !! Energy of the fast neutral
    integer(long)               :: ebi, tii,tei   !! bin postions in arrays
    !! ---- Solution of differential equation  ---- ! 
    real(double),   dimension(nlevs,nlevs)  :: eigvec, eigvec_inv
    real(double),   dimension(nlevs)        :: eigval, coef
    real(double),   dimension(nlevs)        :: exp_eigval_dt 
    real(double),    dimension(nlevs)        :: dens !! Density of neutrals 
    integer(long)                           :: n !! counter 
    photons=0.d0
    !! --------------- Check if inputs are valid for colrad -------------- !!
    if(sum(states).lt.colrad_threshold .and. inputs%npa.eq.0)then
       if(neut_type.eq.2)print*, 'threshold!',ac
       return
    endif
    denp=cell(ac(1),ac(2),ac(3))%plasma%denp
    dene=cell(ac(1),ac(2),ac(3))%plasma%dene
    deni=cell(ac(1),ac(2),ac(3))%plasma%deni  
    ti=cell(ac(1),ac(2),ac(3))%plasma%ti
    te=cell(ac(1),ac(2),ac(3))%plasma%te
    vrot=cell(ac(1),ac(2),ac(3))%plasma%vrot

    !! IF the temperature or density is too low, stop simulation
    !! => particles in the SOL
    !! (stopped by return of photons=0.!)
    if(ti.le.0.05.or.te.le.0.05.or.denp.lt.1.d12)then
       if(  neut_type.eq.nbif_type.or. &  !! Store density for NBI simulation!
            neut_type.eq.nbih_type.or. &
            neut_type.eq.nbit_type)then
          dens(:)=states*dt/nlaunch!![neutrals/(cm^3)]!!
          cell(ac(1),ac(2),ac(3))%neut_dens(neut_type,:)= & 
               cell(ac(1),ac(2),ac(3))%neut_dens(neut_type,:)+real(dens(:),float)
       endif
       !! do not stop simulation in the case of NPA until the particle 
       !! has reached the detector
       if(inputs%npa.eq.1)then !! NPA simulation !!
          dens(:)=states*dt/nlaunch!![neutrals/(cm^3)]!!
          cell(ac(1),ac(2),ac(3))%neut_dens(neut_type,:)= & 
               cell(ac(1),ac(2),ac(3))%neut_dens(neut_type,:) &
               +real(dens(:),float)
          photons=1.d0
          if(npa%at_detector)then
             npa%v(npa%counter,:)=vn(:)
             npa%wght(npa%counter)=sum(dens)*grid%dv/dt !! [neutrals]
             npa%at_detector=.false.
             photons=0.d0
          endif
       endif
       return
    endif
    !! NORMAL START OF COLRAD
    !! ------------------ Get Matrix with rate coeffients ----------------- !!
    !! - DEUTERIUM ions (Impact excitation/ionization,charge exchange)
    vnet_square=dot_product(vn-vrot,vn-vrot)          ![cm/s]
    eb=v_to_E*inputs%ab*vnet_square ![kev]

    !! DEUTERIUM (Impact excitation/ionization,charge exchange)
    ebi= floor(eb/atomic%d_eb_qp)+1   
    if(ebi.ge.atomic%nr_eb_qp)stop 'Eb out of range of qptable!'
    tii= floor(ti/atomic%d_ti_qp)+1
    if(tii.ge.atomic%nr_ti_qp)stop 'Ti out of range of qptable!'
    call table_interp(atomic%qp(:,:,:,:),atomic%d_eb_qp,atomic%d_ti_qp &
         ,ebi,tii,eb,ti,qp)
    qp=qp*denp ![1/s]  

    !! IMPURITIES
    ebi= floor(eb/atomic%d_eb_qi)+1   
    if(ebi.ge.atomic%nr_eb_qi)stop 'Eb out of range of qitable!'
    tii= floor(ti/atomic%d_ti_qi)+1
    if(tii.ge.atomic%nr_ti_qi)stop 'Ti out of range of qitable!'
    call table_interp(atomic%qi(:,:,:,:),atomic%d_eb_qi,atomic%d_ti_qi &
         ,ebi,tii,eb,ti,qi)
    qi=qi*deni ![1/s]  

    !! ELECTRONS
    ebi= floor(eb/atomic%d_eb_qe)+1   
    if(ebi.ge.atomic%nr_eb_qe)stop 'Eb out of range of qetable!'
    tei= floor(te/atomic%d_te_qe)+1
    if(tei.ge.atomic%nr_te_qe)stop 'Te out of range of qetable!'
    call table_interp(atomic%qe(:,:,:,:),atomic%d_eb_qe,atomic%d_te_qe &
         ,ebi,tei,eb,te,qe)
    qe=qe*dene ![1/s]  

    !! - Write off-diagnonal elements (populating transitions) - !!
    matrix=atomic%einstein(1:nlevs,1:nlevs)  &     
         +              qp(1:nlevs,1:nlevs)  &
         +              qi(1:nlevs,1:nlevs)  &
         +              qe(1:nlevs,1:nlevs)
    !! - Write diagonal elements (depopulating transitions) - !!
    do n=1,nlevs         
       matrix(n,n)=&
            - sum(atomic%einstein(:,n)) &
            - sum(qp(:,n)) &
            - sum(qi(:,n)) &
            - sum(qe(:,n))
    enddo
    !! WWH  Add Hutchinson's truncation correction
!    matrix(nlevs,nlevs)=10.*matrix(nlevs,nlevs)


    call eigen(matrix, eigvec, eigval)
    call matinv(eigvec, eigvec_inv)
    coef = matmul(eigvec_inv, states)!coeffs determined from states at t=0
    exp_eigval_dt = exp(eigval*dt)   ! to improve speed (used twice)
    states(:) = matmul(eigvec, coef * exp_eigval_dt)  ![neutrals/cm^3/s]!
    dens(:)   = matmul(eigvec,coef*(exp_eigval_dt-1.d0)/eigval)/nlaunch
    !$OMP CRITICAL(col_rad)
    cell(ac(1),ac(2),ac(3))%neut_dens(neut_type,:)= & 
         cell(ac(1),ac(2),ac(3))%neut_dens(neut_type,:) &
         + real(dens(:),float)![neutrals/cm^3]!
    !$OMP END CRITICAL(col_rad)
    !! -------------------- determine photon flux ------------------------ !!
    photons=dens(3)*atomic%einstein(2,3) !! - [Ph/(s*cm^3)] - !!
  end subroutine colrad
  
  !***************************************************************************
  !-----------spectrum--------------------------------------------------------
  !***************************************************************************
  subroutine spectrum(vi,ac,pos,photons,neut_type,wavout,intout)
    !!spectrum.pro computes the wavelengths of emitted photons
    real(double), dimension(:), intent(in) :: vi!!velocitiy of neutral [cm/s]
    integer(long),dimension(:), intent(in) :: ac  !!actual cell
    integer(long)             , intent(in) :: neut_type!!type of neutral
    real(double)              , intent(in) :: photons !! photons from colrad
    real(double), dimension(3), intent(in) :: pos    !! mean position in cell
    real(double), dimension(9), intent(out), optional :: intout!!intensity
    real(double), dimension(9), intent(out), optional :: wavout !!wavelength 
    real(double)               :: lambda_Doppler , cos_los_Efield, E
    real(double), dimension(9) ::intens!!intensity vector
    real(double), dimension(9) ::wavel !!wavelength vector[A
    real(double), dimension(3) :: vp  !!unit vector of sight line
    real(double), dimension(3) :: vn  ! vi in m/s
    real(double), dimension(3) :: efield  !E-field (static + vxB)
    real(double), dimension(3) :: bfield  !B-field
    integer(long)              :: i, ichan, bin  !counter, wavelengths bins
    if(inputs%npa.eq.1)return
    loop_over_channels: do ichan=1,spec%nchan
       if(cell(ac(1),ac(2),ac(3))%los_wght(ichan).le.0.)cycle loop_over_channels
       !! vector directing towards the optical head
       vp(1)=pos(1)-spec%xyzhead(ichan,1) 
       vp(2)=pos(2)-spec%xyzhead(ichan,2) 
       vp(3)=pos(3)-spec%xyzhead(ichan,3) 
       vp=vp/sqrt(dot_product(vp,vp))
       ! Calculate Doppler shift
       vn=vi*0.01d0 ! [m/s]
       lambda_Doppler = lambda0*(1.d0 + dot_product(vn,vp)/c0)
       !! Calculate Stark Splitting
       ! Calcualate E-field
       bfield(:) = cell(ac(1),ac(2),ac(3))%plasma%B
       efield(:) = cell(ac(1),ac(2),ac(3))%plasma%E
       efield(1) = efield(1) +  vn(2)*bfield(3) - vn(3)*bfield(2)
       efield(2) = efield(2) - (vn(1)*bfield(3) - vn(3)*bfield(1))
       efield(3) = efield(3) +  vn(1)*bfield(2) - vn(2)*bfield(1)
       E=sqrt(dot_product(efield,efield))
       !Stark Splitting
       wavel =  lambda_Doppler + E * stark_wavel ![A]
       !Intensities of stark components
       if (E .eq. 0.d0) then 
          cos_los_Efield = 0.d0 
       else 
          cos_los_Efield = dot_product(vp,efield) / E
       endif
       intens(1:3)=stark_intens(1:3)*(1.d0-cos_los_Efield**2.d0)!Pi(parallel)
       intens(4:6)=stark_intens(4:6)*(1.d0+cos_los_Efield**2.d0)!Sigma(perp)
       intens(7:9)=stark_intens(7:9)*(1.d0-cos_los_Efield**2.d0)!Pi(parallel)
       !! --- E.g. mirrors may change the pi to sigma intensity ratio  --- !!
       intens(4:6) = intens(4:6)*inputs%sigma_pi_ratio
       !! --- normalize and multiply with photon density from colrad --- !!
       intens      = intens/sum(intens)*photons 
       if(present(wavout))then
          wavout=wavel
          intout=intens
          return
       endif
       !! ---------------------- Store spectra ---------------------- !!
       do i=1,9                                                       
          bin=int(((wavel(i)-spec%lambdamin)/spec%dlambda)+.5)
          if (bin.lt.1)              bin = 1                        
          if (bin.gt.spec%nlambda) bin = spec%nlambda  
          !$OMP CRITICAL(spec_trum)   
          spec%spectra(bin,ichan,neut_type)= &
               spec%spectra(bin,ichan,neut_type) &
               + real(intens(i),float)*cell(ac(1),ac(2),ac(3))%los_wght(ichan)
          !$OMP END CRITICAL(spec_trum)
       enddo
    enddo loop_over_channels
  end subroutine spectrum
  
  
  !*****************************************************************************
  !------------track------------------------------------------------------------
  !*****************************************************************************
  subroutine track(vin, rin, tcell, icell,pos, ncell)
    !!track computes the path of a neutral through a the FIDAcode grid 
    real(double), dimension(:)  , intent(in)   :: rin  ! initial position
    real(double), dimension(:)  , intent(in)   :: vin  ! velocitiy
    integer(long)               , intent(out)  :: ncell! number of cells
    real(double), dimension(:)  , intent(out)  :: tcell! time per cell
    integer(long),dimension(:,:), intent(out)  :: icell! cell indices
    real(double), dimension(:,:), intent(out)  :: pos  ! mean position in cell
    integer(long)              :: cc    !!step number along the track
    integer(long),dimension(3) :: p,l    !!indices of the cells
    real(double), dimension(3) :: dt_arr !!time to cell boundary
    real(double)               :: dt     !!min time to cell boundary
    real(double), dimension(3) :: vn !!velocitiy that can be changed
    real(double), dimension(3) :: ri !!position of ray  
    tcell=0.d0 ; icell=0 ;  ; pos=0.d0 ; ncell=0
    vn(:)=vin(:) ;  ri(:)=rin(:)
    !! define actual cell
    minpos=minloc(abs(rin(1)-grid%xxc))
    p(1)=minpos(1)
    minpos=minloc(abs(rin(2)-grid%yyc))
    p(2)=minpos(1)  
    minpos=minloc(abs(rin(3)-grid%zzc))
    p(3)=minpos(1) 
    !! Fudge zero velocity components to avoid overflow error
    where (vn(:).eq.0) vn(:) = 0.001   
    !! Start tracking routine
    icell(:,1)=p(:)
    cc=1 
    !!loop along track of neutral
    tracking: do while(cc.lt.(grid%ntrack))
       l(:)=p(:)
       where(vn(:).gt.0.d0) l(:)=p(:)+1 
       if ( l(1).gt.grid%nx.or.& 
            l(2).gt.grid%ny.or.&
            l(3).gt.grid%nz) exit tracking  
       !time needed to go to next cell
       dt_arr(1)=(grid%xx(l(1))-ri(1))/vn(1)
       dt_arr(2)=(grid%yy(l(2))-ri(2))/vn(2)
       dt_arr(3)=(grid%zz(l(3))-ri(3))/vn(3)
       minpos=minloc(dt_arr) 
       dt=dt_arr(minpos(1))
       pos(:,cc) = ri(:) + vn(:)*dt*0.5  !! mean postion in cell
       ri(:)     = ri(:) + vn(:)*dt
       if (vn(minpos(1)).gt.0.d0)  then 
          p(minpos(1))=p(minpos(1))+1
       else
          p(minpos(1))=p(minpos(1))-1
       endif
       if (any(p.le.0))exit tracking    
       tcell(cc)=dt
       icell(:,cc+1)=p(:)
       cc=cc+1
    enddo tracking
    ncell=cc-1
  end subroutine track
  
  !*****************************************************************************
  !-----------ndmc (NBI)--------------------------------------------------------
  !*****************************************************************************
  subroutine ndmc
    integer(long)                          :: indmc     !! counter for markers
    integer(long)                          :: efrac     !! counter energy frac
    integer(long)                          :: type      !! full half third En
    real(double)                           :: nlaunch   !! nr. of markers
    real(double)                           :: nneutrals !! # NBI particles 
    real(double), dimension(3)             :: vnbi      !! velocities(full..)
    real(double), dimension(3)             :: rnbi      !! initial position
    !!Tracking routine output
    integer(long)                          :: jj     !! counter for track
    integer(long)                          :: ncell  !! number of cells
    real(double), dimension(  grid%ntrack) :: tcell  !! time per cell
    integer(long),dimension(3,grid%ntrack) :: icell  !! index of cells
    real(double), dimension(3,grid%ntrack) :: pos    !! mean position in cell
    integer(long),dimension(3)             :: ac     !! actual cell 
    !!collisional radiative model and spectrum calculation
    real(double), dimension(nlevs)         :: states
    real(double)                           :: photons
    print*,'    # of markers: ',inputs%nr_ndmc
    !! ------------- calculate nr. of injected neutrals ---------------- !!
    !! # of injected neutrals = NBI power/energy_per_particle
    nneutrals=1.d6*nbi%pinj/ (1.d3*nbi%einj*e0 &
         *( nbi%species_mix(1)      &
         +  nbi%species_mix(2)/2.d0 &
            +  nbi%species_mix(3)/3.d0 ) )
    !! ------------------ loop over the markers ------------------------ !!
    nlaunch=dble(inputs%nr_ndmc)
    !$OMP PARALLEL DO private(indmc,efrac,vnbi,rnbi,tcell,icell,pos,ncell,states,ac,photons,type,jj)
    loop_over_markers: do indmc=1,inputs%nr_ndmc
       !if(modulo(indmc,10000).eq.0) print*,indmc
       energy_fractions: do efrac=1,3
          type=efrac+1  !! type is for colrad and spectrum
          !! (type = 2: full energy, =3: half energy, =4: third energy
          call mc_nbi(vnbi(:),efrac,rnbi(:))
          if(rnbi(1).eq.-1)cycle loop_over_markers
          call track(vnbi,rnbi,tcell,icell,pos,ncell)
          if(ncell.eq.0) cycle loop_over_markers
          !! --------- solve collisional radiative model along track ----- !!
          states=0.d0
          states(1)=nneutrals*nbi%species_mix(efrac)/grid%dv 
          loop_along_track: do jj=1,ncell
             ac=icell(:,jj)
             call colrad(ac,vnbi,tcell(jj),states,photons,type,nlaunch)
             if(photons.gt.0.d0) then
                if(inputs%nospec.eq.0) &
                call spectrum(vnbi(:),ac,pos(:,jj),photons,type)
             endif
          enddo loop_along_track
       enddo energy_fractions
    enddo loop_over_markers
    !$OMP END PARALLEL DO
  end subroutine ndmc
  

  !*****************************************************************************
  !-----------Bremsstrahlung ---------------------------------------------------
  !*****************************************************************************
  subroutine bremsstrahlung
    integer(long)     :: i,j,k,ichan   !! indices of cells
    real(double)      :: ne,zeff,te,gaunt
    real(double), dimension(:)  , allocatable :: lambda_arr,brems
    !! ------------------------ calculate wavelength array ------------------ !!
    print*, 'calculate the bremsstrahung!'
    allocate(lambda_arr(spec%nlambda),brems(spec%nlambda))
    do i=1,spec%nlambda 
       lambda_arr(i)=(i-0.5)*spec%dlambda+spec%lambdamin ! [A]
    enddo 
    loop_along_x: do i = 1, grid%Nx
       loop_along_y:    do j = 1,grid%Ny
          loop_along_z:     do k = 1, grid%Nz
             loop_over_channels: do ichan=1,spec%nchan
                if(cell(i,j,k)%los_wght(ichan).le.0.)cycle loop_over_channels
                ne=cell(i,j,k)%plasma%dene     ![cm^3]
                zeff=cell(i,j,k)%plasma%zeff   !       
                te=cell(i,j,k)%plasma%te*1000. ! [eV]
                if(te .le. 0.)cycle loop_over_channels
                gaunt=5.542-(3.108-log(te/1000.))*(0.6905-0.1323/zeff)
                brems(:)=7.57d-9*gaunt*ne**2*zeff/(lambda_arr(:)*sqrt(te)) &
                     *exp(-h_planck*c0/(lambda_arr(:)*te))
                spec%spectra(:,ichan,brems_type) =  &
                     spec%spectra(:,ichan,brems_type)  &
                     +brems(:)*10. &        !!from [1/A] to 1/[nm]
                     *cell(i,j,k)%los_wght(ichan)*1.e-2 !! integration
             enddo loop_over_channels
          enddo loop_along_z
       enddo loop_along_y
    enddo loop_along_x
    deallocate(lambda_arr,brems)
  end subroutine bremsstrahlung
  

  !*****************************************************************************
  !-------------- Direct charge exchange calculation---------------------------
  !*****************************************************************************
  subroutine dcx
    integer(long)                          :: i,j,k   !! indices of cells
    real(double), dimension(3)             :: randomu   
    integer(long)                          :: idcx    !! counter
    real(double), dimension(3)             :: ri      !! start position
    real(double), dimension(3)             :: vhalo   !! velocity bulk plasma 
    integer(long),dimension(3)             :: ac      !! actual cell 
    !! Determination of the CX probability
    real(double), dimension(nlevs)         :: denn    !!  neutral dens (n=1-4)
    real(double), dimension(nlevs)         :: prob    !!  Prob. for CX 
    real(double), dimension(3)             :: vnbi    !! Velocity of NBIneutrals
    real(double), dimension(nlevs)         :: rates   !! Rate coefficiants forCX
    !! Collisiional radiative model along track
    real(double), dimension(nlevs)         :: states  !! Density of n-states
    integer(long)                          :: ncell
    real(double), dimension(  grid%ntrack) :: tcell   !! time per cell
    integer(long),dimension(3,grid%ntrack) :: icell   !! index of cells
    real(double), dimension(3,grid%ntrack) :: pos     !! mean position in cell
    integer(long)                          :: jj      !! counter along track
    real(double)                           :: photons !! photon flux
    real(double), dimension(grid%nx,grid%ny,grid%nz)::papprox !! approx. density
    real(double)                           :: papprox_tot
    real(double)                           :: nlaunch !! Number of markers
    real(double)                          :: launch  !! float number of markers
    real(double)                          :: remainder !! used in nlaunch
    real(double),dimension(1)             :: randomu1 !! used in nlaunch
    papprox=0.d0
    papprox_tot=0.d0
    !! ------------- calculate papprox needed for guess of nlaunch --------!!
    do i = 1, grid%Nx
       do j = 1, grid%Ny
          do k = 1, grid%Nz
             papprox(i,j,k)=dble(sum(cell(i,j,k)%neut_dens(nbif_type,:))  &
                  +              sum(cell(i,j,k)%neut_dens(nbih_type,:))  &
                  +              sum(cell(i,j,k)%neut_dens(nbit_type,:))) &
                  *(cell(i,j,k)%plasma%denp-cell(i,j,k)%plasma%denf)
             if(cell(i,j,k)%rho.lt.1.1)papprox_tot=papprox_tot+papprox(i,j,k)
          enddo
       enddo
    enddo
    ! Loop through all of the cells
    print*,'    # of markers: ',inputs%nr_dcx
    !$OMP PARALLEL DO private(i,j,k,idcx,randomu,ac,vhalo,ri,photons,rates,prob, &
    !$OMP& jj,states,vnbi,tcell,icell,pos,ncell,nlaunch,denn)
    loop_along_x: do i = 1, grid%Nx
       loop_along_y:    do j = 1,grid%Ny
          loop_along_z:     do k = 1, grid%Nz
             !if(papprox(i,j,k).le.0.d0) cycle loop_along_z 
             !! -----Decide how many markers to follow in this cell--------!!
             !! WWH --- Treat fractions more accurately
             call randu(randomu1)
             launch=inputs%nr_dcx*papprox(i,j,k)/papprox_tot
             nlaunch=int(launch)
             remainder=launch-nlaunch
             if (randomu1(1) .le. remainder) then
               nlaunch=nlaunch+1
             endif
             if (nlaunch .lt. 1) cycle loop_along_z
             !nlaunch = int(inputs%nr_dcx * papprox(i,j,k)/papprox_tot + 1.)
             !! ------------- loop over the markers ---------------------- !!
             loop_over_dcx: do idcx=1,int(nlaunch)
                ac=(/i,j,k/)
                !! ---------------- calculate ri,vi and track -------------!!   
                call mc_halo(ac,vhalo(:))
                call randu(randomu)  
                ri(1)=grid%xx(ac(1))+ grid%dr(1)*randomu(1)
                ri(2)=grid%yy(ac(2))+ grid%dr(2)*randomu(2) 
                ri(3)=grid%zz(ac(3))+ grid%dr(3)*randomu(3)  
                call track(vhalo(:), ri(:), tcell, icell,pos, ncell)
                if(ncell.eq.0) cycle loop_over_dcx 
                !! ---------------- calculate CX probability ------------- !!
                ac=icell(:,1) 
                prob=0.d0
                vnbi=ri(:)-nbi%xyz_pos(:)
                vnbi=vnbi/sqrt(dot_product(vnbi,vnbi))*nbi%vinj
                ! CX with full energetic NBI neutrals ------ !!
                denn(:)=dble(cell(ac(1),ac(2),ac(3))%neut_dens(nbif_type,:))
                call neut_rates(denn,vhalo,vnbi,rates)
                prob=prob + rates
                ! CX with half energetic NBI neutrals ------ !!
                denn(:)=dble(cell(ac(1),ac(2),ac(3))%neut_dens(nbih_type,:))
                call neut_rates(denn,vhalo,vnbi/sqrt(2.d0),rates)
                prob=prob + rates
                ! CX with third energetic NBI neutrals ------ !!
                denn(:)=dble(cell(ac(1),ac(2),ac(3))%neut_dens(nbit_type,:))
                call neut_rates(denn,vhalo,vnbi/sqrt(3.d0),rates)
                prob=prob + rates
                if(sum(prob).le.0.)cycle loop_over_dcx
                !! --------- solve collisional radiative model along track-!!
                states=prob*(cell(i,j,k)%plasma%denp-cell(i,j,k)%plasma%denf)
                loop_along_track: do jj=1,ncell
                   ac=icell(:,jj)
                   call colrad(ac,vhalo(:),tcell(jj),states &
                        ,photons,halo_type,nlaunch)
                   if(photons.le.0.d0)cycle loop_over_dcx 
                   if(inputs%nospec.eq.0) call spectrum(vhalo(:),ac(:),pos(:,jj),photons,halo_type)
                enddo loop_along_track
             enddo loop_over_dcx
          enddo loop_along_z
       enddo loop_along_y
    enddo loop_along_x
    !$OMP END PARALLEL DO
  end subroutine dcx
  
  !*****************************************************************************
  !-------------------------- halo -------------------------------------------
  !*****************************************************************************
  subroutine halo
    integer(long)                          :: i,j,k !indices of cells   
    integer(long)                          :: ihalo !! counter
    real(double), dimension(3)             :: randomu   
    real(double), dimension(3)             :: ri    !! start position
    real(double), dimension(3)             :: vihalo!! velocity bulk plasma ion
    integer(long),dimension(3)             :: ac    !! actual cell 
    !! Determination of the CX probability
    real(double), dimension(nlevs)         :: denn    !!  neutral dens (n=1-4)
    real(double), dimension(nlevs)         :: prob    !!  Prob. for CX 
    real(double), dimension(nlevs)         :: rates   !! Rate coefficiants forC
    real(double), dimension(3)             :: vnhalo  !! v of halo neutral
    integer(long)                          :: in      !! index over halo neutral
    !! Collisiional radiative model along track
    real(double), dimension(nlevs)         :: states  ! Density of n-states
    integer(long)                          :: ncell
    real(double), dimension(  grid%ntrack) :: tcell  !! time per cell
    integer(long),dimension(3,grid%ntrack) :: icell    !! index of cells
    real(double), dimension(3,grid%ntrack) :: pos    !! mean position in cell
    integer(long)                          :: jj       !! counter along track
    real(double)                           :: photons  !! photon flux
    real(double), dimension(grid%nx,grid%ny,grid%nz)::papprox !! approx. density
    real(double)                           :: papprox_tot 
    !! Halo iteration
    integer(long)                           :: hh !! counters
    real(double)                             :: dcx_dens, halo_iteration_dens
    real(double)                            :: nlaunch !! Number of markers
    real(double)                          :: launch  !! float number of markers
    real(double)                          :: remainder !! used in nlaunch
    real(double),dimension(1)             :: randomu1 !! used in nlaunch
    dcx_dens=0.d0
    do i=1,grid%Nx 
       do j=1,grid%Ny 
          do k=1,grid%Nz 
             if(sum(cell(i,j,k)%neut_dens(halo_type,:)).gt.0.)then
                dcx_dens=dcx_dens+sum(dble(cell(i,j,k)%neut_dens(halo_type,:)))
                cell(i,j,k)%neut_dens(s1type,:) =  &
                     cell(i,j,k)%neut_dens(halo_type,:)
             endif
          enddo
       enddo
    enddo
    if(dcx_dens.eq.0)stop 'the denisty of DCX-neutrals is too small!'
    iterations: do hh=1,20

       !! ------------- calculate papprox needed for guess of nlaunch --------!!
       papprox=0.d0
       papprox_tot=0.d0
       do i = 1, grid%Nx
          do j = 1, grid%Ny
             do k = 1, grid%Nz
                papprox(i,j,k)=dble(sum(cell(i,j,k)%neut_dens(s1type,:))) &
                     *(cell(i,j,k)%plasma%denp-cell(i,j,k)%plasma%denf)
                if(cell(i,j,k)%rho.lt.1.1)papprox_tot=papprox_tot+papprox(i,j,k)
             enddo
          enddo
       enddo

  
       print*, '    # of markers: ' ,inputs%nr_halo
       !$OMP PARALLEL DO private(i,j,k,nlaunch,ihalo,ac,vihalo,randomu,ri,tcell,icell, &
       !$OMP& pos,ncell,prob,denn,in,vnhalo,rates,states,jj,photons)
       loop_along_x: do i = 1, grid%Nx
          loop_along_y:    do j = 1, grid%Ny
             loop_along_z:     do k = 1, grid%Nz
                ! if(papprox(i,j,k)*inputs%nr_halo.le.1.d23) cycle loop_along_z 
                if(papprox(i,j,k).le.0.d0) cycle loop_along_z
                !! ------Decide how many markers to follow in this cell-------!!
               !! WWH --- Treat fractions more accurately
               call randu(randomu1)
               launch=inputs%nr_halo*papprox(i,j,k)/papprox_tot
               nlaunch=int(launch)
               remainder=launch-nlaunch
               if (randomu1(1) .le. remainder) then
                 nlaunch=nlaunch+1
               endif
               if (nlaunch .lt. 1) cycle loop_along_z
                !nlaunch= int(inputs%nr_halo * papprox(i,j,k)/papprox_tot + 1.)
                !! ------------- loop over the markers ---------------------- !!
                loop_over_halos: do ihalo=1,int(nlaunch)
                   ac=(/i,j,k/)
                   !! ---------------- calculate ri,vhalo and track ----------!!
                   call mc_halo( ac, vihalo(:))
                   call randu(randomu)  
                   ri(1)=grid%xx(ac(1))+ grid%dr(1)*randomu(1)
                   ri(2)=grid%yy(ac(2))+ grid%dr(2)*randomu(2) 
                   ri(3)=grid%zz(ac(3))+ grid%dr(3)*randomu(3)  
                   call track(vihalo(:), ri(:),tcell,icell,pos,ncell)
                   if(ncell.eq.0)cycle loop_over_halos
                   !! ---------------- calculate CX probability --------------!!
                   ac=icell(:,1) 
                   prob=0.d0
                   !CX with HALO neutrals
                   denn(:)=dble(cell(ac(1),ac(2),ac(3))%neut_dens(s1type,:))
                   do in=1,int(nr_halo_neutrate)
                      call mc_halo( ac(:), vnhalo(:))
                      call neut_rates(denn,vihalo,vnhalo,rates)
                      prob=prob+rates/nr_halo_neutrate
                   enddo
                   if(sum(prob).le.0.)cycle loop_over_halos
                   !! --------- solve collisional radiative model along track-!!
                   states=prob*(cell(i,j,k)%plasma%denp-cell(i,j,k)%plasma%denf)
                   loop_along_track: do jj=1,ncell
                      ac=icell(:,jj)
                      call colrad(ac,vihalo(:),tcell(jj),states &
                           ,photons,s2type,nlaunch)
                      if(photons.le.0.d0)cycle loop_over_halos 
                      if(inputs%nospec.eq.0) call spectrum(vihalo(:),ac(:),pos(:,jj),photons,halo_type)
                   enddo loop_along_track
                enddo loop_over_halos
             enddo loop_along_z
          enddo loop_along_y
       enddo loop_along_x
       !$OMP END PARALLEL DO
       halo_iteration_dens=0.d0
       do i=1,grid%Nx 
          do j=1,grid%Ny 
             do k=1,grid%Nz
                halo_iteration_dens=halo_iteration_dens &
                     +sum(dble(cell(i,j,k)%neut_dens(s2type,:)))
                cell(i,j,k)%neut_dens(halo_type,:) = &
                      cell(i,j,k)%neut_dens(halo_type,:) &
                      + cell(i,j,k)%neut_dens(s2type,:)
                cell(i,j,k)%neut_dens(s1type,:)= &
                     cell(i,j,k)%neut_dens(s2type,:)
                cell(i,j,k)%neut_dens(s2type,:)= 0.
             enddo
          enddo
       enddo
       if(halo_iteration_dens/dcx_dens.gt.1)exit iterations
       inputs%nr_halo=floor(inputs%nr_dcx*halo_iteration_dens/dcx_dens)
       if(inputs%nr_halo.lt.inputs%nr_dcx*0.01)exit iterations
    enddo iterations    
  end subroutine halo
  !*****************************************************************************
  !-----------FIDA simulation---------------------------------------------------
  !*****************************************************************************
  subroutine fida      
    integer(long)                         :: i,j,k   !! indices  x,y,z  of cells
    integer(long)                         :: iion
    real(double), dimension(3)            :: ri      !! start position
    real(double), dimension(3)            :: vi      !! velocity of fast ions
    integer(long),dimension(3)            :: ac      !! new actual cell 
    !! Determination of the CX probability
    real(double), dimension(nlevs)        :: denn    !!  neutral dens (n=1-4)
    real(double), dimension(nlevs)        :: prob    !! Prob. for CX 
    real(double), dimension(3)            :: vnbi    !! Velocity of NBI neutrals
    real(double), dimension(3)            :: vnhalo  !! v of halo neutral
    integer(long)                         :: in      !! index of neut rates
    real(double), dimension(nlevs)        :: rates   !! Rate coefficiants for CX
    !! Collisiional radiative model along track
    real(double), dimension(nlevs)        :: states  ! Density of n-states
    integer(long)                         :: ncell
    real(double), dimension(  grid%ntrack):: tcell   !! time per cell
    integer(long),dimension(3,grid%ntrack):: icell   !! index of cells
    real(double), dimension(3,grid%ntrack):: pos     !! mean position in cell
    integer(long)                         :: jj     !! counter along track
    real(double)                          :: photons !! photon flux 
    real(double), dimension(grid%nx,grid%ny,grid%nz)::papprox !! approx. density
    real(double)                          :: dray !! (for NPA)
    real(double), dimension(3)            :: ray     !! ray towards NPA
    integer(long),dimension(3)            :: npa_cell! cell where npa is reached
    real(double)                          :: papprox_tot 
    real(double)                          :: nlaunch !! Number of markers
    real(double)                          :: launch  !! float number of markers
    real(double)                          :: remainder !! used in nlaunch
    real(double),dimension(1)             :: randomu !! used in nlaunch
    !! ------------- calculate papprox needed for guess of nlaunch --------!!
    papprox=0.d0
    papprox_tot=0.d0
    do i = 1, grid%Nx
       do j = 1, grid%Ny
          do k = 1,grid%Nz
             cell(i,j,k)%neut_dens(s1type,:)=cell(i,j,k)%neut_dens(halo_type,:)
             cell(i,j,k)%neut_dens(s2type,:)=0.d0
             if(inputs%npa.eq.1)then
                if(.not.(any(cell(i,j,k)%los_wght(:).gt.0)))cycle
             endif
             papprox(i,j,k)=dble(sum(cell(i,j,k)%neut_dens(nbif_type,:))  &
                  +          sum(cell(i,j,k)%neut_dens(nbih_type,:))  &
                  +          sum(cell(i,j,k)%neut_dens(nbit_type,:))  &
                  +          sum(cell(i,j,k)%neut_dens(halo_type,:))) &
                  *          cell(i,j,k)%plasma%denf
             if(cell(i,j,k)%rho.lt.1.1)papprox_tot=papprox_tot+papprox(i,j,k)
          enddo
       enddo
    enddo
    !$OMP PARALLEL DO private(i,j,k,nlaunch,iion,ac,vi,ri,dray,ray,minpos,npa_cell, &
    !$OMP& tcell,icell,pos,ncell,jj,prob,denn,rates,vnbi,in,vnhalo,states,photons)
    loop_along_x: do i = 1, grid%Nx
       !print*,i, 'of ', grid%Nx
       loop_along_y:    do j = 1, grid%Ny
          loop_along_z:     do k = 1, grid%Nz
             if(papprox(i,j,k).le.0.d0) cycle loop_along_z 
             !! -----Decide how many markers to follow in this cell------- !!
             !! WWH --- Treat fractions more accurately
             call randu(randomu)
             launch=inputs%nr_fida*papprox(i,j,k)/papprox_tot
             nlaunch=int(launch)
             remainder=launch-nlaunch
             if (randomu(1) .le. remainder) then
               nlaunch=nlaunch+1
             endif
             if (nlaunch .lt. 1) cycle loop_along_z
             !nlaunch= int(inputs%nr_fida * papprox(i,j,k)/papprox_tot + 1.)
             !! ------------- loop over the markers ---------------------- !!
             loop_over_fast_ions: do iion=1,int(nlaunch)
                ac=(/i,j,k/)
                !! ---------------- calculate vi, ri and track --------- !!
                call mc_fastion(ac, vi(:)) 
                if(sum(vi).eq.0)cycle loop_over_fast_ions
                call mc_start  (ac, vi(:),  ri(:))
                !! -------- check if track ends at the NPA detector ---- !!
                if(inputs%npa.eq.1)then 
                   npa%at_detector=.false.
                   !! check if track ends at the NPA detector
                   !! uses only the first detector of the array !!
                   dray=sqrt(dot_product(ri-spec%xyzhead(1,:) &
                        ,ri-spec%xyzhead(1,:)))
                   ray=ri(:)+vi(:)/sqrt(dot_product(vi,vi))*dray
                   if(sqrt(dot_product(ray-spec%xyzhead(1,:) &
                        ,ray-spec%xyzhead(1,:)))&
                        .lt.npa%size(1)) then
                      !! define cell where NPA is reached
                      minpos=minloc(abs(ray(1)-(grid%xxc)))
                      npa_cell(1)=minpos(1)
                      minpos=minloc(abs(ray(2)-(grid%yyc)))
                      npa_cell(2)=minpos(1)  
                      minpos=minloc(abs(ray(3)-(grid%zzc)))
                      npa_cell(3)=minpos(1)
                   else
                      cycle loop_over_fast_ions   
                   endif
                endif
                call track(vi(:), ri(:), tcell, icell,pos, ncell)
                if(ncell.eq.0)cycle loop_over_fast_ions
                !! ---------------- calculate CX probability --------------!!
                ac=icell(:,1) !! new actual cell maybe due to gyro orbit!
                prob=0.d0
                vnbi=ri(:)-nbi%xyz_pos(:)
                vnbi=vnbi/sqrt(dot_product(vnbi,vnbi))*nbi%vinj
                ! CX with full energetic NBI neutrals
                denn(:)=dble(cell(ac(1),ac(2),ac(3))%neut_dens(nbif_type,:))
                call neut_rates(denn,vi,vnbi,rates)
                prob=prob + rates
                ! CX with half energetic NBI neutrals
                denn(:)=dble(cell(ac(1),ac(2),ac(3))%neut_dens(nbih_type,:))
                call neut_rates(denn,vi,vnbi/sqrt(2.d0),rates)
                prob=prob + rates
                ! CX with third energetic NBI neutrals
                denn(:)=dble(cell(ac(1),ac(2),ac(3))%neut_dens(nbit_type,:))
                call neut_rates(denn,vi,vnbi/sqrt(3.d0),rates)
                prob=prob + rates
                ! CX with HALO neutrals
                denn(:)=dble(cell(ac(1),ac(2),ac(3))%neut_dens(halo_type,:))
                do in=1,int(nr_halo_neutrate)
                   call mc_halo( ac(:), vnhalo(:))
                   call neut_rates(denn,vi,vnhalo,rates)
                   prob=prob + rates/nr_halo_neutrate
                enddo
                
                if(sum(prob).le.0.)cycle loop_over_fast_ions
                !! --------- solve collisional radiative model along track-!!
                states=prob*cell(i,j,k)%plasma%denf
                
                loop_along_track: do jj=1,ncell        
                   ac=icell(:,jj)
                   if(inputs%npa.eq.1)then
                      if(  ac(1).eq.npa_cell(1).and. &
                           ac(2).eq.npa_cell(2).and. &
                           ac(3).eq.npa_cell(3))then
                         npa%at_detector=.true.
                         npa%counter=npa%counter+1
                         if(npa%counter.gt.inputs%nr_npa)stop'too many neutrals'
                         npa%fpos(npa%counter,:)=pos(:,1)
                         npa%ipos(npa%counter,:)=pos(:,jj)
                      endif
                   endif
                   call colrad(ac(:),vi(:),tcell(jj) &
                        ,states,photons,s2type,nlaunch)
                   if(photons.le.0.d0)cycle loop_over_fast_ions
                   if(inputs%nospec.eq.0) call spectrum(vi(:),ac(:),pos(:,jj),photons,fida_type)
                enddo loop_along_track
             enddo loop_over_fast_ions
          enddo loop_along_z
       enddo loop_along_y
    enddo loop_along_x
    !$OMP END PARALLEL DO 
  end subroutine fida

 
  !*****************************************************************************
  !----------- Calculation of weight functions----------------------------------
  !*****************************************************************************
  subroutine weight_function
      real(double)                   :: radius
    real(double)                   :: photons !! photon flux 
    real(double), dimension(nlevs) :: fdens,hdens,tdens,halodens
    real(double), dimension(3)     :: bvec,avec,cvec,evec,vrot,los_vec
    real(double)                   :: theta
    real(double)                   :: ti,te,dene,denp,deni,length
    real(double)                   :: rad,max_wght
    integer(long)                  :: nwav
    real(double),dimension(:)    ,allocatable  :: wav_arr,central_wavel
    integer(long)                  :: ii,i,j,k,l   !! indices wavel,vx,vy,vz
    real(double),dimension(:,:,:),allocatable :: wfunct
    real(double),dimension(:)    ,allocatable :: ebarr,ptcharr,phiarr
    real(double)                   :: sinus
    real(double),dimension(3)      :: vi,vi_norm
    real(double)                   :: vabs
    real(double),dimension(9)      :: intens !!intensity vector
    real(double),dimension(9)      :: wavel  !!wavelength vector [A)
   !! Determination of the CX probability
    real(double),dimension(3)      :: vnbi_f,vnbi_h,vnbi_t !! Velocity of NBI neutrals 
    real(double),dimension(3)      :: vhalo  !! v of halo neutral
    integer(long)                  :: in      !! index of neut rates
    real(double),dimension(nlevs)  :: rates   !! Rate coefficiants for CX
    real(double),dimension(nlevs)  :: states  ! Density of n-states
    !! COLRAD
    real(double)                   :: dt  !!time interval in cell
    !! ---- Solution of differential equation  ---- ! 
    integer(long),dimension(3)           :: ac  !!actual cell
    real(double), dimension(3)           :: pos !! position of mean cell
    integer(long)                        :: cc 
    real(double),dimension(  grid%nx*grid%ny*grid%nz) :: wght   !! radiation wght per cell 
    real(double),dimension(  grid%nx*grid%ny*grid%nz) :: los_wght !! los wght 
    real(double),dimension(grid%nx,grid%ny,grid%nz,spec%nchan) :: los_weight !! los wght
    integer(long)                          :: ichan
    character(120)                         :: filename
    !! length through cloud of neutrals
    real(double), dimension(3,grid%ntrack) :: pos_out
    real(double), dimension(3)             :: pos_edge
    integer(long)                          :: ic,jc,kc,jj
    integer(long)                          :: ncell  !! number of cells
    real(double), dimension(  grid%ntrack) :: tcell  !! time per cell
    integer(long),dimension(3,grid%ntrack) :: icell  !! index of cells
    real(double)                           :: wght2
    !! DEFINE wavelength array
    nwav=floor((inputs%wavel_end_wght-inputs%wavel_start_wght)/inputs%dwav_wght)
    allocate(wav_arr(nwav+1))
    allocate(central_wavel(nwav+1))
    print*,nwav,' wavelengths to be simulated!'
    do i=1,nwav+1
       wav_arr(i)=dble(i-1.)*inputs%dwav_wght+inputs%wavel_start_wght
    enddo
    central_wavel=wav_arr+0.5*inputs%dwav_wght
    wav_arr=wav_arr*10. !![A]


    !! define pitch, energy and gyro angle arrays
    !! define energy - array
    print*, 'nr of energies, pitches and gyro angles', inputs%nr_wght
    print*, 'maximal energy: ', inputs%emax_wght
    allocate(ebarr(inputs%nr_wght))  
    do i=1,inputs%nr_wght
       ebarr(i)=dble(i-0.5)*inputs%emax_wght/dble(inputs%nr_wght)
    enddo
    !! define pitch - array
    allocate(ptcharr(inputs%nr_wght))
    do i=1,inputs%nr_wght
       ptcharr(i)=dble(i-0.5)*2./dble(inputs%nr_wght)-1.
    enddo
    !! define gyro - array
    allocate(phiarr(inputs%nr_wght))
    do i=1,inputs%nr_wght
       phiarr(i)=dble(i-0.5)*2.d0*pi/dble(inputs%nr_wght)
    enddo
    !! define storage arrays
    allocate(wfunct(nwav,inputs%nr_wght,inputs%nr_wght))
  
    !! Open file for the outputs
    filename=trim(adjustl(result_dir))//"/weight_function.bin" 
    open (66,form='unformatted', file =filename,access='stream')
    write(66)real(inputs%shot_number,float)
    write(66)real(inputs%time,float)
    write(66)real(inputs%ichan_wght,float) 
    write(66)real(inputs%nr_wght,float)  !!Nr of energies
    write(66)real(ebarr(2)-ebarr(1),float)   
    write(66)real(inputs%nr_wght,float)  !!Nr. of pitches
    write(66)real(abs(ptcharr(2)-ptcharr(1)),float)     
    write(66)real(nwav,float) 
    write(66)real(inputs%dwav_wght,float) 
    write(66)real(inputs%wavel_start_wght,float) 
    write(66)real(inputs%wavel_end_wght,float)
    if(inputs%ichan_wght.gt.0) then
       write(66)real(1.,float)
    else
       write(66)real(spec%nchan,float)
    endif
   
    !!save the los-weights into an array
    !! because the structure is over-written
    do i=1,grid%nx 
       do j=1,grid%ny 
          do k=1,grid%nz
             los_weight(i,j,k,:)=cell(i,j,k)%los_wght(:)
          enddo
       enddo
    enddo

    loop_over_channels: do ichan=1,spec%nchan
       if(inputs%ichan_wght.gt.0) then
          if(ichan.ne.inputs%ichan_wght)cycle loop_over_channels
       endif
       write(66)real(ichan,float)
       print*,'channel:',ichan
       radius=sqrt(spec%xyzlos(ichan,1)**2 &
            +spec%xyzlos(ichan,2)**2)
       print*,'Radius:',radius
       !! Calcullate mean kinetic profiles...
       cc=0       ; max_wght=0.d0 ; los_wght=0.d0 ; wght=0.d0
       fdens=0.d0 ; hdens=0.d0    ; tdens=0.d0    ; halodens=0.d0
       bvec=0.d0  ; evec=0.d0
       ti=0.d0    ; te=0.d0
       dene=0.d0  ; denp=0.d0     ; deni=0.d0
       vrot=0.d0  ; pos=0.d0
       do i=1,grid%nx 
          do j=1,grid%ny 
             do k=1,grid%nz
                if(los_weight(i,j,k,ichan).gt.0.)then
                   cc=cc+1
                   los_wght(cc)=los_weight(i,j,k,ichan)
                   !! determine mean values like the halo density along LOS
                   wght(cc)=dble(cell(i,j,k)%neut_dens(nbif_type,3)   &
                        + cell(i,j,k)%neut_dens(nbih_type,3) &
                        + cell(i,j,k)%neut_dens(nbit_type,3)  &
                        + cell(i,j,k)%neut_dens(halo_type,3))*los_wght(cc)
                   if (wght(cc).gt.max_wght)max_wght=wght(cc)
                   fdens=fdens &
                        +dble(cell(i,j,k)%neut_dens(nbif_type,:))*los_wght(cc) 
                   hdens=hdens &
                        +dble(cell(i,j,k)%neut_dens(nbih_type,:))*los_wght(cc) 
                   tdens=tdens &
                        +dble(cell(i,j,k)%neut_dens(nbit_type,:))*los_wght(cc)
                   halodens=halodens &
                        +dble(cell(i,j,k)%neut_dens(halo_type,:))*los_wght(cc)
                   bvec(:)=bvec(:)+cell(i,j,k)%plasma%B(:)  * wght(cc)
                   evec(:)=evec(:)+cell(i,j,k)%plasma%E(:)  * wght(cc)
                   ti     =ti     +cell(i,j,k)%plasma%ti    * wght(cc)
                   te     =te     +cell(i,j,k)%plasma%te    * wght(cc)
                   dene   =dene   +cell(i,j,k)%plasma%dene  * wght(cc)
                   denp   =denp   +cell(i,j,k)%plasma%denp  * wght(cc)
                   deni   =deni   +cell(i,j,k)%plasma%deni  * wght(cc)
                   vrot(:)=vrot(:)+cell(i,j,k)%plasma%vrot(:)*wght(cc)
                   pos(:)=pos(:)+(/grid%xxc(i),grid%yyc(j),grid%zzc(k)/) &
                        *wght(cc)
                endif
             enddo
          enddo
       enddo
       length=sum(los_wght(:)*wght(:))/max_wght ! (FWHM)
       print*,'intersection length of NBI and LOS: ', length
       rad=sum(wght)
       pos =pos     / rad
       vnbi_f=pos(:)-nbi%xyz_pos(:)
       vnbi_f=vnbi_f/sqrt(dot_product(vnbi_f,vnbi_f))*nbi%vinj
       vnbi_h=vnbi_f/sqrt(2.d0)
       vnbi_t=vnbi_f/sqrt(3.d0)     
       !! use only cell 111 for the calculation!
       ac=(/1,1,1/)
       !! normalize quantities
       denp=denp / rad
       cell(ac(1),ac(2),ac(3))%plasma%denp=denp
       dene=dene / rad
       cell(ac(1),ac(2),ac(3))%plasma%dene=dene
       deni=deni / rad  
       cell(ac(1),ac(2),ac(3))%plasma%deni=deni
       ti=ti     / rad
       cell(ac(1),ac(2),ac(3))%plasma%ti=ti
       te=te     / rad
       cell(ac(1),ac(2),ac(3))%plasma%te=te
       vrot=vrot / rad
       cell(ac(1),ac(2),ac(3))%plasma%vrot=vrot
       bvec=bvec    / rad
       print*, '|B|: ',sqrt(dot_product(bvec,bvec)), ' T'
       cell(ac(1),ac(2),ac(3))%plasma%B=bvec   
       evec=evec    / rad
       cell(ac(1),ac(2),ac(3))%plasma%E=evec 
       !! set los_wght to 1 only for one channel (needed for spectrum routine)
       cell(ac(1),ac(2),ac(3))%los_wght(:)=0.
       cell(ac(1),ac(2),ac(3))%los_wght(ichan)=1.
       !! Determine the angle between the B-field and the Line of Sight
       los_vec(1)=pos(1)-spec%xyzhead(ichan,1) 
       los_vec(2)=pos(2)-spec%xyzhead(ichan,2) 
       los_vec(3)=pos(3)-spec%xyzhead(ichan,3) 
       !! normalize los_vec and bvec and determine angle between B and LOS
       los_vec=los_vec/sqrt(dot_product(los_vec,los_vec)) 
       bvec=bvec/sqrt(dot_product(bvec,bvec))
       theta=180.-acos(dot_product(bvec,los_vec))*180./pi
       print*,'Angle between B and LOS [deg]:', theta
       !! write angle and radius into the output file
       write(66)real(theta,float)
       write(66)real(radius,float)
       !! --- calculate vectors avec,cvec that are perpendicular to bvec -- !!
       if (abs(bvec(3)).eq.1) then
          avec=(/1.d0,0.d0,0.d0/)
          cvec=(/0.d0,1.d0,0.d0/)
       else 
          if (bvec(3).eq.0.) then
             avec=(/0.d0,0.d0,1.d0/)
             cvec=(/bvec(2),-bvec(1), 0.d0/)/sqrt(bvec(1)**2+bvec(2)**2)
          else
             avec=(/bvec(2),-bvec(1),0.d0/)/sqrt(bvec(1)**2+bvec(2)**2)
             cvec=(/avec(2),-avec(1),(avec(1)*bvec(2)-avec(2)*bvec(1))/bvec(3)/)
             cvec=cvec/sqrt(dot_product(cvec,cvec))
          endif
       endif
       !! START calculation of weight functions
       print*, 'nwav: ' ,nwav
       wfunct       = 0.d0
       !! do the main simulation  !! 
       !$OMP PARALLEL DO private(i,j,k,vabs,sinus,vi,states,    &
       !$OMP& rates,in,vhalo,dt,photons,wavel,intens,l,ii, &
       !$OMP& tcell,icell,pos_out,ncell,pos_edge,cc,max_wght,   &
       !$OMP& los_wght,wght,jj,ic,jc,kc,wght2,length,vi_norm)
       do i = 1, inputs%nr_wght
          vabs = sqrt(ebarr(i)/(v_to_E*inputs%ab))
          do j = 1, inputs%nr_wght
             sinus = sqrt(1.d0-ptcharr(j)**2)
             do k = 1, inputs%nr_wght
                vi_norm(:)=sinus*cos(phiarr(k))*avec+ptcharr(j) &
                     *bvec+sinus*sin(phiarr(k))*cvec
                call track(-vi_norm,pos,tcell,icell,pos_out,ncell)
                pos_edge=pos_out(:,ncell-2)
                !! now determine the track length throught the grid!
                call track(vi_norm,pos_edge,tcell,icell,pos_out,ncell)
                cc=0 ; max_wght=0.d0 ; los_wght=0.d0 ;wght=0.d0
                loop_along_track: do jj=1,ncell
                   ic=icell(1,jj)
                   jc=icell(2,jj)
                   kc=icell(3,jj)
                   wght2=dble(cell(ic,jc,kc)%neut_dens(nbif_type,3) &
                        + cell(ic,jc,kc)%neut_dens(nbih_type,3) &
                        + cell(ic,jc,kc)%neut_dens(nbit_type,3) &
                        + cell(ic,jc,kc)%neut_dens(halo_type,3))
                   if (wght2.gt.0)then
                      cc=cc+1
                      los_wght(cc)=    tcell(jj)
                      wght(cc)=wght2 * tcell(jj)
                      if (wght(cc).gt.max_wght)max_wght=wght(cc)
                   endif
                enddo loop_along_track
                length=sum(los_wght(:)*wght(:))/max_wght ! (FWHM)
                !! determine time by length and velocity 
                dt=length/vabs
                !! -------------- calculate CX probability -------!!
                ! CX with full energetic NBI neutrals
                states=0.d0
                vi(:) = vi_norm(:)*vabs
                call neut_rates(fdens,vi,vnbi_f,rates)
                states=states + rates
                ! CX with half energetic NBI neutrals
                call neut_rates(hdens,vi,vnbi_h,rates)
                states=states + rates
                ! CX with third energetic NBI neutrals
                call neut_rates(tdens,vi,vnbi_t,rates)
                states=states + rates
                ! CX with HALO neutrals
                do in=1,int(nr_halo_neutrate)
                   call mc_halo( ac(:),vhalo(:))
                   call neut_rates(halodens,vi,vhalo,rates)
                   states=states + rates/nr_halo_neutrate
                enddo
                call colrad(ac,vi,dt,states,photons,nbif_type,1.d0)
                !! photons: [Ph*cm/s/fast-ion]-!!
                call spectrum(vi,ac,pos,1.d0,nbif_type,wavel,intens)
                stark_components: do l=1,9 
                   wavelength_ranges: do ii=1,nwav
                      if (wavel(l).ge.wav_arr(ii).and. &
                           wavel(l).lt.wav_arr(ii+1)) then 
                         wfunct(ii,i,j) = wfunct(ii,i,j) &
                              + intens(l)*photons/dble(inputs%nr_wght)
                      endif
                   enddo wavelength_ranges
                enddo stark_components
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
       !!wfunct:[Ph*cm/s] !!
       do ii=1,nwav
          do i=1,inputs%nr_wght
             do j=1,inputs%nr_wght
                write(66)real(wfunct(ii,i,j),float)
                !![Ph*cm/s/fast-ion] 
             enddo
          enddo
       enddo
    enddo loop_over_channels
    close (66)
    print*, 'weight function written to: ',filename
    deallocate(ebarr)  
    deallocate(ptcharr)
    deallocate(phiarr)
    deallocate(wav_arr)
    deallocate(central_wavel)
    deallocate(wfunct)    
  end subroutine weight_function
end module application
!*****************************************************************************
!-----------Main program ---------------------------------------------------
!*****************************************************************************
program fidasim
  use application
  implicit none 
  integer(long), dimension(8)              :: time_arr,time_start,time_end !Time array
  integer(4)                               :: seed     !seed random generators
  integer(long)                            :: i,j,k
  integer(long)                            :: hour,minu,sec

  !! measure time
  call date_and_time (values=time_start)
  !! get filename of input
  call getarg(1, result_dir)
  seed = 4
  colrad_threshold=1.d6 ![cm^3/s/mc_marker]

  !! 1.e6 is a very low density rate per marker (typially ~1.e12)!
  !! it should not affect the calculation
  call randseed(seed)  
  call read_inputs
  call read_atomic
  if(inputs%nospec.eq.0) call read_los
  call read_plasma
  !! set the neutral density array to zero !!
  do i=1,grid%nx 
     do j=1,grid%ny 
        do k=1,grid%nz
           cell(i,j,k)%neut_dens(:,:)=real(0.,float)
        enddo
     enddo
  enddo

  !! calculate level of bremsstrahlung
  if(inputs%f90brems.eq.1) then
    call bremsstrahlung
  else
    call read_bremsstrahlung
  endif

  if(inputs%load_neutrals.eq.1)then
     call read_neutrals()
  else
     !! -------------------------- ndmc (NBI)----------------------------- !! 
     call date_and_time (values=time_arr)
     print*, 'ndmc:    ' ,time_arr(5:7)
     call ndmc
     !! do the HALO calcualtion only if enough markers are defined!
     if(inputs%nr_halo.gt.100)then
        !! -------------------------- DCX (Direct charge exchange) ---------- !!
        call date_and_time (values=time_arr)
        print*, 'dcx:    ' ,time_arr(5:7)
        call dcx
        !! ------------------------- HALO ----------------------------------- !!
        call date_and_time (values=time_arr)
        print*, 'halo:'   ,time_arr(5:7)
        call halo
     endif
     !! ---------- write output        ----------------------------------- !!   
     call write_neutrals()
      if(inputs%nospec.eq.0) call write_nbi_halo_spectra()
  endif
  !! ---------------- FIDA (main FIDA/NPA simultaion)--------------------- !!
  if(inputs%nofida.eq.0)then    
     call date_and_time (values=time_arr)
     print*, 'da main:    ' ,time_arr(5:7)
     call read_fbm
     print*,'start fida'
     call fida
     !! ------- Store Spectra and neutral densities in binary files ------ !!
     call write_fida_spectra()
     !! ---------------- Store NPA simulation ------------ !!
     if(inputs%npa.eq.1) call write_npa()
  endif
  !! ----------- Calculation of weight functions -----------------------!!
  if(inputs%calc_wght.eq.1) then 
     colrad_threshold=0. !! to speed up simulation!
     call date_and_time (values=time_arr)
     print*, 'weight function:    ' ,time_arr(5:7)
     call weight_function()
  endif
  call date_and_time (values=time_arr)
  call date_and_time (values=time_end)
  print*, 'END: hour, minute, second: ' ,time_arr(5:7)
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
    
  print*, 'duration:                  ',hour,minu,sec
  print*,  inputs%shot_number , inputs%time

  !! --------------- Finally, deallocate allocated arrays ---------------- !!
  !! atomic structure
  deallocate(atomic%qp)
  deallocate(atomic%neut)
  deallocate(atomic%qi)
  deallocate(atomic%qe)
  deallocate(atomic%einstein) 
  !! distribution structure
  if(allocated(distri%energy))then
     deallocate(distri%energy)
     deallocate(distri%pitch)
     deallocate(distri%fbm)
  endif
  !! grid and cell structure
  if(inputs%nospec.eq.0)then
     do i = 1, grid%Nx
        do j = 1, grid%Ny
           do k = 1, grid%Nz 
              deallocate(cell(i,j,k)%los_wght) 
           enddo
        enddo
     enddo
     deallocate(spec%xyzlos)
     deallocate(spec%spectra)
     deallocate(spec%xyzhead)
  endif
  deallocate(cell)
  deallocate(grid%xx)
  deallocate(grid%yy)
  deallocate(grid%zz)
  !!deallocate npa arrays
  if(inputs%npa .eq. 1)then
     deallocate(npa%v) 
     deallocate(npa%ipos) 
     deallocate(npa%fpos) 
     deallocate(npa%wght)
  endif

end program fidasim
 
