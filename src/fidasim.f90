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

    allocate(spec_chords%inter(beam_grid%nx,beam_grid%ny,beam_grid%nz))
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

    if(inputs%calc_neutron.ge.1)then
        allocate(neutron%rate(particles%nclass))
        neutron%rate = 0.d0
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
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''

        !! ---------- DCX (Direct charge exchange) ---------- !!
        call date_and_time (values=time_arr)
        if(inputs%verbose.ge.1) then
            write(*,'(A,I2,":",I2.2,":",I2.2)') 'dcx:     ' , &
                  time_arr(5),time_arr(6),time_arr(7)
        endif
        call dcx()
        if(inputs%dump_dcx.eq.1) call write_dcx()
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''

        !! ---------- HALO ---------- !!
        call date_and_time (values=time_arr)
        if(inputs%verbose.ge.1) then
            write(*,'(A,I2,":",I2.2,":",I2.2)') 'halo:    ' , &
                  time_arr(5),time_arr(6),time_arr(7)
        endif
        call halo()
        !! ---------- WRITE NEUTRALS ---------- !!
        call write_neutrals()
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
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
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
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
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
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
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
    endif

    if(inputs%calc_npa.ge.1) then
        call write_npa()
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
    endif

    !! -------------------------------------------------------------------
    !! ------------------- Calculation of neutron flux -------------------
    !! -------------------------------------------------------------------
    if(inputs%calc_neutron.ge.1) then
        call date_and_time (values=time_arr)
        if(inputs%verbose.ge.1) then
            write(*,'(A,I2,":",I2.2,":",I2.2)') 'neutron rate:    ',  &
                  time_arr(5),time_arr(6),time_arr(7)
        endif
        if(inputs%dist_type.eq.1) then
            call neutron_f()
        else
            call neutron_mc()
        endif
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
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
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
    endif

    if(inputs%calc_npa_wght.ge.1) then
        call date_and_time (values=time_arr)
        if(inputs%verbose.ge.1) then
            write(*,'(A,I2,":",I2.2,":",I2.2)') 'npa weight function:    ',  &
                  time_arr(5),time_arr(6),time_arr(7)
        endif
        call npa_weights()
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
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
