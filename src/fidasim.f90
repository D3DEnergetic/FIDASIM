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
        endif
        if(inputs%calc_dcx.ge.1) then
            allocate(spec%dcx(n_stark,inputs%nlambda,spec_chords%nchan,n_thermal))
            spec%dcx = 0.d0
        endif
        if(inputs%calc_halo.ge.1) then
            allocate(spec%halo(n_stark,inputs%nlambda,spec_chords%nchan,n_thermal))
            spec%halo = 0.d0
        endif
        if(inputs%calc_cold.ge.1) then
            allocate(spec%cold(n_stark,inputs%nlambda,spec_chords%nchan,n_thermal))
            spec%cold = 0.d0
        endif
        if(inputs%calc_fida.ge.1) then
            allocate(spec%fida(n_stark,inputs%nlambda,spec_chords%nchan,particles%nclass))
            spec%fida = 0.d0
        endif
        if(inputs%calc_pfida.ge.1) then
            allocate(spec%pfida(n_stark,inputs%nlambda,spec_chords%nchan,particles%nclass))
            spec%pfida = 0.d0
        endif
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
                write(*,*) 'nbi:     ' , time(time_start)
            endif
            call nbi_spec()
            if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
        endif

        if(inputs%calc_dcx.ge.1) then
            if(inputs%verbose.ge.1) then
                write(*,*) 'dcx:     ' , time(time_start)
            endif
            call dcx_spec()
            if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
        endif

        if(inputs%calc_halo.ge.1) then
            if(inputs%verbose.ge.1) then
                write(*,*) 'halo:    ' , time(time_start)
            endif
            call halo_spec()
            if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
        endif
    else
        if(inputs%calc_beam.ge.1) then
            !! ----------- BEAM NEUTRALS ---------- !!
            if(inputs%calc_nbi_dens.ge.1) then
                if(inputs%verbose.ge.1) then
                    write(*,*) 'nbi:     ' , time(time_start)
                endif
                call ndmc()
                if(inputs%verbose.ge.1) write(*,'(30X,a)') ''

                if(inputs%calc_birth.eq.1)then
                    if(inputs%verbose.ge.1) then
                        write(*,*) 'write birth:    ' , time(time_start)
                    endif
                    call write_birth_profile()
                    if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
                endif
            endif

            !! ---------- DCX (Direct charge exchange) ---------- !!
            if(inputs%calc_dcx_dens.ge.1) then
                if(inputs%verbose.ge.1) then
                    write(*,*) 'dcx:     ' , time(time_start)
                endif
                call dcx()
                if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
            endif

            !! ---------- HALO ---------- !!
            if(inputs%calc_halo_dens.ge.1) then
                if(inputs%verbose.ge.1) then
                    write(*,*) 'halo:    ' , time(time_start)
                endif
                call halo()
                if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
            endif

            !! ---------- WRITE NEUTRALS ---------- !!
            if(inputs%verbose.ge.1) then
                write(*,*) 'write neutrals:    ' , time(time_start)
            endif
#ifdef _MPI
            if(my_rank().eq.0) call write_neutrals()
#else
            call write_neutrals()
#endif
            if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
        endif
    endif

    !! -----------------------------------------------------------------------
    !!------------------------------ COLD D-ALPHA ----------------------------
    !! -----------------------------------------------------------------------
    if(inputs%calc_cold.ge.1) then
        if(inputs%verbose.ge.1) then
            write(*,*) 'cold:    ' ,time(time_start)
        endif
        call cold_spec()
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
    endif

    !! -----------------------------------------------------------------------
    !!----------------------------- BREMSSTRAHLUNG ---------------------------
    !! -----------------------------------------------------------------------
    if(inputs%calc_brems.ge.1) then
        if(inputs%verbose.ge.1) then
            write(*,*) 'bremsstrahlung:    ' ,time(time_start)
        endif
        call bremsstrahlung()
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
    endif

    !! -----------------------------------------------------------------------
    !! --------------------- CALCULATE the FIDA RADIATION --------------------
    !! -----------------------------------------------------------------------
    if(inputs%calc_fida.ge.1)then
        if(inputs%verbose.ge.1) then
            write(*,*) 'fida:    ' ,time(time_start)
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
            write(*,*) 'pfida:   ' ,time(time_start)
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
            write(*,*) 'write spectra:    ' , time(time_start)
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
            write(*,*) 'npa:     ' ,time(time_start)
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
            write(*,*) 'pnpa:     ' ,time(time_start)
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
            write(*,*) 'write npa:    ' , time(time_start)
        endif
        call write_npa()
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
    endif

    !! -------------------------------------------------------------------
    !! ------------------- Calculation of neutron flux -------------------
    !! -------------------------------------------------------------------
    if(inputs%calc_neutron.ge.1) then
        if(inputs%verbose.ge.1) then
            write(*,*) 'neutron rate:    ', time(time_start)
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
        if(inputs%verbose.ge.1) then
            write(*,*) 'fida weight function:    ', time(time_start)
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
            write(*,*) 'npa weight function:    ', time(time_start)
        endif
        call npa_weights()
        if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
    endif

#ifdef _MPI
    call cleanup_mpi()
#endif

    if(inputs%verbose.ge.1) then
        write(*,*) 'END: hour:minute:second ', time(time_start)
    endif

end program fidasim
