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

    ! >>> [JFCM, 2025-03-14] >>>
    integer :: time_0, time_1, count_rate, ss, rr
    real :: elapsed_time
    ! <<< [JFCM, 2025-03-14] <<<
    ! >>> [JFCM, 2025-09-02] >>>
    integer :: n_birth_wall
    integer :: n_birth_nbi
    ! <<< [JFCM, 2025-09-02] <<<

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
    !! >>> [JFCM, 2025-10-27] >>>
    ! Update reservour size:
    reservoir_size = inputs%reservoir_size
    write(*,*) "RESERVOIR_SIZE modified and set to ",reservoir_size
    write(*,*) ""
    !! <<< [JFCM, 2025-10-27] <<<

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
    ! >>> [JFCM, 2025_07_23] >>>
    call define_vacuum_vessel()
    ! <<< [JFCM, 2025_07_23] <<<
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
        ! allocate(birth%dens(3, &
        !                     beam_grid%nx, &
        !                     beam_grid%ny, &
        !                     beam_grid%nz))
        !! >>> [jfcm, 2024-11-23] >>>
        allocate(birth%dens(5, &
                            beam_grid%nx, &
                            beam_grid%ny, &
                            beam_grid%nz))
        !! >>> [jfcm, 2024-11-23] >>>

        !! >>> [JFCM, 2025-09-02] >>>
        ! TODO: need to allow using births_per_marker and more then one birth per marker for NBI
        n_birth_wall = 0
        do ss = 1,size(vessel%surface)
          do rr = 1,size(vessel%surface(ss)%region)
            if (vessel%surface(ss)%region(rr)%enable_source) then
              n_birth_wall = n_birth_wall + vessel%surface(ss)%region(rr)%source%num_markers
            endif
          end do
        end do
        !! <<< [JFCM, 2025-09-02] <<<

        !! >>> [JFCM, 2025-09-02] >>>
        ! allocate(birth%part(int(3*inputs%n_birth*inputs%n_nbi)))
        n_birth_nbi = int(3*inputs%n_birth*inputs%n_nbi)
        allocate(birth%part(n_birth_nbi + n_birth_wall))
        !! <<< [JFCM, 2025-09-02] <<<

    endif

    if (inputs%calc_sink.ge.1) then
      ! allocate(sink%dens(n_thermal, &
      !                     beam_grid%nx, &
      !                     beam_grid%ny, &
      !                     beam_grid%nz))
      ! allocate(sink%part(inputs%n_dcx))
      ! allocate(sink%part(size(birth%part)))
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

                ! >>> [JFCM, 2025-03-14] >>>
                call SYSTEM_CLOCK(time_0,count_rate)
                ! <<< [JFCM, 2025-03-14] <<<

                ! NBI calculation:
                call ndmc()

                ! >>> [JFCM, 2025-03-14] >>>
                call SYSTEM_CLOCK(time_1)
                elapsed_time = real(time_1 - time_0) / real(count_rate)  ! Convert to seconds
                print *, "NDMC elapsed time: ", elapsed_time, " seconds"
                ! <<< [JFCM, 2025-03-14] <<<

                ! >>> [JFCM, 2025-09-02] >>>
                ! Wall source calculation:
                do ss = 1,size(vessel%surface)
                  do rr = 1,size(vessel%surface(ss)%region)
                    if (vessel%surface(ss)%region(rr)%enable_source) then
                      write(*,*) "calculate wall source"
                      ! call calculate_wall_source_process(ss,rr)
                    endif
                end do
              end do
                ! <<< [JFCM, 2025-09-02] <<<

                if(inputs%verbose.ge.1) write(*,'(30X,a)') ''

                if(inputs%calc_birth.eq.1)then
                    if(inputs%verbose.ge.1) then
                        write(*,*) 'write 0th gen birth:    ' , time_string(time_start)
                    endif
                    call write_birth_profile()
                    if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
                endif

                if (inputs%calc_birth.eq.1 .and. inputs%calc_sink.ge.1) then
                  if(inputs%verbose.ge.1) then
                      write(*,*) 'CALCULATE_DCX_PROCESS:    ' , time_string(time_start)
                  endif
                  call calculate_dcx_process

                  if(inputs%verbose.ge.1) then
                      write(*,*) 'write 0th gen sink and 1st gen birth:    ' , time_string(time_start)
                  endif
                  call write_sink_profile()
                  call write_birth_profile(gen=1)
                  if(inputs%verbose.ge.1) write(*,'(30X,a)') ''
                endif

                !! >>> [JFCM, 2025-07-02] >>>
                !! TODO: need switches for user to enable this
                !! TODO: consider using n_halo_gen = 0 for nbi + dcx only, n_halo_gen > 0 for nbi + dcx + halo
                !! TODO: Incorporate write_<src>_profile into calculate_dcx_process
                !! TODO: Fix the photon accumulation in dcx and halo process subroutines to mimic orignal dcx and halo subroutines, then add switch to skip or enable calc_spec so what we can move dcx and halo process subroutines to where dcx and halo are located.
                !! TODO: Bring the wall condition inputs to the user interface
                if (inputs%calc_birth.eq.1 .and. inputs%calc_sink.ge.1 .and. inputs%enable_halo.ge.1 ) then
                  if(inputs%verbose.ge.1) then
                      write(*,*) 'CALCULATE_HALO_PROCESS:    ' , time_string(time_start)
                  endif
                  call calculate_halo_process

                endif
                !! <<< [JFCM, 2025-07-02] <<<

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
            if(inputs%verbose.ge.1) then
                write(*,*) 'write neutrals:    ' , time_string(time_start)
            endif
#ifdef _MPI
            if(my_rank().eq.0) call write_neutrals()
#else
            call write_neutrals()
#endif
            if(inputs%verbose.ge.1) write(*,'(30X,a)') ''

            !! Deallocating birth and sink:
            ! >>>>>>>>>>> [jfcm, 2024-11-23] >>>>>>>>>>>>>>>
            ! if (inputs%calc_birth.ge.1) deallocate(birth%dens,birth%part)
            ! if (inputs%calc_sink.ge.1) deallocate(sink%dens,sink%part)
            ! Since we disabled the deallocation step in write_birth_profile(), we need to deallocate here
            ! <<<<<<<<<<< [jfcm, 2024-11-23] <<<<<<<<<<<<<<<
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
