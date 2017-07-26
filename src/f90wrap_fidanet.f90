! Module fidanet defined in file fidanet.f90

subroutine f90wrap_settables
    use fidanet, only: settables
    implicit none
    
    call settables()
end subroutine f90wrap_settables

subroutine f90wrap_setplasma(dene, denp, denimp, te, ti)
    use fidanet, only: setplasma
    implicit none
    
    real(8), intent(in) :: dene
    real(8), intent(in) :: denp
    real(8), intent(in) :: denimp
    real(8), intent(in) :: te
    real(8), intent(in) :: ti
    call setplasma(dene=dene, denp=denp, denimp=denimp, te=te, ti=ti)
end subroutine f90wrap_setplasma

subroutine f90wrap_setinputs(ai, ab, impq)
    use fidanet, only: setinputs
    implicit none
    
    real(8), intent(in) :: ai
    real(8), intent(in) :: ab
    integer, intent(in) :: impq
    call setinputs(ai=ai, ab=ab, impq=impq)
end subroutine f90wrap_setinputs

subroutine f90wrap_calcvn(i_type, eb, vn)
    use fidanet, only: calcvn
    implicit none
    
    integer, intent(in) :: i_type
    real(8), intent(in) :: eb
    real(8), dimension(3), intent(inout) :: vn
    call calcvn(i_type=i_type, eb=eb, vn=vn)
end subroutine f90wrap_calcvn

subroutine f90wrap_setstates(newstates, states)
    use fidanet, only: setstates
    implicit none
    
    real(8), dimension(6), intent(in) :: newstates
    real(8), dimension(6), intent(inout) :: states
    call setstates(newstates=newstates, states=states)
end subroutine f90wrap_setstates

subroutine f90wrap_testcol(i_type, eb, dt, states, dens, photons)
    use fidanet, only: testcol
    implicit none
    
    integer, intent(in) :: i_type
    real(8), intent(in) :: eb
    real(8), intent(in) :: dt
    real(8), dimension(6), intent(inout) :: states
    real(8), dimension(6), intent(inout) :: dens
    real(8), intent(out) :: photons
    call testcol(i_type=i_type, eb=eb, dt=dt, states=states, dens=dens, &
        photons=photons)
end subroutine f90wrap_testcol

subroutine f90wrap_fidanet__get__ai(f90wrap_ai)
    use fidanet, only: fidanet_ai => ai
    implicit none
    real(8), intent(out) :: f90wrap_ai
    
    f90wrap_ai = fidanet_ai
end subroutine f90wrap_fidanet__get__ai

subroutine f90wrap_fidanet__set__ai(f90wrap_ai)
    use fidanet, only: fidanet_ai => ai
    implicit none
    real(8), intent(in) :: f90wrap_ai
    
    fidanet_ai = f90wrap_ai
end subroutine f90wrap_fidanet__set__ai

subroutine f90wrap_fidanet__get__ab(f90wrap_ab)
    use fidanet, only: fidanet_ab => ab
    implicit none
    real(8), intent(out) :: f90wrap_ab
    
    f90wrap_ab = fidanet_ab
end subroutine f90wrap_fidanet__get__ab

subroutine f90wrap_fidanet__set__ab(f90wrap_ab)
    use fidanet, only: fidanet_ab => ab
    implicit none
    real(8), intent(in) :: f90wrap_ab
    
    fidanet_ab = f90wrap_ab
end subroutine f90wrap_fidanet__set__ab

subroutine f90wrap_fidanet__array__states(dummy_this, nd, dtype, dshape, dloc)
    use libfida
    use fidanet, only: fidanet_states => states
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(fidanet_states)
    dloc = loc(fidanet_states)
end subroutine f90wrap_fidanet__array__states

subroutine f90wrap_fidanet__array__dens(dummy_this, nd, dtype, dshape, dloc)
    use libfida
    use fidanet, only: fidanet_dens => dens
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(fidanet_dens)
    dloc = loc(fidanet_dens)
end subroutine f90wrap_fidanet__array__dens

subroutine f90wrap_fidanet__get__photons(f90wrap_photons)
    use fidanet, only: fidanet_photons => photons
    implicit none
    real(8), intent(out) :: f90wrap_photons
    
    f90wrap_photons = fidanet_photons
end subroutine f90wrap_fidanet__get__photons

subroutine f90wrap_fidanet__set__photons(f90wrap_photons)
    use fidanet, only: fidanet_photons => photons
    implicit none
    real(8), intent(in) :: f90wrap_photons
    
    fidanet_photons = f90wrap_photons
end subroutine f90wrap_fidanet__set__photons

! End of module fidanet defined in file fidanet.f90

