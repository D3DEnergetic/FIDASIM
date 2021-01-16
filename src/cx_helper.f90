module cx_helper
    implicit none
    private
    public :: cx_matrix
    contains
    subroutine cx_matrix(velocity,nlevels,rates) bind(C,name="cx_matrix")!{{{
        !! Calculate the charge-exchange cross-section matrix for a given relative collision velocity in [cm/s].
        !! The matrix gives rates of the incoming (n) & outgoing (n') energy levels of the electron, up to n=6
        !! The rates are given in [cm^-3*s^-1]
        use libfida, only : SimulationInputs,bb_cx_rates,read_tables, nlevs,inputs,impurity_charge
        use, intrinsic :: iso_c_binding, only: c_int, c_double
        implicit none
        real(c_double), intent(in) :: velocity
        integer, intent(in) :: nlevels
        real(c_double), dimension(nlevels,nlevels), intent(inout) :: rates !! catches the output from bb_cx_rates
        integer :: n_in !! incoming neutral energy level (principle quantum number)
        real(c_double), dimension(nlevels) :: denn  !! density [cm^-3] of each energy state of neutrals
        real(c_double), dimension(3) :: vn = [0.0,0.0,0.0] !! cm/s x,y,z components of velocity of neutral
        real(c_double), parameter, dimension(3) :: vi = [0.0,0.0,0.0]!! velocity of ion
        if (nlevels .gt. nlevs) then
            print *,"nlevels must be less than ",nlevs
            stop
        endif 
        vn(1) = velocity
        impurity_charge = 6

        inputs%verbose = 1
        inputs%tables_file = '/Users/lmorton/Code/FIDASIM_lib/tables/atomic_tables.h5'
        inputs%calc_neutron = 0
        call read_tables
        do n_in=1,6
            denn = 0.0
            denn(n_in) = 1.0
            call bb_cx_rates(denn,vn,vi,rates(n_in,:))
        enddo
        !! CX reaction rate in units [1/s] (ie, reactions per incident fast ion per second, use the fast ion velocity to get a decay length, or the fast ion density to get a volumetric rate)
        !! The rate has 'nlev' values, for each of the energy levels that the post-CX fast neutral could have
    end subroutine cx_matrix!}}}
end module cx_helper