module cx_helper
    implicit none
    contains
    subroutine cx_matrix(velocity,rates)
        !! Calculate the charge-exchange cross-section matrix for a given relative collision velocity in [cm/s].
        !! The matrix gives rates of the incoming (n) & outgoing (n') energy levels of the electron, up to n=6
        !! The rates are given in [cm^-3*s^-1]
        use libfida, only : SimulationInputs,bb_cx_rates,read_tables, nlevs,inputs,impurity_charge
        implicit none
        integer, parameter   :: Float64 = 8 !!Use 8-byte floats (ie, Float64)
        real(Float64), intent(in) :: velocity
        real(Float64), dimension(nlevs,nlevs), intent(out) :: rates !! catches the output from bb_cx_rates
        integer :: n_in !! incoming neutral energy level (principle quantum number)
        real(Float64), dimension(nlevs) :: denn  !! density [cm^-3] of each energy state of neutrals
        real(Float64), dimension(3) :: vn = [0.0,0.0,0.0] !! cm/s x,y,z components of velocity of neutral
        real(Float64), parameter, dimension(3) :: vi = [0.0,0.0,0.0]!! velocity of ion
       
        vn(1) = velocity
        impurity_charge = 6
        
        !!type(SimulationInputs) :: inputs !! needed to specify a few key things used in the rate calculator
        inputs%verbose = 1
        inputs%tables_file = '/Users/lmorton/Code/FIDASIM_lib/tables/atomic_tables.h5'
        inputs%calc_neutron = 0
        call read_tables
        do n_in=1,6
            denn = 0.0
            denn(n_in) = 1.0
            call bb_cx_rates(denn,vn,vi,rates(n_in,:))
            print *,rates(n_in,:)
        enddo
        !! CX reaction rate in units [1/s] (ie, reactions per incident fast ion per second, use the fast ion velocity to get a decay length, or the fast ion density to get a volumetric rate)
        !! The rate has 'nlev' values, for each of the energy levels that the post-CX fast neutral could have
    end subroutine cx_matrix
end module cx_helper