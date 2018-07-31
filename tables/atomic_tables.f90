!+This file contains the routines for calculating atomic cross sections and reaction rates for FIDASIM
module atomic_tables
!+Library for calculating atomic cross sections and reaction rate coefficients for Hydrogen interactions
!+
!+###References
!+
!+1. [W.L. Wiese, M.W. Smith, and B.M. Glennon. *Atomic Transition Probabilities. Volume 1. Hydrogen through Neon*.
!+National Bureau of Standards Washington DC Institute for Basic Standards, 1966.](http://www.dtic.mil/dtic/tr/fulltext/u2/634145.pdf)
!+2. [R.K. Janev, D. Reiter, and  U. Samm. *Collision processes in low-temperature hydrogen plasmas*.
!+Forschungszentrum Jülich, Zentralbibliothek, 2003.](http://www.eirene.de/report_4105.pdf)
!+3. [M. O'Mullane. *Review of proton impact driven ionisation from the excited levels in neutral hydrogen beams*.
!+ADAS note, 2009.](http://www.adas.ac.uk/notes/adas_c09-01.pdf)
!+4. [ADAS: Atomic Data and Analysis Structure](http://www.adas.ac.uk/)
!+5. [R.K. Janev and J.J. Smith. *Cross sections for collision processes of hydrogen atoms
!+with electrons, protons and multiply charged ions.* Atomic and Plasma-Material Interaction Data for Fusion:
!+Volume 4, 1993.](http://www-pub.iaea.org/books/IAEABooks/1839/Atomic-and-Plasma-Material-Interaction-Data-for-Fusion)
!+6. [Reinhold, C. O., R. E. Olson, and W. Fritsch. *Excitation of atomic hydrogen by fully stripped ions.*
!+Physical Review A 41.9 1990.](http://journals.aps.org/pra/abstract/10.1103/PhysRevA.41.4837)
!+7. [Bosch, H-S., and G. M. Hale. *Improved formulas for fusion cross-sections and thermal reactivities.*
!+ Nuclear fusion 32.4 1992.](http://iopscience.iop.org/article/10.1088/0029-5515/32/4/I07/meta)
!+8. [Aladdin Database: R.K. Janev, W.D. Langer, K. Evans Jr., D.E. Post Jr. H-HE-PLASMA (1987)](https://www-amdis.iaea.org/ALADDIN/collision.html)
use H5LT
use HDF5
use hdf5_utils
#ifdef _MPI
use mpi_utils
#endif

IMPLICIT NONE

interface bt_maxwellian
    !+Calculates the reaction rate coefficients given beam energy `eb` and target temperature `T`
    !+where the velocity distribution of the target is a Maxwellian
    module procedure bt_maxwellian_eb
    module procedure bt_maxwellian_n, bt_maxwellian_n_m
    module procedure bt_maxwellian_q_n, bt_maxwellian_q_n_m
end interface

integer, parameter, private   :: Int32   = 4
    !+ Defines a 32 bit integer
integer, parameter, private   :: Int64   = 8
    !+ Defines a 64 bit integer
integer, parameter, private   :: Float32 = 4
    !+ Defines a 32 bit floating point real
integer, parameter, private   :: Float64 = 8
    !+ Defines a 64 bit floating point real

real(Float64), parameter :: PI = 3.14159265d0
real(Float64), parameter :: e_amu = 5.485799093287202d-4
    !+ Atomic mass of an electron [amu]
real(Float64), parameter :: H1_amu = 1.00782504d0
    !+ Atomic mass of Hydrogen-1 (protium) [amu]
real(Float64), parameter :: H2_amu = 2.0141017778d0
    !+ Atomic mass of Hydrogen-2 (deuterium) [amu]
real(Float64), parameter :: H3_amu = 3.0160492d0
    !+ Atomic mass of Hydrogen-3 (tritium) [amu]
real(Float64), parameter :: He3_amu = 3.0160293d0
    !+ Atomic mass of Helium-3 [amu]
real(Float64), parameter :: B_amu = 10.81d0
    !+ Atomic mass of Boron [amu]
real(Float64), parameter :: C_amu = 12.011d0
    !+ Atomic mass of Carbon [amu]

integer, parameter :: B_q = 5
    !+ Proton number of Boron
integer, parameter :: C_q = 6
    !+ Proton number of Carbon

real(Float64), dimension(15,15), parameter :: EINSTEIN = reshape([ &                                                  !(n,m)
0.d0,4.699d8,5.575d7,1.278d7,4.125d6,1.644d6,7.568d5,3.869d5,2.143d5,1.263d5,7.834d4,5.066d4,3.393d4,2.341d4,1.657d4,&!(:,1)
0.d0,0.d0   ,4.410d7,8.419d6,2.530d6,9.732d5,4.389d5,2.215d5,1.216d5,7.122d4,4.397d4,2.834d4,1.893d4,1.303d4,9.210d3,&!(:,2)
0.d0,0.d0   ,0.d0   ,8.986d6,2.201d6,7.783d5,3.358d5,1.651d5,8.905d4,5.156d4,3.156d4,2.021d4,1.343d4,9.211d3,6.490d3,&!(:,3)
0.d0,0.d0   ,0.d0   ,0.d0   ,2.699d6,7.711d5,3.041d5,1.424d5,7.459d4,4.235d4,2.556d4,1.620d4,1.069d4,7.288d3,5.110d3,&!(:,4)
0.d0,0.d0   ,0.d0   ,0.d0   ,0.d0   ,1.025d6,3.253d5,1.388d5,6.908d4,3.800d4,2.246d4,1.402d4,9.148d3,6.185d3,4.308d3,&!(:,5)
0.d0,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,4.561d5,1.561d5,7.065d4,3.688d4,2.110d4,1.288d4,8.271d3,5.526d3,3.815d3,&!(:,6)
0.d0,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,2.272d5,8.237d4,3.905d4,2.117d4,1.250d4,7.845d3,5.156d3,3.516d3,&!(:,7)
0.d0,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,1.233d5,4.676d4,2.301d4,1.287d4,7.804d3,5.010d3,3.359d3,&!(:,8)
0.d0,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,7.141d4,2.812d4,1.427d4,8.192d3,5.080d3,3.325d3,&!(:,9)
0.d0,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,4.377d4,1.774d4,9.231d3,5.417d3,3.324d3,&!(:,10)
0.d0,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,2.799d4,1.163d4,6.186d3,3.699d3,&!(:,11)
0.d0,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,1.857d4,7.884d3,4.271d3,&!(:,12)
0.d0,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,1.271d4,5.496d3,&!(:,13)
0.d0,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,8.933d3,&!(:,14)
0.d0,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ]&!(:,15)
, [15,15])
    !+ Einstein coefficients for spontaneous emission from state initial state `n` to final state `m`
    !+
    !+References:
    !+
    !+* H - Table A in Ref. 1 [[atomic_tables(module)]]

!!Loop Parallization Settings
integer :: istart = 1
    !+ Starting loop counter (1 if OpenMP, processor number if MPI)
integer :: istep = 1
    !+ Loop step size (1 if OpenMP, number of processes if MPI)
logical :: verbose = .True.
    !+ Indicates whether process is verbose

contains

function p_cx_1_janev(Erel) result(sigma)
    !+Calculates total cross section for proton-Hydrogen charge exchange interactions from the \(n=1\) state at energy `Erel`
    !+
    !+###Equation
    !+ $$H^+ + H(1) \rightarrow H(\forall m) + H^+$$
    !+###References
    !+* Eq. 44 and Table 9 in Ref. 2 [[atomic_tables(module)]]
    real(Float64), intent(in) :: Erel
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]
    real(Float64), dimension(6), parameter :: a = [3.2345d0,  2.3588d2,  2.3713d0, &
                                                   3.8371d-2, 3.8068d-6, 1.1832d-10 ]
        !+ Fitting Parameters from Table 9 in Ref. 2
    real(Float64), parameter :: n = 1.d0
    real(Float64) :: Ehat

    Ehat = Erel * n**2.0

    sigma = (1.d-16*a(1)*(n**4))*log(a(2)/Ehat + a(3)) / &
            (1.d0+a(4)*Ehat + a(5)*Ehat**(3.5) + a(6)*Ehat**(5.4))

end function p_cx_1_janev

function p_cx_2_janev(Erel) result(sigma)
    !+Calculates total cross section for proton-Hydrogen charge exchange interactions from the \(n=2\) state at energy `Erel`
    !+
    !+###Equation
    !+ $$H^+ + H(2) \rightarrow H(\forall m) + H^+$$
    !+###References
    !+* Eq. 44 and Table 9 in Ref. 2 [[atomic_tables(module)]]
    real(Float64), intent(in) :: Erel
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]
    real(Float64), dimension(6), parameter :: a = [9.2750d-1, 6.5040d3, 2.0699d1, &
                                                   1.3405d-2, 3.0842d-6, 1.1832d-10 ]
        !+ Fitting Parameters from Table 9 in Ref. 2
    real(Float64), parameter :: n = 2.d0
    real(Float64) :: Ehat

    Ehat = Erel * n**2.0

    sigma = (1.d-16*a(1)*(n**4))*log(a(2)/Ehat + a(3)) / &
            (1.d0+a(4)*Ehat + a(5)*Ehat**(3.5) + a(6)*Ehat**(5.4))

end function p_cx_2_janev

function p_cx_3_janev(Erel) result(sigma)
    !+Calculates total cross section for proton-Hydrogen charge exchange interactions from the \(n=3\) state at energy `Erel`
    !+
    !+###Equation
    !+ $$H^+ + H(3) \rightarrow H(\forall m) + H^+$$
    !+###References
    !+* Eq. 44 and Table 9 in Ref. 2 [[atomic_tables(module)]]
    real(Float64), intent(in) :: Erel
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]
    real(Float64), dimension(6), parameter :: a = [3.7271d-1, 2.7645d6,  1.4857d3, &
                                                   1.5720d-3, 3.0842d-6, 1.1832d-10 ]
        !+ Fitting Parameters from Table 9 in Ref. 2


    real(Float64), parameter :: n = 3.d0
    real(Float64) :: Ehat

    Ehat = Erel * n**2.0

    sigma = (1.d-16*a(1)*(n**4))*log(a(2)/Ehat + a(3)) / &
            (1.d0+a(4)*Ehat + a(5)*Ehat**(3.5) + a(6)*Ehat**(5.4))

end function p_cx_3_janev

function p_cx_n_janev(Erel, n) result(sigma)
    !+Calculates cross section for proton-Hydrogen charge exchange interactions from the \(n \geq 4\) state at energy `Erel`
    !+
    !+###Equation
    !+ $$H^+ + H(n \geq 4) \rightarrow H(\forall m) + H^+$$
    !+###References
    !+* Eq. 44 and Table 9 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: Erel
        !+ Relative collision energy [keV/amu]
    integer, intent(in)       :: n
        !+ Initial atomic energy level/state
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]
    real(Float64), dimension(6), parameter :: a = [2.1336d-1, 1.0000d10, 1.3426d6, &
                                                   1.8184d-3, 3.0842d-6, 1.1832d-10 ]
        !+ Fitting Parameters from Table 9 in Ref. 2

    real(Float64) :: Ehat

    if(n.lt.4) then
        write(*,'(a)') "P_CX_N_JANEV: n cannot be less than 4"
        stop
    endif

    Ehat = Erel * n**2.0
    sigma = (1.d-16*a(1)*(n**4))*log(a(2)/Ehat + a(3)) / &
            (1.d0+a(4)*Ehat + a(5)*Ehat**(3.5) + a(6)*Ehat**(5.4))

end function p_cx_n_janev

function p_cx_janev(Erel,n) result(sigma)
    !+Calculates total cross section for proton-Hydrogen charge exchange interactions from the `n` state at energy `Erel`
    !+
    !+###Equation
    !+ $$H^+ + H(n) \rightarrow H(m) + H^+$$
    !+###References
    !+* Eq. 44 and Table 9 in Ref. 2 [[atomic_tables(module)]]
    real(Float64), intent(in) :: Erel
        !+ Relative collision energy [keV/amu]
    integer, intent(in)       :: n
        !+ Initial atomic energy level/state
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    integer :: i

    i = min(n,4)

    select case (i)
        case (0)
            stop
        case (1)
            sigma = p_cx_1_janev(Erel)
        case (2)
            sigma = p_cx_2_janev(Erel)
        case (3)
            sigma = p_cx_3_janev(Erel)
        case DEFAULT
            sigma = p_cx_n_janev(Erel, n)
    end select

end function p_cx_janev

function aljan1(energy, pcf) result(sigma)
    !+ Fit function for [[p_cx_1_2_janev]]
    real(Float64), intent(in) :: energy
        !+ Relative collision energy [eV]
    real(Float64), dimension(:), intent(in) :: pcf
        !+ Fit Coefficients
    real(Float64) :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64) :: aloge1, aloge, xjan, xcon, emin, emax
    integer :: i

    emin = pcf(2)
    emax = pcf(3)
    if((energy.lt.emin).or.(energy.gt.emax)) then
        sigma = 0
        return
    endif

    aloge1 = log(energy)
    aloge = aloge1
    xjan = pcf(4)
    do i = 5, 12
        xcon = pcf(i)*aloge
        xjan = xjan + xcon
        aloge = aloge * aloge1
    enddo
    sigma = exp(xjan)
end

function p_cx_1_2_janev(Erel) result(sigma)
    !+Calculates cross section for a proton-Hydrogen charge exchange interaction
    !+from the \(n=1\) state to the \(m=2\) state at energy `Erel`
    !+
    !+###Equation
    !+ $$H^+ + H(1) \rightarrow H(2) + H^+$$
    !+###References
    !+* Ref. 8 [[atomic_tables(module)]]
    real(Float64), intent(in) :: Erel
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(12), parameter :: c2s = [0.d0, 262.d0, 20000.d0, &
                                                     -13273.250877640001d0,   &
                                                      13175.766145199999d0,   &
                                                     -5683.9321578580002d0,   &
                                                      1386.3097801490001d0,   &
                                                     -208.97945613069999d0,   &
                                                      19.92976245274d0,       &
                                                     -1.1738005761570001d0,   &
                                                      0.039024228107669999d0, &
                                                     -0.00056062403399319998d0]

    real(Float64), dimension(12), parameter :: c2p = [0.d0, 19.d0 ,20000.d0,   &
                                                     -21.975719499349999d0,   &
                                                     -47.425022512600002d0,   &
                                                      36.280131405959999d0,   &
                                                     -14.23003075866d0,       &
                                                      3.2730902401440001d0,   &
                                                     -0.45579289122599997d0,  &
                                                      0.037735883474579998d0, &
                                                     -0.001707904867106d0,    &
                                                      3.251203344615d-5]

    real(Float64) :: e, sigma2s, sigma2p,logsig1,logsig2,slope

    e = Erel*1.d3

    sigma2s = aljan1(e,c2s)
    sigma2p = aljan1(e,c2p)

    if(e.gt.c2p(2)) then
        sigma = sigma2s + sigma2p
    else
        ! linearly extrapolate in log-space
        logsig2 = log(aljan1(c2p(2)+1,c2s) + aljan1(c2p(2)+1,c2p))
        logsig1 = log(aljan1(c2p(2),c2s) + aljan1(c2p(2),c2p))
        slope = (logsig2 - logsig1)/(log(c2p(2)+1) - log(c2p(2)))
        sigma = exp(slope*(log(e) - log(c2p(2))) + logsig1)
    endif

end

function p_cx_1_1_adas(Erel) result(sigma)
    !+Calculates cross section for a proton-Hydrogen charge exchange interaction
    !+from the \(n=1\) state to the \(m=1\) state at energy `Erel`
    !+
    !+###Equation
    !+ $$H^+ + H(1) \rightarrow H(1) + H^+$$
    !+###References
    !+* Ref. 4 [[atomic_tables(module)]]
    real(Float64), intent(in) :: Erel
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(7), parameter :: a = [-3.496092687d2, 4.724931484d2, &
                                                   -2.720493064d2, 8.158564625d1, &
                                                   -1.339790721d1, 1.138706949d0, &
                                                   -3.914774156d-2 ]
    real(Float64) :: e, ee, fac, l, p

    e = Erel*1.d3
    if(e.ge.1.d3) then
        ee = max(e,1.0)
        fac = 1.d0
    else
        ee = 1.0d3
        fac = Erel**(-0.2)
    endif

    l = log10(ee)

    p = a(1) + a(2)*l + a(3)*l**2.0 + a(4)*l**3.0 + &
        a(5)*l**4.0   + a(6)*l**5.0 + a(7)*l**6.0

    sigma = fac*(10.d0**p)

end function p_cx_1_1_adas

function p_cx_1_2_adas(Erel) result(sigma)
    !+Calculates cross section for a proton-Hydrogen charge exchange interaction
    !+from the \(n=1\) state to the \(m=2\) state at energy `Erel`
    !+
    !+###Equation
    !+ $$H^+ + H(1) \rightarrow H(2) + H^+$$
    !+###References
    !+* Ref. 4 [[atomic_tables(module)]]
    real(Float64), intent(in) :: Erel
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(9), parameter :: a = [-4.036239511d3, 6.941235312d3, &
                                                   -5.186974866d3, 2.194885201d3, &
                                                   -5.765960509d2, 9.653534186d1, &
                                                   -1.008066138d1, 6.010731909d-1,&
                                                   -1.567417031d-2 ]
    real(Float64) :: e, ee, fac, l, p

    e = Erel*1.d3
    if(e.ge.1.d3) then
        ee = max(e,1.0)
        fac = 1.d0
    else
        ee = 1.0d3
        fac = Erel**(0.5)
    endif

    l = log10(ee)

    p = a(1) + a(2)*l + a(3)*l**2.0 + a(4)*l**3.0 + &
        a(5)*l**4.0   + a(6)*l**5.0 + a(7)*l**6.0 + &
        a(8)*l**7.0   + a(9)*l**8.0

    sigma = fac*(10.d0**p)

end function p_cx_1_2_adas

function p_cx_1_3_adas(Erel) result(sigma)
    !+Calculates cross section for a proton-Hydrogen charge exchange interaction
    !+from the \(n=1\) state to the \(m=3\) state at energy `Erel`
    !+
    !+###Equation
    !+ $$H^+ + H(1) \rightarrow H(3) + H^+$$
    !+###References
    !+* Ref. 4 [[atomic_tables(module)]]
    real(Float64), intent(in) :: Erel
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(10), parameter :: a = [7.037287586d4, -1.479161477d5, &
                                                    1.370120708d5, -7.343180122d4, &
                                                    2.509832081d4, -5.674317075d3, &
                                                    8.487767749d2, -8.102284612d1, &
                                                    4.480007503d0, -1.093512342d-1 ]
    real(Float64) :: e, ee, fac, l, p

    e = Erel*1.d3
    if(e.ge.2.d3) then
        ee = max(e,1.0)
        fac = 1.d0
    else
        ee = 2.0d3
        fac = (Erel**(1.4))/2.8
    endif

    l = log10(ee)

    p = a(1) + a(2)*l + a(3)*l**2.0 + a(4)*l**3.0 + &
        a(5)*l**4.0   + a(6)*l**5.0 + a(7)*l**6.0 + &
        a(8)*l**7.0   + a(9)*l**8.0 + a(10)*l**9.0

    sigma = fac*(10.d0**p)

end function p_cx_1_3_adas

function p_cx_1_4_adas(Erel) result(sigma)
    !+Calculates cross section for a proton-Hydrogen charge exchange interaction
    !+from the \(n=1\) state to the \(m=4\) state at energy `Erel`
    !+
    !+###Equation
    !+ $$H^+ + H(1) \rightarrow H(4) + H^+$$
    !+###References
    !+* Ref. 4 [[atomic_tables(module)]]
    real(Float64), intent(in) :: Erel
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(10), parameter :: a = [6.826447557d4, -1.431980004d5, &
                                                    1.323968679d5, -7.083995050d4, &
                                                    2.417608863d4, -5.458418789d3, &
                                                    8.154875237d2, -7.776012846d1, &
                                                    4.295431731d0, -1.047567211d-1 ]
    real(Float64) :: e, ee, fac, l, p

    e = Erel*1.d3
    if(e.ge.2.d3) then
        ee = max(e,1.0)
        fac = 1.d0
    else
        ee = 2.0d3
        fac = (Erel**(2.0))/4.0
    endif

    l = log10(ee)

    p = a(1) + a(2)*l + a(3)*l**2.0 + a(4)*l**3.0 + &
        a(5)*l**4.0   + a(6)*l**5.0 + a(7)*l**6.0 + &
        a(8)*l**7.0   + a(9)*l**8.0 + a(10)*l**9.0

    sigma = fac*(10.d0**p)

end function p_cx_1_4_adas

function p_cx_1(Erel,m_max) result(sigma)
    !+Calculates an array of cross section for proton-Hydrogen charge exchange interactions
    !+from the \(n=1\) state to m = 1..`m_max` states at energy `Erel`
    !+
    !+@note Cross sections are normalized to the total cross sections calculated by
    !+[[p_cx_janev(proc)]]
    !+###Equation
    !+ $$H^+ + H(1) \rightarrow H(m=1..m_{max}) + H^+$$
    !+###References
    !+* Ref. 2 [[atomic_tables(module)]]
    !+* Ref. 4 [[atomic_tables(module)]]
    !+* Ref. 8 [[atomic_tables(module)]]
    real(Float64), intent(in)       :: Erel
        !+ Relative collision energy [keV/amu]
    integer, intent(in)             :: m_max
        !+ Number of `m` states to calculate
    real(Float64), dimension(m_max) :: sigma
        !+ Array of cross sections where the index refers to the `m`'th state [\(cm^2\)]

    integer :: i
    real(Float64) :: norm_fac

    sigma = 0.d0
    do i=1,m_max
        select case (i)
            case (1)
                if(Erel.le.2.0) then
                    sigma(1) = p_cx_janev(Erel, 1)
                else
                    sigma(1) = p_cx_1_1_adas(Erel)
                endif
            case (2)
                if(Erel.le.2.0) then
                    sigma(2) = p_cx_1_2_janev(Erel)
                else
                    sigma(2) = p_cx_1_2_adas(Erel)
                endif
            case (3)
                sigma(3) = p_cx_1_3_adas(Erel)
            case (4)
                sigma(4) = p_cx_1_4_adas(Erel)
            case DEFAULT
                sigma(i) = 0.d0
        end select
    enddo

    !Normalize to Janev to be consistent with other n levels (p_cx_2/3/...)
    norm_fac = p_cx_janev(Erel, 1)/sum(sigma)
    sigma = norm_fac*sigma

end function p_cx_1

function p_cx_2_2_adas(Erel) result(sigma)
    !+Calculates cross section for a proton-Hydrogen charge exchange interaction
    !+from the \(n=2\) state to the \(m=2\) state at energy `Erel`
    !+
    !+###Equation
    !+ $$H^+ + H(2) \rightarrow H(2) + H^+$$
    !+###References
    !+* Ref. 4 [[atomic_tables(module)]]
    real(Float64), intent(in) :: Erel
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(11), parameter :: a2s = [-1.896015167d6, 4.431727330d6, &
                                                      -4.627815357d6, 2.843068107d6, &
                                                      -1.137952956d6, 3.100801094d5, &
                                                      -5.825744660d4, 7.452319142d3, &
                                                      -6.212350647d2, 3.047712749d1, &
                                                      -6.682658463d-1 ]

    real(Float64), dimension(11), parameter :: a2p = [-1.614213508d5, 3.772469288d5, &
                                                      -3.924736424d5, 2.393127027d5, &
                                                      -9.470300966d4, 2.541276100d4, &
                                                      -4.682860453d3, 5.851219013d2, &
                                                      -4.744504549d1, 2.254460913d0, &
                                                      -4.767235839d-2 ]
    real(Float64), parameter :: n = 2.d0

    real(Float64) :: e, ee, fac, l, sigma2s, sigma2p

    e = Erel * 1.d3 * n**2.0
    if(Erel.le.1.5d2) then
        ee = max(e, 1.d3)
        fac = 1.d0
    else
        ee = 1.5e5 * n**2.d0
        fac = 2.d15 * ((e*1.d-3)**(-5.5))
    endif

    l = log10(ee)

    sigma2s = a2s(1) + a2s(2)*l + a2s(3)*l**2.0 + a2s(4)*l**3.0 + &
              a2s(5)*l**4.0     + a2s(6)*l**5.0 + a2s(7)*l**6.0 + &
              a2s(8)*l**7.0     + a2s(9)*l**8.0 + a2s(10)*l**9.0 + a2s(11)*l**10.0
    sigma2s = 10.d0**(sigma2s)

    sigma2p = a2p(1) + a2p(2)*l + a2p(3)*l**2.0 + a2p(4)*l**3.0 + &
              a2p(5)*l**4.0     + a2p(6)*l**5.0 + a2p(7)*l**6.0 + &
              a2p(8)*l**7.0     + a2p(9)*l**8.0 + a2p(10)*l**9.0 + a2p(11)*l**10.0
    sigma2p = 10.d0**(sigma2p)

    sigma = fac*(0.25*sigma2s + 0.75*sigma2p)

end function p_cx_2_2_adas

function p_cx_2_3_adas(Erel) result(sigma)
    !+Calculates cross section for a proton-Hydrogen charge exchange interaction
    !+from the \(n=2\) state to the \(m=3\) state at energy `Erel`
    !+
    !+###Equation
    !+ $$H^+ + H(2) \rightarrow H(3) + H^+$$
    !+###References
    !+* Ref. 4 [[atomic_tables(module)]]
    real(Float64), intent(in) :: Erel
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(11), parameter :: a2s = [-3.513030327d5, 9.281116596d5, &
                                                      -1.086843398d6, 7.437325055d5, &
                                                      -3.296609685d5, 9.897503768d4, &
                                                      -2.039707143d4, 2.850670244d3, &
                                                      -2.587092857d2, 1.377382945d1, &
                                                      -3.268306303d-1 ]

    real(Float64), dimension(11), parameter :: a2p = [-1.901264631d5, 5.124716103d5, &
                                                      -6.101921504d5, 4.234717934d5, &
                                                      -1.899866398d5, 5.764464326d4, &
                                                      -1.199087959d4, 1.689900512d3, &
                                                      -1.545334374d2, 8.285001228d0, &
                                                      -1.978656474d-1 ]

    real(Float64), parameter :: n = 2.d0

    real(Float64) :: ee, l, sigma2s, sigma2p

    ee = max(Erel * 1.d3 * n**2.d0, 1.d3)

    l = log10(ee)

    sigma2s = a2s(1) + a2s(2)*l + a2s(3)*l**2.0 + a2s(4)*l**3.0 + &
              a2s(5)*l**4.0     + a2s(6)*l**5.0 + a2s(7)*l**6.0 + &
              a2s(8)*l**7.0     + a2s(9)*l**8.0 + a2s(10)*l**9.0 + a2s(11)*l**10.0
    sigma2s = 10.d0**(sigma2s)

    sigma2p = a2p(1) + a2p(2)*l + a2p(3)*l**2.0 + a2p(4)*l**3.0 + &
              a2p(5)*l**4.0     + a2p(6)*l**5.0 + a2p(7)*l**6.0 + &
              a2p(8)*l**7.0     + a2p(9)*l**8.0 + a2p(10)*l**9.0 + a2p(11)*l**10.0
    sigma2p = 10.d0**(sigma2p)

    sigma = (0.25*sigma2s + 0.75*sigma2p)

end function p_cx_2_3_adas

subroutine m_spread(n, m_max, sigma_tot, sigma)
    !+ Spreads the total charge exchange cross section, `sigma_tot`,
    !+ among the non-filled m states of `sigma` according to an exponential
    integer, intent(in)                            :: n
        !+ Initial atomic energy level/state
    integer, intent(in)                            :: m_max
        !+ Number of m states in `sigma`
    real(Float64), intent(in)                      :: sigma_tot
        !+ Amount of "cross section" to spread about the non-filled m state of sigma
    real(Float64), dimension(m_max), intent(inout) :: sigma
        !+ Array of cross sections from the `n` state to m=1..`m_max` [\(cm^2\)]
    real(Float64) :: En, Em
    real(Float64) :: norm_fac
    real(Float64), dimension(m_max) :: sigma_m
    integer :: m

    sigma_m = 0.d0
    En = 13.6/(real(n)**2.0)
    do m=1,m_max
        Em = 13.6/(real(m)**2.0)
        if(sigma(m).eq.0.d0) then
            sigma_m(m) = (sigma_tot/sqrt(2.0*PI))*exp(-0.5*(En-Em)**2.0)
        endif
    enddo

    norm_fac = sigma_tot/sum(sigma_m)
    do m=1,m_max
        if(sigma(m).eq.0.d0) sigma(m) = sigma_m(m)*norm_fac
        if(sigma(m).ne.sigma(m)) sigma(m) = 0.d0
    enddo

end subroutine m_spread

function p_cx_2(Erel,m_max) result(sigma)
    !+Calculates an array of cross sections for proton-Hydrogen charge exchange interactions
    !+from the \(n=2\) state to m = 1..`m_max` states at energy `Erel`
    !+
    !+@note
    !+Cross sections are normalized to the total cross sections calculated by
    !+[[p_cx_janev(proc)]].
    !+
    !+@note
    !+Cross sections for the \(n=2 \rightarrow m=1\) states are calculated via
    !+equivalence principle using [[p_cx_1_2_adas(proc)]].
    !+
    !+@note
    !+Cross Sections for \(m \geq 4\) are calculated by "spreading" their
    !+expected total cross sections among the \(m \geq 4\) states.
    !+
    !+###Equation
    !+ $$H^+ + H(2) \rightarrow H(m=1..m_{max}) + H^+$$
    !+###References
    !+* Ref. 2 [[atomic_tables(module)]]
    !+* Ref. 4 [[atomic_tables(module)]]
    real(Float64), intent(in)       :: Erel
        !+ Relative collision energy [keV/amu]
    integer, intent(in)             :: m_max
        !+ Number of `m` states to calculate
    real(Float64), dimension(m_max) :: sigma
        !+ Array of cross sections where the index refers to the `m`'th state [\(cm^2\)]

    real(Float64), parameter :: n2 = 4.d0
    integer :: i
    real(Float64) :: En, Em, sigma_n, norm_fac

    sigma = 0.d0
    do i=1,min(m_max,3)
        select case (i)
            case (1)
                sigma(1) = p_cx_1_2_adas(Erel*n2)/n2
            case (2)
                sigma(2) = p_cx_2_2_adas(Erel)
            case (3)
                sigma(3) = p_cx_2_3_adas(Erel)
        end select
    enddo
    sigma_n = max(p_cx_janev(Erel, 2) - sum(sigma), 0.d0)

    call m_spread(2,m_max,sigma_n,sigma)

    norm_fac = p_cx_janev(Erel, 2)/sum(sigma)
    sigma = sigma*norm_fac

end function p_cx_2

function p_cx_3_2_adas(Erel) result(sigma)
    !+Calculates cross section for a proton-Hydrogen charge exchange interaction
    !+from the \(n=3\) state to the \(m=2\) state at energy `Erel`
    !+
    !+###Equation
    !+ $$H^+ + H(3) \rightarrow H(2) + H^+$$
    !+###References
    !+* Ref. 4 [[atomic_tables(module)]]
    real(Float64), intent(in) :: Erel
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(11), parameter :: a = [-1.149224555d6, 2.750368877d6, &
                                                    -2.942222842d6, 1.852584954d6, &
                                                    -7.603284323d5, 2.125284465d5, &
                                                    -4.097580431d4, 5.380901722d3, &
                                                    -4.606297192d2, 2.321345254d1, &
                                                    -5.230186707d-1 ]

    real(Float64), parameter :: n = 3.0
    real(Float64) :: ee, fac, l, p

    if(Erel.lt.90.0) then
        ee = max(Erel * 1.d3 * n**2.0, 1.d3) !keV to eV
        fac = 1.d0
    else
        ee = 90.0 * 1.d3 * n**2.0
        fac = 1.d16 * (Erel*n**2.0)**(-5.5)
    endif

    l = log10(ee)

    p = a(1) + a(2)*l + a(3)*l**2.0 + a(4)*l**3.0 + &
        a(5)*l**4.0   + a(6)*l**5.0 + a(7)*l**6.0 + &
        a(8)*l**7.0   + a(9)*l**8.0 + a(10)*l**9.0 + a(11)*l**10.d0

    sigma = fac*(10.d0**p)

end function p_cx_3_2_adas

function p_cx_3_3_adas(Erel) result(sigma)
    !+Calculates cross section for a proton-Hydrogen charge exchange interaction
    !+from the \(n=3\) state to the \(m=3\) state at energy `Erel`
    !+
    !+###Equation
    !+ $$H^+ + H(3) \rightarrow H(3) + H^+$$
    !+###References
    !+* Ref. 4 [[atomic_tables(module)]]
    real(Float64), intent(in) :: Erel
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(10), parameter :: a = [-4.302808608d4, 9.499298161d4, &
                                                    -9.264698488d4, 5.236947172d4, &
                                                    -1.890479538d4, 4.519068626d3, &
                                                    -7.152485009d2, 7.227063167d1, &
                                                    -4.230036444d0, 1.092702525d-1 ]

    real(Float64), parameter :: n = 3.0
    real(Float64) :: ee, fac, l, p

    if(Erel.lt.90.0) then
        ee = max(Erel * 1.d3 * n**2.0, 1.d3) !keV to eV
        fac = 1.d0
    else
        ee = 90.0 * 1.d3 * n**2
        fac = 0.85d16 *(Erel*n**2.0)**(-5.5)
    endif

    l = log10(ee)

    p = a(1) + a(2)*l + a(3)*l**2.0 + a(4)*l**3.0 + &
        a(5)*l**4.0   + a(6)*l**5.0 + a(7)*l**6.0 + &
        a(8)*l**7.0   + a(9)*l**8.0 + a(10)*l**9.0

    sigma = fac*(10.d0**p)

end function p_cx_3_3_adas

function p_cx_3_4_adas(Erel) result(sigma)
    !+Calculates cross section for a proton-Hydrogen charge exchange interaction
    !+from the \(n=3\) state to the \(m=4\) state at energy `Erel`
    !+
    !+###Equation
    !+ $$H^+ + H(3) \rightarrow H(4) + H^+$$
    !+###References
    !+* Ref. 4 [[atomic_tables(module)]]
    real(Float64), intent(in) :: Erel
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(9), parameter :: a = [ 1.705303425d4,-3.316878090d4, &
                                                    2.792556433d4,-1.330264490d4, &
                                                    3.921666688d3,-7.327555138d2, &
                                                    8.476342861d1,-5.551987930d0, &
                                                    1.577120745d-1 ]
    real(Float64), parameter :: n = 3.0
    real(Float64) :: ee, fac, l, p

    if(Erel.lt.90.0) then
        ee = max(Erel * 1.d3 * n**2.0, 1.d3) !keV to eV
        fac = 1.d0
    else
        ee = 90.0 * 1.d3 * n**2.0
        fac = 0.82d16 *(Erel*n**2.0)**(-5.5)
    endif

    l = log10(ee)

    p = a(1) + a(2)*l + a(3)*l**2.0 + a(4)*l**3.0 + &
        a(5)*l**4.0   + a(6)*l**5.0 + a(7)*l**6.0 + &
        a(8)*l**7.0   + a(9)*l**8.0

    sigma = fac*(10.d0**p)

end function p_cx_3_4_adas

function p_cx_3_5_adas(Erel) result(sigma)
    !+Calculates cross section for a proton-Hydrogen charge exchange interaction
    !+from the \(n=3\) state to the \(m=5\) state at energy `Erel`
    !+
    !+###Equation
    !+ $$H^+ + H(3) \rightarrow H(5) + H^+$$
    !+###References
    !+* Ref. 4 [[atomic_tables(module)]]
    real(Float64), intent(in) :: Erel
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(11), parameter :: a = [-2.786268232d2, 4.269683825d4, &
                                                    -8.973561028d4, 8.365732310d4, &
                                                    -4.524587937d4, 1.563630402d4, &
                                                    -3.580391824d3, 5.432527332d2, &
                                                    -5.267599631d1, 2.962329657d0, &
                                                    -7.362649692d-2 ]
    real(Float64), parameter :: n = 3.0
    real(Float64) :: ee, fac, l, p

    ee = max(Erel * 1.d3 * n**2.0, 1.d3) !keV to eV
    fac = 1.d0

    l = log10(ee)

    p = a(1) + a(2)*l + a(3)*l**2.0 + a(4)*l**3.0 + &
        a(5)*l**4.0   + a(6)*l**5.0 + a(7)*l**6.0 + &
        a(8)*l**7.0   + a(9)*l**8.0 + a(10)*l**9.0 + a(11)*l**10.0

    sigma = fac*(10.d0**p)

end function p_cx_3_5_adas

function p_cx_3_6inf_adas(Erel) result(sigma)
    !+Calculates total cross section for a proton-Hydrogen charge exchange interaction
    !+from the \(n=3\) state to \(\forall \; m \geq 6\) states at energy `Erel`
    !+
    !+###Equation
    !+ $$H^+ + H(3) \rightarrow H(\forall \; m \geq 6) + H^+$$
    !+###References
    !+* Ref. 4 [[atomic_tables(module)]]
    real(Float64), intent(in) :: Erel
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(11), parameter :: a = [ 7.146969470d5,-1.665413326d6, &
                                                     1.735840441d6,-1.065792786d6, &
                                                     4.269334710d5,-1.165954977d5, &
                                                     2.198700496d4,-2.827160468d3, &
                                                     2.372409350d2,-1.173264972d1, &
                                                     2.596865877d-1 ]
    real(Float64), parameter :: n = 3.0
    real(Float64) :: ee, fac, l, p

    if(Erel.lt.90.0) then
        ee = max(Erel * 1.d3 * n**2.0, 1.d3) !keV to eV
        fac = 1.d0
    else
        ee = 90.0 * 1.d3 * n**2.0
        fac = 2.d20 *(Erel*n**2.0)**(-7.0)
    endif

    l = log10(ee)

    p = a(1) + a(2)*l + a(3)*l**2.0 + a(4)*l**3.0 + &
        a(5)*l**4.0   + a(6)*l**5.0 + a(7)*l**6.0 + &
        a(8)*l**7.0   + a(9)*l**8.0 + a(10)*l**9.0 + a(11)*l**10.0

    sigma = fac*(10.d0**p)

end function p_cx_3_6inf_adas

function p_cx_3(Erel,m_max) result(sigma)
    !+Calculates an array of cross sections for proton-Hydrogen charge exchange interactions
    !+from the \(n=3\) state to m = 1..`m_max` states at energy `Erel`
    !+
    !+@note
    !+Cross sections are normalized to the total cross sections calculated by
    !+[[p_cx_janev(proc)]].
    !+
    !+@note
    !+Cross sections for the \(n=3 \rightarrow m=1\) states are calculated via
    !+equivalence principle using [[p_cx_1_3_adas(proc)]].
    !+
    !+@note
    !+Cross Sections for \(m \geq 6\) are calculated by "spreading" their
    !+expected total cross sections among the \( m \geq 6\) states.
    !+
    !+###Equation
    !+ $$H^+ + H(3) \rightarrow H(m=1..m_{max}) + H^+$$
    !+###References
    !+* Ref. 4 [[atomic_tables(module)]]
    real(Float64), intent(in)       :: Erel
        !+ Relative collision energy [keV/amu]
    integer, intent(in)             :: m_max
        !+ Number of `m` states to calculate
    real(Float64), dimension(m_max) :: sigma
        !+ Array of cross sections where the index refers to the `m`'th state [\(cm^2\)]

    real(Float64), parameter :: n2 = 9.d0
    real(Float64) :: eb, En, Em, sigma_m6, norm_fac
    real(Float64), dimension(m_max) :: sigma1

    sigma = 0.d0
    sigma1 = 0.d0
    sigma1 = p_cx_1(Erel*n2,m_max)
    sigma(1) = p_cx_1_3_adas(Erel*n2)/n2

    sigma(2) = p_cx_3_2_adas(Erel)
    sigma(3) = p_cx_3_3_adas(Erel)
    sigma(4) = p_cx_3_4_adas(Erel)

    if(m_max.ge.5) then
        sigma(5) = p_cx_3_5_adas(Erel)
    endif

    if(m_max.ge.6) then
        sigma_m6 = p_cx_3_6inf_adas(Erel)
        call m_spread(3, m_max, sigma_m6, sigma)
    endif

    norm_fac = p_cx_janev(Erel, 3)/sum(sigma)
    sigma = sigma*norm_fac

end function p_cx_3

function p_cx_n(Erel, n, m_max) result(sigma)
    !+Calculates an array of cross sections for proton-Hydrogen charge exchange interactions
    !+from the `n` state to m = 1..`m_max` states at energy `Erel`
    !+
    !+@note
    !+Cross sections are normalized to the total cross sections calculated by
    !+[[p_cx_janev(proc)]].
    !+
    !+@note
    !+Cross sections for some transitions are calculated via the equivalence principle or
    !+by "spreading" their expected total cross sections among the non-filled m states.
    !+
    !+###Equation
    !+ $$H^+ + H(n) \rightarrow H(m=1..m_{max}) + H^+$$
    !+###References
    !+* Ref. 2 [[atomic_tables(module)]]
    !+* Ref. 4 [[atomic_tables(module)]]
    !+* Ref. 8 [[atomic_tables(module)]]
    real(Float64), intent(in)       :: Erel
        !+ Relative collision energy [keV/amu]
    integer, intent(in)             :: n
        !+ Initial atomic energy level/state
    integer, intent(in)             :: m_max
        !+ Number of `m` states to calculate
    real(Float64), dimension(m_max) :: sigma
        !+ Array of cross sections where the index refers to the `m`'th state [\(cm^2\)]

    real(Float64), dimension(m_max) :: sigma2,sigma3
    real(Float64) :: sigma_n,e,norm_fac

    sigma = 0.d0
    select case (n)
        case (0)
            stop
        case (1)
            sigma = p_cx_1(Erel,m_max)
            return
        case (2)
            sigma = p_cx_2(Erel,m_max)
            return
        case (3)
            sigma = p_cx_3(Erel,m_max)
            return
        case (4)
            e = Erel*n**2.0
            sigma2 = p_cx_2(e/(2.0**2.0),m_max)
            sigma(1) = p_cx_1_4_adas(e/(1.0**2.0))*(1.d0/n)**2.0
            sigma(2) = sigma2(4)*(2.d0/n)**2.0
            sigma(3) = p_cx_3_4_adas(e/(3.0**2.0))*(3.d0/n)**2.0
        case (5)
            e = Erel*n**2.0
            sigma2 = p_cx_2(e/(2.0**2.0),m_max)
            sigma(2) = sigma2(5)*(2.d0/n)**2.0
            sigma(3) = p_cx_3_5_adas(e/(3.0**2.0))*(3.d0/n)**2.0
        case (6)
            e = Erel*n**2.0
            sigma2 = p_cx_2(e/(2.0**2.0),m_max)
            sigma(2) = sigma2(6)*(2.d0/n)**2.0
            sigma3 = p_cx_3(e/(3.0**2.0),m_max)*(3.d0/n)**2.0
            sigma(3) = sigma3(6)
        case DEFAULT
    end select

    sigma_n = max(p_cx_janev(Erel,n) - sum(sigma),0.0)
    call m_spread(n, m_max, sigma_n, sigma)

    norm_fac = p_cx_janev(Erel, n)/sum(sigma)
    sigma = norm_fac*sigma

end function p_cx_n

function p_cx_n_m(Erel, n, m) result(sigma)
    !+Calculates cross section for a proton-Hydrogen charge exchange interaction
    !+from the `n` state to the `m` state at energy `Erel`
    !+
    !+###Equation
    !+ $$H^+ + H(n) \rightarrow H(m) + H^+$$
    !+###References
    !+* Ref. 2 [[atomic_tables(module)]]
    !+* Ref. 4 [[atomic_tables(module)]]
    !+* Ref. 8 [[atomic_tables(module)]]
    real(Float64), intent(in) :: Erel
        !+ Relative collision energy [keV/amu]
    integer, intent(in)       :: n
        !+ Initial atomic energy level/state
    integer, intent(in)       :: m
        !+ Final atomic energy level/state
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    integer :: m_max = 12
    real(Float64), dimension(12) :: sigma_m

    sigma_m = p_cx_n(Erel, n, m_max)
    if(m.le.0) then
        sigma = sum(sigma_m)
    else
        sigma = sigma_m(m)
    endif

end function p_cx_n_m

function p_cx(Erel, n_max, m_max) result(sigma)
    !+Calculates a matrix of cross sections for proton-Hydrogen charge exchange interactions
    !+from the \(n=1..n_{max} \rightarrow m=1..m_{max}\) states at energy `Erel`
    !+
    !+###Equation
    !+ $$H^+ + H(n=1..n_{max}) \rightarrow H(m=1..m_{max}) + H^+$$
    !+###References
    !+* Ref. 2 [[atomic_tables(module)]]
    !+* Ref. 4 [[atomic_tables(module)]]
    !+* Ref. 8 [[atomic_tables(module)]]
    real(Float64), intent(in)             :: Erel
        !+ Relative collision energy [keV/amu]
    integer, intent(in)                   :: n_max
        !+ Number of initial atomic energy levels/states
    integer, intent(in)                   :: m_max
        !+ Number of final atomic energy levels/states
    real(Float64), dimension(n_max,m_max) :: sigma
        !+ Matrix of cross sections where the subscripts correspond
        !+ to the \(n \rightarrow m\) transitions: p_cx[n,m] [\(cm^2\)]

    real(Float64), dimension(12,12) :: sigma_full

    integer :: n, m

    do n=1,12
        sigma_full(n,:) = p_cx_n(Erel, n, 12)
    enddo

    sigma = sigma_full(1:n_max,1:m_max)

end function p_cx

!proton-Hydrogen impact ionization
function p_ioniz_1_janev(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact ionization interaction
    !+from the \(n=1\) state at energy `eb`
    !+
    !+###Equation
    !+ $$H^+ + H(1) \rightarrow H^+ + H^+ + e$$
    !+###References
    !+* Eq. 40 and Table 8 in Ref. 2 [[atomic_tables(module)]]
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(8), parameter :: b = [2.0160d-3, 3.7154d0,  &
                                                   3.9890d-2, 3.1413d-1, &
                                                   2.1254d0,  6.3990d3,  &
                                                   6.1897d1,  9.2731d3 ]
        !+ Fitting Parameters from Table 8 in Ref. 2

    real(Float64), parameter :: n2 = 1.d0
    real(Float64) :: Ehat
    real(Float64) :: p1, p2, p3

    Ehat = eb*n2
    p1 = b(1)*(n2)**2.0
    p2 = Ehat**b(2) * exp(-b(3)*Ehat) / (1.d0 + b(4)*Ehat**b(5))
    p3 = (b(6)* exp(-b(7)/Ehat) *log(1.d0  +b(8)*Ehat) ) /Ehat
    sigma = 1.0d-16 * p1 * (p2 + p3)

end function p_ioniz_1_janev

function p_ioniz_2_omullane(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact ionization interaction
    !+from the \(n=2\) state at energy `eb`
    !+
    !+###Equation
    !+ $$H^+ + H(2) \rightarrow H^+ + H^+ + e$$
    !+###References
    !+* Eq. 5 and Table 1 in Ref. 3 [[atomic_tables(module)]]
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]
    real(Float64), dimension(8), parameter :: b = [3.9330d-3, 1.8188d0,  &
                                                   1.8870d-2, 6.7489d-3, &
                                                   1.3768d0, 6.8852d2,   &
                                                   9.6435d1, 5.6515d23 ]
        !+ Fitting Parameters from Table 1 in Ref. 3

    real(Float64), parameter :: n2 = 4.d0
    real(Float64) :: Ehat
    real(Float64) :: p1, p2, p3

    Ehat = eb*n2
    p1 = b(1)*(n2)**2.0
    p2 = Ehat**b(2) * exp(-b(3)*Ehat) / (1.d0 + b(4)*Ehat**b(5))
    p3 = (b(6)* exp(-b(7)/Ehat) *log(1.d0  +b(8)*Ehat) ) /Ehat
    sigma = 1.0d-16 * p1 * (p2 + p3)

end function p_ioniz_2_omullane

function p_ioniz_3_omullane(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact ionization interaction
    !+from the \(n=3\) state at energy `eb`
    !+
    !+###Equation
    !+ $$H^+ + H(3) \rightarrow H^+ + H^+ + e$$
    !+###References
    !+* Eq. 5 and Table 1 in Ref. 3 [[atomic_tables(module)]]
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]
    real(Float64), dimension(8), parameter :: b = [1.1076d-2, 1.6197d0,  &
                                                   6.7154d-3, 5.1188d-3, &
                                                   1.8549d0,  2.3696d2,  &
                                                   7.8286d1,  1.0926d23 ]
        !+ Fitting Parameters from Table 1 in Ref. 3

    real(Float64), parameter :: n2 = 9.d0
    real(Float64) :: Ehat
    real(Float64) :: p1, p2, p3

    Ehat = eb*n2
    p1 = b(1)*(n2)**2.0
    p2 = Ehat**b(2) * exp(-b(3)*Ehat) / (1.d0 + b(4)*Ehat**b(5))
    p3 = (b(6)* exp(-b(7)/Ehat) *log(1.d0  +b(8)*Ehat) ) /Ehat
    sigma = 1.0d-16 * p1 * (p2 + p3)

end function p_ioniz_3_omullane

function p_ioniz_4_omullane(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact ionization interaction
    !+from the \(n=4\) state at energy `eb`
    !+
    !+###Equation
    !+ $$H^+ + H(4) \rightarrow H^+ + H^+ + e$$
    !+###References
    !+* Eq. 5 and Table 1 in Ref. 3 [[atomic_tables(module)]]
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]
    real(Float64), dimension(8), parameter :: b = [1.1033d-2, 1.6281d0,  &
                                                   5.5955d-3, 7.2023d-3, &
                                                   1.7358d0,  2.2755d2,  &
                                                   8.6339d1,  3.9151d29 ]
        !+ Fitting Parameters from Table 1 in Ref. 3

    real(Float64), parameter :: n2 = 16.d0
    real(Float64) :: Ehat
    real(Float64) :: p1, p2, p3

    Ehat = eb*n2
    p1 = b(1)*(n2)**2.0
    p2 = Ehat**b(2) * exp(-b(3)*Ehat) / (1.d0 + b(4)*Ehat**b(5))
    p3 = (b(6)* exp(-b(7)/Ehat) *log(1.d0  +b(8)*Ehat) ) /Ehat
    sigma = 1.0d-16 * p1 * (p2 + p3)

end function p_ioniz_4_omullane

function p_ioniz_5_omullane(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact ionization interaction
    !+from the \(n=5\) state at energy `eb`
    !+
    !+###Equation
    !+ $$H^+ + H(5) \rightarrow H^+ + H^+ + e$$
    !+###References
    !+* Eq. 5 and Table 1 in Ref. 3 [[atomic_tables(module)]]
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]
    real(Float64), dimension(8), parameter :: b = [1.1297d-2, 1.8685d0,  &
                                                   1.5038d-2, 1.1195d-1, &
                                                   1.0538d0,  8.6096d2,  &
                                                   8.9939d1,  1.9249d4 ]
        !+ Fitting Parameters from Table 1 in Ref. 3

    real(Float64), parameter :: n2 = 25.d0
    real(Float64) :: Ehat
    real(Float64) :: p1, p2, p3

    Ehat = eb*n2
    p1 = b(1)*(n2)**2.0
    p2 = Ehat**b(2) * exp(-b(3)*Ehat) / (1.d0 + b(4)*Ehat**b(5))
    p3 = (b(6)* exp(-b(7)/Ehat) *log(1.d0  +b(8)*Ehat) ) /Ehat
    sigma = 1.0d-16 * p1 * (p2 + p3)

end function p_ioniz_5_omullane

function p_ioniz_n(eb,n) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact ionization interaction
    !+from the `n`th state at energy `eb`
    !+
    !+###Equation
    !+ $$H^+ + H(n) \rightarrow H^+ + H^+ + e$$
    !+###References
    !+* Eq. 40 and Table 8 in Ref. 2 for \(n=1\) [[atomic_tables(module)]]
    !+* Eq. 5 and Table 1 in Ref. 3 for \(n \geq 2\) [[atomic_tables(module)]]
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    integer, intent(in)       :: n
        !+ Initial atomic energy level/state
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    select case (n)
        case (0)
            stop
        case (1)
            sigma = p_ioniz_1_janev(eb)
        case (2)
            sigma = p_ioniz_2_omullane(eb)
        case (3)
            sigma = p_ioniz_3_omullane(eb)
        case (4)
            sigma = p_ioniz_4_omullane(eb)
        case DEFAULT
            sigma = p_ioniz_5_omullane(eb)*(n/5.d0)**4
    end select

end function p_ioniz_n

function p_ioniz(eb,n_max) result(sigma)
    !+Calculates an array of cross sections for proton-Hydrogen impact ionization interactions
    !+from the n = 1..`n_max` state at energy `eb`
    !+
    !+###Equation
    !+ $$H^+ + H(n=1..n_{max}) \rightarrow H^+ + H^+ + e$$
    !+###References
    !+* Eq. 40 and Table 8 in Ref. 2 for \(n=1\) [[atomic_tables(module)]]
    !+* Eq. 5 and Table 1 in Ref. 3 for \(n \geq 2\) [[atomic_tables(module)]]
    real(Float64), intent(in)       :: eb
        !+ Relative collision energy [keV/amu]
    integer, intent(in)             :: n_max
        !+ Number of initial atomic energy level/state
    real(Float64), dimension(n_max) :: sigma
        !+ Array of cross sections where the index refers to the `n`'th state [\(cm^2\)]

    integer :: i

    do i=1,n_max
        sigma(i) = p_ioniz_n(eb,i)
    enddo

end function p_ioniz

!! proton-Hydrogen impact excitation
function p_excit_1_2_janev(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact excitation transition from
    !+the \(n=1\) state to the \(m=2\) state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(1) \rightarrow H^+ + H(2) $$
    !+
    !+###References
    !+* Eq. 29.b and Table 4 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]
    real(Float64), dimension(10), parameter :: a = [34.433d0, 8.5476d0,  &
                                                    7.8501d0, -9.2217d0, &
                                                    1.8020d-2, 1.6931d0, &
                                                    1.9422d-3, 2.9068d0, &
                                                    44.507d0, 0.56870d0 ]
        !+ Fitting parameters from Table 4 in Ref. 2

    sigma = 1.d-16 * a(1) * ( a(2)*exp(-a(3)*eb)/(eb**a(4)) + &
                     a(5)*exp(-a(6)/eb)/(1.+a(7)*eb**a(8))  + &
                     exp(-a(9)/eb)*log(1.+a(10)*eb)/eb )

end function p_excit_1_2_janev

function p_excit_1_3_janev(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact excitation transition from
    !+the \(n=1\) state to the \(m=3\) state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(1) \rightarrow H^+ + H(3) $$
    !+
    !+###References
    !+* Eq. 30 and Table 5 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]
    real(Float64), dimension(8), parameter :: b = [ 6.1950d0, 5.5162d-3,  &
                                                    0.29114d0, -4.5264d0, &
                                                    6.0311d0, -2.0679d0,  &
                                                    35.773d0, 0.54818d0 ]
        !+ Fitting parameters from Table 5 in Ref. 2

    sigma = 1.d-16 * b(1) * (b(2)*exp(-b(3)*eb)/ &
                     (eb**b(4)+b(5)*eb**b(6)) + &
                     exp(-b(7)/eb)*log(1.+b(8)*eb)/eb)

end function p_excit_1_3_janev

function p_excit_1_4_janev(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact excitation transition from
    !+the \(n=1\) state to the \(m=4\) state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(1) \rightarrow H^+ + H(4) $$
    !+
    !+###References
    !+* Eq. 30 and Table 5 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]
    real(Float64), dimension(8), parameter :: b = [2.0661d0, 5.1335d-4,  &
                                                   0.28953d0, -2.2849d0, &
                                                   0.11528d0, -4.8970d0, &
                                                   34.975d0, 0.91213d0 ]
        !+ Fitting parameters from Table 5 in Ref. 2

    sigma = 1.d-16 * b(1) * (b(2)*exp(-b(3)*eb)/ &
                     (eb**b(4)+b(5)*eb**b(6)) + &
                     exp(-b(7)/eb)*log(1.+b(8)*eb)/eb)

end function p_excit_1_4_janev

function p_excit_1_5_janev(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact excitation transition from
    !+the \(n=1\) state to the \(m=5\) state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(1) \rightarrow H^+ + H(5) $$
    !+
    !+###References
    !+* Eq. 30 and Table 5 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]
    real(Float64), dimension(8), parameter :: b = [1.2449d0, 3.0826d-4,   &
                                                   0.31063d0, -2.4161d0,  &
                                                   0.024664d0, -6.3726d0, &
                                                   32.291d0, 0.21176d0 ]
        !+ Fitting parameters from Table 5 in Ref. 2

    sigma = 1.d-16 * b(1) * (b(2)*exp(-b(3)*eb)/ &
                     (eb**b(4)+b(5)*eb**b(6)) + &
                     exp(-b(7)/eb)*log(1.+b(8)*eb)/eb)

end function p_excit_1_5_janev

function p_excit_1_6_janev(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact excitation transition from
    !+the \(n=1\) state to the \(m=6\) state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(1) \rightarrow H^+ + H(6) $$
    !+
    !+###References
    !+* Eq. 30 and Table 5 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]
    real(Float64), dimension(8), parameter :: b = [0.63771d0, 3.2949d-4,  &
                                                   0.25757d0, -2.2950d0,  &
                                                   0.050796d0, -5.5986d0, &
                                                   37.174d0, 0.39265d0 ]
        !+ Fitting parameters from Table 5 in Ref. 2

    sigma = 1.d-16 * b(1) * (b(2)*exp(-b(3)*eb)/ &
                     (eb**b(4)+b(5)*eb**b(6)) + &
                     exp(-b(7)/eb)*log(1.+b(8)*eb)/eb)

end function p_excit_1_6_janev

function p_excit_1_janev(eb, m_max) result(sigma)
    !+Calculates an array of cross sections for a proton-Hydrogen impact excitation transitions from
    !+the \(n=1\) state to the \(m=1..{m_max}\) state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(1) \rightarrow H^+ + H(m=1..m_{max}), m \gt 1 $$
    !+
    !+###References
    !+* Eq. 29.b and Table 4 in Ref. 2 for \(m = 2\) [[atomic_tables(module)]]
    !+* Eq. 30 and Table 5 in Ref. 2 for \(m = 3-6\) [[atomic_tables(module)]]
    !+* Eq. 31 and Table 5 in Ref. 2 for \(m \gt 6\) [[atomic_tables(module)]]
    !+
    real(Float64), intent(in)       :: eb
        !+ Relative collision energy [keV/amu]
    integer, intent(in)             :: m_max
        !+ Number of m states to calculate
    real(Float64), dimension(m_max) :: sigma
        !+ Array of cross sections where the m'th index refer to the transition
        !+ from \(n=1\) to m [\(cm^2\)]

    integer :: m
    sigma = 0.d0

    do m=1,m_max
        select case (m)
            case (1)
                sigma(1) = 0.d0
            case (2)
                sigma(2) = p_excit_1_2_janev(eb)
            case (3)
                sigma(3) = p_excit_1_3_janev(eb)
            case (4)
                sigma(4) = p_excit_1_4_janev(eb)
            case (5)
                sigma(5) = p_excit_1_5_janev(eb)
            case (6)
                sigma(6) = p_excit_1_6_janev(eb)
            case DEFAULT
                sigma(m) = p_excit_1_6_janev(eb)*(6.0/real(m))**3.0
        end select
    enddo

end function p_excit_1_janev

function p_excit_2_3_janev(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact excitation transition from
    !+the \(n=2\) state to the \(m=3\) state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(2) \rightarrow H^+ + H(3) $$
    !+
    !+###References
    !+* Eq. 32 and Table 6 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]
    real(Float64), dimension(6), parameter :: c = [394.51d0, 0.013597d0, &
                                                   0.16565d0, -0.8949d0, &
                                                   21.606d0,  0.62426d0  ]
        !+ Fitting parameters from Table 6 in Ref. 2

    sigma = 1.d-16 * c(1)*(c(2)*exp(-c(3)*eb)/(eb**c(4)) + &
                     exp(-c(5)/eb)*log(1.+c(6)*eb)/eb)

end function p_excit_2_3_janev

function p_excit_2_4_janev(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact excitation transition from
    !+the \(n=2\) state to the \(m=4\) state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(2) \rightarrow H^+ + H(4) $$
    !+
    !+###References
    !+* Eq. 32 and Table 6 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]
    real(Float64), dimension(6), parameter :: c = [50.744d0, 0.014398d0, &
                                                   0.31584d0, -1.4799d0, &
                                                   19.416d0,   4.0262d0  ]
        !+ Fitting parameters from Table 6 in Ref. 2

    sigma = 1.d-16 * c(1)*(c(2)*exp(-c(3)*eb)/(eb**c(4)) + &
                     exp(-c(5)/eb)*log(1.+c(6)*eb)/eb)

end function p_excit_2_4_janev

function p_excit_2_5_janev(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact excitation transition from
    !+the \(n=2\) state to the \(m=5\) state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(2) \rightarrow H^+ + H(5) $$
    !+
    !+###References
    !+* Eq. 32 and Table 6 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]
    real(Float64), dimension(6), parameter :: c = [18.264d0, 0.013701d0, &
                                                   0.31711d0, -1.4775d0, &
                                                   18.973d0,   2.9056d0  ]
        !+ Fitting parameters from Table 6 in Ref. 2

    sigma = 1.d-16 * c(1)*(c(2)*exp(-c(3)*eb)/(eb**c(4)) + &
                     exp(-c(5)/eb)*log(1.+c(6)*eb)/eb)

end function p_excit_2_5_janev

function p_excit_2_6_janev(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact excitation transition from
    !+the \(n=2\) state to the \(m=6\) state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(2) \rightarrow H^+ + H(6) $$
    !+
    !+###References
    !+* Eq. 33 and Table 6 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: A = 4.61d-1

    sigma = A*p_excit_2_5_janev(eb)

end function p_excit_2_6_janev

function p_excit_2_7_janev(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact excitation transition from
    !+the \(n=2\) state to the \(m=7\) state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(2) \rightarrow H^+ + H(7) $$
    !+
    !+###References
    !+* Eq. 33 and Table 6 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: A = 2.475d-1

    sigma = A*p_excit_2_5_janev(eb)

end function p_excit_2_7_janev

function p_excit_2_8_janev(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact excitation transition from
    !+the \(n=2\) state to the \(m=8\) state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(2) \rightarrow H^+ + H(8) $$
    !+
    !+###References
    !+* Eq. 33 and Table 6 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: A = 1.465d-1

    sigma = A*p_excit_2_5_janev(eb)

end function p_excit_2_8_janev

function p_excit_2_9_janev(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact excitation transition from
    !+the \(n=2\) state to the \(m=9\) state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(2) \rightarrow H^+ + H(9) $$
    !+
    !+###References
    !+* Eq. 33 and Table 6 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: A = 9.2d-2

    sigma = A*p_excit_2_5_janev(eb)

end function p_excit_2_9_janev

function p_excit_2_10_janev(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact excitation transition from
    !+the \(n=2\) state to the \(m=10\) state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(2) \rightarrow H^+ + H(10) $$
    !+
    !+###References
    !+* Eq. 33 and Table 6 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: A = 6.05d-2

    sigma = A*p_excit_2_5_janev(eb)

end function p_excit_2_10_janev

function p_excit_2_janev(eb, m_max) result(sigma)
    !+Calculates an array of cross sections for a proton-Hydrogen impact excitation transitions from
    !+the \(n=2\) state to the \(m=1..{m_max}\) state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(2) \rightarrow H^+ + H(m=1..m_{max}), m \gt 2$$
    !+
    !+###References
    !+* Eq. 32 and Table 6 in Ref. 2 for \(m \le 5\) [[atomic_tables(module)]]
    !+* Eq. 33 and Table 6 in Ref. 2 for \(m = 6-10\) [[atomic_tables(module)]]
    !+* Eq. 34 and Table 6 in Ref. 2 for \(m \gt 10\) [[atomic_tables(module)]]
    !+
    real(Float64), intent(in)       :: eb
        !+ Relative collision energy [keV/amu]
    integer, intent(in)             :: m_max
        !+ Number of m states to calculate
    real(Float64), dimension(m_max) :: sigma
        !+ Array of cross sections where the m'th index refer to the transition
        !+ from \(n=2\) to m [\(cm^2\)]

    integer :: m
    sigma = 0.d0

    do m=1,m_max
        select case (m)
            case (1)
                sigma(1) = 0.d0
            case (2)
                sigma(2) = 0.d0
            case (3)
                sigma(3) = p_excit_2_3_janev(eb)
            case (4)
                sigma(4) = p_excit_2_4_janev(eb)
            case (5)
                sigma(5) = p_excit_2_5_janev(eb)
            case (6)
                sigma(6) = p_excit_2_6_janev(eb)
            case (7)
                sigma(7) = p_excit_2_7_janev(eb)
            case (8)
                sigma(8) = p_excit_2_8_janev(eb)
            case (9)
                sigma(9) = p_excit_2_9_janev(eb)
            case (10)
                sigma(10) = p_excit_2_10_janev(eb)
            case DEFAULT
                sigma(m) = sigma(10)*(10.0/real(m))**3.0
        end select
    enddo

end function p_excit_2_janev

function p_excit_3_4_janev(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact excitation transition from
    !+the \(n=3\) state to the \(m=4\) state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(3) \rightarrow H^+ + H(4) $$
    !+
    !+###References
    !+* Eq. 35 and Table 7 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]
    real(Float64), dimension(6), parameter :: c = [1247.5d0,  0.068781d0, &
                                                   0.521176d0, -1.2722d0, &
                                                   11.319d0,    2.6235d0  ]
        !+ Fitting parameters from Table 7 in Ref. 2

    sigma = 1.d-16 * c(1)*(c(2)*exp(-c(3)*eb)/(eb**c(4)) + &
                     exp(-c(5)/eb)*log(1.+c(6)*eb)/eb)

end function p_excit_3_4_janev

function p_excit_3_5_janev(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact excitation transition from
    !+the \(n=3\) state to the \(m=5\) state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(3) \rightarrow H^+ + H(5) $$
    !+
    !+###References
    !+* Eq. 35 and Table 7 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]
    real(Float64), dimension(6), parameter :: c = [190.59d0, 0.073307d0, &
                                                   0.54177d0, -1.2894d0, &
                                                   11.096d0,   2.9098d0  ]
        !+ Fitting parameters from Table 7 in Ref. 2

    sigma = 1.d-16 * c(1)*(c(2)*exp(-c(3)*eb)/(eb**c(4)) + &
                     exp(-c(5)/eb)*log(1.+c(6)*eb)/eb)

end function p_excit_3_5_janev

function p_excit_3_6_janev(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact excitation transition from
    !+the \(n=3\) state to the \(m=6\) state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(3) \rightarrow H^+ + H(6) $$
    !+
    !+###References
    !+* Eq. 35 and Table 7 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]
    real(Float64), dimension(6), parameter :: c = [63.494d0, 0.077953d0, &
                                                   0.53461d0, -1.2881d0, &
                                                   11.507d0,   4.3417d0  ]
        !+ Fitting parameters from Table 7 in Ref. 2

    sigma = 1.d-16 * c(1)*(c(2)*exp(-c(3)*eb)/(eb**c(4)) + &
                     exp(-c(5)/eb)*log(1.+c(6)*eb)/eb)

end function p_excit_3_6_janev

function p_excit_3_7_janev(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact excitation transition from
    !+the \(n=3\) state to the \(m=7\) state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(3) \rightarrow H^+ + H(7) $$
    !+
    !+###References
    !+* Eq. 36 and Table 7 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: A = 4.67d-1

    sigma = A*p_excit_3_6_janev(eb)

end function p_excit_3_7_janev

function p_excit_3_8_janev(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact excitation transition from
    !+the \(n=3\) state to the \(m=8\) state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(3) \rightarrow H^+ + H(8) $$
    !+
    !+###References
    !+* Eq. 36 and Table 7 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: A = 2.545d-1

    sigma = A*p_excit_3_6_janev(eb)

end function p_excit_3_8_janev

function p_excit_3_9_janev(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact excitation transition from
    !+the \(n=3\) state to the \(m=9\) state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(3) \rightarrow H^+ + H(9) $$
    !+
    !+###References
    !+* Eq. 36 and Table 7 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: A = 1.54d-1

    sigma = A*p_excit_3_6_janev(eb)

end function p_excit_3_9_janev

function p_excit_3_10_janev(eb) result(sigma)
    !+Calculates cross section for a proton-Hydrogen impact excitation transition from
    !+the \(n=3\) state to the \(m=10\) state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(3) \rightarrow H^+ + H(10) $$
    !+
    !+###References
    !+* Eq. 36 and Table 7 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: A = 1.0d-1

    sigma = A*p_excit_3_6_janev(eb)

end function p_excit_3_10_janev

function p_excit_3_janev(eb, m_max) result(sigma)
    !+Calculates an array of cross sections for proton-Hydrogen impact excitation transitions from
    !+the \(n=3\) state to the \(m=1..{m_max}\) state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(3) \rightarrow H^+ + H(m=1..m_{max}), m \gt 3 $$
    !+
    !+###References
    !+* Eq. 35 and Table 7 in Ref. 2 for \(m \le 6\) [[atomic_tables(module)]]
    !+* Eq. 36 and Table 7 in Ref. 2 for \(m = 7-10\) [[atomic_tables(module)]]
    !+* Eq. 37 and Table 7 in Ref. 2 for \(m \gt 10\) [[atomic_tables(module)]]
    !+
    real(Float64), intent(in)       :: eb
        !+ Relative collision energy [keV/amu]
    integer, intent(in)             :: m_max
        !+ Number of m states to calculate
    real(Float64), dimension(m_max) :: sigma
        !+ Array of cross sections where the m'th index refer to the transition
        !+ from \(n=3\) to m [\(cm^2\)]

    integer :: m
    sigma = 0.d0

    do m=1,m_max
        select case (m)
            case (1)
                sigma(1) = 0.d0
            case (2)
                sigma(2) = 0.d0
            case (3)
                sigma(3) = 0.d0
            case (4)
                sigma(4) = p_excit_3_4_janev(eb)
            case (5)
                sigma(5) = p_excit_3_5_janev(eb)
            case (6)
                sigma(6) = p_excit_3_6_janev(eb)
            case (7)
                sigma(7) = p_excit_3_7_janev(eb)
            case (8)
                sigma(8) = p_excit_3_8_janev(eb)
            case (9)
                sigma(9) = p_excit_3_9_janev(eb)
            case (10)
                sigma(10) = p_excit_3_10_janev(eb)
            case DEFAULT
                sigma(m) = sigma(10)*(10.0/real(m))**3.0
        end select
    enddo

end function p_excit_3_janev

function p_excit_n(eb, n, m_max) result(sigma)
    !+Calculates an array of cross sections for a proton-Hydrogen impact excitation transitions from
    !+the `n` state to the \(m=1..{m_max}\) state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(n) \rightarrow H^+ + H(m=1..m_{max}), m \gt n $$
    !+
    !+###References
    !+* Eq. 29.b and Table 4 in Ref. 2 for \(n = 1\) and \(m = 2\) [[atomic_tables(module)]]
    !+* Eq. 30 and Table 5 in Ref. 2 for \(n = 1\) and \(m = 3-6\) [[atomic_tables(module)]]
    !+* Eq. 31 and Table 5 in Ref. 2 for \(n = 1\) and \(m \gt 6\) [[atomic_tables(module)]]
    !+* Eq. 32 and Table 6 in Ref. 2 for \(n = 2\) and \(m \le 5\) [[atomic_tables(module)]]
    !+* Eq. 33 and Table 6 in Ref. 2 for \(n = 2\) and \(m = 6-10\) [[atomic_tables(module)]]
    !+* Eq. 34 and Table 6 in Ref. 2 for \(n = 2\) and \(m \gt 10\) [[atomic_tables(module)]]
    !+* Eq. 35 and Table 7 in Ref. 2 for \(n = 3\) and \(m \le 6\) [[atomic_tables(module)]]
    !+* Eq. 36 and Table 7 in Ref. 2 for \(n = 3\) and \(m = 7-10\) [[atomic_tables(module)]]
    !+* Eq. 37 and Table 7 in Ref. 2 for \(n = 3\) and \(m \gt 10\) [[atomic_tables(module)]]
    !+* Eq. 38-39 in Ref. 2 for \(n \gt 3\) and \(m \gt 4\) [[atomic_tables(module)]]
    !+
    real(Float64), intent(in)       :: eb
        !+ Relative collision energy [keV/amu]
    integer, intent(in)             :: n
        !+ Initial atomic energy level/state
    integer, intent(in)             :: m_max
        !+ Number of m states to calculate
    real(Float64), dimension(m_max) :: sigma
        !+ Array of cross sections where the m'th index refer to the transition
        !+ from `n` to m [\(cm^2\)]

    integer :: m
    real(Float64) :: nf, mf, Etil, s, D, A, G, L, F
    real(Float64) :: y, zpl, zmi, C2pl, C2mi, H
    sigma = 0.d0

    select case (n)
        case (0)
            stop
        case (1)
            sigma = p_excit_1_janev(eb,m_max)
        case (2)
            sigma = p_excit_2_janev(eb,m_max)
        case (3)
            sigma = p_excit_3_janev(eb,m_max)
        case DEFAULT
            nf = real(n)
            m_loop: do m=1,m_max
                if(n.ge.m) then
                    sigma(m) = 0.d0
                    cycle m_loop
                endif
                mf = real(m)
                Etil = Eb/25.0
                s = (mf-nf)

                D = exp(-1.0/(nf*mf*Etil**2.0))
                A = 8.0/(3.0*s)*(mf/(s*nf))**3*(0.184-0.04/s**(2.0/3.0)) * &
                    (1.0 - 0.2*s/(nf*mf))**(1.0 + 2.0*s)
                G = 0.5*( Etil*nf**2.0/(mf - 1.0/mf) )**3.
                L = log(1.0 + 0.53*Etil**2.0 * nf*(mf - 2.0/mf)/(1.0 + 0.4*Etil))
                F = ( 1.0 - 0.3*s*D/(nf*mf) )**(1.0 + 2.0*s)

                y = 1.0/( 1.0 - D*log(18*s)/(4.0*s) )
                zpl = 2.0/(Etil * nf**2 * ( (2.0 - (nf/mf)**2)**0.5 + 1.0))
                zmi = 2.0/(Etil * nf**2 * ( (2.0 - (nf/mf)**2)**0.5 - 1.0))
                C2pl = zpl**2 * log(1.0 + 2.0*zpl/3.0)/(2.0*y + 3.0*zpl/2.0)
                C2mi = zmi**2 * log(1.0 + 2.0*zmi/3.0)/(2.0*y + 3.0*zmi/2.0)

                H = C2mi - C2pl

                sigma(m) = ((8.8d-17*n**4)/Etil)*(A*L*D + F*G*H)
            enddo m_loop
    end select

end function p_excit_n

function p_excit_n_m(eb, n, m) result(sigma)
    !+Calculates the cross section for a proton-Hydrogen impact excitation transition from
    !+the `n` state to the `m` state at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(n) \rightarrow H^+ + H(m), m \gt n $$
    !+
    !+###References
    !+* Eq. 29.b and Table 4 in Ref. 2 for \(n = 1\) and \(m = 2\) [[atomic_tables(module)]]
    !+* Eq. 30 and Table 5 in Ref. 2 for \(n = 1\) and \(m = 3-6\) [[atomic_tables(module)]]
    !+* Eq. 31 and Table 5 in Ref. 2 for \(n = 1\) and \(m \gt 6\) [[atomic_tables(module)]]
    !+* Eq. 32 and Table 6 in Ref. 2 for \(n = 2\) and \(m \le 5\) [[atomic_tables(module)]]
    !+* Eq. 33 and Table 6 in Ref. 2 for \(n = 2\) and \(m = 6-10\) [[atomic_tables(module)]]
    !+* Eq. 34 and Table 6 in Ref. 2 for \(n = 2\) and \(m \gt 10\) [[atomic_tables(module)]]
    !+* Eq. 35 and Table 7 in Ref. 2 for \(n = 3\) and \(m \le 6\) [[atomic_tables(module)]]
    !+* Eq. 36 and Table 7 in Ref. 2 for \(n = 3\) and \(m = 7-10\) [[atomic_tables(module)]]
    !+* Eq. 37 and Table 7 in Ref. 2 for \(n = 3\) and \(m \gt 10\) [[atomic_tables(module)]]
    !+* Eq. 38-39 in Ref. 2 for \(n \gt 3\) and \(m \gt 4\) [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    integer, intent(in)       :: n
        !+ Initial atomic energy level/state
    integer, intent(in)       :: m
        !+ Final atomic energy level/state
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(12) :: sigma_m

    sigma_m = p_excit_n(eb, n, 12)
    if(m.le.0) then
        sigma = sum(sigma_m)
    else
        sigma = sigma_m(m)
    endif

end function p_excit_n_m

function p_excit(eb, n_max, m_max) result(sigma)
    !+Calculates a matrix of cross sections for a proton-Hydrogen impact excitation transitions
    !+from the \(n=1..n_{max} \rightarrow m=1..m_{max}\) states at energy `eb`
    !+
    !+###Equation
    !+$$ H^+ + H(n=1..n_{max}) \rightarrow H^+ + H(m=1..m_{max}), m \gt n $$
    !+
    !+###References
    !+* Eq. 29.b and Table 4 in Ref. 2 for \(n = 1\) and \(m = 2\) [[atomic_tables(module)]]
    !+* Eq. 30 and Table 5 in Ref. 2 for \(n = 1\) and \(m = 3-6\) [[atomic_tables(module)]]
    !+* Eq. 31 and Table 5 in Ref. 2 for \(n = 1\) and \(m \gt 6\) [[atomic_tables(module)]]
    !+* Eq. 32 and Table 6 in Ref. 2 for \(n = 2\) and \(m \le 5\) [[atomic_tables(module)]]
    !+* Eq. 33 and Table 6 in Ref. 2 for \(n = 2\) and \(m = 6-10\) [[atomic_tables(module)]]
    !+* Eq. 34 and Table 6 in Ref. 2 for \(n = 2\) and \(m \gt 10\) [[atomic_tables(module)]]
    !+* Eq. 35 and Table 7 in Ref. 2 for \(n = 3\) and \(m \le 6\) [[atomic_tables(module)]]
    !+* Eq. 36 and Table 7 in Ref. 2 for \(n = 3\) and \(m = 7-10\) [[atomic_tables(module)]]
    !+* Eq. 37 and Table 7 in Ref. 2 for \(n = 3\) and \(m \gt 10\) [[atomic_tables(module)]]
    !+* Eq. 38-39 in Ref. 2 for \(n \gt 3\) and \(m \gt 4\) [[atomic_tables(module)]]
    !+
    real(Float64), intent(in)             :: eb
        !+ Relative collision energy [keV/amu]
    integer, intent(in)                   :: m_max
        !+ Number of initial atomic energy levels/states
    integer, intent(in)                   :: n_max
        !+ Number of final atomic energy levels/states
    real(Float64), dimension(n_max,m_max) :: sigma
        !+ Matrix of cross sections where the subscripts correspond
        !+ to the \(n \rightarrow m\) transitions: p_excit[n,m] [\(cm^2\)]

    real(Float64), dimension(12,12) :: sigma_full

    integer :: n, m

    do n=1,12
        sigma_full(n,:) = p_excit_n(eb, n, 12)
    enddo

    sigma = sigma_full(1:n_max,1:m_max)

end function p_excit

function e_ioniz_1_janev(eb) result(sigma)
    !+Calculates cross section for a electron-Hydrogen impact ionization from
    !+the \(n=1\) state at energy `eb`
    !+
    !+###Equation
    !+$$ e + H(1) \rightarrow e + H^+ + e$$
    !+
    !+###References
    !+* Eq. 14 and Table 3 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]
    integer, parameter :: n = 1
        !+ Initial atomic energy level/state
    real(Float64), dimension(6), parameter :: A = [ 0.18450d0, -0.032226d0, &
                                                   -0.034539d0, 1.4003d0,   &
                                                   -2.8115d0, 2.2986d0 ]
        !+ Fitting parameters from Table 3 in Ref. 2
    real(Float64) :: Edn2
    real(Float64) :: e, x
    real(Float64) :: s

    Edn2 = 13.6/real(n)**2
    e = eb * 1.d3 !keV to eV

    x = (1.0 - Edn2/e)
    s = A(2)*x + A(3)*(x**2.0) + A(4)*(x**3.0) + A(5)*(x**4.0) + A(6)*(x**5.0)
    sigma = ((1.d-13)/(Edn2*e))*(A(1)*log(e/Edn2) + s)
    sigma = max(sigma,0.d0)

end function e_ioniz_1_janev

function e_ioniz_2_janev(eb) result(sigma)
    !+Calculates cross section for a electron-Hydrogen impact ionization from
    !+the \(n=2\) state at energy `eb`
    !+
    !+###Equation
    !+$$ e + H(2) \rightarrow e + H^+ + e$$
    !+
    !+###References
    !+* Eq. 14 and Table 3 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]
    integer, parameter :: n = 2
        !+ Initial atomic energy level/state
    real(Float64), dimension(6), parameter :: A = [ 0.14784d0, 0.0080871d0, &
                                                   -0.062270d0, 1.9414d0,   &
                                                   -2.1980d0, 0.95894d0 ]
        !+ Fitting parameters from Table 3 in Ref. 2

    real(Float64) :: Edn2
    real(Float64) :: e, x
    real(Float64) :: s

    Edn2 = 13.6/real(n)**2
    e = eb * 1.d3 !keV to eV

    x = (1.0 - Edn2/e)
    s = A(2)*x + A(3)*(x**2.0) + A(4)*(x**3.0) + A(5)*(x**4.0) + A(6)*(x**5.0)
    sigma = ((1.d-13)/(Edn2*e))*(A(1)*log(e/Edn2) + s)
    sigma = max(sigma, 0.d0)

end function e_ioniz_2_janev

function e_ioniz_3_janev(eb) result(sigma)
    !+Calculates cross section for a electron-Hydrogen impact ionization from
    !+the \(n=3\) state at energy `eb`
    !+
    !+###Equation
    !+$$ e + H(3) \rightarrow e + H^+ + e$$
    !+
    !+###References
    !+* Eq. 14 and Table 3 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]
    integer, parameter :: n = 3
        !+ Initial atomic energy level/state
    real(Float64), dimension(6), parameter :: A = [0.058463d0, -0.051272d0, &
                                                   0.85310d0, -0.57014d0,   &
                                                   0.76684d0, 0.00d0 ]
        !+ Fitting parameters from Table 3 in Ref. 2

    real(Float64) :: Edn2
    real(Float64) :: e, x
    real(Float64) :: s

    Edn2 = 13.6/real(n)**2
    e = eb * 1.d3 !keV to eV

    if(e.ge.1.5) then
        x = (1.0 - Edn2/e)
        s = A(2)*x + A(3)*(x**2.0) + A(4)*(x**3.0) + A(5)*(x**4.0) + A(6)*(x**5.0)
        sigma = ((1.d-13)/(Edn2*e))*(A(1)*log(e/Edn2) + s)
    else
        sigma = 0.d0
    endif

end function e_ioniz_3_janev

function e_ioniz_n(eb, n) result(sigma)
    !+Calculates cross section for a electron-Hydrogen impact ionization from
    !+the `n` state at energy `eb`
    !+
    !+###Equation
    !+$$ e + H(n) \rightarrow e + H^+ + e$$
    !+
    !+###References
    !+* Eq. 14 and Table 3 in Ref. 2 [[atomic_tables(module)]]
    !+* Eq. 15-16 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    integer, intent(in)       :: n
        !+ Initial atomic energy level/state
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64) :: rn, xn, Edn2
    real(Float64) :: g0, g1, g2, An, b, Bn

    select case (n)
        case (0)
            stop
        case (1)
            sigma = e_ioniz_1_janev(eb)
        case (2)
            sigma = e_ioniz_2_janev(eb)
        case (3)
            sigma = e_ioniz_3_janev(eb)
        case DEFAULT
            rn = 1.94/n**1.57
            Edn2 = 13.6/n**2.0
            xn = (eb*1.d3)/Edn2

            g0 = 0.9935 + 0.2328/n - 0.1296/n**2.0
            g1 = -(1.0/n)*(0.6282 - 0.5598/n + 0.5299/n**2.0)
            g2 =  (1.0/n**2.0)*(0.3887 - 1.181/n + 1.47/n**2.0)

            An = 32.0*n/(3.0*sqrt(3.0)*PI)*(g0/3.0 + g1/4.0 + g2/5.0)
            b  = (1.0/n)*(4.0 - 18.63/n + 36.24/n**2.0 - 28.09/n**3.0)
            Bn = (2.0/3.0)*(n**2.0)*(5.0 + b)

            if(xn.gt.1) then
                sigma = 1.76*n**2/xn*(1.0 - exp(-rn*xn)) * &
                        (An*log(xn) + (Bn - An*log(2.0*n**2)) * &
                        (1.0 - 1.0/xn)**2)*1.e-16
            else
                sigma = 0.d0
            endif
    end select

end function e_ioniz_n

function e_ioniz(eb, n_max) result(sigma)
    !+Calculates an array of cross sections for a electron-Hydrogen impact ionization from
    !+the \(n=1..n_{max}\) states at energy `eb`
    !+
    !+###Equation
    !+$$ e + H(n=1..n_{max}) \rightarrow e + H^+ + e$$
    !+
    !+###References
    !+* Eq. 14 and Table 3 in Ref. 2 [[atomic_tables(module)]]
    !+* Eq. 15-16 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in)       :: eb
        !+ Collision energy [keV]
    integer, intent(in)             :: n_max
        !+ Number of initial atomic energy levels/states to calculate
    real(Float64), dimension(n_max) :: sigma
        !+ Array of cross sections where the n'th index refers to a ionization from the n'th state [\(cm^2\)]

    integer :: i

    do i=1,n_max
        sigma(i) = e_ioniz_n(eb,i)
    enddo

end function e_ioniz

function e_excit_1_2_janev(eb) result(sigma)
    !+Calculates cross section for a electron-Hydrogen impact excitation transition from
    !+the \(n=1\) state to the \(m=2\) state at energy `eb`
    !+
    !+###Equation
    !+$$ e + H(1) \rightarrow e + H(2) $$
    !+
    !+###References
    !+* Eq. 4 and Table 1 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: sigma0=5.984d0
    real(Float64), parameter :: deltaE=10.2d0
    real(Float64), parameter :: a=0.228d0
        !+ Fitting parameter from Table 2 in Ref. 2
    real(Float64), parameter :: b=0.1865d0
        !+ Fitting paramter from Table 2 in Ref. 2
    real(Float64), parameter :: c=0.5025d0
        !+ Fitting paramter from Table 2 in Ref. 2
    real(Float64), dimension(6), parameter :: An = [ 4.4979d0, 1.4182d0, &
                                                    -20.877d0, 49.735d0, &
                                                    -46.249d0, 17.442d0 ]
        !+ Fitting parameters from Table 2 in Ref. 2

    real(Float64) :: ecoll, x, s

    ecoll = eb*1.d3
    x = (ecoll)/deltaE

    if((ecoll.gt.10.2).and.(ecoll.le.11.56)) then
        sigma = 1.d-16 * (a + b*(ecoll - deltaE))
        return
    endif

    if((ecoll.ge.11.56).and.(ecoll.le.12.23)) then
        sigma = 1.d-16 * c
        return
    endif

    if(ecoll.ge.12.23) then
        s = An(2) + An(3)/x + An(4)/x**2.0 + An(5)/x**3.0 + An(6)/x**4.0
        sigma = 1.d-16 * sigma0/(deltaE*x) * (An(1)*log(x) + s)
        return
    endif

    if(x.le.1.0) then
        sigma = 0.0
        return
    endif

end function e_excit_1_2_janev

function e_excit_1_3_janev(eb) result(sigma)
    !+Calculates cross section for a electron-Hydrogen impact excitation transition from
    !+the \(n=1\) state to the \(m=3\) state at energy `eb`
    !+
    !+###Equation
    !+$$ e + H(1) \rightarrow e + H(3) $$
    !+
    !+###References
    !+* Eq. 5 and Table 2 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: sigma0 = 5.984d0
    real(Float64), parameter :: deltaE = 12.09d0
        !+ Fitting parameter from Table 2 in Ref. 2
    real(Float64), parameter :: alpha = 0.38277d0
        !+ Fitting parameter from Table 2 in Ref. 2
    real(Float64), dimension(5), parameter :: A = [ 0.75448d0, 0.42956d0, &
                                                   -0.58288d0, 1.0693d0,  &
                                                    0.d0 ]
        !+ Fitting parameters from Table 2 in Ref. 2

    real(Float64) :: ecoll, x, s

    ecoll = eb*1.d3
    x=ecoll/deltaE
    s = A(2) + A(3)/x + A(4)/x**2.0 + A(5)/x**3.0

    sigma = 1.d-16 * sigma0/(deltaE*x) * (1.0 - 1.0/x)**alpha * (A(1)*log(x) + s)

    if(x.le.1.0) sigma = 0.d0

end function e_excit_1_3_janev

function e_excit_1_4_janev(eb) result(sigma)
    !+Calculates cross section for a electron-Hydrogen impact excitation transition from
    !+the \(n=1\) state to the \(m=4\) state at energy `eb`
    !+
    !+###Equation
    !+$$ e + H(1) \rightarrow e + H(4) $$
    !+
    !+###References
    !+* Eq. 5 and Table 2 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: sigma0 = 5.984d0
    real(Float64), parameter :: deltaE = 12.75d0
        !+ Fitting parameter from Table 2 in Ref. 2
    real(Float64), parameter :: alpha = 0.41844d0
        !+ Fitting parameter from Table 2 in Ref. 2
    real(Float64), dimension(5), parameter :: A = [ 0.24300d0, 0.24846d0, &
                                                    0.19701d0, 0.d0,      &
                                                    0.d0 ]
        !+ Fitting parameters from Table 2 in Ref. 2

    real(Float64) :: ecoll, x, s

    ecoll = eb*1.d3
    x=ecoll/deltaE
    s = A(2) + A(3)/x + A(4)/x**2.0 + A(5)/x**3.0

    sigma = 1.d-16 * sigma0/(deltaE*x) * (1.0 - 1.0/x)**alpha * (A(1)*log(x) + s)

    if(x.le.1.0) sigma = 0.d0

end function e_excit_1_4_janev

function e_excit_1_5_janev(eb) result(sigma)
    !+Calculates cross section for a electron-Hydrogen impact excitation transition from
    !+the \(n=1\) state to the \(m=5\) state at energy `eb`
    !+
    !+###Equation
    !+$$ e + H(1) \rightarrow e + H(5) $$
    !+
    !+###References
    !+* Eq. 5 and Table 2 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: sigma0 = 5.984d0
    real(Float64), parameter :: deltaE = 13.06d0
        !+ Fitting parameter from Table 2 in Ref. 2
    real(Float64), parameter :: alpha = 0.45929d0
        !+ Fitting parameter from Table 2 in Ref. 2
    real(Float64), dimension(5), parameter :: A = [ 0.11508d0, 0.13092d0, &
                                                    0.23581d0, 0.d0,      &
                                                    0.d0 ]
        !+ Fitting parameters from Table 2 in Ref. 2

    real(Float64) :: ecoll, x, s

    ecoll = eb*1.d3
    x=ecoll/deltaE
    s = A(2) + A(3)/x + A(4)/x**2.0 + A(5)/x**3.0

    sigma = 1.d-16 * sigma0/(deltaE*x) * (1.0 - 1.0/x)**alpha * (A(1)*log(x) + s)

    if(x.le.1.0) sigma = 0.d0

end function e_excit_1_5_janev

function e_excit_f(n, m) result(fnm)
    !+ Oscillator strength for a `n`\(\rightarrow\)`m` transition due to electron-Hydrogen impact excitation
    !+
    !+###References
    !+* Eqs. 11-13 in Ref. 2 [[atomic_tables(module)]]
    integer, intent(in) :: n
        !+ Initial atomic energy level/state
    integer, intent(in) :: m
        !+ Final atomic energy level/state
    real(Float64)       :: fnm
        !+ Oscillator strength

    real(Float64), dimension(3) :: g
    real(Float64) :: x, nf, mf, gs

    nf = real(n)
    mf = real(m)
    x = 1.0 - (nf/mf)**2.0

    select case (n)
        case (1)
            g = [1.133,-0.4059,0.0714]
        case (2)
            g = [1.0785,-0.2319,0.02947]
        case DEFAULT
            g(1) = 0.9935 + 0.2328/nf - 0.1296/nf**2
            g(2) =-1.0/nf * (0.6282 - 0.5598/nf + 0.5299/nf**2)
            g(3) = 1.0/nf**2.0 * (0.3887 - 1.1810/nf + 1.4700/nf**2)
    end select

    gs = g(1) + g(2)/x + g(3)/x**2
    fnm = 32.0/(3.0*sqrt(3.0)*PI) * nf/mf**3 * 1/x**3 * gs

end function e_excit_f

function e_excit_1_janev(eb, m_max) result(sigma)
    !+Calculates an array of cross sections for a electron-Hydrogen impact excitation transition from
    !+the \(n=1\) state to the \(m=1..m_{max}\) state at energy `eb`
    !+
    !+###Equation
    !+$$ e + H(1) \rightarrow e + H(m=1..m_{max}) $$
    !+
    !+###References
    !+* Eq. 5 and Table 2 in Ref. 2 [[atomic_tables(module)]]
    !+* Eqs. 6-7 in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in)       :: eb
        !+ Collision energy [keV]
    integer, intent(in)             :: m_max
        !+ Number of m states to calculate
    real(Float64), dimension(m_max) :: sigma
        !+ Array of cross sections where the m'th index refers to
        !+an excitation from the \(n=1\) state to the m'th state [\(cm^2\)]

    integer :: m
    real(Float64) :: x, y, A, B, deltaE

    do m=1,m_max
        select case (m)
            case (1)
                sigma(1) = 0.d0
            case (2)
                sigma(2) = e_excit_1_2_janev(eb)
            case (3)
                sigma(3) = e_excit_1_3_janev(eb)
            case (4)
                sigma(4) = e_excit_1_4_janev(eb)
            case (5)
                sigma(5) = e_excit_1_5_janev(eb)
            case DEFAULT
                y = 1.0 - (1.d0/m)**2.0
                deltaE = 13.6*y
                x = (eb*1.d3)/deltaE

                A = 2.0 * e_excit_f(1,m)/y
                B = 4.0/(m**3.0 * y)*(1.0 + 4.0/(3.0*y) - 0.603/y**2.0)

                sigma(m) = 1.76e-16/(y*x)*(1.0 - exp(-0.45*y*x))* &
                           (A*(log(x) + 1.0/(2.0*x)) + (B - A*log(2.0/y))* &
                           (1.0 - 1.0/x))

                if(x.le.1.0) sigma(m) = 0.d0
        end select
    enddo

end function e_excit_1_janev

function e_excit_2_3_janev(eb) result(sigma)
    !+Calculates cross section for a electron-Hydrogen impact excitation transition from
    !+the \(n=2\) state to the \(m=3\) state at energy `eb`
    !+
    !+###Equation
    !+$$ e + H(2) \rightarrow e + H(3) $$
    !+
    !+###References
    !+* Eq. 5 in Ref. 2 [[atomic_tables(module)]]
    !+* Section 2.1.1 B in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: sigma0 = 5.984d0
    real(Float64), parameter :: deltaE = 1.8888d0
        !+ Energy difference between \(n=1\) and \(n=3\):
        !+ \(\Delta E = 13.6\left (\frac{1}{2^2} - \frac{1}{3^2}\right )\)
    real(Float64), parameter :: alpha = 1.3196d0
        !+ Fitting parameter from Section 2.1.1 B in Ref. 2
    real(Float64), dimension(5), parameter :: A = [ 38.906d0, 5.2373d0, 119.25d0, &
                                                   -595.39d0, 816.71d0]
        !+ Fitting parameters from Section 2.1.1 B in Ref. 2

    real(Float64) :: ecoll, x, s

    ecoll = eb*1.d3
    x = ecoll/deltaE
    s = A(2) + A(3)/x + A(4)/x**2.0 + A(5)/x**3.0

    sigma = 1.d-16 * sigma0/(deltaE*x) * (1.0 - 1.0/x)**alpha * (A(1)*log(x) + s)

    if(x.le.1.0) sigma = 0.d0

end function e_excit_2_3_janev

function e_excit_n(eb, n, m_max) result(sigma)
    !+Calculates an array of cross sections for a electron-Hydrogen impact excitation transition from
    !+the `n` state to the \(m=1..m_{max}\) state at energy `eb`
    !+
    !+###Equation
    !+$$ e + H(n) \rightarrow e + H(m=1..m_{max}), m \gt n $$
    !+
    !+###References
    !+* Eq. 5 and Table 2 in Ref. 2 [[atomic_tables(module)]]
    !+* Eqs. 6-7 in Ref. 2 [[atomic_tables(module)]]
    !+* Section 2.1.1 B in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in)       :: eb
        !+ Collision energy [keV]
    integer, intent(in)             :: n
        !+ Initial atomic energy level/state
    integer, intent(in)             :: m_max
        !+ Number of m states to calculate
    real(Float64), dimension(m_max) :: sigma
        !+ Array of cross sections where the m'th index refers to
        !+an excitation from the `n`\(\rightarrow\)`m` state [\(cm^2\)]

    integer :: m
    real(Float64) :: nf, mf
    real(Float64) :: x, y, A, B, bn, r, deltaE

    nf = real(n)

    if(n.eq.1) then
        sigma = e_excit_1_janev(eb, m_max)
    else
        m_loop: do m=1,m_max
            mf = real(m)
            if(n.ge.m) then
                sigma(m) = 0.d0
                cycle m_loop
            endif
            if((n.eq.2).and.(m.eq.3)) then
                sigma(m) = e_excit_2_3_janev(eb)
            else
                deltaE=13.6*(1.0/nf**2 - 1.0/mf**2)
                x = (eb*1.d3)/deltaE
                y = 1.0 - (nf/mf)**2
                r = 1.94/nf**1.57

                A = 2.0 * nf**2 * e_excit_f(n,m)/y
                bn = 1.0/nf*(4.0 - 18.63/nf + 36.24/nf**2 - 28.09/nf**3)
                B = 4.0 * nf**4/(mf**3*y**2)*(1.0 + 4.0/(3.0*y) + bn/y**2.0)

                sigma(m) = 1.76e-16*nf**2/(y*x)*(1.0 - exp(-r*y*x))* &
                           (A*(log(x) + 1.0/(2.0*x)) + (B - A*log(2.0*n**2.0/y))* &
                           (1.0 - 1.0/x))

                if(x.le.1.0) sigma(m) = 0.d0
            endif
        enddo m_loop
    endif

end function e_excit_n

function e_excit_n_m(eb, n, m) result(sigma)
    !+Calculates an array of cross sections for a electron-Hydrogen impact excitation transition from
    !+the `n` \(\rightarrow\) `m` state at energy `eb`
    !+
    !+###Equation
    !+$$ e + H(n) \rightarrow e + H(m), m \gt n $$
    !+
    !+###References
    !+* Eq. 5 and Table 2 in Ref. 2 [[atomic_tables(module)]]
    !+* Eqs. 6-7 in Ref. 2 [[atomic_tables(module)]]
    !+* Section 2.1.1 B in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in)       :: eb
        !+ Collision energy [keV]
    integer, intent(in)             :: n
        !+ Initial atomic energy level/state
    integer, intent(in)       :: m
        !+ Final atomic energy level/state
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(12) :: sigma_m

    sigma_m = e_excit_n(eb, n, 12)
    if(m.le.0) then
        sigma = sum(sigma_m)
    else
        sigma = sigma_m(m)
    endif

end function e_excit_n_m

function e_excit(eb, n_max, m_max) result(sigma)
    !+Calculates a matrix of cross section for a proton-Hydrogen impact excitation transition
    !+from the \(n=1..n_{max} \rightarrow m=1..m_{max}\) states at energy `eb`
    !+
    !+###Equation
    !+$$ e + H(n=1..n_{max}) \rightarrow e + H(m=1..m_{max}), m \gt n $$
    !+
    !+###References
    !+* Eq. 5 and Table 2 in Ref. 2 [[atomic_tables(module)]]
    !+* Eqs. 6-7 in Ref. 2 [[atomic_tables(module)]]
    !+* Section 2.1.1 B in Ref. 2 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in)             :: eb
        !+ Relative collision energy [keV/amu]
    integer, intent(in)                   :: m_max
        !+ Number of initial atomic energy levels/states
    integer, intent(in)                   :: n_max
        !+ Number of final atomic energy levels/states
    real(Float64), dimension(n_max,m_max) :: sigma
        !+ Matrix of cross sections where the subscripts correspond
        !+ to the \(n \rightarrow m\) transitions: e_excit[n,m] [\(cm^2\)]

    real(Float64), dimension(12,12) :: sigma_full

    integer :: n

    do n=1,12
        sigma_full(n,:) = e_excit_n(eb, n, 12)
    enddo

    sigma = sigma_full(1:n_max,1:m_max)

end function e_excit

!Impurities
!A[q]_cx_[n]_[source]
function B5_cx_1_adas(eb) result(sigma)
    !+ Calculates the total charge exchange cross section for a Neutral Hydrogen atom
    !+in the \(n=1\) state colliding with a fully stripped Boron ion at energy `eb`
    !+
    !+###Equation
    !+$$ B^{5+} + H(1) \rightarrow B^{4+} + H^+ $$
    !+
    !+###References
    !+* Ref. 4 [[atomic_tables(module)]]
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(7), parameter :: A = [ 1.174052518d3, -1.793561728d3, &
                                                    1.117522436d3, -3.679435571d2, &
                                                    6.750816878d1, -6.542029074d0, &
                                                    2.614113716d-1 ]
    real(Float64) :: e, l, p

    e = max(eb*1.d3,10.0)
    l = log10(e)

    if(e.le.4.d5) then
        p = A(1) + A(2)*l + A(3)*l**2 + A(4)*l**3 + &
            A(5)*l**4 + A(6)*l**5 + A(7)*l**6
        sigma = 10.d0**p
    else
        sigma = 0.d0
    endif

end function B5_cx_1_adas

function B5_cx_2_adas(eb) result(sigma)
    !+ Calculates the total charge exchange cross section for a Neutral Hydrogen atom
    !+in the \(n=2\) state colliding with a fully stripped Boron ion at energy `eb`
    !+
    !+###Equation
    !+$$ B^{5+} + H(2) \rightarrow B^{4+} + H^+ $$
    !+
    !+###References
    !+* Ref. 4 [[atomic_tables(module)]]
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(10), parameter :: A = [6.603246818d1, -3.072575676d2, &
                                                    5.030801019d2, -4.585636345d2, &
                                                    2.568666393d2, -9.185150382d1, &
                                                    2.100012584d1, -2.964174788d0, &
                                                    2.346396110d-1, -7.943766873d-3]
    real(Float64) :: e, l, p

    e = max(eb*1.d3,10.0)
    l = log10(e)

    if(e.le.1.d5) then
        p = A(1) + A(2)*l + A(3)*l**2 + A(4)*l**3 + &
            A(5)*l**4 + A(6)*l**5 + A(7)*l**6 +     &
            A(8)*l**7 + A(9)*l**8 + A(10)*l**9
        sigma = 10.d0**p
    else
        sigma = 0.d0
    endif

end function B5_cx_2_adas

function C6_cx_1_adas(eb) result(sigma)
    !+ Calculates the total charge exchange cross section for a Neutral Hydrogen atom
    !+in the \(n=1\) state colliding with a fully stripped Carbon ion at energy `eb`
    !+
    !+###Equation
    !+$$ C^{6+} + H(1) \rightarrow C^{5+} + H^+ $$
    !+
    !+###References
    !+* Ref. 4 [[atomic_tables(module)]]
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(7), parameter :: A = [2.007882674d2, -3.546893286d2, &
                                                   2.381542403d2, -8.355431742d1, &
                                                   1.617519888d1, -1.638152470d0, &
                                                   6.768953863d-2 ]

    real(Float64) :: e, l, p

    e = max(eb*1.d3,1.5d3)
    l = log10(e)

    p = A(1) + A(2)*l + A(3)*l**2 + A(4)*l**3 + &
        A(5)*l**4 + A(6)*l**5 + A(7)*l**6
    sigma = 10.d0**p

end function C6_cx_1_adas

function C6_cx_2_adas(eb) result(sigma)
    !+ Calculates the total charge exchange cross section for a Neutral Hydrogen atom
    !+in the \(n=2\) state colliding with a fully stripped Carbon ion at energy `eb`
    !+
    !+###Equation
    !+$$ C^{6+} + H(2) \rightarrow C^{5+} + H^+ $$
    !+
    !+###References
    !+* Ref. 4 [[atomic_tables(module)]]
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(11), parameter :: A = [9.151879441d5, -2.134573133d6, &
                                                    2.223792624d6, -1.362648703d6, &
                                                    5.438401343d5, -1.477110500d5, &
                                                    2.764972254d4, -3.522105245d3, &
                                                    2.921934171d2, -1.425552507d1, &
                                                    3.106007048d-1 ]
    real(Float64) :: e, l, p

    e = max(eb*1.d3,1.5d3)*2.0**2
    l = log10(e)

    p = A(1) + A(2)*l + A(3)*l**2 + A(4)*l**3 + &
        A(5)*l**4 + A(6)*l**5 + A(7)*l**6 +     &
        A(8)*l**7 + A(9)*l**8 + A(10)*l**9 + A(11)*l**10
    sigma = 10.d0**p

end function C6_cx_2_adas

function C6_cx_3_adas(eb) result(sigma)
    !+ Calculates the total charge exchange cross section for a Neutral Hydrogen atom
    !+in the \(n=3\) state colliding with a fully stripped Carbon ion at energy `eb`
    !+
    !+###Equation
    !+$$ C^{6+} + H(3) \rightarrow C^{5+} + H^+ $$
    !+
    !+###References
    !+* Ref. 4 [[atomic_tables(module)]]
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(11), parameter :: A = [9.208877916d5, -2.147294379d6, &
                                                    2.236451628d6, -1.370042347d6, &
                                                    5.466461899d5, -1.484338816d5, &
                                                    2.777765778d4, -3.537459450d3, &
                                                    2.933884362d2, -1.430994136d1, &
                                                    3.117002878d-1 ]
    real(Float64) :: e, l, p

    e = max(eb*1.d3,1.5d3)*3.0**2
    l = log10(e)

    p = A(1) + A(2)*l + A(3)*l**2 + A(4)*l**3 + &
        A(5)*l**4 + A(6)*l**5 + A(7)*l**6 +     &
        A(8)*l**7 + A(9)*l**8 + A(10)*l**9 + A(11)*l**10
    sigma = 10.d0**p

end function C6_cx_3_adas

function Aq_cx_n_adas(eb, q, n) result(sigma)
    !+ Calculates the total charge exchange cross section for a Neutral Hydrogen atom
    !+in the `n` state colliding with a ion with charge `q` at energy `eb`
    !+
    !+@note Returns 0 if ADAS cross sections are not available for given inputs
    !+
    !+###Equation
    !+$$ A^{q+} + H(n) \rightarrow A^{(q-1)+} + H^+ $$
    !+
    !+###References
    !+* Ref. 4 [[atomic_tables(module)]]
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    integer, intent(in)       :: q
        !+ Ion charge
    integer, intent(in)       :: n
        !+ Initial atomic energy level/state
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    sigma = 0.d0
    select case (q)
        case (5)
            if(n.eq.1) sigma = B5_cx_1_adas(eb)
            if(n.eq.2) sigma = B5_cx_2_adas(eb)
        case (6)
            if(n.eq.1) sigma = C6_cx_1_adas(eb)
            if(n.eq.2) sigma = C6_cx_2_adas(eb)
            if(n.eq.3) sigma = C6_cx_3_adas(eb)
        case DEFAULT
            sigma = 0.d0
    end select

end function Aq_cx_n_adas

function B5_cx_1_janev(eb) result(sigma)
    !+ Calculates the total charge exchange cross section for a Neutral Hydrogen atom
    !+in the \(n=1\) state colliding with a fully stripped Boron ion at energy `eb`
    !+
    !+###Equation
    !+$$ B^{5+} + H(1) \rightarrow B^{4+} + H^+ $$
    !+
    !+###References
    !+* Page 166 in Ref. 5 [[atomic_tables(module)]]
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(11), parameter :: A = [31.226d0, 1.1442d0,    &
                                                    4.8372d-8, 3.0961d-10, &
                                                    4.7205d0, 6.2844d-7,   &
                                                    3.1297d0, 0.12556d0,   &
                                                    0.30098d0, 5.9607d-2,  &
                                                   -0.57923d0 ]

    sigma = 1.d-16*A(1)*(exp(-A(2)/eb**A(8)) /  &
           (1.0 + A(3)*eb**2 + A(4)*eb**A(5) +  &
            A(6)*eb**A(7)) + A(9)*exp(-A(10)*eb) /eb**A(11))

end function B5_cx_1_janev

function C6_cx_1_janev(eb) result(sigma)
    !+ Calculates the total charge exchange cross section for a Neutral Hydrogen atom
    !+in the \(n=1\) state colliding with a fully stripped Carbon ion at energy `eb`
    !+
    !+###Equation
    !+$$ C^{6+} + H(1) \rightarrow C^{5+} + H^+ $$
    !+
    !+###References
    !+* Page 168 in Ref. 5 [[atomic_tables(module)]]
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(11), parameter :: A = [418.18d0, 2.1585d0,   &
                                                    3.4808d-4, 5.3333d-9, &
                                                    4.6556d0, 0.33755d0,  &
                                                    0.81736d0, 0.27874d0, &
                                                    1.8003d-6, 7.1033d-2, &
                                                    0.53261d0 ]

    sigma = 1.d-16*A(1)*(exp(-A(2)/eb**A(8)) /  &
           (1.0 + A(3)*eb**2 + A(4)*eb**A(5) +  &
            A(6)*eb**A(7)) + A(9)*exp(-A(10)*eb) /eb**A(11))

end function C6_cx_1_janev

function Aq_cx_n_janev(eb, q, n) result(sigma)
    !+ Calculates the total charge exchange cross section for a Neutral Hydrogen atom
    !+in the `n` state colliding with a ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(n) \rightarrow A^{(q-1)+} + H^+, q \gt 3 $$
    !+
    !+###References
    !+* Page 166 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 168 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 174 in Ref. 5 [[atomic_tables(module)]]
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    integer, intent(in)       :: q
        !+ Ion charge
    integer, intent(in)       :: n
        !+ Initial atomic energy level/state
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: A = 1.507d5
    real(Float64), parameter :: B = 1.974d-5

    real(Float64) :: etil, nf, qf

    nf = real(n)
    qf = real(q)

    if((n.eq.1).and.(q.eq.5)) then
        sigma = B5_cx_1_janev(eb)
        return
    endif

    if((n.eq.1).and.(q.eq.6)) then
        sigma = C6_cx_1_janev(eb)
        return
    endif

    if(n.le.1) then
        sigma = 0.d0
        return
    endif

    etil = eb*(nf**2.0)/(qf**0.5)

    sigma = qf*nf**4 * 7.04d-16 * A/(etil**3.5 * (1.0 + B*etil**2)) * &
            (1.0 - exp(-2.0*etil**3.5 * (1.0 + B*etil**2)/(3.0*A)))

end function Aq_cx_n_janev

function Aq_cx_n(eb, q, n) result(sigma)
    !+ Calculates the total charge exchange cross section for a Neutral Hydrogen atom
    !+in the `n` state colliding with a ion with charge `q` at energy `eb`
    !+
    !+@note Uses ADAS(Ref. 4) cross sections if available else uses Janev (Ref. 5) cross sections
    !+
    !+###Equation
    !+$$ A^{q+} + H(n) \rightarrow A^{(q-1)+} + H^+, q \gt 3 $$
    !+
    !+###References
    !+* Ref. 4 [[atomic_tables(module)]]
    !+* Page 166 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 168 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 174 in Ref. 5 [[atomic_tables(module)]]
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    integer, intent(in)       :: q
        !+ Ion charge
    integer, intent(in)       :: n
        !+ Initial atomic energy level/state
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    sigma = Aq_cx_n_adas(eb, q, n)
    if(sigma.eq.0.d0) then
        sigma = Aq_cx_n_janev(eb, q, n)
    endif

end function Aq_cx_n

function Aq_cx(eb, q, n_max) result(sigma)
    !+ Calculates an array of total charge exchange cross sections for a Neutral Hydrogen atom
    !+in the n=1...n_max states colliding with a ion with charge `q` at energy `eb`
    !+
    !+@note Uses ADAS(Ref. 4) cross sections if available else uses Janev (Ref. 5) cross sections
    !+
    !+###Equation
    !+$$ A^{q+} + H(n=1..n_{max}) \rightarrow A^{(q-1)+} + H^+, q \gt 3 $$
    !+
    !+###References
    !+* Ref. 4 [[atomic_tables(module)]]
    !+* Page 174 in Ref. 5 [[atomic_tables(module)]]
    real(Float64), intent(in)       :: eb
        !+ Relative collision energy [keV/amu]
    integer, intent(in)             :: q
        !+ Ion charge
    integer, intent(in)             :: n_max
        !+ Number of initial atomic energy levels/states
    real(Float64), dimension(n_max) :: sigma
        !+ Array of cross sections where the n'th index refers to a charge exchange from the n'th state [\(cm^2\)]

    integer :: n

    do n=1,n_max
        sigma(n) = Aq_cx_n(eb, q, n)
    enddo

end function Aq_cx

!Impurity impact ionization
function B5_ioniz_1_janev(eb) result(sigma)
    !+ Calculates the total ionization cross section for a Neutral Hydrogen atom
    !+in the \(n=1\) state colliding with a fully stripped Boron ion at energy `eb`
    !+
    !+###Equation
    !+$$ B^{5+} + H(1) \rightarrow B^{5+} + H^+ + e $$
    !+
    !+###References
    !+* Page 152 in Ref. 5 [[atomic_tables(module)]]
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(8), parameter :: A = [351.52d0, 233.63d0,   &
                                                   3.2952d3, 5.3787d-6,  &
                                                   1.8834d-2, -2.2064d0, &
                                                   7.2074d0, -3.78664d0 ]

    sigma = 1.d-16*A(1)*(exp(-A(2)/eb)*log(1 + A(3)*eb)/eb &
          + A(4)*exp(-A(5)*eb)/((eb**A(6)) + A(7)*(eb**A(8))))

end function B5_ioniz_1_janev

function C6_ioniz_1_janev(eb) result(sigma)
    !+ Calculates the total ionization cross section for a Neutral Hydrogen atom
    !+in the \(n=1\) state colliding with a fully stripped Carbon ion at energy `eb`
    !+
    !+###Equation
    !+$$ C^{6+} + H(1) \rightarrow C^{6+} + H^+ + e $$
    !+
    !+###References
    !+* Page 154 in Ref. 5 [[atomic_tables(module)]]
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma

    real(Float64), dimension(8), parameter :: A = [ 438.36d0, 327.10d0,    &
                                                    1.4444d5, 3.5212d-3,   &
                                                    8.3031d-3, -0.63731d0, &
                                                    1.9116d4, -3.1003d0 ]

    sigma = 1.d-16*A(1)*(exp(-A(2)/eb)*log(1 + A(3)*eb)/eb &
          + A(4)*exp(-A(5)*eb)/((eb**A(6)) + A(7)*(eb**A(8))))

end function C6_ioniz_1_janev

function Aq_ioniz_n_janev(eb, q, n) result(sigma)
    !+ Calculates the generic total ionization cross section for a Neutral Hydrogen atom
    !+in the `n` state colliding with a ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(n) \rightarrow A^{q+} + H^+ + e, n \gt 1, q \gt 3 $$
    !+
    !+###References
    !+* Page 160 in Ref. 5 [[atomic_tables(module)]]
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    integer, intent(in)       :: q
        !+ Ion charge
    integer, intent(in)       :: n
        !+ Initial atomic energy level/state
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: M = 0.283d0
    real(Float64), parameter :: B = 4.04d0
    real(Float64), parameter :: c = 137.d0
    real(Float64), parameter :: g = 0.662d0
    real(Float64), parameter :: lambda = 0.76d0

    real(Float64) :: nf, qf, u, v, sigma_b

    nf = real(n)
    qf = real(q)
    v = sqrt(eb/25.)
    u = nf*v

    sigma_b = 3.52d-16 * (nf**4) * (qf**2)/(u**2) * &
              (M * (log((u**2)/(c**2 - u**2)) - (u**2)/(c**2)) + B - g/u**2)
    sigma_b = max(sigma_b,0.d0)

    sigma = exp(-lambda*qf/u**2)*sigma_b

end function Aq_ioniz_n_janev

function Aq_ioniz_n(eb, q, n) result(sigma)
    !+ Calculates the total ionization cross section for a Neutral Hydrogen atom
    !+in the `n` state colliding with a ion with charge `q` at energy `eb`
    !+
    !+@note Uses specialized cross sections if available else uses generic cross sections
    !+
    !+###Equation
    !+$$ A^{q+} + H(n) \rightarrow A^{(q-1)+} + H^+, q \gt 3 $$
    !+
    !+###References
    !+* Page 152 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 154 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 160 in Ref. 5 [[atomic_tables(module)]]
    real(Float64), intent(in) :: eb
        !+ Relative collision energy [keV/amu]
    integer, intent(in)       :: q
        !+ Ion charge
    integer, intent(in)       :: n
        !+ Initial atomic energy level/state
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    if((q.eq.5).and.(n.eq.1)) then
        sigma = B5_ioniz_1_janev(eb)
        return
    endif

    if((q.eq.6).and.(n.eq.1)) then
        sigma = C6_ioniz_1_janev(eb)
        return
    endif

    sigma = Aq_ioniz_n_janev(eb, q, n)

end function Aq_ioniz_n

function Aq_ioniz(eb, q, n_max) result(sigma)
    !+ Calculates an array of total ionization cross sections for a Neutral Hydrogen atom
    !+in the n=1...n_max states colliding with a ion with charge `q` at energy `eb`
    !+
    !+@note Uses specialized cross sections if available else uses generic cross sections
    !+
    !+###Equation
    !+$$ A^{q+} + H(n=1..n_{max}) \rightarrow A^{(q-1)+} + H^+, q \gt 3 $$
    !+
    !+###References
    !+* Page 152 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 154 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 160 in Ref. 5 [[atomic_tables(module)]]
    real(Float64), intent(in)       :: eb
        !+ Relative collision energy [keV/amu]
    integer, intent(in)             :: q
        !+ Ion charge
    integer, intent(in)             :: n_max
        !+ Number of initial states n to calculate
    real(Float64), dimension(n_max) :: sigma
        !+ Array of cross sections where the n'th index refers to a ionization from the n'th state [\(cm^2\)]

    integer :: n

    do n=1,n_max
        sigma(n) = Aq_ioniz_n(eb, q, n)
    enddo

end function Aq_ioniz

!Impurity impact excitation
function Aq_excit_1_2_janev(eb, q) result(sigma)
    !+Calculates the excitation cross section for a neutral Hydrogen atom transitioning from
    !+the \(n=1\) state to the \(m=2\) state due to a collision an ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(1) \rightarrow A^{q+} + H(2), q \gt 4 $$
    !+
    !+###References
    !+* Page 132 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    integer, intent(in)       :: q
        !+ Ion charge
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(6), parameter :: A = [38.738d0, 37.033d0,   &
                                                   0.39862d0, 7.7582d-5, &
                                                   0.25402d0, -2.7418d0 ]
    real(Float64) :: Etil, xsi, qf

    qf = real(q)
    etil = eb/qf
    xsi = 2.**(0.5238*(1 - sqrt(2.0/qf)))
    sigma = qf*1.d-16*xsi*A(1)*(exp(-A(2)/etil)*log(1 + A(3)*etil)/etil &
          + A(4)*exp(-A(5)*etil)/etil**A(6))

end function Aq_excit_1_2_janev

function Aq_excit_1_3_janev(eb, q) result(sigma)
    !+Calculates the excitation cross section for a neutral Hydrogen atom transitioning from
    !+the \(n=1\) state to the \(m=3\) state due to a collision an ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(1) \rightarrow A^{q+} + H(3), q \gt 4 $$
    !+
    !+###References
    !+* Page 134 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    integer, intent(in)       :: q
        !+ Ion charge
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(6), parameter :: A = [4.3619d0, 57.451d0,  &
                                                   21.001d0, 2.3292d-4, &
                                                   0.083130d0, -2.2364d0 ]
    real(Float64) :: Etil, xsi, qf

    qf = real(q)
    etil = eb/qf
    xsi = 2.**(0.5238*(1 - sqrt(2.0/qf)))
    sigma = qf*1.d-16*xsi*A(1)*(exp(-A(2)/etil)*log(1 + A(3)*etil)/etil &
          + A(4)*exp(-A(5)*etil)/etil**A(6))

end function Aq_excit_1_3_janev

function Aq_excit_1_4_janev(eb, q) result(sigma)
    !+Calculates the excitation cross section for a neutral Hydrogen atom transitioning from
    !+the \(n=1\) state to the \(m=4\) state due to a collision an ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(1) \rightarrow A^{q+} + H(4), q \gt 4 $$
    !+
    !+###References
    !+* Page 134 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    integer, intent(in)       :: q
        !+ Ion charge
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(6), parameter :: A = [1.3730d0, 60.710d0,  &
                                                   31.797d0, 2.0207d-4, &
                                                   0.082513d0, -2.3055d0 ]
    real(Float64) :: Etil, xsi, qf

    qf = real(q)
    etil = eb/qf
    xsi = 2.**(0.5238*(1 - sqrt(2.0/qf)))
    sigma = qf*1.d-16*xsi*A(1)*(exp(-A(2)/etil)*log(1 + A(3)*etil)/etil &
          + A(4)*exp(-A(5)*etil)/etil**A(6))

end function Aq_excit_1_4_janev

function Aq_excit_1_5_janev(eb, q) result(sigma)
    !+Calculates the excitation cross section for a neutral Hydrogen atom transitioning from
    !+the \(n=1\) state to the \(m=5\) state due to a collision an ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(1) \rightarrow A^{q+} + H(5), q \gt 4 $$
    !+
    !+###References
    !+* Page 136 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    integer, intent(in)       :: q
        !+ Ion charge
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(6), parameter :: A = [0.56565d0, 67.333d0, &
                                                   55.290d0, 2.1595d-4, &
                                                   0.081624d0, -2.1971d0 ]
    real(Float64) :: Etil, xsi, qf

    qf = real(q)
    etil = eb/qf
    xsi = 2.**(0.5238*(1 - sqrt(2.0/qf)))
    sigma = qf*1.d-16*xsi*A(1)*(exp(-A(2)/etil)*log(1 + A(3)*etil)/etil &
          + A(4)*exp(-A(5)*etil)/etil**A(6))

end function Aq_excit_1_5_janev

function Aq_excit_1_janev(eb, q, m_max) result(sigma)
    !+Calculates an array of the excitation cross sections for a neutral Hydrogen atom transitioning from
    !+the \(n=1\) state to the m=1..`m_max` states due to a collision an ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(1) \rightarrow A^{q+} + H(m=1..m_{max}), q \gt 4 $$
    !+
    !+###References
    !+* Page 132 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 134 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 136 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in)       :: eb
        !+ Collision energy [keV]
    integer, intent(in)             :: q
        !+ Ion charge
    integer, intent(in)             :: m_max
        !+ Number of m states to calculate
    real(Float64), dimension(m_max) :: sigma
        !+ Array of cross sections where the m'th index refers to
        !+an excitation from the \(n=1\) state to the m'th state [\(cm^2\)]

    integer :: m

    sigma = 0.d0
    do m=1,m_max
        select case (m)
            case (1)
                sigma(1) = 0.d0
            case (2)
                sigma(2) = Aq_excit_1_2_janev(eb, q)
            case (3)
                sigma(3) = Aq_excit_1_3_janev(eb, q)
            case (4)
                sigma(4) = Aq_excit_1_4_janev(eb, q)
            case (5)
                sigma(5) = Aq_excit_1_5_janev(eb, q)
            case DEFAULT
                sigma(m) = sigma(5)*(5.0/m)**3.0
        end select
    enddo

end function Aq_excit_1_janev

function Aq_excit_2_3_janev(eb, q) result(sigma)
    !+Calculates the excitation cross section for a neutral Hydrogen atom transitioning from
    !+the \(n=2\) state to the \(m=3\) state due to a collision an ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(2) \rightarrow A^{q+} + H(3), q \gt 3 $$
    !+
    !+###References
    !+* Page 138 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    integer, intent(in)       :: q
        !+ Ion charge
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(6), parameter :: A = [358.03d0, 25.283d0,   &
                                                   1.4726d0, 0.014398d0, &
                                                   0.12207d0, -0.86210d0 ]
    real(Float64) :: etil, qf, xsi

    qf = real(q)
    etil = eb/qf
    xsi = 2.0**(0.5238*(1 - sqrt(2.0/qf)))
    sigma = qf*1.e-16*xsi*A(1)*(exp(-A(2)/etil)*log(1 + A(3)*etil)/etil &
          + A(4)*exp(-A(5)*etil)/etil**A(6))

end function Aq_excit_2_3_janev

function Aq_excit_2_4_janev(eb, q) result(sigma)
    !+Calculates the excitation cross section for a neutral Hydrogen atom transitioning from
    !+the \(n=2\) state to the \(m=4\) state due to a collision an ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(2) \rightarrow A^{q+} + H(4), q \gt 3 $$
    !+
    !+###References
    !+* Page 138 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    integer, intent(in)       :: q
        !+ Ion charge
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(6), parameter :: A = [50.744d0, 19.416d0,   &
                                                   4.0262d0, 0.014398d0, &
                                                   0.31584d0, -1.4799d0 ]
    real(Float64) :: etil, qf, xsi

    qf = real(q)
    etil = eb/qf
    xsi = 2.0**(0.5238*(1 - sqrt(2.0/qf)))
    sigma = qf*1.e-16*xsi*A(1)*(exp(-A(2)/etil)*log(1 + A(3)*etil)/etil &
          + A(4)*exp(-A(5)*etil)/etil**A(6))

end function Aq_excit_2_4_janev

function Aq_excit_2_5_janev(eb, q) result(sigma)
    !+Calculates the excitation cross section for a neutral Hydrogen atom transitioning from
    !+the \(n=2\) state to the \(m=5\) state due to a collision an ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(2) \rightarrow A^{q+} + H(5), q \gt 3 $$
    !+
    !+###References
    !+* Page 138 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    integer, intent(in)       :: q
        !+ Ion charge
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(6), parameter :: A = [18.264d0, 18.973d0,   &
                                                   2.9056d0, 0.013701d0, &
                                                   0.31711d0, -1.4775d0 ]
    real(Float64) :: etil, qf, xsi

    qf = real(q)
    etil = eb/qf
    xsi = 2.0**(0.5238*(1 - sqrt(2.0/qf)))
    sigma = qf*1.e-16*xsi*A(1)*(exp(-A(2)/etil)*log(1 + A(3)*etil)/etil &
          + A(4)*exp(-A(5)*etil)/etil**A(6))

end function Aq_excit_2_5_janev

function Aq_excit_2_6_janev(eb, q) result(sigma)
    !+Calculates the excitation cross section for a neutral Hydrogen atom transitioning from
    !+the \(n=2\) state to the \(m=6\) state due to a collision an ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(2) \rightarrow A^{q+} + H(6), q \gt 3 $$
    !+
    !+###References
    !+* Page 140 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    integer, intent(in)       :: q
        !+ Ion charge
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: A = 0.4610d0

    real(Float64) :: hi

    hi = 2.0**(0.397*(1.0 - sqrt(2.0/q)))

    sigma = A*hi*Aq_excit_2_5_janev(eb, q)

end function Aq_excit_2_6_janev

function Aq_excit_2_7_janev(eb, q) result(sigma)
    !+Calculates the excitation cross section for a neutral Hydrogen atom transitioning from
    !+the \(n=2\) state to the \(m=7\) state due to a collision an ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(2) \rightarrow A^{q+} + H(7), q \gt 3 $$
    !+
    !+###References
    !+* Page 140 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    integer, intent(in)       :: q
        !+ Ion charge
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: A = 0.2475d0

    real(Float64) :: hi

    hi = 2.0**(0.397*(1.0 - sqrt(2.0/q)))

    sigma = A*hi*Aq_excit_2_5_janev(eb, q)

end function Aq_excit_2_7_janev

function Aq_excit_2_8_janev(eb, q) result(sigma)
    !+Calculates the excitation cross section for a neutral Hydrogen atom transitioning from
    !+the \(n=2\) state to the \(m=8\) state due to a collision an ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(2) \rightarrow A^{q+} + H(8), q \gt 3 $$
    !+
    !+###References
    !+* Page 140 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    integer, intent(in)       :: q
        !+ Ion charge
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: A = 0.1465d0

    real(Float64) :: hi

    hi = 2.0**(0.397*(1.0 - sqrt(2.0/q)))

    sigma = A*hi*Aq_excit_2_5_janev(eb, q)

end function Aq_excit_2_8_janev

function Aq_excit_2_9_janev(eb, q) result(sigma)
    !+Calculates the excitation cross section for a neutral Hydrogen atom transitioning from
    !+the \(n=2\) state to the \(m=9\) state due to a collision an ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(2) \rightarrow A^{q+} + H(9), q \gt 3 $$
    !+
    !+###References
    !+* Page 140 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    integer, intent(in)       :: q
        !+ Ion charge
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: A = 0.092d0

    real(Float64) :: hi

    hi = 2.0**(0.397*(1.0 - sqrt(2.0/q)))

    sigma = A*hi*Aq_excit_2_5_janev(eb, q)

end function Aq_excit_2_9_janev

function Aq_excit_2_10_janev(eb, q) result(sigma)
    !+Calculates the excitation cross section for a neutral Hydrogen atom transitioning from
    !+the \(n=2\) state to the \(m=10\) state due to a collision an ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(2) \rightarrow A^{q+} + H(10), q \gt 3 $$
    !+
    !+###References
    !+* Page 140 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    integer, intent(in)       :: q
        !+ Ion charge
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: A = 0.0605d0

    real(Float64) :: hi

    hi = 2.0**(0.397*(1.0 - sqrt(2.0/q)))

    sigma = A*hi*Aq_excit_2_5_janev(eb, q)

end function Aq_excit_2_10_janev

function Aq_excit_2_janev(eb, q, m_max) result(sigma)
    !+Calculates an array of the excitation cross sections for a neutral Hydrogen atom transitioning from
    !+the \(n=2\) state to the m=1..`m_max` states due to a collision an ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(2) \rightarrow A^{q+} + H(m=1..m_{max}), q \gt 4, m \gt n $$
    !+
    !+###References
    !+* Page 138 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 140 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in)       :: eb
        !+ Collision energy [keV]
    integer, intent(in)             :: q
        !+ Ion charge
    integer, intent(in)             :: m_max
        !+ Number of m states to calculate
    real(Float64), dimension(m_max) :: sigma
        !+ Array of cross sections where the m'th index refers to
        !+an excitation from the \(n=2\) state to the m'th state [\(cm^2\)]

    integer :: m

    sigma = 0.d0
    do m=1,m_max
        select case (m)
            case (1)
                sigma(1) = 0.d0
            case (2)
                sigma(2) = 0.d0
            case (3)
                sigma(3) = Aq_excit_2_3_janev(eb, q)
            case (4)
                sigma(4) = Aq_excit_2_4_janev(eb, q)
            case (5)
                sigma(5) = Aq_excit_2_5_janev(eb, q)
            case (6)
                sigma(6) = Aq_excit_2_6_janev(eb, q)
            case (7)
                sigma(7) = Aq_excit_2_7_janev(eb, q)
            case (8)
                sigma(8) = Aq_excit_2_8_janev(eb, q)
            case (9)
                sigma(9) = Aq_excit_2_9_janev(eb, q)
            case (10)
                sigma(10) = Aq_excit_2_10_janev(eb, q)
            case DEFAULT
                sigma(m) = sigma(10)*(10.0/m)**3.0
        end select
    enddo

end function Aq_excit_2_janev

function Aq_excit_3_4_janev(eb, q) result(sigma)
    !+Calculates the excitation cross section for a neutral Hydrogen atom transitioning from
    !+the \(n=3\) state to the \(m=4\) state due to a collision an ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(3) \rightarrow A^{q+} + H(4), q \gt 3 $$
    !+
    !+###References
    !+* Page 142 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    integer, intent(in)       :: q
        !+ Ion charge
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(6), parameter :: A = [1247.5d0, 11.319d0,   &
                                                   2.6235d0, 0.068781d0, &
                                                   0.521176d0, -1.2722d0 ]

    real(Float64) :: Etil, qf, xsi

    qf = real(q)
    etil = eb/qf
    xsi = 2.0**(0.397*(1 - sqrt(2.0/qf)))
    sigma = qf*1.e-16*xsi*A(1)*(exp(-A(2)/etil)*log(1+A(3)*etil)/etil &
          + A(4)*exp(-A(5)*etil)/etil**A(6))

end function Aq_excit_3_4_janev

function Aq_excit_3_5_janev(eb, q) result(sigma)
    !+Calculates the excitation cross section for a neutral Hydrogen atom transitioning from
    !+the \(n=3\) state to the \(m=5\) state due to a collision an ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(3) \rightarrow A^{q+} + H(5), q \gt 3 $$
    !+
    !+###References
    !+* Page 142 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    integer, intent(in)       :: q
        !+ Ion charge
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(6), parameter :: A = [190.59d0, 11.096d0,   &
                                                   2.9098d0, 0.073307d0, &
                                                   0.54177d0, -1.2894d0 ]

    real(Float64) :: Etil, qf, xsi

    qf = real(q)
    etil = eb/qf
    xsi = 2.0**(0.397*(1 - sqrt(2.0/qf)))
    sigma = qf*1.e-16*xsi*A(1)*(exp(-A(2)/etil)*log(1+A(3)*etil)/etil &
          + A(4)*exp(-A(5)*etil)/etil**A(6))

end function Aq_excit_3_5_janev

function Aq_excit_3_6_janev(eb, q) result(sigma)
    !+Calculates the excitation cross section for a neutral Hydrogen atom transitioning from
    !+the \(n=3\) state to the \(m=6\) state due to a collision an ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(3) \rightarrow A^{q+} + H(6), q \gt 3 $$
    !+
    !+###References
    !+* Page 142 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    integer, intent(in)       :: q
        !+ Ion charge
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(6), parameter :: A = [63.494d0, 11.507d0,   &
                                                   4.3417d0, 0.077953d0, &
                                                   0.53461d0, -1.2881d0 ]

    real(Float64) :: Etil, qf, xsi

    qf = real(q)
    etil = eb/qf
    xsi = 2.0**(0.397*(1 - sqrt(2.0/qf)))
    sigma = qf*1.e-16*xsi*A(1)*(exp(-A(2)/etil)*log(1+A(3)*etil)/etil &
          + A(4)*exp(-A(5)*etil)/etil**A(6))

end function Aq_excit_3_6_janev

function Aq_excit_3_7_janev(eb, q) result(sigma)
    !+Calculates the excitation cross section for a neutral Hydrogen atom transitioning from
    !+the \(n=3\) state to the \(m=7\) state due to a collision an ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(3) \rightarrow A^{q+} + H(7), q \gt 3 $$
    !+
    !+###References
    !+* Page 144 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    integer, intent(in)       :: q
        !+ Ion charge
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: A = 0.4670d0

    real(Float64) :: hi

    hi=2.0**(0.397*(1.0 - sqrt(2.0/q)))
    sigma=hi*A*Aq_excit_3_6_janev(eb, q)

end function Aq_excit_3_7_janev

function Aq_excit_3_8_janev(eb, q) result(sigma)
    !+Calculates the excitation cross section for a neutral Hydrogen atom transitioning from
    !+the \(n=3\) state to the \(m=8\) state due to a collision an ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(3) \rightarrow A^{q+} + H(8), q \gt 3 $$
    !+
    !+###References
    !+* Page 144 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    integer, intent(in)       :: q
        !+ Ion charge
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: A = 0.2545d0

    real(Float64) :: hi

    hi=2.0**(0.397*(1.0 - sqrt(2.0/q)))
    sigma=hi*A*Aq_excit_3_6_janev(eb, q)

end function Aq_excit_3_8_janev

function Aq_excit_3_9_janev(eb, q) result(sigma)
    !+Calculates the excitation cross section for a neutral Hydrogen atom transitioning from
    !+the \(n=3\) state to the \(m=9\) state due to a collision an ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(3) \rightarrow A^{q+} + H(9), q \gt 3 $$
    !+
    !+###References
    !+* Page 144 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    integer, intent(in)       :: q
        !+ Ion charge
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: A = 0.1540d0

    real(Float64) :: hi

    hi=2.0**(0.397*(1.0 - sqrt(2.0/q)))
    sigma=hi*A*Aq_excit_3_6_janev(eb, q)

end function Aq_excit_3_9_janev

function Aq_excit_3_10_janev(eb, q) result(sigma)
    !+Calculates the excitation cross section for a neutral Hydrogen atom transitioning from
    !+the \(n=3\) state to the \(m=10\) state due to a collision an ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(3) \rightarrow A^{q+} + H(10), q \gt 3 $$
    !+
    !+###References
    !+* Page 144 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in) :: eb
        !+ Collision energy [keV]
    integer, intent(in)       :: q
        !+ Ion charge
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), parameter :: A = 0.1d0

    real(Float64) :: hi

    hi=2.0**(0.397*(1.0 - sqrt(2.0/q)))
    sigma=hi*A*Aq_excit_3_6_janev(eb, q)

end function Aq_excit_3_10_janev

function Aq_excit_3_janev(eb, q, m_max) result(sigma)
    !+Calculates an array of the excitation cross sections for a neutral Hydrogen atom transitioning from
    !+the \(n=3\) state to the m=1..`m_max` states due to a collision an ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(3) \rightarrow A^{q+} + H(m=1..m_{max}), q \gt 4, m \gt n $$
    !+
    !+###References
    !+* Page 142 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 144 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in)       :: eb
        !+ Collision energy [keV]
    integer, intent(in)             :: q
        !+ Ion charge
    integer, intent(in)             :: m_max
        !+ Number of m states to calculate
    real(Float64), dimension(m_max) :: sigma
        !+ Array of cross sections where the m'th index refers to
        !+an excitation from the \(n=3\) state to the m'th state [\(cm^2\)]

    integer :: m

    sigma = 0.d0
    do m=1,m_max
        select case (m)
            case (1)
                sigma(1) = 0.d0
            case (2)
                sigma(2) = 0.d0
            case (3)
                sigma(3) = 0.d0
            case (4)
                sigma(4) = Aq_excit_3_4_janev(eb, q)
            case (5)
                sigma(5) = Aq_excit_3_5_janev(eb, q)
            case (6)
                sigma(6) = Aq_excit_3_6_janev(eb, q)
            case (7)
                sigma(7) = Aq_excit_3_7_janev(eb, q)
            case (8)
                sigma(8) = Aq_excit_3_8_janev(eb, q)
            case (9)
                sigma(9) = Aq_excit_3_9_janev(eb, q)
            case (10)
                sigma(10) = Aq_excit_3_10_janev(eb, q)
            case DEFAULT
                sigma(m) = sigma(10)*(10.0/m)**3.0
        end select
    enddo

end function Aq_excit_3_janev

function Aq_excit_n_janev(eb, q, n, m_max) result(sigma)
    !+Calculates an array of the generic excitation cross sections for a neutral Hydrogen atom transitioning from
    !+the `n` state to the m=1..`m_max` states due to a collision an ion with charge `q` at energy `eb`
    !+
    !+###Equation
    !+$$ A^{q+} + H(n) \rightarrow A^{q+} + H(m=1..m_{max}), q \gt 3, m \gt n, n \gt 3 $$
    !+
    !+###References
    !+* Page 146 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in)       :: eb
        !+ Collision energy [keV]
    integer, intent(in)             :: q
        !+ Ion charge
    integer, intent(in)             :: n
        !+ Initial atomic energy level/state
    integer, intent(in)             :: m_max
        !+ Number of m states to calculate
    real(Float64), dimension(m_max) :: sigma
        !+ Array of cross sections where the m'th index refers to
        !+an excitation from the `n` state to the m'th state [\(cm^2\)]

    integer :: m
    real(Float64) :: nf, mf, qf, etil, hi, s
    real(Float64) :: D, A, G, L, F, H, y, zpl, zmi, C2pl, C2mi

    nf = real(n)
    qf = real(q)
    sigma = 0.d0
    m_loop: do m=1,m_max
        mf = real(m)
        if(n.ge.m) then
            sigma(m) = 0.d0
            cycle m_loop
        endif
        etil = eb/(25.0*qf)
        hi = 2.0**(0.322*(1.0 - sqrt(2.0/qf)))
        s = (mf - nf)

        D = exp(-1.0/(nf*mf*etil**2))
        A = 8.0/(3.0*s) * (mf/(s*nf))**3 * (0.184 - 0.04/s**(2.0/3.0)) * &
           (1.0 - 0.2*s/(nf*mf))**(1.0 + 2.0*s)
        G = 0.5*(etil*nf**2.0 / (mf - 1.0/mf))**3.0
        L = log(1.0 + 0.53*etil**2.0 * nf*(mf - 2.0/mf )/(1.0 + 0.4*etil))
        F = (1.0 - 0.3*s*D/(nf*mf))**(1.0 + 2.0*s)

        y = 1.0/(1.0 - D*log(18*s)/(4.0*s))
        zpl = 2.0/(etil*nf**2*(sqrt(2.0 - nf**2/mf**2) + 1.0))
        zmi = 2.0/(etil*nf**2*(sqrt(2.0 - nf**2/mf**2) - 1.0))
        C2pl = zpl**2*log(1.0 + 2.0*zpl/3.0)/(2.0*y + 3.0*zpl/2.0)
        C2mi = zmi**2*log(1.0 + 2.0*zmi/3.0)/(2.0*y + 3.0*zmi/2.0)
        H = C2mi - C2pl

        sigma(m) = q*hi*8.86e-17*nf**4/etil*(A*D*L+F*G*H)
    enddo m_loop

end function Aq_excit_n_janev

function Aq_excit_n(eb, q, n, m_max) result(sigma)
    !+Calculates an array of the excitation cross sections for a neutral Hydrogen atom transitioning from
    !+the `n` state to the m=1..`m_max` states due to a collision an ion with charge `q` at energy `eb`
    !+
    !+@note Uses specialized cross sections if available else uses generic cross sections
    !+
    !+###Equation
    !+$$ A^{q+} + H(n) \rightarrow A^{q+} + H(m=1..m_{max}), q \gt 3, m \gt n, n \gt 3 $$
    !+
    !+###References
    !+* Page 132 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 134 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 136 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 138 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 140 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 142 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 142 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 144 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 146 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in)       :: eb
        !+ Collision energy [keV]
    integer, intent(in)             :: q
        !+ Ion charge
    integer, intent(in)             :: n
        !+ Initial atomic energy level/state
    integer, intent(in)             :: m_max
        !+ Number of m states to calculate
    real(Float64), dimension(m_max) :: sigma
        !+ Array of cross sections where the m'th index refers to
        !+an excitation from the `n` state to the m'th state [\(cm^2\)]

    select case (n)
        case (0)
            stop
        case (1)
            sigma = Aq_excit_1_janev(eb, q, m_max)
        case (2)
            sigma = Aq_excit_2_janev(eb, q, m_max)
        case (3)
            sigma = Aq_excit_3_janev(eb, q, m_max)
        case DEFAULT
            sigma = Aq_excit_n_janev(eb, q, n, m_max)
    end select

end function Aq_excit_n

function Aq_excit_n_m(eb ,q, n, m) result(sigma)
    !+Calculates the excitation cross section for a neutral Hydrogen atom transitioning from
    !+the `n`\(\rightarrow\)`m` state due to a collision an ion with charge `q` at energy `eb`
    !+
    !+@note Uses specialized cross sections if available else uses generic cross sections
    !+
    !+###Equation
    !+$$ A^{q+} + H(n) \rightarrow A^{q+} + H(m), q \gt 3, m \gt n $$
    !+
    !+###References
    !+* Page 132 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 134 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 136 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 138 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 140 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 142 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 142 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 144 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 146 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in)       :: eb
        !+ Collision energy [keV]
    integer, intent(in)             :: q
        !+ Ion charge
    integer, intent(in)             :: n
        !+ Initial atomic energy level/state
    integer, intent(in)             :: m
        !+ Final atomic energy level/state
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(12) :: sigma_m

    sigma_m = Aq_excit_n(eb, q, n, 12)
    if(m.le.0) then
        sigma = sum(sigma_m)
    else
        sigma = sigma_m(m)
    endif

end function Aq_excit_n_m

function Aq_excit(eb, q, n_max, m_max) result(sigma)
    !+Calculates an matrix of the excitation cross sections for a neutral Hydrogen atom transitioning from
    !+the n=1..`n_max` state to the m=1..`m_max` states due to a collision an ion with charge `q` at energy `eb`
    !+
    !+@note Uses specialized cross sections if available else uses generic cross sections
    !+
    !+###Equation
    !+$$ A^{q+} + H(n=1..n_{max}) \rightarrow A^{q+} + H(m=1..m_{max}), q \gt 3, m \gt n$$
    !+
    !+###References
    !+* Page 132 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 134 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 136 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 138 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 140 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 142 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 142 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 144 in Ref. 5 [[atomic_tables(module)]]
    !+* Page 146 in Ref. 5 [[atomic_tables(module)]]
    !+
    real(Float64), intent(in)              :: eb
        !+ Collision energy [keV]
    integer, intent(in)                    :: q
        !+ Ion charge
    integer, intent(in)                    :: n_max
        !+ Number of n states to calculate
    integer, intent(in)                    :: m_max
        !+ Number of m states to calculate
    real(Float64), dimension(n_max, m_max) :: sigma
        !+ Matrix of cross sections where the subscripts refers to
        !+an excitation from the `n` state to the m'th state: Aq_excit[n,m] [\(cm^2\)]

    real(Float64), dimension(12,12) :: sigma_full
    integer :: n, m

    do n=1,12
        sigma_full(n,:) = Aq_excit_n(eb, q, n, 12)
    enddo

    sigma = sigma_full(1:n_max,1:m_max)

end function Aq_excit

function d_d_fusion_t(eb) result(sigma)
    !+Calculates total cross section at a given Deuterium energy, `eb`,
    !+for the Tritium branch of Deuterium-Deutrium nuclear reactions
    !+
    !+###Equation
    !+$$ D + D \rightarrow T(1.01 MeV) + p(3.02 MeV) (50%)$$
    !+
    !+###References
    !+* Equations 8-9
    !+* Table IV in Ref. 7 [[atomic_tables(module)]]
    real(Float64), intent(in) :: eb
        !+ Deuterium energy [keV]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(5), parameter :: A = [ 5.5576d4,  2.1054d2,  &
                                                   -3.2638d-2, 1.4987d-6, &
                                                    1.8181d-10 ]
    real(Float64), dimension(4), parameter :: B = [0.d0,0.d0,0.d0,0.d0]
    real(Float64), parameter :: Bg = 31.3970
    real(Float64) :: S, E

    E = min(max(eb,0.5),5000.0)

    S = (A(1) + E*(A(2) + E*(A(3) + E*(A(4) + E*A(5))))) / &
        (1    + E*(B(1) + E*(B(2) + E*(B(3) + E*B(4)))))

    sigma = (1.0d-27)*(S/(E*exp(Bg/sqrt(E))))

end function d_d_fusion_t

function d_d_fusion_he(eb) result(sigma)
    !+Calculates total cross section at a given deuterium energy, `eb`,
    !+for the Helium-3 branch of Deuterium-Deutrium nuclear reactions
    !+
    !+###Equation
    !+$$ D + D \rightarrow He^3(0.82 MeV) + n(2.45 MeV) (50%)$$
    !+
    !+###References
    !+* Equations 8-9
    !+* Table IV in Ref. 7 [[atomic_tables(module)]]
    real(Float64), intent(in) :: eb
        !+ Deuterium energy [keV]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(5), parameter :: A = [ 5.3701d4,  3.3027d2,  &
                                                   -1.2706d-1, 2.9327d-5, &
                                                   -2.5151d-9 ]
    real(Float64), dimension(4), parameter :: B = [0.d0,0.d0,0.d0,0.d0]
    real(Float64), parameter :: Bg = 31.3970
    real(Float64) :: S, E

    E = min(max(eb,0.5),4900.0)

    S = (A(1) + E*(A(2) + E*(A(3) + E*(A(4) + E*A(5))))) / &
        (1    + E*(B(1) + E*(B(2) + E*(B(3) + E*B(4)))))

    sigma = (1.0d-27)*(S/(E*exp(Bg/sqrt(E))))

end function d_d_fusion_he

function d_t_fusion(eb) result(sigma)
    !+Calculates total cross section at a given deuterium energy, `eb`,
    !+for Deuterium-Tritium nuclear reactions in the range [0.5-550 keV]
    !+
    !+###Equation
    !+$$ D + T \rightarrow He^4(3.5 MeV) + n(14.1 MeV)$$
    !+
    !+###References
    !+* Equations 8-9
    !+* Table IV, VI in Ref. 7 [[atomic_tables(module)]]
    real(Float64), intent(in) :: eb
        !+ Deuterium energy [keV]
    real(Float64)             :: sigma
        !+ Cross Section [\(cm^2\)]

    real(Float64), dimension(5), parameter :: A1 = [ 6.927d4,  7.454d8, &
                                                    2.050d6, 5.2002d4,  &
                                                    0.d0 ]
    real(Float64), dimension(4), parameter :: B1 = [ 6.38d1,  -9.95d-1, &
                                                    6.981d-5, 1.728d-4  ]

    real(Float64), dimension(5), parameter :: A2 = [-1.4714d6, 0.d0,  &
                                                     0.d0, 0.d0, 0.d0 ]
    real(Float64), dimension(4), parameter :: B2 = [-8.4127d-3, 4.7983d-6, &
                                                    -1.0748d-9, 8.5184d-14 ]
    real(Float64), parameter :: Bg = 34.3827

    real(Float64), dimension(5) :: A
    real(Float64), dimension(4) :: B
    real(Float64) :: S, E

    E = min(max(eb,0.5),4700.0)

    if(E.le.530.0) then
        A = A1
        B = B1
    else
        A = A2
        B = B2
    endif

    S = (A(1) + E*(A(2) + E*(A(3) + E*(A(4) + E*A(5))))) / &
        (1    + E*(B(1) + E*(B(2) + E*(B(3) + E*B(4)))))

    sigma = (1.0d-27)*(S/(E*exp(Bg/sqrt(E))))

end function d_t_fusion

function simpsons_rule(f, dx) result(I)
    !+ Performs 1D integration using Simpsons rule
    !+
    !+ ###References
    !+* [Simpson's rule](http://mathworld.wolfram.com/SimpsonsRule.html)
    real(Float64), dimension(:), intent(in) :: f
        !+ Array of equally spaced \(f(x)\) values
    real(Float64), intent(in)               :: dx
        !+ Spacing between x values
    real(Float64)                           :: I

    integer :: s, ii

    s = size(f)
    I = 0.d0
    if(mod(s,2).eq.1) then
        write(*,'(a)') "Length of array must be even"
        return
    endif

    I = f(1)
    do ii=2,s-1
        if(mod(ii,2).eq.1) then
            I = I + 4.0*f(ii)
        else
            I = I + 2.0*f(ii)
        endif
    enddo
    I = I + f(s)
    I = (dx/3.0)*I

end function simpsons_rule

subroutine bt_maxwellian_eb(fn, T, eb, am, ab, rate)
    !+ Calculates Maxwellian reaction rate for a beam with atomic mass `ab` and energy `eb`
    !+firing into a target with atomic mass `am` and temperature `T` which has a cross section given by the function `fn`
    interface
        function fn(a)
            !+Cross section function
            real(8)              :: fn !sigma
            real(8), intent(in)  :: a !eb
        end function fn
    end interface
    real(Float64), intent(in)  :: T
        !+Target temperature [keV]
    real(Float64), intent(in)  :: eb
        !+Beam energy [keV]
    real(Float64), intent(in)  :: am
        !+Target atomic mass [amu]
    real(Float64), intent(in)  :: ab
        !+Beam atomic mass [amu]
    real(Float64), intent(out) :: rate
        !+Reaction Rate [\(cm^3/s\)]

    integer :: n_vr
    real(Float64) :: vr_max, dvr
    real(Float64), dimension(32) :: vr
    real(Float64), dimension(32) :: fr
    integer :: n_vz
    real(Float64) :: vz_max,dvz
    real(Float64), dimension(62) :: vz
    real(Float64), dimension(62) :: fz
    real(Float64) :: T_per_amu, eb_per_amu, ared, sig, sig_eff
    real(Float64) :: zb, u2_to_erel, u2, erel, v_therm, dE

    integer :: i, j

    n_vr = 32
    vr_max = 4.d0
    dvr = vr_max/(n_vr - 1.d0)
    do i=1,n_vr
        vr(i) = (i-1)*dvr
    enddo

    n_vz = 62
    vz_max = 4.d0
    dvz = 2.0*vz_max/(n_vz - 1.d0)
    do i=1,n_vz
        vz(i) = (i-1)*dvz - vz_max
    enddo

    T_per_amu = max(T, 1.d-6)/am
    eb_per_amu = eb/ab
    ared = am*ab/(am + ab)

    v_therm = 1.384d6 * sqrt(T_per_amu*1.d3)
    zb = sqrt(eb_per_amu/T_per_amu)
    u2_to_erel = ared*T_per_amu

    fz = 0.d0
    fr = 0.d0

    do i=1,n_vz
        do j=1,n_vr
            u2 = (zb - vz(i))**2.0 + vr(j)**2.0
            erel = u2_to_erel*u2
            sig = fn(erel)
            fr(j) = sig*sqrt(u2)*exp(-(vz(i)**2.0 + vr(j)**2.0))*vr(j)
        enddo
        fz(i) = simpsons_rule(fr, dvr)
    enddo

    sig_eff = (2.0/sqrt(PI))*simpsons_rule(fz, dvz)
    rate = sig_eff*v_therm

end subroutine bt_maxwellian_eb

subroutine bt_maxwellian_n(fn, T, eb, am, ab, n, rate)
    !+ Calculates Maxwellian reaction rate for a beam with atomic mass `ab`, energy `eb`, and energy level `n`
    !+firing into a target with atomic mass `am` and temperature `T` which has a cross section given by the function `fn`
    interface
        function fn(a, b)
            !+Cross section function
            real(8)              :: fn !sigma
            real(8), intent(in)  :: a !eb
            integer, intent(in)  :: b !n
        end function fn
    end interface
    real(Float64), intent(in)  :: T
        !+Target temperature [keV]
    real(Float64), intent(in)  :: eb
        !+Beam energy [keV]
    real(Float64), intent(in)  :: am
        !+Target atomic mass [amu]
    real(Float64), intent(in)  :: ab
        !+Beam atomic mass [amu]
    integer, intent(in)        :: n
        !+Initial atomic energy level/state
    real(Float64), intent(out) :: rate
        !+Reaction Rate [\(cm^3/s\)]

    logical :: dxc
    integer :: n_vr
    real(Float64) :: vr_max, dvr
    real(Float64), dimension(32) :: vr
    real(Float64), dimension(32) :: fr
    integer :: n_vz
    real(Float64) :: vz_max,dvz
    real(Float64), dimension(62) :: vz
    real(Float64), dimension(62) :: fz
    real(Float64) :: T_per_amu, eb_per_amu, ared, sig, sig_eff
    real(Float64) :: zb, u2_to_erel, u2, erel, v_therm, dE

    integer :: i, j

    n_vr = 32
    vr_max = 4.d0
    dvr = vr_max/(n_vr - 1.d0)
    do i=1,n_vr
        vr(i) = (i-1)*dvr
    enddo

    n_vz = 62
    vz_max = 4.d0
    dvz = 2.0*vz_max/(n_vz - 1.d0)
    do i=1,n_vz
        vz(i) = (i-1)*dvz - vz_max
    enddo

    T_per_amu = max(T, 1.d-6)/am
    eb_per_amu = eb/ab
    ared = am*ab/(am + ab)
    dE = (13.6d-3)/(n**2.0)

    v_therm = 1.384d6 * sqrt(T_per_amu*1.d3)
    zb = sqrt(eb_per_amu/T_per_amu)
    u2_to_erel = ared*T_per_amu
    if(ared.lt.0.5) ared = 1.0 !for electron interactions

    fz = 0.d0
    fr = 0.d0

    do i=1,n_vz
        do j=1,n_vr
            u2 = (zb - vz(i))**2.0 + vr(j)**2.0
            erel = u2_to_erel*u2
            if(erel.ge.dE) then
                sig = fn(erel/ared,n)
            else
                sig = 0.d0
            endif
            fr(j) = sig*sqrt(u2)*exp(-(vz(i)**2.0 + vr(j)**2.0))*vr(j)
        enddo
        fz(i) = simpsons_rule(fr, dvr)
    enddo

    sig_eff = (2.0/sqrt(PI))*simpsons_rule(fz, dvz)
    rate = sig_eff*v_therm

end subroutine bt_maxwellian_n

subroutine bt_maxwellian_q_n(fqn, q, T, eb, am, ab, n, rate)
    !+ Calculates Maxwellian reaction rate for a beam with atomic mass `ab`, energy `eb`, and energy level `n`
    !+firing into a target with atomic mass `am`, temperature `T`, and charge `q`  which has a cross section given by the function `fqn`
    interface
        function fqn(a, b, c)
            !+Cross section function
            real(8)             :: fqn !sigma
            real(8), intent(in) :: a !eb
            integer, intent(in) :: b !q
            integer, intent(in) :: c !n
        end function fqn
    end interface
    integer, intent(in)        :: q
        !+Target charge
    real(Float64), intent(in)  :: T
        !+Target temperature [keV]
    real(Float64), intent(in)  :: eb
        !+Beam energy [keV]
    real(Float64), intent(in)  :: am
        !+Target atomic mass [amu]
    real(Float64), intent(in)  :: ab
        !+Beam atomic mass [amu]
    integer, intent(in)        :: n
        !+Initial atomic energy level/state
    real(Float64), intent(out) :: rate
        !+Reaction Rate [\(cm^3/s\)]

    integer :: n_vr
    real(Float64) :: vr_max, dvr
    real(Float64), dimension(32) :: vr
    real(Float64), dimension(32) :: fr
    integer :: n_vz
    real(Float64) :: vz_max, dvz
    real(Float64), dimension(62) :: vz
    real(Float64), dimension(62) :: fz
    real(Float64) :: T_per_amu, eb_per_amu, ared, sig, sig_eff
    real(Float64) :: zb, u2_to_erel, u2, erel, v_therm, dE

    integer :: i, j

    n_vr = 32
    vr_max = 4.d0
    dvr = vr_max/(n_vr - 1.d0)
    do i=1,n_vr
        vr(i) = (i-1)*dvr
    enddo

    n_vz = 62
    vz_max = 4.d0
    dvz = 2.0*vz_max/(n_vz - 1.d0)
    do i=1,n_vz
        vz(i) = (i-1)*dvz - vz_max
    enddo

    T_per_amu = max(T, 1.d-6)/am
    eb_per_amu = eb/ab
    ared = am*ab/(am + ab)
    dE = (13.6d-3)/(n**2.0)

    v_therm = 1.384d6 * sqrt(T_per_amu*1.d3)
    zb = sqrt(eb_per_amu/T_per_amu)
    u2_to_erel = ared*T_per_amu
    if(ared.lt.0.5) ared = 1.0

    fz = 0.d0
    fr = 0.d0

    do i=1,n_vz
        do j=1,n_vr
            u2 = (zb - vz(i))**2.0 + vr(j)**2.0
            erel = u2_to_erel*u2
            if(erel.ge.dE) then
                sig = fqn(erel/ared, q, n)
            else
                sig = 0.d0
            endif
            fr(j) = sig*sqrt(u2)*exp(-(vz(i)**2.0 + vr(j)**2.0))*vr(j)
        enddo
        fz(i) = simpsons_rule(fr, dvr)
    enddo

    sig_eff = (2.0/sqrt(PI))*simpsons_rule(fz, dvz)
    rate = sig_eff*v_therm

end subroutine bt_maxwellian_q_n

subroutine bt_maxwellian_n_m(fnm, T, eb, am, ab, n, m, rate, deexcit)
    !+ Calculates Maxwellian reaction rate for a `n`\(\rightarrow)`m` transition due to a beam with atomic mass `ab` and energy `eb`
    !+firing into a target with atomic mass `am` and temperature `T` which has a cross section given by the function `fnm`
    interface
        function fnm(a, b, c)
            !+Cross section function
            real(8)             :: fnm !sigma
            real(8), intent(in) :: a !eb
            integer, intent(in) :: b !n
            integer, intent(in) :: c !m
        end function fnm
    end interface
    real(Float64), intent(in)     :: T
        !+Target temperature [keV]
    real(Float64), intent(in)     :: eb
        !+Beam energy [keV]
    real(Float64), intent(in)     :: am
        !+Target atomic mass [amu]
    real(Float64), intent(in)     :: ab
        !+Beam atomic mass [amu]
    integer, intent(in)           :: n
        !+Initial atomic energy level/state
    integer, intent(in)           :: m
        !+Final atomic energy level/state
    real(Float64), intent(out)    :: rate
        !+Reaction Rate [\(cm^3/s\)]
    logical, intent(in), optional :: deexcit
        !+Calculate de-excitation reaction rate

    logical :: dxc
    integer :: n_vr
    real(Float64) :: vr_max, dvr
    real(Float64), dimension(32) :: vr
    real(Float64), dimension(32) :: fr
    integer :: n_vz
    real(Float64) :: vz_max, dvz
    real(Float64), dimension(62) :: vz
    real(Float64), dimension(62) :: fz
    real(Float64) :: T_per_amu, eb_per_amu, ared, sig, sig_eff
    real(Float64) :: zb, u2_to_erel, u2, erel, dE, factor, En, Em, v_therm

    integer :: i, j

    if(present(deexcit)) then
        dxc = deexcit
    else
        dxc = .False.
    endif

    n_vr = 32
    vr_max = 4.d0
    dvr = vr_max/(n_vr - 1.d0)
    do i=1,n_vr
        vr(i) = (i-1)*dvr
    enddo

    n_vz = 62
    vz_max = 4.d0
    dvz = 2.0*vz_max/(n_vz - 1.d0)
    do i=1,n_vz
        vz(i) = (i-1)*dvz - vz_max
    enddo

    En = (13.6d-3)*(1.0 - (1.d0/n)**2.0)
    Em = (13.6d-3)*(1.0 - (1.d0/m)**2.0)
    dE = Em - En

    T_per_amu = max(T, 1.d-6)/am
    eb_per_amu = eb/ab
    ared = am*ab/(am + ab)

    v_therm = 1.384d6 * sqrt(T_per_amu*1.d3)
    zb = sqrt(eb_per_amu/T_per_amu)
    u2_to_erel = ared*T_per_amu
    if(ared.lt.0.5) ared = 1.0

    fz = 0.d0
    fr = 0.d0

    do i=1,n_vz
        do j=1,n_vr
            u2 = (zb - vz(i))**2.0 + vr(j)**2.0
            erel = u2_to_erel*u2
            if(dxc) then
                factor = (erel + dE)/erel
                erel = erel + dE
            else
                factor = 1.0
            endif
            if(erel.ge.dE) then
                sig = fnm(erel/ared, n, m)
            else
                sig = 0.d0
            endif
            fr(j) = factor*sig*sqrt(u2)*exp(-(vz(i)**2.0 + vr(j)**2.0))*vr(j)
        enddo
        fz(i) = simpsons_rule(fr, dvr)
    enddo

    sig_eff = (2.0/sqrt(PI))*simpsons_rule(fz, dvz)
    rate = sig_eff*v_therm
    if(dxc) rate = rate*(real(n)/real(m))**2.0

end subroutine bt_maxwellian_n_m

subroutine bt_maxwellian_q_n_m(fqnm, q, T, eb, am, ab, n, m, rate, deexcit)
    !+ Calculates Maxwellian reaction rate for a `n`\(\rightarrow)`m` transition due to a beam with atomic mass `ab` and energy `eb`
    !+firing into a target with atomic mass `am`, temperature `T`, and charge `q` which has a cross section given by the function `fqnm`
    interface
        function fqnm(a, b, c, d)
            !+Cross section function
            real(8)             :: fqnm !sigma
            real(8), intent(in) :: a !eb
            integer, intent(in) :: b !q
            integer, intent(in) :: c !n
            integer, intent(in) :: d !m
        end function fqnm
    end interface
    integer, intent(in)           :: q
        !+Target charge
    real(Float64), intent(in)     :: T
        !+Target temperature [keV]
    real(Float64), intent(in)     :: eb
        !+Beam energy [keV]
    real(Float64), intent(in)     :: am
        !+Target atomic mass [amu]
    real(Float64), intent(in)     :: ab
        !+Beam atomic mass [amu]
    integer, intent(in)           :: n
        !+Initial atomic energy level/state
    integer, intent(in)           :: m
        !+Final atomic energy level/state
    real(Float64), intent(out)    :: rate
        !+Reaction Rate [\(cm^3/s\)]
    logical, intent(in), optional :: deexcit
        !+Calculate de-excitation reaction rate

    logical :: dxc
    integer :: n_vr
    real(Float64) :: vr_max, dvr
    real(Float64), dimension(32) :: vr
    real(Float64), dimension(32) :: fr
    integer :: n_vz
    real(Float64) :: vz_max, dvz
    real(Float64), dimension(62) :: vz
    real(Float64), dimension(62) :: fz
    real(Float64) :: T_per_amu, eb_per_amu, ared, sig, sig_eff
    real(Float64) :: zb, u2_to_erel, u2, erel, dE, factor, En, Em, v_therm

    integer :: i, j

    if(present(deexcit)) then
        dxc = deexcit
    else
        dxc = .False.
    endif

    n_vr = 32
    vr_max = 4.d0
    dvr = vr_max/(n_vr - 1.d0)
    do i=1,n_vr
        vr(i) = (i-1)*dvr
    enddo

    n_vz = 62
    vz_max = 4.d0
    dvz = 2.0*vz_max/(n_vz - 1.d0)
    do i=1,n_vz
        vz(i) = (i-1)*dvz - vz_max
    enddo

    En = (13.6d-3)*(1.0 - (1.d0/n)**2.0)
    Em = (13.6d-3)*(1.0 - (1.d0/m)**2.0)
    dE = Em - En

    T_per_amu = max(T, 1.d-6)/am
    eb_per_amu = eb/ab
    ared = am*ab/(am + ab)

    v_therm = 1.384d6 * sqrt(T_per_amu*1.d3)
    zb = sqrt(eb_per_amu/T_per_amu)
    u2_to_erel = ared*T_per_amu
    if(ared.lt.0.5) ared = 1.0

    fz = 0.d0
    fr = 0.d0

    do i=1,n_vz
        do j=1,n_vr
            u2 = (zb - vz(i))**2.0 + vr(j)**2.0
            erel = u2_to_erel*u2
            if(dxc) then
                factor = (erel + dE)/erel
                erel = erel + dE
            else
                factor = 1.0
            endif
            if(erel.ge.dE) then
                sig = fqnm(erel/ared, q, n, m)
            else
                sig = 0.d0
            endif
            fr(j) = factor*sig*sqrt(u2)*exp(-(vz(i)**2.0 + vr(j)**2.0))*vr(j)
        enddo
        fz(i) = simpsons_rule(fr, dvr)
    enddo

    sig_eff = (2.0/sqrt(PI))*simpsons_rule(fz, dvz)
    rate = sig_eff*v_therm
    if(dxc) rate = rate*(real(n)/real(m))**2.0

end subroutine bt_maxwellian_q_n_m

subroutine write_einstein(id, n_max, m_max)
    !+ Write Einstein coefficients to HDF5 file
    integer(HID_T), intent(inout) :: id
        !+ HDF5 file or group object id
    integer, intent(in)           :: n_max
        !+ Number of initial atomic energy states to calculate
    integer, intent(in)           :: m_max
        !+ Number of final atomic energy states to calculate

    real(Float64), dimension(n_max,m_max) :: ein

    integer(HID_T) :: gid
    integer(HSIZE_T), dimension(1) :: dim1
    integer(HSIZE_T), dimension(2) :: dim2
    integer :: error

    if(verbose) then
        ein(:,:) = EINSTEIN(1:n_max,1:m_max)

        call h5gcreate_f(id, "spontaneous", gid, error)

        dim1 = [1]
        dim2 = [n_max, m_max]

        call h5ltmake_dataset_int_f(gid, "n_max", 0, dim1, [n_max], error)
        call h5ltmake_dataset_int_f(gid, "m_max", 0, dim1, [m_max], error)

        call h5ltmake_compressed_dataset_double_f(gid, "einstein", 2, dim2, ein, error)

        call h5ltset_attribute_string_f(id, "spontaneous", "description", &
             "Atomic rates for spontaneous emission/deexcitation", error)
        call h5ltset_attribute_string_f(gid, "n_max", "description", &
             "Number of initial energy levels", error)
        call h5ltset_attribute_string_f(gid, "m_max", "description", &
             "Number of final energy levels", error)

        call h5ltset_attribute_string_f(gid, "einstein", "description", &
             "n/m resolved einstein coefficients: einstein(n,m)", error)
        call h5ltset_attribute_string_f(gid, "einstein", "units", "1/s", error)
        call h5ltset_attribute_string_f(gid, "einstein", "reaction", &
             "H(n) -> H(m) + ph, n > m", error)

        call h5gclose_f(gid, error)
    endif

end subroutine write_einstein

subroutine write_bb_H_H(id, namelist_file, n_max, m_max)
    !+ Write Hydrogen-Hydrogen interaction cross sections to a HDF5 file
    integer(HID_T), intent(inout) :: id
        !+ HDF5 file or group object id
    character(len=*), intent(in)  :: namelist_file
        !+ Namelist file that contains settings
    integer, intent(in)           :: n_max
        !+ Number of initial atomic energy states to calculate
    integer, intent(in)           :: m_max
        !+ Number of final atomic energy states to calculate

    logical :: calculate
    real(Float64) :: emin
    real(Float64) :: emax
    integer :: nenergy

    real(Float64) :: eb
    real(Float64) :: dlogE
    real(Float64), dimension(:), allocatable :: ebarr

    real(Float64), dimension(:,:), allocatable :: ioniz
    real(Float64), dimension(:,:,:), allocatable :: cx, excit

    integer(HID_T) :: gid
    integer(HSIZE_T), dimension(1) :: dim1
    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(3) :: dim3

    integer :: i, cnt, error
    logical :: exis

    NAMELIST /H_H_cross/ calculate, nenergy, emin, emax

    calculate = .True.; nenergy = 200; emin = 1.d-3 ; emax = 8.d2

    inquire(file=namelist_file,exist=exis)
    if(.not.exis) then
        write(*,'(a,a)') 'WRITE_BB_H_H: Input file does not exist: ',trim(namelist_file)
        write(*,'(a)') 'Continuing with default settings...'
    else
        open(13,file=namelist_file)
        read(13,NML=H_H_cross)
        close(13)
    endif

    if(.not.calculate) return

    allocate(ebarr(nenergy))
    allocate(ioniz(n_max,nenergy))
    allocate(cx(n_max,m_max,nenergy))
    allocate(excit(n_max,m_max,nenergy))

    ebarr = 0.d0
    ioniz = 0.d0
    cx = 0.d0
    excit = 0.d0

    if(verbose) then
        write(*,'(a)') "---- H-H cross sections settings ----"
        write(*,'(T2,"Emin = ",e9.2, " keV")') emin
        write(*,'(T2,"Emax = ",e9.2, " keV")') emax
        write(*,'(T2,"Nenergy = ", i4)') nenergy
        write(*,*) ''
    endif

    cnt = 0
    dlogE = (log10(emax) - log10(emin))/(nenergy - 1)
    !$OMP PARALLEL DO private(i, eb)
    do i=istart, nenergy, istep
        eb = 10.d0**(log10(emin) + (i-1)*dlogE)
        ebarr(i) = eb

        cx(:,:,i) = p_cx(eb, n_max, m_max)
        excit(:,:,i) = p_excit(eb, n_max, m_max)
        ioniz(:,i) = p_ioniz(eb, n_max)
        cnt = cnt + 1
        if(verbose) WRITE(*,'(f7.2,"%",a,$)') 100*cnt*istep/real(nenergy),char(13)
    enddo
    !$OMP END PARALLEL DO

#ifdef _MPI
    call parallel_sum(ebarr)
    call parallel_sum(cx)
    call parallel_sum(excit)
    call parallel_sum(ioniz)
#endif

    if(verbose) then
        call h5gcreate_f(id, "H_H", gid, error)

        dim1 = [1]
        dim2 = [n_max, nenergy]
        dim3 = [n_max, m_max, nenergy]

        call h5ltmake_dataset_int_f(gid, "nenergy", 0, dim1, [nenergy], error)
        call h5ltmake_dataset_int_f(gid, "n_max", 0, dim1, [n_max], error)
        call h5ltmake_dataset_int_f(gid, "m_max", 0, dim1, [m_max], error)
        call h5ltmake_dataset_double_f(gid, "dlogE", 0, dim1, [dlogE], error)
        call h5ltmake_dataset_double_f(gid, "emin", 0, dim1, [emin], error)
        call h5ltmake_dataset_double_f(gid, "emax", 0, dim1, [emax], error)

        dim1 = [nenergy]
        call h5ltmake_compressed_dataset_double_f(gid, "energy", 1, dim1, ebarr, error)
        call h5ltmake_compressed_dataset_double_f(gid, "cx", 3, dim3, cx, error)
        call h5ltmake_compressed_dataset_double_f(gid, "ionization", 2, dim2, ioniz, error)
        call h5ltmake_compressed_dataset_double_f(gid, "excitation", 3, dim3, excit, error)

        call h5ltset_attribute_string_f(id, "H_H", "description", &
             "Cross sections for Hydrogen-Hydrogen interactions", error)
        call h5ltset_attribute_string_f(gid, "nenergy", "description", &
             "Number of nucleon energy values", error)
        call h5ltset_attribute_string_f(gid, "n_max", "description", &
             "Number of initial energy levels", error)
        call h5ltset_attribute_string_f(gid, "m_max", "description", &
             "Number of final energy levels", error)
        call h5ltset_attribute_string_f(gid, "energy", "description", &
             "Nucleon energy values", error)
        call h5ltset_attribute_string_f(gid, "energy", "units", "keV/amu", error)
        call h5ltset_attribute_string_f(gid, "dlogE", "description", &
             "Energy spacing in log-10", error)
        call h5ltset_attribute_string_f(gid, "dlogE", "units", "log10(keV/amu)", error)
        call h5ltset_attribute_string_f(gid, "emin","description", &
             "Minimum energy", error)
        call h5ltset_attribute_string_f(gid, "emin", "units", "keV/amu", error)
        call h5ltset_attribute_string_f(gid, "emax","description", &
             "Maximum energy", error)
        call h5ltset_attribute_string_f(gid, "emax", "units", "keV/amu", error)

        call h5ltset_attribute_string_f(gid, "cx", "description", &
             "n/m resolved charge exchange cross sections: cx(n,m,energy)", error)
        call h5ltset_attribute_string_f(gid, "cx", "units", "cm^2", error)
        call h5ltset_attribute_string_f(gid, "cx", "reaction", &
             "H(+) + H(n) -> H(m) + H(+)", error)

        call h5ltset_attribute_string_f(gid, "excitation", "description", &
             "n/m resolved excitation cross sections: excitation(n,m,energy)", error)
        call h5ltset_attribute_string_f(gid, "excitation", "units", "cm^2", error)
        call h5ltset_attribute_string_f(gid, "excitation", "reaction", &
             "H(+) + H(n) -> H(+) + H(m), m > n", error)

        call h5ltset_attribute_string_f(gid, "ionization", "description", &
             "n resolved ionization cross sections: ionization(n,energy)", error)
        call h5ltset_attribute_string_f(gid, "ionization", "units", "cm^2", error)
        call h5ltset_attribute_string_f(gid, "ionization", "reaction", &
             "H(+) + H(n) -> H(+) + H(+) + e-", error)

        call h5gclose_f(gid, error)
    endif

    deallocate(ebarr, cx, excit, ioniz)

end subroutine write_bb_H_H

subroutine write_bb_H_e(id, namelist_file, n_max, m_max)
    !+ Write Hydrogen-Electron interaction cross sections to a HDF5 file
    integer(HID_T), intent(inout) :: id
        !+ HDF5 file or group object id
    character(len=*), intent(in)  :: namelist_file
        !+ Namelist file that contains settings
    integer, intent(in)           :: n_max
        !+ Number of initial atomic energy states to calculate
    integer, intent(in)           :: m_max
        !+ Number of final atomic energy states to calculate

    logical :: calculate
    real(Float64) :: emin
    real(Float64) :: emax
    integer :: nenergy

    real(Float64) :: eb
    real(Float64) :: dlogE
    real(Float64), dimension(:), allocatable :: ebarr

    real(Float64), dimension(:,:), allocatable :: ioniz
    real(Float64), dimension(:,:,:), allocatable :: excit

    integer(HID_T) :: gid
    integer(HSIZE_T), dimension(1) :: dim1
    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(3) :: dim3

    integer :: i, cnt, error
    logical :: exis

    NAMELIST /H_e_cross/ calculate, nenergy, emin, emax

    calculate = .True.; nenergy = 200; emin = 1.d-3 ; emax = 8.d2

    inquire(file=namelist_file,exist=exis)
    if(.not.exis) then
        write(*,'(a,a)') 'WRITE_BB_H_E: Input file does not exist: ',trim(namelist_file)
        write(*,'(a)') 'Continuing with default settings...'
    else
        open(13,file=namelist_file)
        read(13,NML=H_e_cross)
        close(13)
    endif

    if(.not.calculate) return

    allocate(ebarr(nenergy))
    allocate(ioniz(n_max,nenergy))
    allocate(excit(n_max,m_max,nenergy))

    ebarr = 0.d0
    ioniz = 0.d0
    excit = 0.d0

    if(verbose) then
        write(*,'(a)') "---- H-e cross sections settings ----"
        write(*,'(T2,"Emin = ",e9.2, " keV")') emin
        write(*,'(T2,"Emax = ",e9.2, " keV")') emax
        write(*,'(T2,"Nenergy = ", i4)') nenergy
        write(*,*) ''
    endif

    cnt = 0
    dlogE = (log10(emax) - log10(emin))/(nenergy - 1)
    !$OMP PARALLEL DO private(i, eb)
    do i=istart, nenergy, istep
        eb = 10.d0**(log10(emin) + (i-1)*dlogE)
        ebarr(i) = eb

        excit(:,:,i) = e_excit(eb, n_max, m_max)
        ioniz(:,i) = e_ioniz(eb, n_max)

        cnt = cnt + 1
        if(verbose) WRITE(*,'(f7.2,"%",a,$)') 100*cnt*istep/real(nenergy),char(13)
    enddo
    !$OMP END PARALLEL DO

#ifdef _MPI
    call parallel_sum(ebarr)
    call parallel_sum(excit)
    call parallel_sum(ioniz)
#endif

    if(verbose) then
        call h5gcreate_f(id, "H_e", gid, error)

        dim1 = [1]
        dim2 = [n_max, nenergy]
        dim3 = [n_max, m_max, nenergy]

        call h5ltmake_dataset_int_f(gid, "nenergy", 0, dim1, [nenergy], error)
        call h5ltmake_dataset_int_f(gid, "n_max", 0, dim1, [n_max], error)
        call h5ltmake_dataset_int_f(gid, "m_max", 0, dim1, [m_max], error)
        call h5ltmake_dataset_double_f(gid, "dlogE", 0, dim1, [dlogE], error)
        call h5ltmake_dataset_double_f(gid, "emin", 0, dim1, [emin], error)
        call h5ltmake_dataset_double_f(gid, "emax", 0, dim1, [emax], error)


        dim1 = [nenergy]
        call h5ltmake_compressed_dataset_double_f(gid, "energy", 1, dim1, ebarr, error)
        call h5ltmake_compressed_dataset_double_f(gid, "ionization", 2, dim2, ioniz, error)
        call h5ltmake_compressed_dataset_double_f(gid, "excitation", 3, dim3, excit, error)

        call h5ltset_attribute_string_f(id, "H_e", "description", &
             "Cross sections for Hydrogen-Electron interactions", error)
        call h5ltset_attribute_string_f(gid, "nenergy", "description", &
             "Number of nucleon energy values", error)
        call h5ltset_attribute_string_f(gid, "n_max", "description", &
             "Number of initial energy levels", error)
        call h5ltset_attribute_string_f(gid, "m_max", "description", &
             "Number of final energy levels", error)
        call h5ltset_attribute_string_f(gid, "energy", "description", &
             "Nucleon energy values", error)
        call h5ltset_attribute_string_f(gid, "energy", "units", "keV/amu", error)
        call h5ltset_attribute_string_f(gid, "dlogE", "description", &
             "Energy spacing in log-10", error)
        call h5ltset_attribute_string_f(gid, "dlogE", "units", "log10(keV/amu)", error)
        call h5ltset_attribute_string_f(gid, "emin","description", &
             "Minimum Energy", error)
        call h5ltset_attribute_string_f(gid, "emin", "units", "keV/amu", error)
        call h5ltset_attribute_string_f(gid, "emax","description", &
             "Maximum Energy", error)
        call h5ltset_attribute_string_f(gid, "emax", "units", "keV/amu", error)

        call h5ltset_attribute_string_f(gid, "excitation", "description", &
             "n/m resolved excitation cross sections: excitation(n,m,energy)", error)
        call h5ltset_attribute_string_f(gid, "excitation", "units", "cm^2", error)
        call h5ltset_attribute_string_f(gid, "excitation", "reaction", &
             "e- + H(n) -> e- + H(m), m > n", error)

        call h5ltset_attribute_string_f(gid, "ionization", "description", &
             "n resolved ionization cross sections: ionization(n,energy)", error)
        call h5ltset_attribute_string_f(gid, "ionization", "units", "cm^2", error)
        call h5ltset_attribute_string_f(gid, "ionization", "reaction", &
             "e- + H(n) -> e- + H(+) + e-", error)

        call h5gclose_f(gid, error)
    endif

    deallocate(ebarr, ioniz, excit)

end subroutine write_bb_H_e

subroutine write_bb_H_Aq(id, namelist_file, n_max, m_max)
    !+ Write Hydrogen-Impurity interaction cross sections to a HDF5 file
    integer(HID_T), intent(inout) :: id
        !+ HDF5 file or group object id
    character(len=*), intent(in)  :: namelist_file
        !+ Namelist file that contains settings
    integer, intent(in)           :: n_max
        !+ Number of initial atomic energy states to calculate
    integer, intent(in)           :: m_max
        !+ Number of final atomic energy states to calculate

    logical :: calculate
    integer :: q(10) = 0
    real(Float64) :: emin
    real(Float64) :: emax
    integer :: nenergy

    real(Float64) :: eb
    real(Float64) :: dlogE
    real(Float64), dimension(:), allocatable :: ebarr

    real(Float64), dimension(:,:), allocatable :: cx, ioniz
    real(Float64), dimension(:,:,:), allocatable :: excit

    integer(HID_T) :: gid
    integer(HSIZE_T), dimension(1) :: dim1
    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(3) :: dim3

    character(len=10) :: aname
    character(len=5) :: asym
    integer :: i, iq, cnt, error
    logical :: exis

    NAMELIST /H_Aq_cross/ calculate, q, nenergy, emin, emax

    nenergy = 200; emin = 1.d-3 ; emax = 8.d2
    calculate = .True.; q(1) = 6

    inquire(file=namelist_file,exist=exis)
    if(.not.exis) then
        write(*,'(a,a)') 'WRITE_BB_H_Aq: Input file does not exist: ',trim(namelist_file)
        write(*,'(a)') 'Continuing with default settings...'
    else
        open(13,file=namelist_file)
        read(13,NML=H_Aq_cross)
        close(13)
    endif

    if(.not.calculate) return

    q_loop: do iq=1, size(q)
        if(q(iq).eq.0) cycle q_loop

        allocate(ebarr(nenergy))
        allocate(ioniz(n_max,nenergy))
        allocate(cx(n_max,nenergy))
        allocate(excit(n_max,m_max,nenergy))

        ebarr = 0.d0
        ioniz = 0.d0
        cx = 0.d0
        excit = 0.d0

        select case (q(iq))
            case (5)
                aname = "Boron"
                asym = "H_B5"
            case (6)
                aname = "Carbon"
                asym = "H_C6"
            case DEFAULT
                write(aname,'("Impurity-",i1)') q(iq)
                write(asym,'("H_A",i1)') q(iq)
        end select

        if(verbose) then
            write(*,'(a)') "---- H-"//trim(adjustl(aname))//" cross sections settings ----"
            write(*,'(T2,"q = ", i2)') q(iq)
            write(*,'(T2,"Emin = ",e9.2, " keV")') emin
            write(*,'(T2,"Emax = ",e9.2, " keV")') emax
            write(*,'(T2,"Nenergy = ", i4)') nenergy
            write(*,*) ''
        endif

        cnt = 0
        dlogE = (log10(emax) - log10(emin))/(nenergy - 1)
        !$OMP PARALLEL DO private(i, eb)
        do i=istart, nenergy, istep
            eb = 10.d0**(log10(emin) + (i-1)*dlogE)
            ebarr(i) = eb

            cx(:,i) = Aq_cx(eb, q(iq), n_max)
            ioniz(:,i) = Aq_ioniz(eb, q(iq), n_max)
            excit(:,:,i) = Aq_excit(eb, q(iq), n_max, m_max)
            cnt = cnt + 1
            if(verbose) WRITE(*,'(f7.2,"%",a,$)') 100*cnt*istep/real(nenergy),char(13)
        enddo

#ifdef _MPI
        call parallel_sum(ebarr)
        call parallel_sum(cx)
        call parallel_sum(excit)
        call parallel_sum(ioniz)
#endif

        if(verbose) then
            call h5gcreate_f(id, trim(adjustl(asym)), gid, error)

            dim1 = [1]
            dim2 = [n_max, nenergy]
            dim3 = [n_max, m_max, nenergy]

            call h5ltmake_dataset_int_f(gid, "nenergy", 0, dim1, [nenergy], error)
            call h5ltmake_dataset_int_f(gid, "n_max", 0, dim1, [n_max], error)
            call h5ltmake_dataset_int_f(gid, "m_max", 0, dim1, [m_max], error)
            call h5ltmake_dataset_double_f(gid, "dlogE", 0, dim1, [dlogE], error)
            call h5ltmake_dataset_double_f(gid, "emin", 0, dim1, [emin], error)
            call h5ltmake_dataset_double_f(gid, "emax", 0, dim1, [emax], error)

            dim1 = [nenergy]
            call h5ltmake_compressed_dataset_double_f(gid, "energy", 1, dim1, ebarr, error)
            call h5ltmake_compressed_dataset_double_f(gid, "cx", 2, dim2, cx, error)
            call h5ltmake_compressed_dataset_double_f(gid, "ionization", 2, dim2, ioniz, error)
            call h5ltmake_compressed_dataset_double_f(gid, "excitation", 3, dim3, excit, error)

            call h5ltset_attribute_string_f(id, trim(adjustl(asym)), "description", &
                 "Cross sections for Hydrogen-"//trim(adjustl(aname))//" interactions", error)
            call h5ltset_attribute_string_f(gid, "nenergy", "description", &
                 "Number of nucleon energy values", error)
            call h5ltset_attribute_string_f(gid, "n_max", "description", &
                 "Number of initial energy levels", error)
            call h5ltset_attribute_string_f(gid, "m_max", "description", &
                 "Number of final energy levels", error)
            call h5ltset_attribute_string_f(gid, "energy", "description", &
                 "Nucleon energy values", error)
            call h5ltset_attribute_string_f(gid, "energy", "units", "keV/amu", error)
            call h5ltset_attribute_string_f(gid, "dlogE", "description", &
                 "Energy spacing in log-10", error)
            call h5ltset_attribute_string_f(gid, "dlogE", "units", "log10(keV/amu)", error)
            call h5ltset_attribute_string_f(gid, "emin","description", &
                 "Minimum energy", error)
            call h5ltset_attribute_string_f(gid, "emin", "units", "keV/amu", error)
            call h5ltset_attribute_string_f(gid, "emax","description", &
                 "Maximum energy", error)
            call h5ltset_attribute_string_f(gid, "emax", "units", "keV/amu", error)

            call h5ltset_attribute_string_f(gid, "cx", "description", &
                 "n resolved charge exchange / electron capture cross sections: cx(n,energy)", error)
            call h5ltset_attribute_string_f(gid, "cx", "units", "cm^2", error)
            call h5ltset_attribute_string_f(gid, "cx", "reaction", &
                 "A(q+) + H(n) -> A((q-1)+) + H(+)", error)

            call h5ltset_attribute_string_f(gid, "excitation", "description", &
                 "n/m resolved excitation cross sections: excitation(n,m,energy)", error)
            call h5ltset_attribute_string_f(gid, "excitation", "units", "cm^2", error)
            call h5ltset_attribute_string_f(gid, "excitation", "reaction", &
                 "A(q+) + H(n) -> A(q+) + H(m), m > n", error)

            call h5ltset_attribute_string_f(gid, "ionization", "description", &
                 "n resolved ionization cross sections: ionization(n,energy)", error)
            call h5ltset_attribute_string_f(gid, "ionization", "units", "cm^2", error)
            call h5ltset_attribute_string_f(gid, "ionization", "reaction", &
                 "A(q+) + H(n) -> A(q+) + H(+) + e-", error)

            call h5gclose_f(gid, error)
        endif

        deallocate(ebarr, ioniz, cx, excit)
    enddo q_loop

end subroutine write_bb_H_Aq

subroutine write_bb_D_D(id, namelist_file)
    !+ Write Deuterium-Deuterium interaction cross sections to a HDF5 file
    integer(HID_T), intent(inout) :: id
        !+ HDF5 file or group object id
    character(len=*), intent(in)  :: namelist_file
        !+ Namelist file that contains settings

    logical :: calculate
    integer :: nbranch = 2
    real(Float64) :: emin
    real(Float64) :: emax
    integer :: nenergy

    real(Float64) :: eb
    real(Float64) :: dlogE
    real(Float64), dimension(:), allocatable :: ebarr

    real(Float64), dimension(:,:), allocatable :: fusion

    integer(HID_T) :: gid
    integer(HSIZE_T), dimension(1) :: dim1
    integer(HSIZE_T), dimension(2) :: dim2

    integer :: i, cnt, error
    logical :: exis

    NAMELIST /D_D_cross/ calculate, nenergy, emin, emax

    calculate = .True.; nenergy = 200; emin = 1.d-3 ; emax = 8.d2

    inquire(file=namelist_file,exist=exis)
    if(.not.exis) then
        write(*,'(a,a)') 'WRITE_BB_D_D: Input file does not exist: ',trim(namelist_file)
        write(*,'(a)') 'Continuing with default settings...'
    else
        open(13,file=namelist_file)
        read(13,NML=D_D_cross)
        close(13)
    endif

    if(.not.calculate) return

    allocate(ebarr(nenergy))
    allocate(fusion(nenergy,nbranch))

    ebarr = 0.d0
    fusion = 0.d0

    if(verbose) then
        write(*,'(a)') "---- D-D cross sections settings ----"
        write(*,'(T2,"Emin = ",e9.2, " keV")') emin
        write(*,'(T2,"Emax = ",e9.2, " keV")') emax
        write(*,'(T2,"Nenergy = ", i4)') nenergy
        write(*,*) ''
    endif

    cnt = 0
    dlogE = (log10(emax) - log10(emin))/(nenergy - 1)
    !$OMP PARALLEL DO private(i, eb)
    do i=istart, nenergy, istep
        eb = 10.d0**(log10(emin) + (i-1)*dlogE)
        ebarr(i) = eb

        fusion(i,1) = d_d_fusion_t(eb)
        fusion(i,2) = d_d_fusion_he(eb)
        cnt = cnt + 1
        if(verbose) WRITE(*,'(f7.2,"%",a,$)') 100*cnt*istep/real(nenergy),char(13)
    enddo
    !$OMP END PARALLEL DO

#ifdef _MPI
    call parallel_sum(ebarr)
    call parallel_sum(fusion)
#endif

    if(verbose) then
        call h5gcreate_f(id, "D_D", gid, error)

        dim1 = [1]
        call h5ltmake_dataset_int_f(gid, "nbranch", 0, dim1, [nbranch], error)
        call h5ltmake_dataset_int_f(gid, "nenergy", 0, dim1, [nenergy], error)
        call h5ltmake_dataset_double_f(gid, "dlogE", 0, dim1, [dlogE], error)
        call h5ltmake_dataset_double_f(gid, "emin", 0, dim1, [emin], error)
        call h5ltmake_dataset_double_f(gid, "emax", 0, dim1, [emax], error)

        dim1 = [nenergy]
        dim2 = [nenergy, nbranch]
        call h5ltmake_compressed_dataset_double_f(gid, "energy", 1, dim1, ebarr, error)
        call h5ltmake_compressed_dataset_double_f(gid, "fusion", 2, dim2, fusion, error)

        call h5ltset_attribute_string_f(id, "D_D", "description", &
             "Cross sections for Deuterium-Deuterium interactions", error)
        call h5ltset_attribute_string_f(gid, "nbranch", "description", &
             "Number of reaction branches", error)
        call h5ltset_attribute_string_f(gid, "nenergy", "description", &
             "Number of energy values", error)
        call h5ltset_attribute_string_f(gid, "energy", "description", &
             "Deuterium energy values", error)
        call h5ltset_attribute_string_f(gid, "energy", "units", "keV", error)
        call h5ltset_attribute_string_f(gid, "dlogE", "description", &
             "Energy spacing in log-10", error)
        call h5ltset_attribute_string_f(gid, "dlogE", "units", "log10(keV)", error)
        call h5ltset_attribute_string_f(gid, "emin","description", &
             "Minimum energy", error)
        call h5ltset_attribute_string_f(gid, "emin", "units", "keV", error)
        call h5ltset_attribute_string_f(gid, "emax","description", &
             "Maximum energy", error)
        call h5ltset_attribute_string_f(gid, "emax", "units", "keV", error)

        call h5ltset_attribute_string_f(gid, "fusion", "description", &
             "Cross sections for the Tritium[1] and He3[2] branches of D-D nuclear reactions: fusion(energy, branch)", error)
        call h5ltset_attribute_string_f(gid, "fusion", "units", "cm^2", error)
        call h5ltset_attribute_string_f(gid, "fusion", "reaction", &
             "D + D -> [1] T(1.01 MeV) + p(3.02 MeV) (50%); [2] He3(0.82 MeV) + n(2.45 MeV) (50%)", error)

        call h5gclose_f(gid, error)
    endif
    deallocate(ebarr, fusion)

end subroutine write_bb_D_D

subroutine write_bb_D_T(id, namelist_file)
    !+ Write Deuterium-Tritium interaction cross sections to a HDF5 file
    integer(HID_T), intent(inout) :: id
        !+ HDF5 file or group object id
    character(len=*), intent(in)  :: namelist_file
        !+ Namelist file that contains settings

    logical :: calculate
    integer :: nbranch = 1
    real(Float64) :: emin
    real(Float64) :: emax
    integer :: nenergy

    real(Float64) :: eb
    real(Float64) :: dlogE
    real(Float64), dimension(:), allocatable :: ebarr

    real(Float64), dimension(:,:), allocatable :: fusion

    integer(HID_T) :: gid
    integer(HSIZE_T), dimension(1) :: dim1
    integer(HSIZE_T), dimension(2) :: dim2

    integer :: i, cnt, error
    logical :: exis

    NAMELIST /D_T_cross/ calculate, nenergy, emin, emax

    calculate = .True.; nenergy = 200; emin = 1.d-3 ; emax = 8.d2

    inquire(file=namelist_file,exist=exis)
    if(.not.exis) then
        write(*,'(a,a)') 'WRITE_BB_D_T: Input file does not exist: ',trim(namelist_file)
        write(*,'(a)') 'Continuing with default settings...'
    else
        open(13,file=namelist_file)
        read(13,NML=D_T_cross)
        close(13)
    endif

    if(.not.calculate) return

    allocate(ebarr(nenergy))
    allocate(fusion(nenergy,nbranch))

    ebarr = 0.d0
    fusion = 0.d0

    if(verbose) then
        write(*,'(a)') "---- D-T cross sections settings ----"
        write(*,'(T2,"Emin = ",e9.2, " keV")') emin
        write(*,'(T2,"Emax = ",e9.2, " keV")') emax
        write(*,'(T2,"Nenergy = ", i4)') nenergy
        write(*,*) ''
    endif

    cnt = 0
    dlogE = (log10(emax) - log10(emin))/(nenergy - 1)
    !$OMP PARALLEL DO private(i, eb)
    do i=istart, nenergy, istep
        eb = 10.d0**(log10(emin) + (i-1)*dlogE)
        ebarr(i) = eb

        fusion(i,1) = d_t_fusion(eb)
        cnt = cnt + 1
        if(verbose) WRITE(*,'(f7.2,"%",a,$)') 100*cnt*istep/real(nenergy),char(13)
    enddo
    !$OMP END PARALLEL DO

#ifdef _MPI
    call parallel_sum(ebarr)
    call parallel_sum(fusion)
#endif

    if(verbose) then
        call h5gcreate_f(id, "D_T", gid, error)

        dim1 = [1]

        call h5ltmake_dataset_int_f(gid, "nbranch", 0, dim1, [nbranch], error)
        call h5ltmake_dataset_int_f(gid, "nenergy", 0, dim1, [nenergy], error)
        call h5ltmake_dataset_double_f(gid, "dlogE", 0, dim1, [dlogE], error)
        call h5ltmake_dataset_double_f(gid, "emin", 0, dim1, [emin], error)
        call h5ltmake_dataset_double_f(gid, "emax", 0, dim1, [emax], error)

        dim1 = [nenergy]
        dim2 = [nenergy, nbranch]
        call h5ltmake_compressed_dataset_double_f(gid, "energy", 1, dim1, ebarr, error)
        call h5ltmake_compressed_dataset_double_f(gid, "fusion", 2, dim2, fusion, error)

        call h5ltset_attribute_string_f(id, "D_T", "description", &
             "Cross sections for Deuterium-Tritium interactions", error)
        call h5ltset_attribute_string_f(gid, "nbranch", "description", &
             "Number of reaction branches", error)
        call h5ltset_attribute_string_f(gid, "nenergy", "description", &
             "Number of energy values", error)
        call h5ltset_attribute_string_f(gid, "energy", "description", &
             "Deuterium energy values", error)
        call h5ltset_attribute_string_f(gid, "energy", "units", "keV", error)
        call h5ltset_attribute_string_f(gid, "dlogE", "description", &
             "Energy spacing in log-10", error)
        call h5ltset_attribute_string_f(gid, "dlogE", "units", "log10(keV)", error)
        call h5ltset_attribute_string_f(gid, "emin","description", &
             "Minimum energy", error)
        call h5ltset_attribute_string_f(gid, "emin", "units", "keV", error)
        call h5ltset_attribute_string_f(gid, "emax","description", &
             "Maximum energy", error)
        call h5ltset_attribute_string_f(gid, "emax", "units", "keV", error)

        call h5ltset_attribute_string_f(gid, "fusion", "description", &
             "Total cross sections for D-T nuclear reactions: fusion(deuterium energy, branch)", error)
        call h5ltset_attribute_string_f(gid, "fusion", "units", "cm^2", error)
        call h5ltset_attribute_string_f(gid, "fusion", "reaction", &
             "D + T -> He4(3.5 MeV) + n(14.1 MeV)", error)

        call h5gclose_f(gid, error)
    endif

    deallocate(ebarr, fusion)

end subroutine write_bb_D_T

subroutine write_bt_H_H(id, namelist_file, n_max, m_max)
    !+ Write Hydrogen-Hydrogen reaction rates to a HDF5 file
    integer(HID_T), intent(inout) :: id
        !+ HDF5 file or group object id
    character(len=*), intent(in)  :: namelist_file
        !+ Namelist file that contains settings
    integer, intent(in)           :: n_max
        !+ Number of initial atomic energy states to calculate
    integer, intent(in)           :: m_max
        !+ Number of final atomic energy states to calculate

    logical :: calculate
    real(Float64) :: emin
    real(Float64) :: emax
    integer :: nenergy

    real(Float64) :: tmin
    real(Float64) :: tmax
    integer :: ntemp

    real(Float64) :: eb
    real(Float64) :: dlogE
    real(Float64), dimension(:), allocatable :: ebarr

    real(Float64) :: ti
    real(Float64) :: dlogT
    real(Float64), dimension(:), allocatable :: tarr

    real(Float64), dimension(:,:,:,:), allocatable :: ioniz
    real(Float64), dimension(:,:,:,:,:), allocatable :: excit, cx

    integer(HID_T) :: gid
    integer(HSIZE_T), dimension(1) :: dim1
    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(3) :: dim3
    integer(HSIZE_T), dimension(4) :: dim4
    integer(HSIZE_T), dimension(5) :: dim5

    integer :: ie, it, ia, n, m, error, cnt
    real(Float64) :: rate
    integer, parameter :: n_bt_amu = 4
    real(Float64), dimension(2,n_bt_amu) :: a
    logical :: exis

    NAMELIST /H_H_rates/ calculate, nenergy, emin, emax, ntemp, tmin, tmax

    nenergy = 100; emin = 1.d-3 ; emax = 4.d2
    ntemp = 100; tmin = 1.d-3 ; tmax = 2.d1
    calculate = .True.

    inquire(file=namelist_file,exist=exis)
    if(.not.exis) then
        write(*,'(a,a)') 'WRITE_BT_H_H: Input file does not exist: ',trim(namelist_file)
        write(*,'(a)') 'Continuing with default settings...'
    else
        open(13,file=namelist_file)
        read(13,NML=H_H_rates)
        close(13)
    endif

    if(.not.calculate) return

    allocate(ebarr(nenergy))
    allocate(tarr(ntemp))
    allocate(ioniz(n_max,nenergy,ntemp,n_bt_amu))
    allocate(cx(n_max,m_max,nenergy,ntemp,n_bt_amu))
    allocate(excit(n_max,m_max,nenergy,ntemp,n_bt_amu))

    ebarr = 0.d0
    tarr = 0.d0
    ioniz = 0.d0
    cx = 0.d0
    excit = 0.d0
    a(:,1) = [H1_amu, H1_amu]
    a(:,2) = [H1_amu, H2_amu]
    a(:,3) = [H2_amu, H1_amu]
    a(:,4) = [H2_amu, H2_amu]

    dlogE = (log10(emax) - log10(emin))/(nenergy - 1)
    do ie=1, nenergy
        ebarr(ie) = 10.d0**(log10(emin) + (ie-1)*dlogE)
    enddo

    dlogT = (log10(tmax) - log10(tmin))/(ntemp - 1)
    do it=1, ntemp
        tarr(it) = 10.d0**(log10(tmin) + (it-1)*dlogT)
    enddo

    if(verbose) then
        write(*,'(a)') "---- H-H reaction rates settings ----"
        write(*,'(T2,"Emin = ",e9.2, " keV")') emin
        write(*,'(T2,"Emax = ",e9.2, " keV")') emax
        write(*,'(T2,"Nenergy = ", i4)') nenergy
        write(*,'(T2,"Tmin = ",e9.2, " keV")') tmin
        write(*,'(T2,"Tmax = ",e9.2, " keV")') tmax
        write(*,'(T2,"Ntemp = ", i4)') ntemp
        write(*,*) ''
    endif

    cnt = 0
    !$OMP PARALLEL DO private(ie, it, ia, n, m, eb, ti, rate)
    do ie=istart, nenergy, istep
        eb = ebarr(ie)
        do it=1, ntemp
            ti = tarr(it)
            do ia=1, n_bt_amu
                do n=1, n_max
                    do m=1, m_max
                        call bt_maxwellian(p_cx_n_m, ti, eb, &
                                           a(2,ia), a(1,ia), n, m, rate)
                        cx(n,m,ie,it,ia) = rate
                        if(m.gt.n) then
                            call bt_maxwellian(p_excit_n_m, ti, eb, &
                                               a(2,ia), a(1,ia), n, m, rate)
                            excit(n,m,ie,it,ia) = rate

                            call bt_maxwellian(p_excit_n_m, ti, eb, &
                                               a(2,ia), a(1,ia), n, m, &
                                               rate, deexcit=.True.)
                            excit(m,n,ie,it,ia) = rate
                        endif
                    enddo
                    call bt_maxwellian(p_ioniz_n, ti, eb, &
                                       a(2,ia), a(1,ia), n, rate)
                    ioniz(n,ie,it,ia) = rate
                enddo
            enddo
            cnt = cnt + 1
            if(verbose) WRITE(*,'(f7.2,"%",a,$)') 100*istep*cnt/real(nenergy*ntemp),char(13)
        enddo
    enddo
    !$OMP END PARALLEL DO

#ifdef _MPI
    call parallel_sum(ebarr)
    call parallel_sum(tarr)
    call parallel_sum(cx)
    call parallel_sum(excit)
    call parallel_sum(ioniz)
#endif

    if(verbose) then
        call h5gcreate_f(id, "H_H", gid, error)

        dim1 = [1]
        dim4 = [n_max, nenergy, ntemp, n_bt_amu]
        dim5 = [n_max, m_max, nenergy, ntemp, n_bt_amu]
        call h5ltmake_dataset_int_f(gid, "nenergy", 0, dim1, [nenergy], error)
        call h5ltmake_dataset_int_f(gid, "ntemp", 0, dim1, [ntemp], error)
        call h5ltmake_dataset_int_f(gid, "n_bt_amu", 0, dim1, [n_bt_amu], error)
        call h5ltmake_dataset_int_f(gid, "n_max", 0, dim1, [n_max], error)
        call h5ltmake_dataset_int_f(gid, "m_max", 0, dim1, [m_max], error)
        call h5ltmake_dataset_double_f(gid, "dlogE", 0, dim1, [dlogE], error)
        call h5ltmake_dataset_double_f(gid, "emin", 0, dim1, [emin], error)
        call h5ltmake_dataset_double_f(gid, "emax", 0, dim1, [emax], error)
        call h5ltmake_dataset_double_f(gid, "dlogT", 0, dim1, [dlogT], error)
        call h5ltmake_dataset_double_f(gid, "tmin", 0, dim1, [tmin], error)
        call h5ltmake_dataset_double_f(gid, "tmax", 0, dim1, [tmax], error)

        dim2 = [2,n_bt_amu]
        call h5ltmake_compressed_dataset_double_f(gid, "bt_amu", 2, dim2, a, error)
        dim1 = [nenergy]
        call h5ltmake_compressed_dataset_double_f(gid, "energy", 1, dim1, ebarr, error)
        dim1 = [ntemp]
        call h5ltmake_compressed_dataset_double_f(gid, "temperature", 1, dim1, tarr, error)
        call h5ltmake_compressed_dataset_double_f(gid, "cx", 5, dim5, cx, error)
        call h5ltmake_compressed_dataset_double_f(gid, "ionization", 4, dim4, ioniz, error)
        call h5ltmake_compressed_dataset_double_f(gid, "excitation", 5, dim5, excit, error)

        call h5ltset_attribute_string_f(id, "H_H", "description", &
             "Beam-Target reaction rates for Hydrogen(beam)-Hydrogen(target) interactions", error)
        call h5ltset_attribute_string_f(gid, "nenergy", "description", &
             "Number of energy values", error)
        call h5ltset_attribute_string_f(gid, "ntemp", "description", &
             "Number of target temperature values", error)
        call h5ltset_attribute_string_f(gid, "n_bt_amu", "description", &
             "Number of beam-target amu combinations", error)
        call h5ltset_attribute_string_f(gid, "n_max", "description", &
             "Number of initial energy levels", error)
        call h5ltset_attribute_string_f(gid, "m_max", "description", &
             "Number of final energy levels", error)

        call h5ltset_attribute_string_f(gid, "bt_amu", "description", &
             "Combinations of beam-target amu's e.g. b_amu, t_amu = bt_amu[:,i]", error)
        call h5ltset_attribute_string_f(gid, "energy", "description", &
             "Energy values", error)
        call h5ltset_attribute_string_f(gid, "energy", "units", "keV", error)
        call h5ltset_attribute_string_f(gid, "dlogE", "description", &
             "Energy spacing in log-10", error)
        call h5ltset_attribute_string_f(gid, "dlogE", "units", "log10(keV)", error)
        call h5ltset_attribute_string_f(gid, "emin","description", &
             "Minimum energy", error)
        call h5ltset_attribute_string_f(gid, "emin", "units", "keV", error)
        call h5ltset_attribute_string_f(gid, "emax","description", &
             "Maximum energy", error)
        call h5ltset_attribute_string_f(gid, "emax", "units", "keV", error)

        call h5ltset_attribute_string_f(gid, "temperature", "description", &
             "Target temperature values", error)
        call h5ltset_attribute_string_f(gid, "temperature", "units", "keV", error)
        call h5ltset_attribute_string_f(gid, "dlogT", "description", &
             "Temperature spacing in log-10", error)
        call h5ltset_attribute_string_f(gid, "dlogT", "units", "log10(keV)", error)
        call h5ltset_attribute_string_f(gid, "tmin","description", &
             "Minimum temperature", error)
        call h5ltset_attribute_string_f(gid, "tmin", "units", "keV", error)
        call h5ltset_attribute_string_f(gid, "tmax","description", &
             "Maximum temperature", error)
        call h5ltset_attribute_string_f(gid, "tmax", "units", "keV", error)

        call h5ltset_attribute_string_f(gid, "cx", "description", &
             "n/m resolved charge exchange reaction rates: cx(n,m,energy,temp,bt_amu)", error)
        call h5ltset_attribute_string_f(gid, "cx", "units", "cm^3/s", error)
        call h5ltset_attribute_string_f(gid, "cx", "reaction", &
             "H(+) + H(n) -> H(m) + H(+)", error)

        call h5ltset_attribute_string_f(gid, "excitation", "description", &
             "n/m resolved (de-)excitation reaction rates: excitation(n,m,energy,temp,bt_amu)", error)
        call h5ltset_attribute_string_f(gid, "excitation", "units", "cm^3/s", error)
        call h5ltset_attribute_string_f(gid, "excitation", "reaction", &
             "H(+) + H(n) -> H(+) + H(m); m > n excitation, m < n de-excitation", error)

        call h5ltset_attribute_string_f(gid, "ionization", "description", &
             "n resolved ionization reaction rates: ionization(n,energy,temp,bt_amu)", error)
        call h5ltset_attribute_string_f(gid, "ionization", "units", "cm^3/s", error)
        call h5ltset_attribute_string_f(gid, "ionization", "reaction", &
             "H(+) + H(n) -> H(+) + H(+) + e-", error)

        call h5gclose_f(gid, error)
    endif

    deallocate(ebarr, tarr, cx, excit, ioniz)

end subroutine write_bt_H_H

subroutine write_bt_H_e(id, namelist_file, n_max, m_max)
    !+ Write Hydrogen-Electron reaction rates to a HDF5 file
    integer(HID_T), intent(inout) :: id
        !+ HDF5 file or group object id
    character(len=*), intent(in)  :: namelist_file
        !+ Namelist file that contains settings
    integer, intent(in)           :: n_max
        !+ Number of initial atomic energy states to calculate
    integer, intent(in)           :: m_max
        !+ Number of final atomic energy states to calculate

    logical :: calculate

    real(Float64) :: emin
    real(Float64) :: emax
    integer :: nenergy

    real(Float64) :: tmin
    real(Float64) :: tmax
    integer :: ntemp

    real(Float64) :: eb
    real(Float64) :: dlogE
    real(Float64), dimension(:), allocatable :: ebarr

    real(Float64) :: ti
    real(Float64) :: dlogT
    real(Float64), dimension(:), allocatable :: tarr

    real(Float64), dimension(:,:,:,:), allocatable :: ioniz
    real(Float64), dimension(:,:,:,:,:), allocatable :: excit

    integer(HID_T) :: gid
    integer(HSIZE_T), dimension(1) :: dim1
    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(3) :: dim3
    integer(HSIZE_T), dimension(4) :: dim4
    integer(HSIZE_T), dimension(5) :: dim5

    integer :: ie, it, ia, n, m, error, cnt
    real(Float64) :: rate
    integer, parameter :: n_bt_amu = 2
    real(Float64), dimension(2,n_bt_amu) :: a
    logical :: exis

    NAMELIST /H_e_rates/ calculate, nenergy, emin, emax, ntemp, tmin, tmax

    nenergy = 100; emin = 1.d-3 ; emax = 4.d2
    ntemp = 100; tmin = 1.d-3 ; tmax = 2.d1
    calculate = .True.

    inquire(file=namelist_file,exist=exis)
    if(.not.exis) then
        write(*,'(a,a)') 'WRITE_BT_H_E: Input file does not exist: ',trim(namelist_file)
        write(*,'(a)') 'Continuing with default settings...'
    else
        open(13,file=namelist_file)
        read(13,NML=H_e_rates)
        close(13)
    endif

    if(.not.calculate) return

    allocate(ebarr(nenergy))
    allocate(tarr(ntemp))
    allocate(ioniz(n_max,nenergy,ntemp,n_bt_amu))
    allocate(excit(n_max,m_max,nenergy,ntemp,n_bt_amu))

    ebarr = 0.d0
    ioniz = 0.d0
    excit = 0.d0
    a(:,1) = [H1_amu, e_amu]
    a(:,2) = [H2_amu, e_amu]

    dlogE = (log10(emax) - log10(emin))/(nenergy - 1)
    do ie=1, nenergy
        ebarr(ie) = 10.d0**(log10(emin) + (ie-1)*dlogE)
    enddo

    dlogT = (log10(tmax) - log10(tmin))/(ntemp - 1)
    do it=1, ntemp
        tarr(it) = 10.d0**(log10(tmin) + (it-1)*dlogT)
    enddo

    if(verbose) then
        write(*,'(a)') "---- H-e reaction rates settings ----"
        write(*,'(T2,"Emin = ",e9.2, " keV")') emin
        write(*,'(T2,"Emax = ",e9.2, " keV")') emax
        write(*,'(T2,"Nenergy = ", i4)') nenergy
        write(*,'(T2,"Tmin = ",e9.2, " keV")') tmin
        write(*,'(T2,"Tmax = ",e9.2, " keV")') tmax
        write(*,'(T2,"Ntemp = ", i4)') ntemp
        write(*,*) ''
    endif

    cnt = 0
    !$OMP PARALLEL DO private(ie, it, ia, n, m, eb, ti, rate)
    do ie=istart, nenergy, istep
        eb = ebarr(ie)
        do it=1, ntemp
            ti = tarr(it)
            do ia=1, n_bt_amu
                do n=1, n_max
                    do m=1, m_max
                        if(m.gt.n) then
                            call bt_maxwellian(e_excit_n_m, ti, eb, &
                                               a(2,ia), a(1,ia), n, m, rate)
                            excit(n,m,ie,it,ia) = rate

                            call bt_maxwellian(e_excit_n_m, ti, eb, &
                                               a(2,ia), a(1,ia), n, m, &
                                               rate, deexcit=.True.)
                            excit(m,n,ie,it,ia) = rate
                        endif
                    enddo
                    call bt_maxwellian(e_ioniz_n, ti, eb, &
                                       a(2,ia), a(1,ia), n, rate)
                    ioniz(n,ie,it,ia) = rate
                enddo
            enddo
            cnt = cnt + 1
            if(verbose) WRITE(*,'(f7.2,"%",a,$)') 100*istep*cnt/real(nenergy*ntemp),char(13)
        enddo
    enddo
    !$OMP END PARALLEL DO

#ifdef _MPI
    call parallel_sum(ebarr)
    call parallel_sum(tarr)
    call parallel_sum(excit)
    call parallel_sum(ioniz)
#endif

    if(verbose) then
        call h5gcreate_f(id, "H_e", gid, error)

        dim1 = [1]
        dim4 = [n_max, nenergy, ntemp, n_bt_amu]
        dim5 = [n_max, m_max, nenergy, ntemp, n_bt_amu]
        call h5ltmake_dataset_int_f(gid, "nenergy", 0, dim1, [nenergy], error)
        call h5ltmake_dataset_int_f(gid, "ntemp", 0, dim1, [ntemp], error)
        call h5ltmake_dataset_int_f(gid, "n_bt_amu", 0, dim1, [n_bt_amu], error)
        call h5ltmake_dataset_int_f(gid, "n_max", 0, dim1, [n_max], error)
        call h5ltmake_dataset_int_f(gid, "m_max", 0, dim1, [m_max], error)
        call h5ltmake_dataset_double_f(gid, "dlogE", 0, dim1, [dlogE], error)
        call h5ltmake_dataset_double_f(gid, "emin", 0, dim1, [emin], error)
        call h5ltmake_dataset_double_f(gid, "emax", 0, dim1, [emax], error)
        call h5ltmake_dataset_double_f(gid, "dlogT", 0, dim1, [dlogT], error)
        call h5ltmake_dataset_double_f(gid, "tmin", 0, dim1, [tmin], error)
        call h5ltmake_dataset_double_f(gid, "tmax", 0, dim1, [tmax], error)

        dim2 = [2,n_bt_amu]
        call h5ltmake_compressed_dataset_double_f(gid, "bt_amu", 2, dim2, a, error)
        dim1 = [nenergy]
        call h5ltmake_compressed_dataset_double_f(gid, "energy", 1, dim1, ebarr, error)
        dim1 = [ntemp]
        call h5ltmake_compressed_dataset_double_f(gid, "temperature", 1, dim1, tarr, error)
        call h5ltmake_compressed_dataset_double_f(gid, "ionization", 4, dim4, ioniz, error)
        call h5ltmake_compressed_dataset_double_f(gid, "excitation", 5, dim5, excit, error)

        call h5ltset_attribute_string_f(id, "H_e", "description", &
             "Beam-Target reaction rates for Hydrogen(beam)-Electron(target) interactions", error)
        call h5ltset_attribute_string_f(gid, "nenergy", "description", &
             "Number of energy values", error)
        call h5ltset_attribute_string_f(gid, "ntemp", "description", &
             "Number of target temperature values", error)
        call h5ltset_attribute_string_f(gid, "n_bt_amu", "description", &
             "Number of beam-target amu combinations", error)
        call h5ltset_attribute_string_f(gid, "n_max", "description", &
             "Number of initial energy levels", error)
        call h5ltset_attribute_string_f(gid, "m_max", "description", &
             "Number of final energy levels", error)

        call h5ltset_attribute_string_f(gid, "bt_amu", "description", &
             "Combinations of beam-target amu's e.g. b_amu, t_amu = bt_amu[:,i]", error)
        call h5ltset_attribute_string_f(gid, "energy", "description", &
             "Energy values", error)
        call h5ltset_attribute_string_f(gid, "energy", "units", "keV", error)
        call h5ltset_attribute_string_f(gid, "dlogE", "description", &
             "Energy spacing in log-10", error)
        call h5ltset_attribute_string_f(gid, "dlogE", "units", "log10(keV)", error)
        call h5ltset_attribute_string_f(gid, "emin","description", &
             "Minimum energy", error)
        call h5ltset_attribute_string_f(gid, "emin", "units", "keV", error)
        call h5ltset_attribute_string_f(gid, "emax","description", &
             "Maximum energy", error)
        call h5ltset_attribute_string_f(gid, "emax", "units", "keV", error)

        call h5ltset_attribute_string_f(gid, "temperature", "description", &
             "Target temperature values", error)
        call h5ltset_attribute_string_f(gid, "temperature", "units", "keV", error)
        call h5ltset_attribute_string_f(gid, "dlogT", "description", &
             "Temperature spacing in log-10", error)
        call h5ltset_attribute_string_f(gid, "dlogT", "units", "log10(keV)", error)
        call h5ltset_attribute_string_f(gid, "tmin","description", &
             "Minimum temperature", error)
        call h5ltset_attribute_string_f(gid, "tmin", "units", "keV", error)
        call h5ltset_attribute_string_f(gid, "tmax","description", &
             "Maximum temperature", error)
        call h5ltset_attribute_string_f(gid, "tmax", "units", "keV", error)

        call h5ltset_attribute_string_f(gid, "excitation", "description", &
             "n/m resolved (de-)excitation reaction rates: excitation(n,m,energy,temp,bt_amu)", error)
        call h5ltset_attribute_string_f(gid, "excitation", "units", "cm^3/s", error)
        call h5ltset_attribute_string_f(gid, "excitation", "reaction", &
             "e- + H(n) -> e- + H(m); m > n excitation, m < n de-excitation", error)

        call h5ltset_attribute_string_f(gid, "ionization", "description", &
             "n resolved ionization reaction rates: ionization(n,energy,temp,bt_amu)", error)
        call h5ltset_attribute_string_f(gid, "ionization", "units", "cm^3/s", error)
        call h5ltset_attribute_string_f(gid, "ionization", "reaction", &
             "e- + H(n) -> e- + H(+) + e-", error)

        call h5gclose_f(gid, error)
    endif

    deallocate(ebarr, tarr, excit, ioniz)

end subroutine write_bt_H_e

subroutine write_bt_H_Aq(id, namelist_file, n_max, m_max)
    !+ Write Hydrogen-Impurity reaction rates to a HDF5 file
    integer(HID_T), intent(inout) :: id
        !+ HDF5 file or group object id
    character(len=*), intent(in)  :: namelist_file
        !+ Namelist file that contains settings
    integer, intent(in)           :: n_max
        !+ Number of initial atomic energy states to calculate
    integer, intent(in)           :: m_max
        !+ Number of final atomic energy states to calculate

    logical :: calculate
    integer :: q(10) = 0
    real(Float64) :: mass

    real(Float64) :: emin
    real(Float64) :: emax
    integer :: nenergy

    real(Float64) :: tmin
    real(Float64) :: tmax
    integer :: ntemp

    real(Float64) :: eb
    real(Float64) :: dlogE
    real(Float64), dimension(:), allocatable :: ebarr

    real(Float64) :: ti
    real(Float64) :: dlogT
    real(Float64), dimension(:), allocatable :: tarr

    real(Float64), dimension(:,:,:,:), allocatable :: ioniz, cx
    real(Float64), dimension(:,:,:,:,:), allocatable :: excit

    integer(HID_T) :: gid
    integer(HSIZE_T), dimension(1) :: dim1
    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(3) :: dim3
    integer(HSIZE_T), dimension(4) :: dim4
    integer(HSIZE_T), dimension(5) :: dim5

    integer :: iq, ie, it, ia, n, m, error, cnt
    real(Float64) :: rate
    integer, parameter :: n_bt_amu = 2
    real(Float64), dimension(2,n_bt_amu) :: a

    character(len=10) :: aname
    character(len=5) :: asym
    logical :: exis

    NAMELIST /H_Aq_rates/ calculate, q, mass, nenergy, emin, emax, ntemp, tmin, tmax

    calculate = .True. ; q(1) = 6 ; mass = C_amu
    nenergy = 100; emin = 1.d-3 ; emax = 4.d2
    ntemp = 100; tmin = 1.d-3 ; tmax = 2.d1
    inquire(file=namelist_file,exist=exis)
    if(.not.exis) then
        write(*,'(a,a)') 'WRITE_BT_H_Aq: Input file does not exist: ',trim(namelist_file)
        write(*,'(a)') 'Continuing with default settings...'
    else
        open(13,file=namelist_file)
        read(13,NML=H_Aq_rates)
        close(13)
    endif

    if(.not.calculate) return

    q_loop: do iq=1, size(q)
        if(q(iq).eq.0) cycle q_loop

        allocate(ebarr(nenergy))
        allocate(tarr(ntemp))
        allocate(ioniz(n_max,nenergy,ntemp,n_bt_amu))
        allocate(cx(n_max,nenergy,ntemp,n_bt_amu))
        allocate(excit(n_max,m_max,nenergy,ntemp,n_bt_amu))


        select case (q(iq))
          case (5)
              aname = "Boron"
              asym = "H_B5"
          case (6)
              aname = "Carbon"
              asym = "H_C6"
          case DEFAULT
              write(aname,'("Impurity-",i1)') q(iq)
              write(asym,'("H_A",i1)') q(iq)
        end select

        ebarr = 0.d0
        ioniz = 0.d0
        cx = 0.d0
        excit = 0.d0
        a(:,1) = [H1_amu, mass]
        a(:,2) = [H2_amu, mass]

        dlogE = (log10(emax) - log10(emin))/(nenergy - 1)
        do ie=1, nenergy
            ebarr(ie) = 10.d0**(log10(emin) + (ie-1)*dlogE)
        enddo

        dlogT = (log10(tmax) - log10(tmin))/(ntemp - 1)
        do it=1, ntemp
            tarr(it) = 10.d0**(log10(tmin) + (it-1)*dlogT)
        enddo

        if(verbose) then
            write(*,'(a)') "---- H-"//trim(adjustl(aname))//" reaction rates settings ----"
            write(*,'(T2,"q = ", i2)') q(iq)
            write(*,'(T2,"mass = ",f7.2, " amu")') mass
            write(*,'(T2,"Emin = ",e9.2, " keV")') emin
            write(*,'(T2,"Emax = ",e9.2, " keV")') emax
            write(*,'(T2,"Nenergy = ", i4)') nenergy
            write(*,'(T2,"Tmin = ",e9.2, " keV")') tmin
            write(*,'(T2,"Tmax = ",e9.2, " keV")') tmax
            write(*,'(T2,"Ntemp = ", i4)') ntemp
            write(*,*) ''
        endif

        cnt = 0
        !$OMP PARALLEL DO private(ie, it, ia, n, m, eb, ti, rate)
        do ie=istart, nenergy, istep
            eb = ebarr(ie)
            do it=1, ntemp
                ti = tarr(it)
                do ia=1, n_bt_amu
                    do n=1, n_max
                        do m=1, m_max
                            if(m.gt.n) then
                                call bt_maxwellian(Aq_excit_n_m, q(iq), ti, eb, &
                                                   a(2,ia), a(1,ia), n, m, rate)
                                excit(n,m,ie,it,ia) = rate

                                call bt_maxwellian(Aq_excit_n_m, q(iq), ti, eb, &
                                                   a(2,ia), a(1,ia), n, m, &
                                                   rate, deexcit=.True.)
                                excit(m,n,ie,it,ia) = rate
                            endif
                        enddo
                        call bt_maxwellian(Aq_cx_n, q(iq), ti, eb, &
                                           a(2,ia), a(1,ia), n, rate)
                        cx(n,ie,it,ia) = rate

                        call bt_maxwellian(Aq_ioniz_n, q(iq), ti, eb, &
                                           a(2,ia), a(1,ia), n, rate)
                        ioniz(n,ie,it,ia) = rate
                    enddo
                enddo
                cnt = cnt + 1
                if(verbose) WRITE(*,'(f7.2,"%",a,$)') 100*istep*cnt/real(nenergy*ntemp),char(13)
            enddo
        enddo
        !$OMP END PARALLEL DO

#ifdef _MPI
        call parallel_sum(ebarr)
        call parallel_sum(tarr)
        call parallel_sum(cx)
        call parallel_sum(excit)
        call parallel_sum(ioniz)
#endif

        if(verbose) then
            call h5gcreate_f(id, trim(adjustl(asym)), gid, error)

            dim1 = [1]
            dim4 = [n_max, nenergy, ntemp, n_bt_amu]
            dim5 = [n_max, m_max, nenergy, ntemp, n_bt_amu]
            call h5ltmake_dataset_int_f(gid, "nenergy", 0, dim1, [nenergy], error)
            call h5ltmake_dataset_int_f(gid, "ntemp", 0, dim1, [ntemp], error)
            call h5ltmake_dataset_int_f(gid, "n_bt_amu", 0, dim1, [n_bt_amu], error)
            call h5ltmake_dataset_int_f(gid, "n_max", 0, dim1, [n_max], error)
            call h5ltmake_dataset_int_f(gid, "m_max", 0, dim1, [m_max], error)
            call h5ltmake_dataset_double_f(gid, "dlogE", 0, dim1, [dlogE], error)
            call h5ltmake_dataset_double_f(gid, "emin", 0, dim1, [emin], error)
            call h5ltmake_dataset_double_f(gid, "emax", 0, dim1, [emax], error)
            call h5ltmake_dataset_double_f(gid, "dlogT", 0, dim1, [dlogT], error)
            call h5ltmake_dataset_double_f(gid, "tmin", 0, dim1, [tmin], error)
            call h5ltmake_dataset_double_f(gid, "tmax", 0, dim1, [tmax], error)

            dim2 = [2,n_bt_amu]
            call h5ltmake_compressed_dataset_double_f(gid, "bt_amu", 2, dim2, a, error)
            dim1 = [nenergy]
            call h5ltmake_compressed_dataset_double_f(gid, "energy", 1, dim1, ebarr, error)
            dim1 = [ntemp]
            call h5ltmake_compressed_dataset_double_f(gid, "temperature", 1, dim1, tarr, error)
            call h5ltmake_compressed_dataset_double_f(gid, "cx", 4, dim4, cx, error)
            call h5ltmake_compressed_dataset_double_f(gid, "ionization", 4, dim4, ioniz, error)
            call h5ltmake_compressed_dataset_double_f(gid, "excitation", 5, dim5, excit, error)

            call h5ltset_attribute_string_f(id, trim(adjustl(asym)), "description", &
                 "Beam-Target reaction rates for Hydrogen(beam)-"//trim(adjustl(aname))// &
                 "(target) interactions", error)
            call h5ltset_attribute_string_f(gid, "nenergy", "description", &
                 "Number of energy values", error)
            call h5ltset_attribute_string_f(gid, "ntemp", "description", &
                 "Number of target temperature values", error)
            call h5ltset_attribute_string_f(gid, "n_bt_amu", "description", &
                 "Number of beam-target amu combinations", error)
            call h5ltset_attribute_string_f(gid, "n_max", "description", &
                 "Number of initial energy levels", error)
            call h5ltset_attribute_string_f(gid, "m_max", "description", &
                 "Number of final energy levels", error)

            call h5ltset_attribute_string_f(gid, "bt_amu", "description", &
                 "Combinations of beam-target amu's e.g. b_amu, t_amu = bt_amu[:,i]", error)
            call h5ltset_attribute_string_f(gid, "energy", "description", &
                 "Energy values", error)
            call h5ltset_attribute_string_f(gid, "energy", "units", "keV", error)
            call h5ltset_attribute_string_f(gid, "dlogE", "description", &
                 "Energy spacing in log-10", error)
            call h5ltset_attribute_string_f(gid, "dlogE", "units", "log10(keV)", error)
            call h5ltset_attribute_string_f(gid, "emin","description", &
                 "Minimum energy", error)
            call h5ltset_attribute_string_f(gid, "emin", "units", "keV", error)
            call h5ltset_attribute_string_f(gid, "emax","description", &
                 "Maximum energy", error)
            call h5ltset_attribute_string_f(gid, "emax", "units", "keV", error)

            call h5ltset_attribute_string_f(gid, "temperature", "description", &
                 "Target temperature values", error)
            call h5ltset_attribute_string_f(gid, "temperature", "units", "keV", error)
            call h5ltset_attribute_string_f(gid, "dlogT", "description", &
                 "Temperature spacing in log-10", error)
            call h5ltset_attribute_string_f(gid, "dlogT", "units", "log10(keV)", error)
            call h5ltset_attribute_string_f(gid, "tmin","description", &
                 "Minimum temperature", error)
            call h5ltset_attribute_string_f(gid, "tmin", "units", "keV", error)
            call h5ltset_attribute_string_f(gid, "tmax","description", &
                 "Maximum temperature", error)
            call h5ltset_attribute_string_f(gid, "tmax", "units", "keV", error)

            call h5ltset_attribute_string_f(gid, "cx", "description", &
                 "n-resolved charge exchange reaction rates: cx(n,energy,temp,bt_amu)", error)
            call h5ltset_attribute_string_f(gid, "cx", "units", "cm^3/s", error)
            call h5ltset_attribute_string_f(gid, "cx", "reaction", &
                 "A(q+) + H(n) -> A((q-1)+) + H(+)", error)

            call h5ltset_attribute_string_f(gid, "excitation", "description", &
                 "n/m resolved (de-)excitation reaction rates: excitation(n,m,energy,temp,bt_amu)", error)
            call h5ltset_attribute_string_f(gid, "excitation", "units", "cm^3/s", error)
            call h5ltset_attribute_string_f(gid, "excitation", "reaction", &
                 "A(q+) + H(n) -> A(q+) + H(m); m > n excitation, m < n de-excitation", error)

            call h5ltset_attribute_string_f(gid, "ionization", "description", &
                 "n resolved ionization reaction rates: ionization(n,energy,temp,bt_amu)", error)
            call h5ltset_attribute_string_f(gid, "ionization", "units", "cm^3/s", error)
            call h5ltset_attribute_string_f(gid, "ionization", "reaction", &
                 "A(q+) + H(n) -> A(q+) + H(+) + e-", error)

            call h5gclose_f(gid, error)
        endif

        deallocate(ebarr, tarr, excit, ioniz, cx)
    enddo q_loop

end subroutine write_bt_H_Aq

subroutine write_bt_D_D(id, namelist_file)
    !+ Write Deuterium-Deuterium reaction rates to a HDF5 file
    integer(HID_T), intent(inout) :: id
        !+ HDF5 file or group object id
    character(len=*), intent(in)  :: namelist_file
        !+ Namelist file that contains settings

    logical :: calculate
    integer :: nbranch = 2
    real(Float64), dimension(2) :: bt_amu = [H2_amu, H2_amu]

    real(Float64) :: emin
    real(Float64) :: emax
    integer :: nenergy

    real(Float64) :: tmin
    real(Float64) :: tmax
    integer :: ntemp

    real(Float64) :: eb
    real(Float64) :: dlogE
    real(Float64), dimension(:), allocatable :: ebarr

    real(Float64) :: ti
    real(Float64) :: dlogT
    real(Float64), dimension(:), allocatable :: tarr

    real(Float64), dimension(:,:,:), allocatable :: fusion

    integer(HID_T) :: gid
    integer(HSIZE_T), dimension(1) :: dim1
    integer(HSIZE_T), dimension(3) :: dim3

    integer :: ie, it, error, cnt
    real(Float64) :: rate_a, rate_b
    logical :: exis

    NAMELIST /D_D_rates/ calculate, nenergy, emin, emax, ntemp, tmin, tmax

    nenergy = 100; emin = 1.d-3 ; emax = 4.d2
    ntemp = 100; tmin = 1.d-3 ; tmax = 2.d1
    calculate = .True.

    inquire(file=namelist_file,exist=exis)
    if(.not.exis) then
        write(*,'(a,a)') 'WRITE_BT_D_D: Input file does not exist: ',trim(namelist_file)
        write(*,'(a)') 'Continuing with default settings...'
    else
        open(13,file=namelist_file)
        read(13,NML=D_D_rates)
        close(13)
    endif

    if(.not.calculate) return

    allocate(ebarr(nenergy))
    allocate(tarr(ntemp))
    allocate(fusion(nenergy,ntemp,nbranch))

    ebarr = 0.d0
    tarr = 0.d0
    fusion = 0.d0

    dlogE = (log10(emax) - log10(emin))/(nenergy - 1)
    do ie=1, nenergy
        ebarr(ie) = 10.d0**(log10(emin) + (ie-1)*dlogE)
    enddo

    dlogT = (log10(tmax) - log10(tmin))/(ntemp - 1)
    do it=1, ntemp
        tarr(it) = 10.d0**(log10(tmin) + (it-1)*dlogT)
    enddo

    if(verbose) then
        write(*,'(a)') "---- D-D reaction rates settings ----"
        write(*,'(T2,"Emin = ",e9.2, " keV")') emin
        write(*,'(T2,"Emax = ",e9.2, " keV")') emax
        write(*,'(T2,"Nenergy = ", i4)') nenergy
        write(*,'(T2,"Tmin = ",e9.2, " keV")') tmin
        write(*,'(T2,"Tmax = ",e9.2, " keV")') tmax
        write(*,'(T2,"Ntemp = ", i4)') ntemp
        write(*,*) ''
    endif

    cnt = 0
    !$OMP PARALLEL DO private(ie, it, eb, ti, rate_a, rate_b)
    do ie=istart, nenergy, istep
        eb = ebarr(ie)
        do it=1, ntemp
            ti = tarr(it)
            call bt_maxwellian(d_d_fusion_t, ti, eb, bt_amu(1), bt_amu(2), rate_a)
            call bt_maxwellian(d_d_fusion_he, ti, eb, bt_amu(2), bt_amu(2), rate_b)
            fusion(ie,it,1) = rate_a
            fusion(ie,it,2) = rate_b
            cnt = cnt + 1
            if(verbose) WRITE(*,'(f7.2,"%",a,$)') 100*istep*cnt/real(nenergy*ntemp),char(13)
        enddo
    enddo
    !$OMP END PARALLEL DO

#ifdef _MPI
    call parallel_sum(ebarr)
    call parallel_sum(tarr)
    call parallel_sum(fusion)
#endif

    if(verbose) then
        call h5gcreate_f(id, "D_D", gid, error)

        dim1 = [1]
        call h5ltmake_dataset_int_f(gid, "nbranch", 0, dim1, [nbranch], error)
        call h5ltmake_dataset_int_f(gid, "nenergy", 0, dim1, [nenergy], error)
        call h5ltmake_dataset_int_f(gid, "ntemp", 0, dim1, [ntemp], error)
        call h5ltmake_dataset_double_f(gid, "dlogE", 0, dim1, [dlogE], error)
        call h5ltmake_dataset_double_f(gid, "emin", 0, dim1, [emin], error)
        call h5ltmake_dataset_double_f(gid, "emax", 0, dim1, [emax], error)
        call h5ltmake_dataset_double_f(gid, "dlogT", 0, dim1, [dlogT], error)
        call h5ltmake_dataset_double_f(gid, "tmin", 0, dim1, [tmin], error)
        call h5ltmake_dataset_double_f(gid, "tmax", 0, dim1, [tmax], error)

        dim1 = [2]
        call h5ltmake_compressed_dataset_double_f(gid, "bt_amu", 1, dim1, bt_amu, error)
        dim1 = [nenergy]
        call h5ltmake_compressed_dataset_double_f(gid, "energy", 1, dim1, ebarr, error)
        dim1 = [ntemp]
        call h5ltmake_compressed_dataset_double_f(gid, "temperature", 1, dim1, tarr, error)
        dim3 = [nenergy, ntemp, nbranch]
        call h5ltmake_compressed_dataset_double_f(gid, "fusion", 3, dim3, fusion, error)

        call h5ltset_attribute_string_f(id, "D_D", "description", &
             "Beam-Target reaction rates for Deuterium(beam)-Deuterium(target) interactions", error)
        call h5ltset_attribute_string_f(gid, "nbranch", "description", &
             "Number of reaction branches", error)
        call h5ltset_attribute_string_f(gid, "nenergy", "description", &
             "Number of energy values", error)
        call h5ltset_attribute_string_f(gid, "ntemp", "description", &
             "Number of target temperature values", error)
        call h5ltset_attribute_string_f(gid, "energy", "description", &
             "Energy values", error)
        call h5ltset_attribute_string_f(gid, "energy", "units", "keV", error)
        call h5ltset_attribute_string_f(gid, "dlogE", "description", &
             "Energy spacing in log-10", error)
        call h5ltset_attribute_string_f(gid, "dlogE", "units", "log10(keV)", error)
        call h5ltset_attribute_string_f(gid, "emin","description", &
             "Minimum energy", error)
        call h5ltset_attribute_string_f(gid, "emin", "units", "keV", error)
        call h5ltset_attribute_string_f(gid, "emax","description", &
             "Maximum energy", error)
        call h5ltset_attribute_string_f(gid, "emax", "units", "keV", error)

        call h5ltset_attribute_string_f(gid, "temperature", "description", &
             "Target temperature values", error)
        call h5ltset_attribute_string_f(gid, "temperature", "units", "keV", error)
        call h5ltset_attribute_string_f(gid, "dlogT", "description", &
             "Temperature spacing in log-10", error)
        call h5ltset_attribute_string_f(gid, "dlogT", "units", "log10(keV)", error)
        call h5ltset_attribute_string_f(gid, "tmin","description", &
             "Minimum temperature", error)
        call h5ltset_attribute_string_f(gid, "tmin", "units", "keV", error)
        call h5ltset_attribute_string_f(gid, "tmax","description", &
             "Maximum temperature", error)
        call h5ltset_attribute_string_f(gid, "tmax", "units", "keV", error)

        call h5ltset_attribute_string_f(gid, "bt_amu", "description", &
             "Isotope mass of the beam and target species respectively", error)
        call h5ltset_attribute_string_f(gid, "bt_amu", "units", "amu", error)

        call h5ltset_attribute_string_f(gid, "fusion", "description", &
             "Beam-Target reaction rates for T/He3 branches of D-D nuclear reactions: fusion(energy, temp, branch)", error)
        call h5ltset_attribute_string_f(gid, "fusion", "units", "cm^3/s", error)
        call h5ltset_attribute_string_f(gid, "fusion", "reaction", &
             "D + D -> [1] T(1.01 MeV) + p(3.02 MeV) (50%); [2] He3(0.82 MeV) + n(2.45 MeV) (50%)", error)

        call h5gclose_f(gid, error)
    endif

    deallocate(ebarr, tarr, fusion)

end subroutine write_bt_D_D

subroutine write_bt_D_T(id, namelist_file)
    !+ Write Deuterium-Tritium reaction rates to a HDF5 file
    integer(HID_T), intent(inout) :: id
        !+ HDF5 file or group object id
    character(len=*), intent(in)  :: namelist_file
        !+ Namelist file that contains settings

    integer :: nbranch = 1
    real(Float64), dimension(2) :: bt_amu = [H2_amu, H3_amu]

    logical :: calculate
    real(Float64) :: emin
    real(Float64) :: emax
    integer :: nenergy

    real(Float64) :: tmin
    real(Float64) :: tmax
    integer :: ntemp

    real(Float64) :: eb
    real(Float64) :: dlogE
    real(Float64), dimension(:), allocatable :: ebarr

    real(Float64) :: ti
    real(Float64) :: dlogT
    real(Float64), dimension(:), allocatable :: tarr

    real(Float64), dimension(:,:,:), allocatable :: fusion

    integer(HID_T) :: gid
    integer(HSIZE_T), dimension(1) :: dim1
    integer(HSIZE_T), dimension(3) :: dim3

    integer :: ie, it, error, cnt
    real(Float64) :: rate
    logical :: exis

    NAMELIST /D_T_rates/ calculate, nenergy, emin, emax, ntemp, tmin, tmax

    nenergy = 100; emin = 1.d-3 ; emax = 4.d2
    ntemp = 100; tmin = 1.d-3 ; tmax = 2.d1
    calculate = .True.

    inquire(file=namelist_file,exist=exis)
    if(.not.exis) then
        write(*,'(a,a)') 'WRITE_BT_D_T: Input file does not exist: ',trim(namelist_file)
        write(*,'(a)') 'Continuing with default settings...'
    else
        open(13,file=namelist_file)
        read(13,NML=D_T_rates)
        close(13)
    endif

    if(.not.calculate) return

    allocate(ebarr(nenergy))
    allocate(tarr(ntemp))
    allocate(fusion(nenergy,ntemp,nbranch))

    ebarr = 0.d0
    tarr = 0.d0
    fusion = 0.d0

    dlogE = (log10(emax) - log10(emin))/(nenergy - 1)
    do ie=1, nenergy
        ebarr(ie) = 10.d0**(log10(emin) + (ie-1)*dlogE)
    enddo

    dlogT = (log10(tmax) - log10(tmin))/(ntemp - 1)
    do it=1, ntemp
        tarr(it) = 10.d0**(log10(tmin) + (it-1)*dlogT)
    enddo

    if(verbose) then
        write(*,'(a)') "---- D-T reaction rates settings ----"
        write(*,'(T2,"Emin = ",e9.2, " keV")') emin
        write(*,'(T2,"Emax = ",e9.2, " keV")') emax
        write(*,'(T2,"Nenergy = ", i4)') nenergy
        write(*,'(T2,"Tmin = ",e9.2, " keV")') tmin
        write(*,'(T2,"Tmax = ",e9.2, " keV")') tmax
        write(*,'(T2,"Ntemp = ", i4)') ntemp
        write(*,*) ''
    endif

    cnt = 0
    !$OMP PARALLEL DO private(ie, it, eb, ti, rate)
    do ie=istart, nenergy, istep
        eb = ebarr(ie)
        do it=1, ntemp
            ti = tarr(it)
            call bt_maxwellian(d_t_fusion, ti, eb, bt_amu(1), bt_amu(2), rate)
            fusion(ie,it,1) = rate
            cnt = cnt + 1
            if(verbose) WRITE(*,'(f7.2,"%",a,$)') 100*istep*cnt/real(nenergy*ntemp),char(13)
        enddo
    enddo
    !$OMP END PARALLEL DO

#ifdef _MPI
    call parallel_sum(ebarr)
    call parallel_sum(tarr)
    call parallel_sum(fusion)
#endif

    if(verbose) then
        call h5gcreate_f(id, "D_T", gid, error)

        dim1 = [1]
        call h5ltmake_dataset_int_f(gid, "nbranch", 0, dim1, [nbranch], error)
        call h5ltmake_dataset_int_f(gid, "nenergy", 0, dim1, [nenergy], error)
        call h5ltmake_dataset_int_f(gid, "ntemp", 0, dim1, [ntemp], error)
        call h5ltmake_dataset_double_f(gid, "dlogE", 0, dim1, [dlogE], error)
        call h5ltmake_dataset_double_f(gid, "emin", 0, dim1, [emin], error)
        call h5ltmake_dataset_double_f(gid, "emax", 0, dim1, [emax], error)
        call h5ltmake_dataset_double_f(gid, "dlogT", 0, dim1, [dlogT], error)
        call h5ltmake_dataset_double_f(gid, "tmin", 0, dim1, [tmin], error)
        call h5ltmake_dataset_double_f(gid, "tmax", 0, dim1, [tmax], error)

        dim1 = [2]
        call h5ltmake_compressed_dataset_double_f(gid, "bt_amu", 1, dim1, bt_amu, error)
        dim1 = [nenergy]
        call h5ltmake_compressed_dataset_double_f(gid, "energy", 1, dim1, ebarr, error)
        dim1 = [ntemp]
        call h5ltmake_compressed_dataset_double_f(gid, "temperature", 1, dim1, tarr, error)
        dim3 = [nenergy, ntemp, nbranch]
        call h5ltmake_compressed_dataset_double_f(gid, "fusion", 3, dim3, fusion, error)

        call h5ltset_attribute_string_f(id, "D_T", "description", &
             "Beam-Target reaction rates for Deuterium(beam)-Tritium(target) interactions", error)
        call h5ltset_attribute_string_f(gid, "nbranch", "description", &
             "Number of reaction branches", error)
        call h5ltset_attribute_string_f(gid, "nenergy", "description", &
             "Number of energy values", error)
        call h5ltset_attribute_string_f(gid, "ntemp", "description", &
             "Number of target temperature values", error)
        call h5ltset_attribute_string_f(gid, "energy", "description", &
             "Energy values", error)
        call h5ltset_attribute_string_f(gid, "energy", "units", "keV", error)
        call h5ltset_attribute_string_f(gid, "dlogE", "description", &
             "Energy spacing in log-10", error)
        call h5ltset_attribute_string_f(gid, "dlogE", "units", "log10(keV)", error)
        call h5ltset_attribute_string_f(gid, "emin","description", &
             "Minimum energy", error)
        call h5ltset_attribute_string_f(gid, "emin", "units", "keV", error)
        call h5ltset_attribute_string_f(gid, "emax","description", &
             "Maximum energy", error)
        call h5ltset_attribute_string_f(gid, "emax", "units", "keV", error)

        call h5ltset_attribute_string_f(gid, "temperature", "description", &
             "Target temperature values", error)
        call h5ltset_attribute_string_f(gid, "temperature", "units", "keV", error)
        call h5ltset_attribute_string_f(gid, "dlogT", "description", &
             "Temperature spacing in log-10", error)
        call h5ltset_attribute_string_f(gid, "dlogT", "units", "log10(keV)", error)
        call h5ltset_attribute_string_f(gid, "tmin","description", &
             "Minimum temperature", error)
        call h5ltset_attribute_string_f(gid, "tmin", "units", "keV", error)
        call h5ltset_attribute_string_f(gid, "tmax","description", &
             "Maximum temperature", error)
        call h5ltset_attribute_string_f(gid, "tmax", "units", "keV", error)

        call h5ltset_attribute_string_f(gid, "bt_amu", "description", &
             "Isotope mass of the beam and target species respectively", error)
        call h5ltset_attribute_string_f(gid, "bt_amu", "units", "amu", error)

        call h5ltset_attribute_string_f(gid, "fusion", "description", &
             "Beam-Target reaction rates for D-T nuclear reactions: fusion(energy, temperature, branch)", error)
        call h5ltset_attribute_string_f(gid, "fusion", "units", "cm^3/s", error)
        call h5ltset_attribute_string_f(gid, "fusion", "reaction", &
             "D + T -> He4(3.5 MeV) + n(14.1 MeV)", error)

        call h5gclose_f(gid, error)
    endif

    deallocate(ebarr, tarr, fusion)

end subroutine write_bt_D_T

subroutine print_default_namelist
    !+ Prints out the default settings as a namelist

    write(*,'(a)') "!Default Atomic Table Settings"
    write(*,'(a)') "&general_settings"
    write(*,'(a)') "n_max = 12,    !Number of initial atomic energy levels"
    write(*,'(a)') "m_max = 12,    !Number of final atomic energy levels"
    write(*,'(a)') "tables_file = './atomic_tables.h5'"
    write(*,'(a)') "/"
    write(*,'(a)') "!Hydrogen-Hydrogen Cross Sections"
    write(*,'(a)') "&H_H_cross"
    write(*,'(a)') "calculate = T, !Calculate Table"
    write(*,'(a)') "nenergy = 200, !Number of energy values"
    write(*,'(a)') "emin = 1.0E-3, !Minimum energy [keV]"
    write(*,'(a)') "emax = 8.0E2   !Maximum energy [keV]"
    write(*,'(a)') "/"
    write(*,'(a)') "!Hydrogen-Electron Cross Sections"
    write(*,'(a)') "&H_e_cross"
    write(*,'(a)') "calculate = T, !Calculate Table"
    write(*,'(a)') "nenergy = 200, !Number of energy values"
    write(*,'(a)') "emin = 1.0E-3, !Minimum energy [keV]"
    write(*,'(a)') "emax = 8.0E2   !Maximum energy [keV]"
    write(*,'(a)') "/"
    write(*,'(a)') "!Hydrogen-Impurity Cross Sections. Up to 10 impurity charges"
    write(*,'(a)') "&H_Aq_cross"
    write(*,'(a)') "calculate = T, !Calculate Table"
    write(*,'(a)') "q(1) = 5,      !Impurity charge: Boron: 5, Carbon: 6, ..."
    write(*,'(a)') "q(2) = 6,      !Impurity charge: Boron: 5, Carbon: 6, ..."
    write(*,'(a)') "nenergy = 200, !Number of energy values"
    write(*,'(a)') "emin = 1.0E-3, !Minimum energy [keV]"
    write(*,'(a)') "emax = 8.0E2   !Maximum energy [keV]"
    write(*,'(a)') "/"
    write(*,'(a)') "!Deuterium-Deuterium Nuclear Cross Sections"
    write(*,'(a)') "&D_D_cross"
    write(*,'(a)') "calculate = T, !Calculate Table"
    write(*,'(a)') "nenergy = 200, !Number of energy values"
    write(*,'(a)') "emin = 1.0E-3, !Minimum energy [keV]"
    write(*,'(a)') "emax = 8.0E2   !Maximum energy [keV]"
    write(*,'(a)') "/"
    write(*,'(a)') "!Hydrogen-Hydrogen Reaction Rates"
    write(*,'(a)') "&H_H_rates"
    write(*,'(a)') "calculate = T, !Calculate Table"
    write(*,'(a)') "nenergy = 100, !Number of energy values"
    write(*,'(a)') "emin = 1.0E-3, !Minimum energy [keV]"
    write(*,'(a)') "emax = 4.0E2,  !Maximum energy [keV]"
    write(*,'(a)') "ntemp = 100,   !Number of temperature values"
    write(*,'(a)') "tmin = 1.0E-3, !Minimum ion temperature [keV]"
    write(*,'(a)') "tmax = 2.0E1   !Maximum ion temperature [keV]"
    write(*,'(a)') "/"
    write(*,'(a)') "!Hydrogen-Electron Reaction Rates"
    write(*,'(a)') "&H_e_rates"
    write(*,'(a)') "calculate = T, !Calculate Table"
    write(*,'(a)') "nenergy = 100, !Number of energy values"
    write(*,'(a)') "emin = 1.0E-3, !Minimum energy [keV]"
    write(*,'(a)') "emax = 4.0E2,  !Maximum energy [keV]"
    write(*,'(a)') "ntemp = 100,   !Number of temperature values"
    write(*,'(a)') "tmin = 1.0E-3, !Minimum electron temperature [keV]"
    write(*,'(a)') "tmax = 2.0E1   !Maximum electron temperature [keV]"
    write(*,'(a)') "/"
    write(*,'(a)') "!Hydrogen-Impurity Reaction Rates. Up to 10 impurity charges"
    write(*,'(a)') "&H_Aq_rates"
    write(*,'(a)') "calculate = T, !Calculate Table"
    write(*,'(a)') "q(1) = 5,      !Impurity charge: Boron: 5, Carbon: 6, ..."
    write(*,'(a)') "q(2) = 6,      !Impurity charge: Boron: 5, Carbon: 6, ..."
    write(*,'(a)') "mass = 12.011, !Impurity mass [amu]"
    write(*,'(a)') "nenergy = 100, !Number of energy values"
    write(*,'(a)') "emin = 1.0E-3, !Minimum energy [keV]"
    write(*,'(a)') "emax = 4.0E2,  !Maximum energy [keV]"
    write(*,'(a)') "ntemp = 100,   !Number of temperature values"
    write(*,'(a)') "tmin = 1.0E-3, !Minimum ion temperature [keV]"
    write(*,'(a)') "tmax = 2.0E1   !Maximum ion temperature [keV]"
    write(*,'(a)') "/"
    write(*,'(a)') "!Deuterium-Deuterium Nuclear Reaction Rates"
    write(*,'(a)') "&D_D_rates"
    write(*,'(a)') "calculate = T, !Calculate Table"
    write(*,'(a)') "nenergy = 100, !Number of energy values"
    write(*,'(a)') "emin = 1.0E-3, !Minimum energy [keV]"
    write(*,'(a)') "emax = 4.0E2,  !Maximum energy [keV]"
    write(*,'(a)') "ntemp = 100,   !Number of temperature values"
    write(*,'(a)') "tmin = 1.0E-3, !Minimum deuterium temperature [keV]"
    write(*,'(a)') "tmax = 2.0E1   !Maximum deuterium temperature [keV]"
    write(*,'(a)') "/"

end subroutine print_default_namelist

end module atomic_tables

program generate_tables
    !+ Tabulates cross sections and reaction rates and writes them to a HDF5 file
    use atomic_tables
    use H5LT
    use HDF5
    use hdf5_utils
    use utilities

#ifdef _OMP
  use omp_lib
#endif

    character(len=200) :: namelist_file
    character(len=3) :: arg

    character(len=200) :: tables_file = ''
    integer :: n_max, m_max

    integer, dimension(8) :: time_start
    integer :: hour, minu, sec
    integer :: argc, max_threads, nthreads
    integer(HID_T) :: fid, gid
    integer :: error
    logical :: exis

    NAMELIST /general_settings/ n_max, m_max, tables_file

    argc = command_argument_count()
    if(argc.eq.0) then
        call print_default_namelist()
        stop
    endif

    if(argc.ge.1) then
        call get_command_argument(1, namelist_file)
    endif

    inquire(file=namelist_file,exist=exis)
    if(.not.exis) then
        write(*,'(a,a)') 'Input file does not exist: ',trim(namelist_file)
        stop
    else
        open(13,file=namelist_file)
        read(13,NML=general_settings)
        close(13)
        n_max = min(n_max,15)
        m_max = min(m_max,15)
    endif

#ifdef _OMP
    max_threads = OMP_get_num_procs()
    if(argc.ge.2) then
        call get_command_argument(2,arg)
        read(arg,'(i3)') nthreads
    else
        nthreads = max_threads
    endif
    max_threads = min(nthreads,max_threads)
    write(*,'(a)') "---- OpenMP settings ----"
    write(*,'(T2,"Number of threads: ",i2)') max_threads
    write(*,*) ''
    call OMP_set_num_threads(max_threads)
#endif

#ifdef _MPI
    call init_mpi()
    istart = my_rank()+1
    istep = num_ranks()
    if(my_rank().ne.0) verbose = .False.
    if(verbose) then
        write(*,'(a)') "---- MPI settings ----"
        write(*,'(T2,"Number of processes: ",i2)') num_ranks()
        write(*,*) ''
    endif
#endif

    if(verbose) then
        write(*,'(a)') "---- General settings ----"
        write(*,'(T2,"n_max = ",i2)') n_max
        write(*,'(T2,"m_max = ",i2)') m_max
        write(*,'(T2,"Tables File: ",a)') trim(tables_file)
        write(*,*) ''
    endif
    !! Check if compression is possible
    call check_compression_availability()

    call date_and_time (values=time_start)

    if(verbose) then
        !! Open HDF5 Interface
        call h5open_f(error)

        !! Create tables file. Overwrites if already exists
        call h5fcreate_f(tables_file, H5F_ACC_TRUNC_F, fid, error)

        !! Create group for cross sections
        call h5gcreate_f(fid, "cross", gid, error)
    endif

    !! Calculate cross sections
    if(verbose) then
        write(*,*) 'Cross Sections:   ',time(time_start)
        write(*,*) ''
    endif
    call write_bb_H_H(gid, namelist_file, n_max, m_max)
    call write_bb_H_e(gid, namelist_file, n_max, m_max)
    call write_bb_H_Aq(gid, namelist_file, n_max, m_max)
    call write_bb_D_D(gid, namelist_file)

    if(verbose) then
        !! Close cross section group
        call h5gclose_f(gid, error)

        !! Create group for reaction rates
        call h5gcreate_f(fid, "rates", gid, error)
    endif

    !! Calculate reaction rates
    if(verbose) then
        write(*,*) 'Reaction Rates:   ',time(time_start)
        write(*,*) ''
    endif
    call write_bt_H_H(gid, namelist_file, n_max, m_max)
    call write_bt_H_e(gid, namelist_file, n_max, m_max)
    call write_bt_H_Aq(gid, namelist_file, n_max, m_max)
    call write_einstein(gid, n_max, m_max)
    call write_bt_D_D(gid, namelist_file)

    if(verbose) then
        !! Close reaction rates group
        call h5gclose_f(gid, error)

        !! Add group attributes
        call h5ltset_attribute_string_f(fid, "/", "description", &
             "Atomic Cross Sections and Rates", error)
        call h5ltset_attribute_string_f(fid, "cross", "description", &
             "Atomic Cross Sections", error)
        call h5ltset_attribute_string_f(fid, "rates", "description", &
             "Atomic Reaction Rates", error)

        !! Close file and HDF5 interface
        call h5fclose_f(fid, error)
        call h5close_f(error)

        write(*,'(a)') "Atomic tables written to "//trim(tables_file)
        write(*,*) ''
    endif

#ifdef _MPI
    call cleanup_mpi()
#endif

    if(verbose) then
        write(*,*) 'END: hour:minute:second ', time(time_start)
    endif

end program
