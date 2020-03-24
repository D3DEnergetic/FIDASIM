[1mdiff --git a/src/fidasim.f90 b/src/fidasim.f90[m
[1mindex e215af3..036041f 100644[m
[1m--- a/src/fidasim.f90[m
[1m+++ b/src/fidasim.f90[m
[36m@@ -8271,21 +8271,29 @@[m [msubroutine bt_cx_rates(plasma, denn, an, vi, rates)[m
 [m
 end subroutine bt_cx_rates[m
 [m
[31m-subroutine get_neutron_rate(plasma, eb, rate)[m
[31m-    !+ Gets neutron rate for a beam with energy `eb` interacting with a target plasma[m
[32m+[m[32msubroutine get_ddpt_rate(plasma, eb, rate, branch)[m
[32m+[m[32m    !+ Gets d(d,p)T rate for a beam with energy `eb` interacting with a target plasma[m
     type(LocalProfiles), intent(in) :: plasma[m
         !+ Plasma Paramters[m
     real(Float64), intent(in)       :: eb[m
         !+ Beam energy [keV][m
     real(Float64), intent(out)      :: rate[m
         !+ Neutron reaction rate [1/s][m
[32m+[m[32m    integer, intent(in), optional   :: branch[m
[32m+[m[32m        !+ Indicates 1 for proton rate and 2 for neutron rate[m
 [m
[31m-    integer :: err_status, neb, nt, ebi, tii, is[m
[32m+[m[32m    integer :: err_status, neb, nt, ebi, tii, is, ib[m
     real(Float64) :: dlogE, dlogT, logEmin, logTmin[m
     real(Float64) :: logeb, logti, lograte[m
     type(InterpolCoeffs2D) :: c[m
     real(Float64) :: b11, b12, b21, b22[m
 [m
[32m+[m[32m    if(present(branch)) then[m
[32m+[m[32m        ib = branch[m
[32m+[m[32m    else[m
[32m+[m[32m        ib = 2[m
[32m+[m[32m    endif[m
[32m+[m
     logeb = log10(eb)[m
     logti = log10(plasma%ti)[m
 [m
[36m@@ -8307,7 +8315,7 @@[m [msubroutine get_neutron_rate(plasma, eb, rate)[m
     b22 = c%b22[m
     if(err_status.eq.1) then[m
         if(inputs%verbose.ge.0) then[m
[31m-            write(*,'(a)') "GET_NEUTRON_RATE: Eb or Ti out of range of D_D table. Setting D_D rates to zero"[m
[32m+[m[32m            write(*,'(a)') "get_ddpt_rate: Eb or Ti out of range of D_D table. Setting D_D rates to zero"[m
             write(*,'("eb = ",ES10.3," [keV]")') eb[m
             write(*,'("ti = ",ES10.3," [keV]")') plasma%ti[m
         endif[m
[36m@@ -8326,12 +8334,13 @@[m [msubroutine get_neutron_rate(plasma, eb, rate)[m
         rate = 0.d0[m
         do is=1,n_thermal[m
             if(thermal_mass(is).eq.H2_amu) then[m
[31m-                rate = rate + plasma%deni(is) * exp(lograte*log_10)[m
[32m+[m[32m                rate = rate + exp(lograte*log_10)[m
[32m+[m[32m             !!!rate = rate + plasma%deni(is) * exp(lograte*log_10)[m
             endif[m
         enddo[m
     endif[m
 [m
[31m-end subroutine get_neutron_rate[m
[32m+[m[32mend subroutine get_ddpt_rate[m
 [m
 subroutine get_proton_rate(plasma, v1, v3, rate)[m
     !+ Gets proton rate for a beam with velocity `v1` interacting with a target plasma[m
[36m@@ -8361,60 +8370,16 @@[m [msubroutine get_proton_rate(plasma, v1, v3, rate)[m
     real(Float64) :: ai, bi, ci, b1, b2[m
 [m
     !! Calculate effective beam energy[m
[31m- !!!vrot = plasma%vrot  ![cm/s][m
[31m-    vrot = [0.0, 0.0, 0.0][m
[32m+[m[32m    vrot = plasma%vrot  ![cm/s][m
     vrel = v1-vrot[m
     vnet_square=dot_product(vrel,vrel)  ![cm/s][m
     eb = v2_to_E_per_amu*fbm%A*vnet_square ![kev][m
     ti = plasma%ti[m
     write(*,'(T2,"ti = ",ES10.3)') ti[m
     write(*,'(T2,"eb = ",ES10.3)') eb[m
[32m+[m[32m    write(*,'(T2,"vr = [",ES10.3,",",ES10.3,",",ES10.3,"]")') vrot(1)/100,vrot(2)/100,vrot(3)/100[m
 [m
[31m-    logeb = log10(eb)[m
[31m-    logti = log10(ti)[m
[31m-[m
[31m-    !!D_D[m
[31m-    err_status = 1[m
[31m-    logEmin = tables%D_D%logemin[m
[31m-    logTmin = tables%D_D%logtmin[m
[31m-    dlogE = tables%D_D%dlogE[m
[31m-    dlogT = tables%D_D%dlogT[m
[31m-    neb = tables%D_D%nenergy[m
[31m-    nt = tables%D_D%ntemp[m
[31m-    call interpol_coeff(logEmin, dlogE, neb, logTmin, dlogT, nt, &[m
[31m-                        logeb, logti, c2D, err_status)[m
[31m-    ebi = c2D%i[m
[31m-    tii = c2D%j[m
[31m-    b11 = c2D%b11[m
[31m-    b12 = c2D%b12[m
[31m-    b21 = c2D%b21[m
[31m-    b22 = c2D%b22[m
[31m-    if(err_status.eq.1) then[m
[31m-        if(inputs%verbose.ge.0) then[m
[31m-            write(*,'(a)') "GET_PROTON_RATE: Eb or Ti out of range of D_D table. Setting D_D rates to zero"[m
[31m-            write(*,'("eb = ",ES10.3," [keV]")') eb[m
[31m-            write(*,'("ti = ",ES10.3," [keV]")') ti[m
[31m-        endif[m
[31m-        rate = 0.d0[m
[31m-        return[m
[31m-    endif[m
[31m-[m
[31m-    lograte = (b11*tables%D_D%log_rate(ebi,tii,1)   + &[m
[31m-               b12*tables%D_D%log_rate(ebi,tii+1,1) + &[m
[31m-               b21*tables%D_D%log_rate(ebi+1,tii,1) + &[m
[31m-               b22*tables%D_D%log_rate(ebi+1,tii+1,1))[m
[31m-[m
[31m-    if (lograte.lt.tables%D_D%minlog_rate) then[m
[31m-        rate = 0.d0[m
[31m-    else[m
[31m-        rate = 0.d0[m
[31m-        do is=1,n_thermal[m
[31m-            if(thermal_mass(is).eq.H2_amu) then[m
[31m-             !!!rate = rate + plasma%deni(is) * exp(lograte*log_10)[m
[31m-                rate = rate + exp(lograte*log_10)[m
[31m-            endif[m
[31m-        enddo[m
[31m-    endif[m
[32m+[m[32m    call get_ddpt_rate(plasma, eb, rate, branch=1)[m
 [m
     !!Calculate anisotropy enhancement/deficit factor[m
     mp = H1_amu*mass_u  ! kg[m
[36m@@ -8455,7 +8420,7 @@[m [msubroutine get_proton_rate(plasma, v1, v3, rate)[m
     ci = b1*abc(3,ei) + b2*abc(3,ei+1)[m
 [m
     k = (ai + bi*cos_theta**2 + ci*cos_theta**4) / (ai+bi/3.+ci/5.)[m
[31m-    k = 1[m
[32m+[m[32m    write(*,'(T2,"k  = ",ES10.3)') k[m
 [m
     rate = k * rate[m
 [m
[36m@@ -11547,7 +11512,7 @@[m [msubroutine neutron_f[m
                             erel = v2_to_E_per_amu*fbm%A*vnet_square ![kev][m
 [m
                             !! Get neutron production rate[m
[31m-                            call get_neutron_rate(plasma, erel, rate)[m
[32m+[m[32m                            call get_ddpt_rate(plasma, erel, rate)[m
                             if(inputs%calc_neutron.ge.2) then[m
                                 neutron%weight(ie,ip,ir,iz,iphi) = neutron%weight(ie,ip,ir,iz,iphi) &[m
                                                                  + rate * factor[m
[36m@@ -11721,7 +11686,7 @@[m [msubroutine neutron_mc[m
                 eb = v2_to_E_per_amu*fast_ion%A*vnet_square ![kev][m
 [m
                 !! Get neutron production rate[m
[31m-                call get_neutron_rate(plasma, eb, rate)[m
[32m+[m[32m                call get_ddpt_rate(plasma, eb, rate)[m
                 rate = rate*fast_ion%weight/ngamma*factor[m
 [m
                 !! Store neutrons[m
[36m@@ -11741,7 +11706,7 @@[m [msubroutine neutron_mc[m
             eb = v2_to_E_per_amu*fast_ion%A*vnet_square ![kev][m
 [m
             !! Get neutron production rate[m
[31m-            call get_neutron_rate(plasma, eb, rate)[m
[32m+[m[32m            call get_ddpt_rate(plasma, eb, rate)[m
             rate = rate*fast_ion%weight*factor[m
 [m
             !! Store neutrons[m
