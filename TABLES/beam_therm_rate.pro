function beam_therm_rate, temp,eb,am,ab,de,sig_name  $
                          ,params=params,deexcit=deexcit,mc=mc
;; Program to calculate beam-plasma reactivity numerically
;; INPUT
;;       temp    Plasma temperature (keV)
;;       eb      Beam energy (keV)
;;       am      Atomic mass of Maxwellian species
;;       ab      Atomic mass of beam species
;;       de      energy difference of final to initial state [keV]
;;       sig_name        name of cross section function (string) 
;;               accepts e in keV/u; returns sigma in cm^2
;; KEYWORD
;;	params	additional parameters needed by function sig_name
;;	deexcit  calculate rate for de-excitation
;; OUTPUT
;;       <sigma*v>       cm^3/s

  if temp eq 0. and eb eq 0. then return,0.
  if keyword_set(mc) then begin
     mass_u=1.6605e-27          ;[kg]   
     ec= 1.602e-19              ;[C]
     nv=5.e6
     ;; define thermal velocity distribution of ions (monte carlo)
     vti=sqrt(temp*1.e3*ec/(mass_u*am)) ;[m/s]
     vi=fltarr(nv,3)
     vi[*,0]=randomn(seed,nv)*vti
     vi[*,1]=randomn(seed,nv)*vti
     vi[*,2]=randomn(seed,nv)*vti
     ;; define beam velocity vector
     vb=sqrt(2.d0*eb*1.d3*ec/(mass_u*ab))       ;[m/s]
     vrel=sqrt((vi[*,0]-vb)^2+vi[*,1]^2+vi[*,2]^2) ; [m/s]
     ared = am*ab/(am+ab) 
     e_coll = 0.5d0*mass_u*ared*vrel^2/ec*1.e-3 ; [kev/amu]  
     factor = replicate(1.,nv)
     if keyword_set(deexcit) then begin
        factor = (e_coll + de)/e_coll*params[0]^2/params[1]^2
        index=where(finite(factor) eq 0,nind)
        if nind gt 0 then factor[index]=0.
        e_coll = e_coll + de
     endif
     if ared ge 0.5 then begin
        ;; ion-ion collision
        e_for_sigma_call = e_coll/ared ;[keV/amu]
     endif else begin
        ;; electron ion collision
        e_for_sigma_call = e_coll ;[keV]
     endelse    
     ;;get the cross sections [cm^2]
     case n_elements(params) of
        1: sigarr=call_function(sig_name,e_for_sigma_call,params[0])
        2: sigarr=call_function(sig_name,e_for_sigma_call,params[0],params[1])
        3: sigarr=call_function(sig_name,e_for_sigma_call,params[0],params[1],params[2])
        else: stop
     endcase 
     return,mean(factor*vrel*100.*sigarr) ; [cm^3/s]
  endif


  ;; Analytic solution
  temp_target_per_am = (temp > 1.d-6)/am
  e_beam_per_amu = eb/ab
  ared = am*ab/(am+ab) 
  ;;thermal velocity of target = sqrt(2*k_b T/m_target)
  v_therm = 1.384d6 * sqrt(temp_target_per_am*1.e3) ;[cm/s] 
  ;;v_beam/v_therm = sqrt(2 E_beam/m_beam) / sqrt(2 k_B T/m_target)
  zb = sqrt(e_beam_per_amu/temp_target_per_am)
  ;;proportionality constant from relative velocity u in units of 
  ;;v_therm to reduced energy (collision energy): 
  ;;Erel = .5 m_red u^2 = .5* m_red * v_therm^2 * (u/v_therm)^2 
  ;;     = .5* m_red * 2 * k_B T /m_target * (u/v_therm)^2 
  ;;     =  m_red * k_B T / m_target   * (u/v_therm)^2
  ;;     =  m_red * temp_target_per_am * (u/v_therm)^2
  u2_to_e_rel =  ared * temp_target_per_am 
  ;;grid in v_R direction (x,y-plane)
  nr = 31
  rmax = 4.d0
  dr = rmax/(nr-1)
  r = dindgen(nr)*dr
  kr = lindgen((nr-1)/2)
  ;;grid in v_Z direction (beam velocity is in this direction)
  nz = 61 
  zmax = 4.d0
  dz = 2.d0*zmax/(nz-1)
  z = dindgen(nz)*dz-zmax
  kz = indgen((nz-1)/2)
  fz = dblarr(nz)
  for iz=0,nz-1 do begin ;; loop along beam axis
     u2 = (zb-z[iz])^2+r[*]^2
     e_coll = u2_to_e_rel * u2[*] ;[keV]   
     factor = replicate(1.,nr)
     if keyword_set(deexcit) then begin
        factor = (e_coll + de)/e_coll
        index=where(finite(factor) eq 0,nind)
        if nind gt 0 then factor[index]=0.
        e_coll = e_coll + de
     endif
     if ared ge 0.5 then begin
        ;; ion-ion collision
        e_for_sigma_call = e_coll/ared ;[keV/amu]
     endif else begin
        e_for_sigma_call = e_coll ;[keV]
     endelse     
     ;;get the cross sections [cm^2] 
     case n_elements(params) of
        1: sigarr=call_function(sig_name,e_for_sigma_call,params[0])
        2: sigarr=call_function(sig_name,e_for_sigma_call,params[0],params[1])
        3: sigarr=call_function(sig_name,e_for_sigma_call,params[0],params[1],params[2])
        else: stop
     endcase 
     ;; sort out NANs and INFINITS
     index=where(finite(sigarr) eq 0,nind)
     if nind gt 0 then sigarr[index]= 1.d-35
     ;; if energy is too low, then set sigma to zero!
     index=where(e_coll lt de,nind)
     if nind gt 0 then sigarr[index]= 0.
  
     f = factor[*]*sigarr[*]*sqrt(u2[*])*exp(-(z[iz]^2+r[*]^2))*r[*]
     ;;Simpson's Rule
     fz[iz] = (-f[0]+f[nr-1]+4.d0*total(f[2*kr+1],/double)+2.*total(f[2*kr],/double))*dr/3.
  endfor
  ;;Simpson's Rule
  sig_eff = (-fz[0]+fz[nz-1]+4.d0*total(fz[2*kz+1],/double)+2.*total(fz[2*kz],/double))*dz/3.
  sig_eff = sig_eff * 2./sqrt(!Dpi)
  rate = sig_eff*v_therm
  if keyword_set(deexcit) then  rate=rate*params[0]^2/params[1]^2
  return, rate



end
