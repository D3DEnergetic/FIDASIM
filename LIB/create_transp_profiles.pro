FUNCTION create_transp_profiles, inputs, $
            no_omega=no_omega,doplot=doplot, sne=sne, ste=ste, sti=sti, simp=simp, srot=srot
;------------------------------------------------------------------------------------------------------
; WWH 9/10/2013
; For f90 fidasim version. 
; A modified version of the code eruskov wrote for IDL fidasim

; PURPOSE:  Write Ne,Te,Ti,Zeff & Omega profile structure from a TRANSP CDF file
;
; INPUTS:   (in inputs structure)
;           transpid     e.g., '34616B02'
;           time         (seconds)
;           profile_dir     directory where input transp CDF file resides 
;
; OUPUTS:   Profile data necessary for FIDA simulations
;
; KEYWORDS: no_omega  set toroidal rotation to zero
;           doplot - plots the five profiles
;           SMOOTHING OPTIONS: each profile can be independently smoothed over
;                              the specified number of TIME points
;

!p.charsize=2   &  	!x.minor=-1   & !y.minor=-1
if keyword_set(sne) OR keyword_set(ste)  OR keyword_set(sti) $
                    OR keyword_set(simp) OR keyword_set(srot) then begin
   	loadct, 13  & !p.color=100  & !p.multi=[0,2,3]  &  !p.psym=1  & !p.thick=1
	device, decompose=0
	window, 0, retain=2, xs=600, ys=800
end

; Find and open netcdf file with transp profile data
cdf_file = STRSPLIT(inputs.cdf_file,'/',/EXTRACT)
cdf_file = STRSPLIT(cdf_file[-1],'_',/EXTRACT)
inputs = create_struct(inputs,"transpid",cdf_file[0])

cdf_file = strupcase(cdf_file[0]) + '.CDF'
slash=strmid(inputs.profile_dir,0,1,/reverse_offset)
if slash ne '/' then begin
   cdfname  = inputs.profile_dir + '/' + cdf_file
endif else begin
   cdfname  = inputs.profile_dir + cdf_file
endelse
test=findfile(cdfname)
if test[0] eq '' then begin
  print," * FATAL ERROR:  Can't find TRANSP file: " + cdfname
  return,{err:1}
endif
print
print, " * Found file: ", cdfname
cdfid=NCDF_Open(cdfname,/nowrite)

;-------------------------------------------------------------------------------
;----  Start reading the profiles and defining profile variables
;----  appropriate for saving in the GA / FIDA profile units format
;-------------------------------------------------------------------------------
z=read_netcdf_data(cdfid, 'ne')
if z.error eq -1 then begin
  print,'ne profile unavailable'
  return,{err:1}
end

tmp=(*z.axis[0]).data
x=reform(tmp(*,0))
t=(*z.axis[1]).data

transp_ne = z.data * 1.e6 ; m-3

z=read_netcdf_data(cdfid, 'te')
if z.error eq -1 then begin
  print,'te profile unavailable'
  return,{err:1}
end
transp_te = z.data  ; ev 

z=read_netcdf_data(cdfid, 'ti')
if z.error eq -1 then begin
  print,'ti profile unavailable'
  return,{err:1}
end
transp_ti = z.data   ; ev 

z=read_netcdf_data(cdfid, 'zeffi')
if z.error eq -1 then begin
  print,'zeffi profile unavailable'
  return,{err:1}
end
transp_zeff = z.data 

if keyword_set(no_omega) then begin
  transp_omega=replicate(0.,n_elements(x),n_elements(t))
endif else begin
  z=read_netcdf_data(cdfid, 'omega')
  if z.error eq -1 then begin
    print,'omega profile unavailable'
    return,{err:1}
  end
  transp_omega = z.data  ; rad/s
endelse

NCDF_Close,cdfid  ; Close the NetCDF file

dummy=min( abs(t-inputs.time), idx)
print, ' * Selecting profiles at :', t[idx], ' s' ;pick the closest timeslice to TOI
print

;-----------------------------------------------------------------------------
;--- Smoothing section
;-----------------------------------------------------------------------------
;------- Ne
if keyword_set(sne) then begin
   z = smooth(transp_ne, [1, sne])
   plot,  x, transp_ne[*,idx], title='Ne' + inputs.transpid
   oplot, x, z[*,idx], psym=0, color=100
   transp_ne = z
end
;------- Zeff
if keyword_set(simp) then begin
   z = smooth(transp_zeff, [1, simp])
   plot,  x, transp_zeff[*,idx], title='Zeff' + inputs.transpid
   oplot, x, z[*,idx], psym=0, color=100
   transp_zeff = z
end
;------- Te
if keyword_set(ste) then begin
   z = smooth(transp_te, [1, ste])
   plot,  x, transp_te[*,idx], title='Te' + inputs.transpid
   oplot, x, z[*,idx], psym=0, color=100
   transp_te = z
end
;------- Ti
if keyword_set(sti) then begin
   z = smooth(transp_ti, [1, sti])
   plot,  x, transp_ti[*,idx], title='Ti' + inputs.transpid
   oplot, x, z[*,idx], psym=0, color=100
   transp_ti = z
end
;------- Omega
if keyword_set(srot) then begin
   z = smooth(transp_omega, [1, srot])
   plot,  x, transp_omega[*,idx], title='Omega' + inputs.transpid
   oplot, x, z[*,idx], psym=0, color=100
   transp_omega = z
end

;======= Optional plotting section ============================================================
if keyword_set(doplot) then begin
	!p.color=220  & !p.multi=[0,2,3]  &  !p.psym=0  & !p.thick=2
	window, 1, retain=2, xs=600, ys=800

	plot, x, transp_ne[*,idx],                  ytitle=' x E13 cm-3', title='Ne  '+inputs.transpid
	plot, x, transp_zeff[*,idx],                ytitle= 'x E13 cm-3', title='Zeff'
	plot, x, transp_te[*,idx],                  ytitle=' keV',        title='Te'
	plot, x, transp_ti[*,idx],    xtitle='rho', ytitle=' keV',        title='Ti'
	plot, x, transp_omega[*,idx], xtitle='rho', ytitle='rad/s',       title='Omega'
endif
;==============================================================================================

		;;SAVE IN PROFILES STRUCTURE
profiles={rho:x,dene:transp_ne[*,idx],te:transp_te[*,idx],ti:transp_ti[*,idx], $
          vtor:transp_omega[*,idx],zeff:transp_zeff[*,idx],err:0} 
return,profiles

END
