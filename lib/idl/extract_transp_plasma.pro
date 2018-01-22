FUNCTION extract_transp_plasma,filename, intime, grid, flux, $
            doplot=doplot, profiles=profiles, $
            sne=sne, ste=ste, sti=sti, simp=simp, srot=srot
    ;+#extract_transp_plasma
    ;+Extracts `plasma` structure from a TRANSP run
    ;+***
    ;+##Arguments
    ;+    **filename**: TRANSP output file e.g. [TRANSP_RUNID].CDF
    ;+
    ;+    **intime**: Time of interest [s]
    ;+
    ;+    **grid**: Interpolation grid
    ;+
    ;+    **flux**: Normalized square root of torodial flux("rho") mapped onto the interpolation grid
    ;+
    ;+##Keyword Arguments
    ;+    **doplot**: Plot profiles
    ;+
    ;+    **profiles**: Set this keyword to a named variable to recieve the plasma profiles as a function of rho
    ;+
    ;+    **s(ne|te|ti|imp|rot)**: Smooth profiles
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> plasma = extract_transp_plasma("./142332H01.CDF", 1.2, grid, flux)
    ;+```

    var_list = ["TRFLX","TFLUX","TIME","NE","TE","TI","ZEFFI","OMEGA","dn0wd","dn0out"]

    zz = read_ncdf(filename,vars = var_list)

    t = zz.time
    dummy=min( abs(t-intime), idx)
    time = double(t[idx])

    transp_ne = zz.ne_[*,idx] ;cm^-3
    transp_te = zz.te[*,idx]*1.d-3  ; kev
    transp_ti = zz.ti[*,idx]*1.d-3   ; kev
    transp_nn = zz.dn0wd[*,idx] ;cm^-3
    transp_zeff = zz.zeffi[*,idx]
    rho_cb = sqrt(zz.trflx[*,idx]/zz.tflux[idx])
    ; center each rho b/c toroidal flux is at cell boundary
    rho = 0.d0*rho_cb
    rho[0] = 0.5*rho_cb[0]
    for i=1,n_elements(rho_cb)-1 do begin
        rho[i] = rho_cb[i] - 0.5*(rho_cb[i] - rho_cb[i-1])
    endfor

    if total(strmatch(tag_names(zz),'OMEGA',/fold_case)) eq 0 then begin
      warn,'OMEGA not found in TRANSP file. Assuming no plasma rotation'
      transp_omega=0.0*rho
    endif else begin
      transp_omega = zz.omega[*,idx]  ; rad/s
    endelse

    print, ' * Selecting profiles at :', time, ' s' ;pick the closest timeslice to TOI

    if keyword_set(doplot) then begin
        !p.charsize=2 & !x.minor=-1 & !y.minor=-1
       	loadct, 13  & !p.color=100  & !p.multi=[0,2,3]  &  !p.psym=1  & !p.thick=1
    	device, decompose=0
    	window, 0, retain=2, xs=600, ys=800
        if keyword_set(sne) then begin
           z = smooth(transp_ne, sne)
           plot,  x, transp_ne, title='Ne' 
           oplot, x, z, psym=0, color=100
           transp_ne = z
        end

        if keyword_set(simp) then begin
           z = smooth(transp_zeff, simp)
           plot,  x, transp_zeff, title='Zeff'
           oplot, x, z, psym=0, color=100
           transp_zeff = z
        end

        if keyword_set(ste) then begin
           z = smooth(transp_te, ste)
           plot,  x, transp_te, title='Te'
           oplot, x, z, psym=0, color=100
           transp_te = z
        end

        if keyword_set(sti) then begin
           z = smooth(transp_ti, sti)
           plot,  x, transp_ti, title='Ti'
           oplot, x, z, psym=0, color=100
           transp_ti = z
        end

        if keyword_set(srot) then begin
           z = smooth(transp_omega, srot)
           plot,  x, transp_omega, title='Omega'
           oplot, x, z, psym=0, color=100
           transp_omega = z
        end
    	!p.color=220  & !p.multi=[0,2,3]  &  !p.psym=0  & !p.thick=2
    	window, 1, retain=2, xs=600, ys=800

    	plot, x, transp_ne,                  ytitle=' x E13 cm-3', title='Ne  '
    	plot, x, transp_zeff,                ytitle= 'x E13 cm-3', title='Zeff'
    	plot, x, transp_te,                  ytitle=' keV',        title='Te'
    	plot, x, transp_ti,    xtitle='rho', ytitle=' keV',        title='Ti'
    	plot, x, transp_omega, xtitle='rho', ytitle='rad/s',       title='Omega'
        plot, x, transp_nn,                  ytitle='cm-3',	   title='Nn'
    endif

    profiles = {rho:rho, $
                dene:transp_ne > 0.0, $
		dennw:transp_nn > 0.0, $
                te:transp_te > 0.0, $
                ti:transp_ti > 0.0, $
                zeff:transp_zeff > 1.0, $
                omega:transp_omega}

    ;; Interpolate onto r-z grid
    dene=interpol(transp_ne,rho,flux) > 0.0
    denn=interpol(transp_nn,rho,flux) > 0.0
    te=interpol(transp_te,rho,flux) > 0.0
    ti=interpol(transp_ti,rho,flux) > 0.0
    zeff=interpol(transp_zeff,rho,flux) > 1.0
    vt = double(grid.r2d*interpol(transp_omega,rho,flux))
    vr = double(replicate(0.0,grid.nr,grid.nz))
    vz = double(replicate(0.0,grid.nr,grid.nz))
    max_flux = max(abs(rho))

    s = size(flux,/dim)
    mask = intarr(s[0],s[1])
    w=where(flux le max_flux) ;where we have profiles
    mask[w] = 1


    ;;SAVE IN PROFILES STRUCTURE
    plasma={data_source:file_expand_path(filename),time:time,mask:mask, $
            dene:dene,denn:denn,te:te,ti:ti,vr:vr,vt:vt,vz:vz,zeff:zeff}

    return,plasma

END
