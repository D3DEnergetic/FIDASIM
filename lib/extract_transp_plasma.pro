FUNCTION extract_transp_plasma,filename, intime, grid, equil, $
            no_omega=no_omega,doplot=doplot, sne=sne, ste=ste, sti=sti, simp=simp, srot=srot

    ;; Plot Settings    
    !p.charsize=2   &  	!x.minor=-1   & !y.minor=-1
    if keyword_set(sne) OR keyword_set(ste)  OR keyword_set(sti) $
                        OR keyword_set(simp) OR keyword_set(srot) then begin
       	loadct, 13  & !p.color=100  & !p.multi=[0,2,3]  &  !p.psym=1  & !p.thick=1
    	device, decompose=0
    	window, 0, retain=2, xs=600, ys=800
    end
    
    var_list = ["X","TIME","NE","TE","TI","ZEFFI"]
    if not keyword_set(no_omega) then begin
        var_list = [var_list, "OMEGA"]
    endif

    zz = read_ncdf(filename,vars = var_list)
    if vars.err eq 1 then begin
        return, 1
    END
    
    transp_ne = zz.ne_ ;cm^-3
    transp_te = zz.te*1.d-3  ; kev
    transp_ti = zz.ti*1.d-3   ; kev
    transp_zeff = zz.zeffi
    x = zz.x
    t = zz.time
    if keyword_set(no_omega) then begin
      transp_omega=replicate(0.,n_elements(x),n_elements(t))
    endif else begin
      transp_omega = zz.omega  ; rad/s
    endelse
    
    dummy=min( abs(t-intime), idx)
    time = t[idx]
    print, ' * Selecting profiles at :', t[idx], ' s' ;pick the closest timeslice to TOI
    print
    
    if keyword_set(sne) then begin
       z = smooth(transp_ne, [1, sne])
       plot,  x, transp_ne[*,idx], title='Ne' + inputs.transpid
       oplot, x, z[*,idx], psym=0, color=100
       transp_ne = z
    end

    if keyword_set(simp) then begin
       z = smooth(transp_zeff, [1, simp])
       plot,  x, transp_zeff[*,idx], title='Zeff' + inputs.transpid
       oplot, x, z[*,idx], psym=0, color=100
       transp_zeff = z
    end

    if keyword_set(ste) then begin
       z = smooth(transp_te, [1, ste])
       plot,  x, transp_te[*,idx], title='Te' + inputs.transpid
       oplot, x, z[*,idx], psym=0, color=100
       transp_te = z
    end

    if keyword_set(sti) then begin
       z = smooth(transp_ti, [1, sti])
       plot,  x, transp_ti[*,idx], title='Ti' + inputs.transpid
       oplot, x, z[*,idx], psym=0, color=100
       transp_ti = z
    end

    if keyword_set(srot) then begin
       z = smooth(transp_omega, [1, srot])
       plot,  x, transp_omega[*,idx], title='Omega' + inputs.transpid
       oplot, x, z[*,idx], psym=0, color=100
       transp_omega = z
    end
    
    if keyword_set(doplot) then begin
    	!p.color=220  & !p.multi=[0,2,3]  &  !p.psym=0  & !p.thick=2
    	window, 1, retain=2, xs=600, ys=800
    
    	plot, x, transp_ne[*,idx],                  ytitle=' x E13 cm-3', title='Ne  '+inputs.transpid
    	plot, x, transp_zeff[*,idx],                ytitle= 'x E13 cm-3', title='Zeff'
    	plot, x, transp_te[*,idx],                  ytitle=' keV',        title='Te'
    	plot, x, transp_ti[*,idx],    xtitle='rho', ytitle=' keV',        title='Ti'
    	plot, x, transp_omega[*,idx], xtitle='rho', ytitle='rad/s',       title='Omega'
    endif

    ;; Interpolate onto r-z grid
    dene=interpol(transp_ne[*,idx],x,equil.flux) > 0.0
    te=interpol(transp_te[*,idx],x,equil.flux) > 0.0
    ti=interpol(transp_ti[*,idx],x,equil.flux) > 0.0
    zeff=interpol(transp_zeff[*,idx],x,equil.flux) > 1.0
    vt = grid.r2d*interpol(trans_omega[*,idx],x,equil.flux)
    vr = replicate(0.0,grid.nr,grid.nz)
    vw = replicate(0.0,grid.nr,grid.nz)
    max_flux = max(abs(rho))

    s = size(flux,/dim)
    mask = intarr(s[0],s[1])
    w=where(flux le max_flux) ;where we have profiles
    mask[w] = 1


    ;;SAVE IN PROFILES STRUCTURE
    profiles={data_source:filename,time:time,max_abs_flux:max(abs(x)),dene:dene,te:te,ti:ti,vr:vr,vt:vt,vw:vw,zeff:zeff} 

    return,profiles

END
