FUNCTION extract_transp_plasma,filename, intime, grid, rhogrid, $
            doplot=doplot, profiles=profiles, dn0out=dn0out, $
            scrapeoff=scrapeoff, rho_scrapeoff=rho_scrapeoff,$
            sne=sne, ste=ste, sti=sti, simp=simp, srot=srot, snn=snn
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
    ;+    **rhogrid**: sqrt(normalized torodial flux) mapped onto the interpolation grid
    ;+
    ;+##Keyword Arguments
    ;+    **doplot**: Plot profiles
    ;+
    ;+    **profiles**: Set this keyword to a named variable to recieve the plasma profiles as a function of rho
    ;+
    ;+    **s(ne|ni|te|ti|imp|rot|nn)**: Smooth profiles
    ;+
    ;+    **dn0out**: Wall Neutral density value `dn0out` variable in transp namelist
    ;+
    ;+    **scrapeoff**: scrapeoff decay length
    ;+
    ;+    **rho_scrapeoff**: scrapeoff length, default = 0.1
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> plasma = extract_transp_plasma("./142332H01.CDF", 1.2, grid, rho)
    ;+```

    var_list = ["X","TRFLX","TFLUX","TIME","NE","NH","ND","NT","NIMP","TE","TI","ZEFFI","OMEGA","DN0WD","XZIMP"]

    zz = read_ncdf(filename,vars = var_list)
    tn = tag_names(zz)

    t = zz.time
    dummy=min( abs(t-intime), idx)
    time = double(t[idx])

    print, ' * Selecting profiles at :', time, ' s' ;pick the closest timeslice to TOI

    impurity_charge = fix(ceil(max(zz.xzimp)))
    ; Densities
    transp_ne = zz.ne_[*,idx] ;cm^-3
    transp_nimp = zz.nimp[*,idx] ;cm^-3
    transp_nn = zz.dn0wd[*,idx] ;cm^-3
    w_thermal = []
    _ = where(tn eq "NH", nmatch)
    if nmatch eq 1 then begin
        transp_nh = zz.nh[*,idx]
        w_thermal = [w_thermal, 0]
    endif else begin
        transp_nh = transp_ne*0
    endelse
    _ = where(tn eq "ND", nmatch)
    if nmatch eq 1 then begin
        transp_nd = zz.nd[*,idx]
        w_thermal = [w_thermal, 1]
    endif else begin
        transp_nd = transp_ne*0
    endelse
    _ = where(tn eq "NT", nmatch)
    if nmatch eq 1 then begin
        transp_nt = zz.nt[*,idx]
        w_thermal = [w_thermal, 2]
    endif else begin
        transp_nt = transp_ne*0
    endelse

    transp_te = zz.te[*,idx]*1.d-3  ; kev
    transp_ti = zz.ti[*,idx]*1.d-3   ; kev
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
      transp_omega=0.0*transp_te
    endif else begin
      transp_omega = zz.omega[*,idx]  ; rad/s
    endelse

    if not keyword_set(dn0out) then dn0out = transp_nn[-1]
    if not keyword_set(scrapeoff) then scrapeoff = 0.0
    if not keyword_set(rho_scrapeoff) then rho_scrapeoff = 0.1

    if scrapeoff gt 0.0 then begin
        drho = abs(rho[-1] - rho[-2])
        rho_sc = rho[-1] + drho*(dindgen(ceil(rho_scrapeoff/drho)) + 1)
        sc = exp(-(rho_sc - rho[-1])/scrapeoff)
        transp_ne   = [transp_ne,transp_ne[-1]*sc]
        transp_nimp = [transp_nimp,transp_nimp[-1]*sc]
        transp_nh   = [transp_nh,transp_nh[-1]*sc]
        transp_nd   = [transp_nd,transp_nd[-1]*sc]
        transp_nt   = [transp_nt,transp_nt[-1]*sc]
        transp_te   = [transp_te,transp_te[-1]*sc]
        transp_ti   = [transp_ti,transp_ti[-1]*sc]
        transp_nn   = [transp_nn,0*sc + dn0out]
        transp_zeff = [transp_zeff, (transp_zeff[-1]-1)*sc + 1]
        transp_omega = [transp_omega,transp_omega[-1]*sc]
        rho = [rho, rho_sc]
    endif

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
           z = smooth(transp_nimp, simp)
           plot,  x, transp_nimp, title='Nimp'
           oplot, x, z, psym=0, color=100
           transp_nimp = z
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

        if keyword_set(snn) then begin
           z = smooth(transp_nn, snn)
           plot,  x, transp_nn, title='Omega'
           oplot, x, z, psym=0, color=100
           transp_nn = z
        end
    	!p.color=220  & !p.multi=[0,2,3]  &  !p.psym=0  & !p.thick=2
    	window, 1, retain=2, xs=600, ys=800

    	plot, x, transp_ne,                  ytitle=' x E13 cm-3', title='Ne  '
    	plot, x, transp_nimp,                ytitle= 'x E13 cm-3', title='Nimp'
    	plot, x, transp_te,                  ytitle=' keV',        title='Te'
    	plot, x, transp_ti,    xtitle='rho', ytitle=' keV',        title='Ti'
    	plot, x, transp_omega, xtitle='rho', ytitle='rad/s',       title='Omega'
        plot, x, transp_nn,                  ytitle='cm-3',	   title='Nn'
    endif

    profiles = {rho:rho, $
                impurity_charge:impurity_charge,$
                dene:transp_ne > 0.0, $
                denimp:transp_nimp > 0.0, $
		denn:transp_nn > 0.0, $
                te:transp_te > 0.0, $
                ti:transp_ti > 0.0, $
                zeff:transp_zeff > 1.0, $
                omega:transp_omega}

    if where(0 eq w_thermal) ge 0 then begin
        profiles = create_struct(profiles,"denh",transp_nh>0.0)
    endif
    if where(1 eq w_thermal) ge 0 then begin
        profiles = create_struct(profiles,"dend",transp_nd>0.0)
    endif
    if where(2 eq w_thermal) ge 0 then begin
        profiles = create_struct(profiles,"dent",transp_nt>0.0)
    endif

    ;; Interpolate onto r-z grid
    dene=interpol(transp_ne,rho,rhogrid) > 0.0
    denimp=interpol(transp_nimp,rho,rhogrid) > 0.0
    denh=interpol(transp_nh,rho,rhogrid) > 0.0
    dend=interpol(transp_nd,rho,rhogrid) > 0.0
    dent=interpol(transp_nt,rho,rhogrid) > 0.0
    denn=(10.d0^interpol(alog10(transp_nn),rho,rhogrid)) > 0.0
    te=interpol(transp_te,rho,rhogrid) > 0.0
    ti=interpol(transp_ti,rho,rhogrid) > 0.0
    zeff= interpol(transp_zeff,rho,rhogrid) > 1.0
    vt = double(grid.r2d*interpol(transp_omega,rho,rhogrid))
    vr = double(replicate(0.0,grid.nr,grid.nz))
    vz = double(replicate(0.0,grid.nr,grid.nz))
    max_rho = max(abs(rho))

    s = size(rhogrid,/dim)
    mask = intarr(s[0],s[1])
    w=where(rhogrid le max_rho) ;where we have profiles
    mask[w] = 1

    deni = dblarr(3,s[0],s[1])
    deni[0,*,*] = denh
    deni[1,*,*] = dend
    deni[2,*,*] = dent
    ai = [1.007276466879d0, 2.013553212745d0,3.01550071632d0]

    ;;SAVE IN PROFILES STRUCTURE
    plasma={data_source:file_expand_path(filename),time:time, mask:mask, impurity_charge:impurity_charge, $
            nthermal:fix(n_elements(w_thermal)), species_mass:ai[w_thermal], dene:dene,deni:deni[w_thermal,*,*], $
            denimp:denimp,denn:denn,te:te,ti:ti,vr:vr,vt:vt,vz:vz,zeff:zeff}

    return,plasma

END
