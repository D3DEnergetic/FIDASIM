FUNCTION proton_table, tname, gname, earray=earray, nrays=nrays, nsteps=nsteps, plot_show=plot_show, flip_bpol=flip_bpol, geometry=geometry
    ;+#proton_table
    ;+Generates charged fusion product table
    ;+***
    ;+##Arguments
    ;+    **tname**: Filename of output .idl table file
    ;+
    ;+    **gname**: GEQDSK file
    ;+
    ;+##Keyword Arguments
    ;+    **earray**: Array of energies (keV) to use
    ;+
    ;+    **nrays**: Number of orbits
    ;+
    ;+    **nsteps**: Number of steps
    ;+
    ;+    **plot_show**: Interactively plot orbit bundles -- (1) on (0) off
    ;+
    ;+    **flip_bpol**: Reverse the direction of `cpasma` in the GEQDSK file -- (1) on (0) off
    ;+
    ;+    **geometry**: Structure containing diagnostic geometry -- tags(rdist,zdist,v,d,rc)
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> tname = 'mast.idl'
    ;+IDL> gname = 'g99999K26'
    ;+IDL> earray = 2730. + 150.*findgen(5)
    ;+IDL> proton_table, tname, gname, earray=earray
    ;+```

    if ~keyword_set(earray) then earray = 2700. + 40*findgen(19)
    if ~keyword_set(nrays) then nrays=50
    if ~keyword_set(nsteps) then nsteps=110
    if ~keyword_set(plot_show) then plot_show=0
    if ~keyword_set(flip_bpol) then flip_bpol=0
    nenergy = n_elements(earray)

    g=readg(gname)  ;Netepenko
    if ~keyword_set(geometry) then begin
        if flip_bpol eq 1 then begin
            g.cpasma = -g.cpasma
            detector_aperture_geometry,g,1,rdist,zdist,v,d,rc
        endif else begin
            detector_aperture_geometry,g,0,rdist,zdist,v,d,rc
        endelse
    endif else begin
        rdist = geometry.rdist
        zdist = geometry.zdist
        v = geometry.v
        d = geometry.d
        rc = geometry.rc
    endelse

    orb0 = orb_proton(g,rdist,zdist,v,d,rc,e0=earray[0],nrays=nrays,nsteps=nsteps,plot_show=plot_show)
    orb = replicate(orb0, nenergy)
    orb[0] = orb0

    for i=1,nenergy-1 do begin
        orb[i] = orb_proton(g,rdist,zdist,v,d,rc,e0=earray[i],nrays=nrays,nsteps=nsteps,plot_show=plot_show)
    endfor

    save,filename=tname,orb,earray

    ;Convert to array for FIDASIM
    size_sight = size(orb0.sightline, /dimension)
    sightline = make_array([[nenergy], size_sight], /double, value=0.0)
    daomega = make_array([[nenergy], size_sight[2:3]], /double, value=0.0)
    nactual = make_array([[nenergy], size_sight[2:3]], /double, value=0.0)
    for i=0,nenergy-1 do begin
        sightline[i,*,*,*,*] = orb[i].sightline
        daomega[i,*,*] = orb[i].daomega
        nactual[i,*,*] = orb[i].nactual
    endfor

    return, {sightline:sightline, daomega:daomega, nactual:nactual, nch:size_sight[3], nrays:nrays, $
             nsteps:nsteps, nenergy:nenergy, earray:earray}

END
