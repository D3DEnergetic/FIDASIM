FUNCTION cfpd_table, tname, gname, earray=earray, amu=amu, z=z, nrays=nrays, step=step, nsteps=nsteps, plot_show=plot_show, flip_bpol=flip_bpol, geometry=geometry
    ;+#cfpd_table
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
    ;+    **amu**: Atomic mass unit
    ;+
    ;+    **z**: Particle charge
    ;+
    ;+    **nrays**: Number of orbits
    ;+
    ;+    **step**: Step size [m]
    ;+
    ;+    **nsteps**: Number of steps
    ;+
    ;+    **plot_show**: Plot orbit bundles -- (1) on (0) off
    ;+
    ;+    **flip_bpol**: Reverse the direction of `cpasma` in the GEQDSK file -- (1) on (0) off
    ;+
    ;+    **geometry**: Structure containing diagnostic geometry -- tags(rdist,zdist,v,d,rc) in [m,rad,s]
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> tname = 'mast.idl'
    ;+IDL> gname = 'g000001.01000'
    ;+IDL> earray = 2460. + 285.*findgen(5)
    ;+IDL> cfpd_table, tname, gname, earray=earray
    ;+```

    if ~keyword_set(earray) then earray = 2700. + 40*findgen(19)
    if ~keyword_set(amu) then amu=1
    if ~keyword_set(z) then z=1.
    if ~keyword_set(nrays) then nrays=50
    if ~keyword_set(nsteps) then nsteps=110
    if ~keyword_set(step) then step=0.01
    if ~keyword_set(plot_show) then plot_show=0
    if ~keyword_set(flip_bpol) then flip_bpol=0
    nenergy = n_elements(earray)

    g=readg(gname)
    if ~keyword_set(geometry) then begin
        print, 'Acquiring diagnostic geometry from definitions in detector_aperture_geometry.pro'
        if flip_bpol eq 1 then begin
            g.cpasma = -g.cpasma
            detector_aperture_geometry,g,1,rdist,zdist,v,d,rc
            print, 'Flipping the direction of B poloidal'
        endif else begin
            detector_aperture_geometry,g,0,rdist,zdist,v,d,rc
        endelse
    endif else begin
        print, 'User provided diagnostic geometry structure'
        rdist = geometry.rdist
        zdist = geometry.zdist
        v = geometry.v
        d = geometry.d
        rc = geometry.rc
    endelse

    orb0 = orb_cfpd(g,rdist,zdist,v,d,rc,e0=earray[0],amu=amu,z=z,nrays=nrays,step=step,nsteps=nsteps,plot_show=plot_show)
    orb = replicate(orb0, nenergy)
    orb[0] = orb0

    for i=1,nenergy-1 do begin
        orb[i] = orb_cfpd(g,rdist,zdist,v,d,rc,e0=earray[i],amu=amu,z=z,nrays=nrays,step=step,nsteps=nsteps,plot_show=plot_show)
    endfor

    save,filename=tname,orb,earray

    ;Convert to array for FIDASIM
    size_sight = size(orb0.sightline, /dimension)
    sightline = make_array([[nenergy], size_sight], /double, value=0.d0)
    daomega = make_array([[nenergy], size_sight[2:3]], /double, value=0.d0)
    nactual = make_array([[nenergy], size_sight[2:3]], /long, value=0)
    for i=0,nenergy-1 do begin
        sightline[i,*,*,*,*] = orb[i].sightline[*,*,*,*]
        daomega[i,*,*] = orb[i].daomega[*,*]
        nactual[i,*,*] = orb[i].nactual[*,*]
    endfor
    sightline *= 100.d0 ;[m] to [cm]
    sightline[*,1,*,*,*] /= 100.d0 ;vphi is ok
    sightline[*,4,*,*,*] /= 100.d0 ;phi is ok
    daomega *= 1.d4 ;[m^2] to [cm^2]

    return, {sightline:sightline, daomega:daomega, nactual:nactual, nrays:long(nrays), $
             nsteps:long(nsteps), nenergy:long(nenergy), earray:double(earray)}

END
