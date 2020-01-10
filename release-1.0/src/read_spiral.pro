FUNCTION read_spiral_header, file

    openr,lun,file,/get_lun

    npart = 0
    weight = 0.0
    gc = 0
    err = 0
    nheader = 0
    line = ''
    header = ''
    readf, lun, line
    while strmid(line,0,1) eq ';' do begin
        nheader = nheader + 1
        header = header + line + string(10B)
        if stregex(line,"particle weight",/fold_case) ne -1 then begin
            weight = float(stregex(line,"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?",/extract))
        endif
        if stregex(line,"number",/fold_case) ne -1 then begin
            npart = long(stregex(line,"[0-9]+",/extract))
        endif
        if stregex(line,"guiding center",/fold_case) ne -1 then begin
            gc = 1
        endif
        if stregex(line,"R_[m]",/fold_case) ne -1 then begin
            if n_elements(strsplit(line)) ne 6 then begin
                err = 1
            endif
        endif
        line = ''
        readf, lun, line
    endwhile
    free_lun,lun
    if gc eq 0 then begin
        error,'Not a guiding center distribution',/halt
    endif

    if npart eq 0 then begin
        error,'Number of particles le 0',/halt
    endif

    return, {nheader:nheader, npart:npart, header:header,weight:weight}

END

FUNCTION finite_struct, s
    n = N_TAGS(s)
    for i=0,n-1 do begin
        if finite(s.(i)) ne 1 then return, 0
    endfor
    return, 1
END

;+ This file contains the procedure to read SPIRAL fast-ion distribution file
FUNCTION read_spiral,file, time=time, ntotal=ntotal, e_range=e_range, p_range=p_range, $
                           particle_weight=particle_weight, btipsign=btipsign
    ;+#read_spiral
    ;+Reads SPIRAL guiding center fast-ion distribution file
    ;+***
    ;+##Arguments
    ;+    **file**: SPIRAL output file
    ;+
    ;+##Keyword Arguments
    ;+    **time**: Time [s]
    ;+
    ;+    **ntotal**: Total number of fast-ions
    ;+
    ;+    **e_range**: Energy range of particles to consider
    ;+
    ;+    **p_range**: Pitch range of particles to consider
    ;+
    ;+    **particle_weight**: Set particle/marker weight such that sum(particle_weights) = ntotal: Defaults to `ntotal`/nparticles
    ;+
    ;+    **btipsign**: Sign of the dot product between the current and magnetic field (Required)
    ;+
    ;+##Return Value
    ;+Distribution structure
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> dist = read_spiral("./mc_159243H06_9",ntotal=1e19,btipsign=-1)
    ;+```

    if not keyword_set(btipsign) then begin
        error,'btipsign is not set.',/halt
    endif

    if not file_test(file) then begin
        error, 'Nonexistent file: '+file,/halt
    endif

    if not keyword_set(time) then begin
        print, 'WARNING: Time is not set for SPIRAL distribution. Setting to 0.0 [s]'
        time = 0.d0
    endif

    header = read_spiral_header(file)
    npart = header.npart
    nhead = header.nheader
    pweight = header.weight

    openr,unit,file, /get_lun

    ; Read header
    head = strarr(nhead)
    readf, unit, head

    ; Read in data
    r = dblarr(npart)
    z = dblarr(npart)
    phi = dblarr(npart)
    energy = dblarr(npart)
    pitch = dblarr(npart)
    cnt = 0L
    for i=0,npart-1 do begin
        s = {r:double(0),phi:double(0),z:double(0),energy:double(0),pitch:double(0)}
        readf,unit,s
        if not finite_struct(s) then begin
            continue
        endif else begin
            r[cnt] = s.r*100
            phi[cnt] = s.phi
            z[cnt] = s.z*100
            energy[cnt] = s.energy
            pitch[cnt] = s.pitch*btipsign ; SPIRAL pitch is defined relative to current
            cnt = cnt+1
        endelse
    endfor
    free_lun, unit

    npart = cnt
    r = r[0:npart-1]
    z = z[0:npart-1]
    phi = phi[0:npart-1]
    energy = energy[0:npart-1]
    pitch = pitch[0:npart-1]
    orbit_class = replicate(1,npart)
    if pweight ne 0.0 then begin
        weight = replicate(pweight,npart)
    endif else begin
        if not keyword_set(particle_weight) then begin
            if not keyword_set(ntotal) then begin
                print, 'WARNING: ntotal is not set. Setting arbitrarily to 1e19'
                ntotal = 1.d19
            endif
            weight = replicate(Ntotal/float(npart), npart)
        endif else begin
            weight = replicate(particle_weight, npart)
        endelse
    endelse

    if not keyword_set(e_range) then begin
        e_range = [min(energy),max(energy)]
    endif
    if not keyword_set(p_range) then begin
        p_range = [min(pitch),max(pitch)]
    endif

    ww = where(energy ge e_range[0] and energy le e_range[1],nw)
    if nw eq 0 then begin
        error,'No particles fall in requested energy range',/halt
    endif
    wwp = where(pitch[ww] ge p_range[0] and pitch[ww] le p_range[1],nwp)
    if nwp eq 0 then begin
        error,'No particles fall in the requested pitch range',/halt
    endif
    ww = ww[wwp]
    nw = n_elements(ww)
    print,'Number of markers: ',npart
    print,'Number of markers in phase space: ',nw

    dist_struct = {type:2,time:time,data_source:file_expand_path(file), $
                   nparticle:long(nw),nclass:1,r:r[ww],z:z[ww],phi:phi[ww],$
                   energy:energy[ww],pitch:pitch[ww],class:orbit_class[ww],$
                   weight:weight[ww]}

    return, dist_struct
END
