FUNCTION read_mc_nubeam,infile,ntotal=ntotal,e_range=e_range,particle_weight = particle_weight,btipsign=btipsign

    if not keyword_set(btipsign) then begin
        print,'ERROR: btipsign is not set.'
        goto, GET_OUT
    endif

    if not keyword_set(ntotal) and not keyword_set(particle_weight) then begin
        print,'WARNING: ntotal not set. Setting arbitrarily to 1e19'
        ntotal = 1.d19
    endif

    zzz=FINDFILE(infile)
    if zzz eq '' then begin
        print, 'ERROR: Nonexistent file:', infile
        goto, GET_OUT
    endif
    
    openr,unit,infile, /get_lun
    line=' '
    readf,unit,line           ; read string
    readf,unit,line
    pos=strpos(line,'N=')
    if pos eq -1 then begin
        print,'ERROR: Second line is missing the number of points'
        goto, GET_OUT
    endif

    parts=str_sep(line, ' ')
    w=where(parts eq 'N=')
    npts=long(parts[w[0]+1])
    if npts lt 5 then begin
        print,'ERROR: Too few points',npts
        goto, GET_OUT
    endif
    
    if not keyword_set(particle_weight) then begin
        particle_weight = ntotal/float(npts)
    endif

    ; Get time
    readf,unit,line
    parts=str_sep(line, ' ')
    w=where(parts eq 'TIME',nw)
    if nw eq 0 then begin
        print,'ERROR: Time not found on 3rd line'
        goto, GET_OUT
    endif
    i=1
    while 1 do begin
        s = parts[w[0]+i]
        if s ne '' and s ne '=' then begin
            time = float(s)
            break
        endif
        i = i+1
    endwhile

    data=fltarr(4,npts)
    
    ready=0
    while not ready do begin
        readf,unit,line  ; Description line
        pos=strpos(line,'R(cm)')
        if pos gt -1 then ready=1
    endwhile
    
    for i=0L,npts-1 do begin
        readf,unit,line           ; read string
        line=strcompress(line)
        parts = str_sep(line, ' ') & np=n_elements(parts)
        if parts(0) eq '' then parts=shift(parts,-1) ; accommodate blank before 1st

        while parts(np-1) eq '' do np=np-1

        if np ne 4 then begin
            print,'ERROR: Wrong number of entries on line:',line
            goto, GET_OUT
        endif else begin
            data[*,i]=parts[0:np-1]
        endelse
        print,format='(f7.2,"%",A,$)',100.0*(i+1)/float(npts),string(13b)
    endfor   
    free_lun, unit

    r=reform(data[0,0:npts-1])
    w=reform(data[1,0:npts-1])
    pitch=reform(data[2,0:npts-1])*btipsign
    energy=reform(data[3,0:npts-1])*1.d-3 ;keV
    weight = replicate(particle_weight,npts)
    orbit_class = replicate(1,npts)
    if not keyword_set(e_range) then begin
        e_range = [min(energy),max(energy)]
    endif

    ww = where(energy ge e_range[0] and energy le e_range[1],nw)
    if nw eq 0 then begin
        print,'ERROR: No particles fall in requested energy range'
        goto, GET_OUT
    endif
    print,'Number of markers: ',npts
    print,'Number of markers in energy range: ',nw

    fbm_struct = {type:2,time:time,data_source:infile, $
                  nparticle:long(nw),nclass:1,r:r[ww],z:w[ww],$
                  energy:energy[ww],pitch:pitch[ww],class:orbit_class[ww],$
                  weight:weight[ww]}

    return, fbm_struct
    GET_OUT:    
END

