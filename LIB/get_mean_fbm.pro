FUNCTION get_mean_fbm,runid,dir=dir,chan=chan

    if not keyword_set(dir) then dir = '.'
	ifile=dir+'/'+runid+'_inputs.cdf'
	nfile=dir+'/'+runid+'_neutrals.cdf'

    inputs = read_ncdf(ifile,vars=['chan_id','los_wght','r_grid','w_grid',$
                                   'FBM','FBM_energy','FBM_pitch','FBM_r2d','FBM_z2d'])

    nen=n_elements(inputs.FBM_energy)
    np =n_elements(inputs.FBM_pitch)
    dE = abs(inputs.FBM_energy[1]-inputs.FBM_energy[0])
    dP = abs(inputs.FBM_pitch[1]-inputs.FBM_pitch[0])

    if n_elements(chan) eq 0 then begin
        nchan=n_elements(inputs.chan_id)
        chan=indgen(nchan)
    endif else nchan=n_elements(chan)

    mean_fbm = dblarr(nen,np,nchan)

    neut=read_ncdf(nfile)
    ndens = total(neut.fdens+neut.hdens+neut.tdens+neut.halodens,4)
    fdens=dblarr(nchan)
    for i=0,nchan-1 do begin
        los = inputs.los_wght[*,*,*,chan[i]]
        wlos = where(los gt 0,nwlos)            
        for j=0,nwlos-1 do begin
            tmp=MIN((inputs.r_grid[wlos[j]]-inputs.fbm_r2d)^2+(inputs.w_grid[wlos[j]]-inputs.fbm_z2d)^2,wh)
            mean_fbm[*,*,i]+=inputs.fbm[*,*,wh]*ndens[wlos[j]]*los[wlos[j]]
        endfor
        mean_fbm[*,*,i] = mean_fbm[*,*,i]/total(ndens[wlos]*los[wlos])
        fdens[i] = total(mean_fbm[*,*,i])*dE*dP
    endfor
    
    return,{fbm:mean_fbm,energy:inputs.FBM_energy,pitch:inputs.FBM_pitch,fdens:fdens,chan_id:inputs.chan_id[chan]}

END
