PRO plot_weights,runid,dir=dir,chan=chan,prod=prod,fida=fida,npa=npa,currentmode=currentmode

	if not keyword_set(dir) then dir=dialog_pickfile(dir='~/FIDASIM/RESULTS/',/directory)
	print,dir
	load_results,runid,results,dir=dir
	
	inputs=results.inputs
	grid=results.grid
	fida_wght=results.fida_weights
	npa_wght=results.npa_weights
	los=results.los
	fbm=results.fbm
	neutrals=results.neutrals

    if not keyword_set(fida) and not keyword_set(npa) then fida=1
	if keyword_set(fida) then begin
		weights=fida_wght
		e_level=2
	endif
	if keyword_set(npa) then begin
		weights=npa_wght
		e_level=indgen(n_elements(neutrals.fdens[0,0,0,*]))
	endif
 	calc_mean_fbm,grid,fbm,weights,los,neutrals,mean_fbm,elevel=e_level

    if fida_wght.ichan gt 0 and not keyword_set(npa) then chan=fida_wght.ichan-1
	if n_elements(chan) ne 0 then begin
		start=chan
		fin=chan
	endif else begin
		start=0
		fin=los.nchan-1
	endelse

	wght=fltarr(n_elements(weights.energy),n_elements(weights.pitch))
	loadct,39,/silent
	for ichan=start,fin do begin
		if keyword_set(fida) then begin
			for ii=0,n_elements(weights.lambda)-1 do begin
				if keyword_set(prod) then begin
					wght=mean_fbm[*,*,ichan]*weights.wfunct[ii,*,*,ichan]
				endif else wght=weights.wfunct[ii,*,*,ichan]
				contour,wght,weights.energy,weights.pitch,nlevels=30,/fill
				wait,1
			endfor
		endif
		if keyword_set(npa)then begin
			if keyword_set(prod) then begin
				wght=mean_fbm[*,*,ichan]*weights.wfunct[*,*,ichan]
			endif else wght=weights.wfunct[*,*,ichan]
            if keyword_set(currentmode) then begin
                for i=0,n_elements(weights.energy)-1 do begin
                    if weights.energy[i] gt 25 then wght[i,*]=wght[i,*]*(weights.energy[i]-25.) else wght[i,*]=0.0
                endfor
            endif
			contour,wght,weights.energy,weights.pitch,nlevels=30,/fill
			wait,1
		endif
	endfor		
	GET_OUT:
END
