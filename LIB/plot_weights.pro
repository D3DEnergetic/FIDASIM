PRO plot_weights,path=path,chan=chan,prod=prod,fida=fida,npa=npa

	if not keyword_set(path) then path=dialog_pickfile(path='~/FIDASIM/RESULTS/',/directory)
	print,path
	load_results,path,results
	
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

	if keyword_set(chan) then begin
		start=chan
		fin=chan
	endif else begin
		start=0
		fin=los.nchan-1
	endelse

	wght=fltarr(weights.nen,weights.npitch)
	loadct,39,/silent
	for ichan=start,fin do begin
		if keyword_set(fida) then begin
			for ii=0,weights.nwav-1 do begin
				if keyword_set(prod) then begin
					wght=mean_fbm[*,*,ichan]*weights.weight_tot[ii,*,*,ichan]
				endif else wght=weights.weight_tot[ii,*,*,ichan]
				contour,wght,weights.energyarr,weights.pitcharr,nlevels=30,/fill
				wait,1
			endfor
		endif
		if keyword_set(npa)then begin
			if keyword_set(prod) then begin
				wght=mean_fbm[*,*,ichan]*weights.weight_tot[*,*,ichan]
			endif else wght=weights.weight_tot[*,*,ichan]
			contour,wght,weights.energyarr,weights.pitcharr,nlevels=30,/fill
			wait,1
		endif
	endfor		
	GET_OUT:
END
