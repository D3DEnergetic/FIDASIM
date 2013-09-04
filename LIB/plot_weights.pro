PRO plot_weights,path=path,chan=chan,prod=prod,fida=fida,npa=npa

	if not keyword_set(path) then path=dialog_pickfile(path='~/FIDASIM/RESULTS/',/directory)
	print,path
	load_results,path,results
	
	inputs=results.inputs
	grid=results.grid
	nbi_halo=results.nbi_halo
	fida_wght=results.fida_weights
	npa_wght=results.npa_weights
	los=results.los
	fbm=results.fbm
	neutrals=results.neutrals


	if keyword_set(fida) then weights=fida_wght
	if keyword_set(npa) then weights=npa_wght

 	calc_mean_fbm,grid,fbm,weights,los,neutrals,mean_fbm

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
                contour,mean_fbm[*,*,ichan],weights.energyarr,weights.pitcharr,nlevels=10
				contour,wght,weights.energyarr,weights.pitcharr,nlevels=30,/fill,/overplot
				wait,1
			endfor
		endif
		if keyword_set(npa)then begin
			if keyword_set(prod) then begin
				wght=mean_fbm[*,*,ichan]*weights.weight_tot[*,*,ichan]
			endif else wght=weights.weight_tot[*,*,ichan]
            contour,mean_fbm[*,*,ichan],weights.energyarr,weights.pitcharr,nlevels=10
			contour,wght,weights.energyarr,weights.pitcharr,nlevels=30,/fill,/overplot
			wait,1
		endif
	endfor		
	GET_OUT:
END
