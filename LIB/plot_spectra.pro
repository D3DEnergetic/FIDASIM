PRO plot_spectra,path=path

	if not keyword_set(path) then path=dialog_pickfile(path='~/FIDASIM/RESULTS/',/directory)
	print,path
	load_results,path,results
	
	grid=results.grid
	fida=results.fida
	nbi_halo=results.nbi_halo
	weights=results.weights
	los=results.los
	fbm=results.fbm
	neutrals=results.neutrals

	if weights.err ne 1 then begin
		if grid.err or fbm.err or los.err or neutrals.err then begin
			print,strupcase('missing files. skipping weight function spectra calculation')
		endif else calc_spectra,grid,fbm,weights,los,neutrals,wlambda,wspec,mean_fbm
	endif 

	if fida.err eq 1 then begin
		print,'MISSING FIDA FILE'
		goto,GET_OUT
	endif
	if nbi_halo.err eq 1 then begin
		print,'MISSING NBI_HALO FILE'
		goto,GET_OUT
	endif

	!p.multi=0
	xran=[647,665]
	loadct,39,/silent
	for ichan=0L,los.nchan-1 do begin
		brems=nbi_halo.brems[*,ichan]
		yran=[1.e15,1.e19]

		plot, [0.],/nodata,xran=xran,yran=yran,/xsty $
		 , ytit='Intensity [Ph/(s m^2 nm sr)]', xtit='lambda [nm]' $
		 ,xthick=linthick,ythick=linthick,/ylog,color=0,background=255
		oplot,fida.lambda,fida.spectra[*,ichan]+brems ,color=254
		oplot,nbi_halo.lambda,nbi_halo.halo[*,ichan]+brems ,color=150
		oplot,nbi_halo.lambda,nbi_halo.full[*,ichan]+brems,color=210
		oplot,nbi_halo.lambda,nbi_halo.half[*,ichan]+brems ,color=200
		oplot,nbi_halo.lambda,nbi_halo.third[*,ichan]+brems ,color=190
		if n_elements(wspec) ne 0 then begin
			oplot,wlambda,wspec[*,ichan]+brems,color=50
		endif
		wait,2
	endfor

	GET_OUT:
END
