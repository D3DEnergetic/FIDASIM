PRO plot_spectra,path=path,fida=fida,nbi=nbi,halo=halo

	if keyword_set(fida) then fida_switch = [1,0,0] else fida_switch=[1,1,1] 
	if keyword_set(nbi) then nbi_switch = [0,1,0] else nbi_switch=[1,1,1]
	if keyword_set(halo) then halo_switch = [0,0,1] else halo_switch=[1,1,1]

	plt=fida_switch+nbi_switch+halo_switch 
	plt-=min(plt)
	if min(plt) eq max(plt) then plt=[1,1,1]
 
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
	loadct,39,/silent

	xran=[647,665]
	minbrems=min(nbi_halo.brems)
	yran=[minbrems,1.e19]

	;;LOOP OVER CHANNELS
	for ichan=0L,los.nchan-1 do begin
		brems=nbi_halo.brems[*,ichan]

		plot, [0.],/nodata,xran=xran,yran=yran,/xsty $
			, ytit='Intensity [Ph/(s m^2 nm sr)]', xtit='lambda [nm]' $
			, xthick=linthick,ythick=linthick,/ylog,color=0,background=255

		if plt[0] ne 0 then begin
			oplot,fida.lambda,fida.spectra[*,ichan]+brems ,color=254
			if n_elements(wspec) ne 0 then begin
				oplot,wlambda,wspec[*,ichan]+brems,color=0,thick=2,linestyle=2
			endif
		endif

		if plt[1] ne 0 then begin
			oplot,nbi_halo.lambda,nbi_halo.full[*,ichan]+brems,color=100
			oplot,nbi_halo.lambda,nbi_halo.half[*,ichan]+brems ,color=150
			oplot,nbi_halo.lambda,nbi_halo.third[*,ichan]+brems ,color=200
		endif

		if plt[2] ne 0 then begin
			oplot,nbi_halo.lambda,nbi_halo.halo[*,ichan]+brems,color=70
		endif
		wait,2
	endfor

	GET_OUT:
END
