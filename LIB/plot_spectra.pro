PRO plot_spectra,runid,dir=dir,chan=chan,fida=fida,nbi=nbi,halo=halo,intens=intens,ps=ps,pretty=pretty

	if keyword_set(fida) then fida_switch = [1,0,0] else fida_switch=[1,1,1] 
	if keyword_set(nbi) then nbi_switch = [0,1,0] else nbi_switch=[1,1,1]
	if keyword_set(halo) then halo_switch = [0,0,1] else halo_switch=[1,1,1]

	plt=fida_switch+nbi_switch+halo_switch 
	plt-=min(plt)
	if min(plt) eq max(plt) then plt=[1,1,1]

	if not keyword_set(dir) then dir=dialog_pickfile(dir='~/',/directory)
	print,dir
	load_results,runid,results,dir=dir
	
	inputs=results.inputs
	grid=results.grid
	spec=results.spectra
	weights=results.fida_weights
    ichan_wght=weights.ichan
	los=results.los
	fbm=results.fbm
	neutrals=results.neutrals

			
	;; LOS IS IN BEAM COORDINATES WE NEED TO CONVERT THEM TO MACHINE COORS.
	alpha=grid.alpha & beta=grid.beta & origin=grid.origin
	xlos=los.xlos
	ylos=los.ylos
	zlos=los.zlos
	
	x =  cos(alpha)*(cos(beta)*xlos + sin(beta)*zlos) $
    - sin(alpha)*ylos + origin[0]

    y =  sin(alpha)*(cos(beta)*xlos + sin(beta)*zlos) $
    + cos(alpha)*ylos + origin[1]

    z = -sin(beta)*xlos + cos(beta)*zlos + origin[2]
	
	rlos=sqrt(x^2.0 + y^2.0)		

	if not keyword_set(intens) then begin
		;;CALCULATE WEIGHT FUNCTIONS
		if weights.err ne 1 then begin
			if grid.err or fbm.err or los.err or neutrals.err then begin
				print,strupcase('missing files. skipping weight function spectra calculation')
			endif else calc_spectra,grid,fbm,weights,los,neutrals,wlambda,wspec,mean_fbm
		endif 

		xran=[647,665]
		minbrems=min(spec.brems)
		print,minbrems
		yran=[minbrems,1.e19]
        if ichan_wght gt 0 then chan=ichan_wght-1
		if n_elements(chan) ne 0 then begin
			startind=chan & endind=chan 
		endif else begin
			startind=0L & endind=spec.nchan-1
		endelse

		;;LOOP OVER CHANNELS
		for ichan=startind,endind do begin
			if spec.err eq 0 then brems=spec.brems[*,ichan] else brems=0
			if keyword_set(pretty) then begin
				plts=plot([0.],/nodata,xrange=xran,yrange=yran $
					, ytitle='Intensity  [Ph/(s $m^2$ nm sr)]', xtitle='$\lambda$ [nm]' $
					, /ylog,min_value=0.1)
				txt=TEXT(100,435,'$R = $'+strtrim(string(rlos[ichan],format='(f10.2)'),1)+' cm',/device)
				if plt[0] ne 0 then begin
					if spec.err eq 0 then begin 
						plt1=plot(spec.lambda,spec.fida[*,ichan]+brems,'r',/overplot,name='Fida',min_value=0.1)
						if n_elements(plt1) ne 0 then targets=[plt1]
					endif
					if n_elements(wspec) ne 0 then begin
						plt2=plot(wlambda,wspec[*,ichan]+brems,'k',thick=2,linestyle=2,/overplot,name='Weight',min_value=0.1)
						if n_elements(plt2) ne 0 and n_elements(plt1) ne 0 then targets=[targets,plt2] else targets=[plt2] 
					endif
				endif
			
				if plt[1] ne 0 and spec.err eq 0 then begin
					plt3=plot(spec.lambda,spec.full[*,ichan]+brems,'c',/overplot,name='Full',min_value=0.1)
					plt4=plot(spec.lambda,spec.half[*,ichan]+brems,'m',/overplot,name='Half',min_value=0.1)
					plt5=plot(spec.lambda,spec.third[*,ichan]+brems,'g',/overplot,name='Third',min_value=0.1)
					if n_elements(targets) eq 0 then targets=[plt3,plt4,plt5] else targets=[targets,plt3,plt4,plt5] 
				endif
			
				if plt[2] ne 0 and spec.err eq 0 then begin
					plt6=plot(spec.lambda,spec.halo[*,ichan]+brems,'b',/overplot,name='Halo',min_value=0.1)
					if n_elements(targets) eq 0 then targets=[plt6] else targets=[targets,plt6] 				
				endif
				if n_elements(targets) ne 0 then leg=legend(target=targets,/device,position=[464,438])
			endif else begin
				!p.multi=0
				loadct,39,/silent
				plot,[0.],/nodata,xrange=xran,yrange=yran $
					, ytitle='Intensity  [Ph/(s m^2 nm sr)]', xtitle='lambda [nm]' $
					, /ylog,min_value=0.1
				xyouts,100,435,'R = '+strtrim(string(rlos[ichan],format='(f10.2)'),1)+' cm',/device
				if plt[0] ne 0 then begin
					if spec.err eq 0 then begin 
						oplot,spec.lambda,spec.fida[*,ichan]+brems,color=253,min_value=0.1
					endif
					if n_elements(wspec) ne 0 then begin
						oplot,wlambda,wspec[*,ichan]+brems,thick=2,linestyle=2,min_value=0.1
					endif
				endif
			
				if plt[1] ne 0 and spec.err eq 0 then begin
					oplot,spec.lambda,spec.full[*,ichan]+brems,color=100,min_value=0.1
					oplot,spec.lambda,spec.half[*,ichan]+brems,color=150,min_value=0.1
					oplot,spec.lambda,spec.third[*,ichan]+brems,color=200,min_value=0.1
				endif
			
				if plt[2] ne 0 and spec.err eq 0 then begin
					oplot,spec.lambda,spec.halo[*,ichan]+brems,color=70,min_value=0.1
				endif
			endelse
			wait,2
		endfor
		if keyword_set(ps) and keyword_set(pretty) then plts.Save,inputs.fidasim_runid+"_spectra.pdf",BORDER=10
	endif else begin
		if inputs.err eq 0 and spec.err eq 0 and spec.err eq 0 then begin
			intensity=dblarr(los.nchan)
			for ichan=0,los.nchan-1 do begin
				intensity[ichan]=total(plt[1]*(spec.full[*,ichan]+spec.half[*,ichan]+spec.third[*,ichan])+$
						  	plt[2]*spec.halo[*,ichan]+plt[0]*spec.fida[*,ichan]+spec.brems[*,ichan])
			endfor
			if keyword_set(pretty) then begin
				plt=plot(rlos,intensity,$
				 	title='Intensity vs. Major Radius ',xtitle='R [cm]',ytitle='Intensity [Ph/(s $m^2$ nm sr)]')
				if keyword_set(ps) then begin
					type=size(ps,/type)
					if type eq 7 then dir=inputs.fidasim_runid+"_intensity.pdf" else dir=inputs.fidasim_runid+"_intensity.pdf" 
					plt.Save,dir,border=10
				endif
			endif else begin
				!p.multi=0
				window,0 & wset,0
				plot,rlos,intensity,title='Intensity vs. Major Radius',xtitle='R [cm]',ytitle='Intensity [Ph/(s m^2 nm sr)]',psym=2
			endelse
		endif else print,'MISSING FILES'
	endelse
	GET_OUT:
END
