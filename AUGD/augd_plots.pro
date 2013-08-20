PRO augd_plots,inputs,grid,nbi,chords,equil,nbgeom,plasma

	ind=long(grid.nz/2.0)
	!p.multi=0
	;;PLOTTING PLANE VIEW BEAMS AND CHORDS
	window,0 & wset,0
	loadct,39,/silent

	;;GET PROPER RANGES
	xmin=min(grid.u) & ymin=min(grid.v) & zmin=min(grid.w)
	xmax=max(grid.u) & ymax=max(grid.v) & zmax=max(grid.w)
	if xmin lt 0 then xmin1=1.2*xmin else xmin1=0.8*xmin
	if xmax lt 0 then xmax1=0.8*xmax else xmax1=1.2*xmax
	if ymin lt 0 then ymin1=1.2*ymin else ymin1=0.8*ymin
	if ymax lt 0 then ymax1=0.8*ymax else ymax1=1.2*ymax
	if zmin lt 0 then zmin1=1.2*zmin else zmin1=0.8*zmin
	if zmax lt 0 then zmax1=0.8*zmax else zmax1=1.2*zmax

	x_range=[xmin1,xmax1] & y_range=[ymin1,ymax1] & z_range=[zmin1,zmax1]

	contour,plasma.denf[*,*,ind],grid.u_grid[*,*,ind],grid.v_grid[*,*,ind],/fill,nlevels=60 $
		,xrange=x_range,yrange=y_range,title='PLANE VIEW',xtitle='U [cm]',ytitle='V [cm]',color=0
	oplot,grid.u_grid,grid.v_grid,psym=3,color=0


	for i=0,chords.nchan-1 do $
		oplot,chords.xlens[i]+[0,2*(chords.xlos[i]-chords.xlens[i])],$
        chords.ylens[i]+[0,2*(chords.ylos[i]-chords.ylens[i])],$
        color=50

	src=nbi.xyz_src
	pos=(nbi.xyz_pos-nbi.xyz_src)*1000+nbi.xyz_src
	oplot,[src[0],pos[0]],[src[1],pos[1]],thick=2,color=230

	;----------------------------------------------
	;;PLOT CROSS SECTION BEAM AND CHORDS 
	window,1 & wset,1
	plot,[0],[0],/nodata,xrange=[0,300.], $
            yrange=100.*[-1,1],$
			color=0,background=255,title='ELEVATION',xtitle='R [cm]',ytitle='Z [cm]'

	oplot,grid.r_grid,grid.w_grid,psym=3,color=0  

	; Lines of sight
	for i=0,chords.nchan-1 do begin
		if chords.zlos[i] ne chords.zlens[i] then begin
			z=(chords.zlos[i]-chords.zlens[i])*findgen(201)/100.+chords.zlens[i]
			x=(chords.xlos[i]-chords.xlens[i])*(z-chords.zlens[i])/ $
			  (chords.zlos[i]-chords.zlens[i]) + chords.xlens[i]
			y=(chords.ylos[i]-chords.ylens[i])*(z-chords.zlens[i])/ $
			  (chords.zlos[i]-chords.zlens[i]) + chords.ylens[i]
			oplot,sqrt(x^2+y^2),z,color=50
		endif else begin 
    		y=(chords.ylos[i]-chords.ylens[i])*findgen(201)/100.+chords.ylens[i]
    		x=(chords.xlos[i]-chords.xlens[i])*(y-chords.ylens[i])/ $
      		  (chords.ylos[i]-chords.ylens[i]) + chords.ylens[i]
		    oplot,sqrt(x^2+y^2),replicate(chords.zlens[i],201),color=50
		endelse
	endfor

	window,2 & wset,2
	!p.multi=[0,2,2,0,1]
	plot,equil.rho_grid,plasma.te,psym=3,color=0,background=255,title='Te (black) Ti (blue)',xtitle='rho',ytitle='keV'
 	oplot,equil.rho_grid,plasma.ti,psym=3,color=50
 	plot,equil.rho_grid,plasma.dene,psym=3,color=0,background=255,title='n_e (black) n_imp (blue) n_ion (green) n_f (red)',xtitle='rho',ytitle='cm^-3'
	oplot,equil.rho_grid,plasma.deni,psym=3,color=50
	oplot,equil.rho_grid,plasma.denp,psym=3,color=150
	oplot,equil.rho_grid,plasma.denf,psym=3,color=250
  	plot,equil.rho_grid,plasma.zeff,psym=3,color=0,background=255,title='zeff',xtitle='rho',ytitle='zeff'
  	plot,equil.rho_grid,plasma.vtor,psym=3,color=0,background=255,title='vtor',xtitle='rho',ytitle='cm/s'
	
	;;PLOT VECTOR FIELDS
	nnx=long(grid.nx/2) & nny=long(grid.ny/2)
	indx=2*lindgen(nnx) & indy=2*lindgen(nny)
    bu=reform(plasma.bu[indx,indy,ind],nnx*nny)
    bv=reform(plasma.bv[indx,indy,ind],nnx*nny)
	vu=reform(plasma.vrotu[indx,indy,ind],nnx*nny)
	vv=reform(plasma.vrotv[indx,indy,ind],nnx*nny)
    uvals1=reform(grid.u_grid[indx,indy,ind],nnx*nny)
    vvals1=reform(grid.v_grid[indx,indy,ind],nnx*nny)
    bfield=vector(bu,bv,uvals1,vvals1,auto_color=1,rgb_table=39,head_angle=20,$
				  title='Magnetic Field',xtitle='U [cm]',ytitle='V [cm]')

	vrotfield=vector(vu,vv,uvals1,vvals1,auto_color=1,rgb_table=39,head_angle=20,$
				     title='Plasma Rotation',xtitle='U [cm]',ytitle='V [cm]')
  	!p.multi=0

END
