PRO d3d_plots,inputs,grid,nbi,fida,equil,nbgeom,plasma

	g=equil.g

    bfieldu=dblarr(grid.nx,grid.ny,grid.nz)
    bfieldv=bfieldu & bfieldw=bfieldu
	denf=bfieldu 
	vrotu=bfieldu & vrotv=bfieldu & vrotw=bfieldu
    uvals=bfieldu & vvals=bfieldu
	
	for i=0L,grid.nx-1 do for j=0L,grid.ny-1 do for k=0L,grid.nz-1 do begin
		l=i+grid.nx*j+grid.nx*grid.ny*k
		;;FAST ION DENSITY
		denf[i,j,k]=plasma.denf[l]

		;;MAGNETIC FIELD
		bfieldu[i,j,k]=equil.b[0,l]
		bfieldv[i,j,k]=equil.b[1,l]
		bfieldw[i,j,k]=equil.b[2,l]

		;;PLASM ROTATION
		vrotu[i,j,k]=plasma.vrot_uvw[0,l]
		vrotv[i,j,k]=plasma.vrot_uvw[1,l]
		vrotw[i,j,k]=plasma.vrot_uvw[2,1]
		uvals[i,j,k]=grid.uc[l]
		vvals[i,j,k]=grid.vc[l]
    endfor
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

	contour,denf[*,*,ind],uvals[*,*,ind],vvals[*,*,ind],/fill,nlevels=60 $
		,xrange=x_range,yrange=y_range,title='PLANE VIEW',xtitle='U [cm]',ytitle='V [cm]',color=0
	oplot,grid.uc,grid.vc,psym=3,color=0


	for i=0,fida.nchan-1 do $
		oplot,fida.xlens[i]+[0,2*(fida.xlos[i]-fida.xlens[i])],$
        fida.ylens[i]+[0,2*(fida.ylos[i]-fida.ylens[i])],$
        color=50

	src=nbi.xyz_src
	pos=(nbi.xyz_pos-nbi.xyz_src)*1000+nbi.xyz_src
	oplot,[src[0],pos[0]],[src[1],pos[1]],thick=2,color=230

	w=where(g.bdry[0,*] gt 0.)
	rmin=100.*min(g.bdry[0,w]) & rmax=100.*max(g.bdry[0,w])
	rmaxis=100.*g.rmaxis
	phi=2.*!pi*findgen(501)/500.
	oplot,rmin*cos(phi),rmin*sin(phi),color=150
	oplot,rmaxis*cos(phi),rmaxis*sin(phi),color=150,linestyle=2
	oplot,rmax*cos(phi),rmax*sin(phi),color=150

	;----------------------------------------------
	;;PLOT CROSS SECTION BEAM AND CHORDS 
	window,1 & wset,1
	plot,[0],[0],/nodata,xrange=[rmin,rmax], $
            yrange=100.*[min(g.lim[1,*]),max(g.lim[1,*])],$
			color=0,background=255,title='ELEVATION',xtitle='R [cm]',ytitle='Z [cm]'

	oplot,grid.r_grid,grid.wc,psym=3,color=0  

	; Lines of sight
	for i=0,fida.nchan-1 do begin
		if fida.zlos[i] ne fida.zlens[i] then begin
			z=(fida.zlos[i]-fida.zlens[i])*findgen(201)/100.+fida.zlens[i]
			x=(fida.xlos[i]-fida.xlens[i])*(z-fida.zlens[i])/ $
			  (fida.zlos[i]-fida.zlens[i]) + fida.xlens[i]
			y=(fida.ylos[i]-fida.ylens[i])*(z-fida.zlens[i])/ $
			  (fida.zlos[i]-fida.zlens[i]) + fida.ylens[i]
			oplot,sqrt(x^2+y^2),z,color=50
		endif else begin 
    		y=(fida.ylos[i]-fida.ylens[i])*findgen(201)/100.+fida.ylens[i]
    		x=(fida.xlos[i]-fida.xlens[i])*(y-fida.ylens[i])/ $
      		  (fida.ylos[i]-fida.ylens[i]) + fida.ylens[i]
		    oplot,sqrt(x^2+y^2),replicate(fida.zlens[i],201),color=50
		endelse
	endfor

	; Equilibrium
	oplot,100.*g.bdry[0,*],100.*g.bdry[1,*],color=150
	oplot,100.*g.lim[0,*],100.*g.lim[1,*],color=0

	window,2 & wset,2
	!p.multi=[0,2,2,0,1]
	plot,equil.rho_grid,plasma.te,psym=3,color=0,background=255,title='Te and Ti',xtitle='rho',ytitle='keV'
 	oplot,equil.rho_grid,plasma.ti,psym=3,color=50
 	plot,equil.rho_grid,plasma.dene,psym=3,color=0,background=255,title='electron density',xtitle='rho',ytitle='cm^-3'
  	plot,equil.rho_grid,plasma.zeff,psym=3,color=0,background=255,title='zeff',xtitle='rho',ytitle='zeff'
	
  	vrotx=transpose(plasma.vrot[0,*])
  	vroty=transpose(plasma.vrot[1,*])
  	vrotz=transpose(plasma.vrot[2,*])
  	plot,equil.rho_grid,sqrt(vrotx^2.0 + vroty^2.0 + vrotz^2.0),psym=3,color=0,background=255,title='vtor',xtitle='rho',ytitle='cm/s'
	
	;;PLOT VECTOR FIELDS
	nnx=long(grid.nx/2) & nny=long(grid.ny/2)
	indx=2*lindgen(nnx) & indy=2*lindgen(nny)
    bu=reform(bfieldu[indx,indy,ind],nnx*nny)
    bv=reform(bfieldv[indx,indy,ind],nnx*nny)
	vu=reform(vrotu[indx,indy,ind],nnx*nny)
	vv=reform(vrotv[indx,indy,ind],nnx*nny)
    uvals1=reform(uvals[indx,indy,ind],nnx*nny)
    vvals1=reform(vvals[indx,indy,ind],nnx*nny)
    bfield=vector(bu,bv,uvals1,vvals1,auto_color=1,rgb_table=39,head_angle=20,$
				  title='Magnetic Field',xtitle='U [cm]',ytitle='V [cm]')

	vrotfield=vector(vu,vv,uvals1,vvals1,auto_color=1,rgb_table=39,head_angle=20,$
				     title='Plasma Rotation',xtitle='U [cm]',ytitle='V [cm]')
  	!p.multi=0

END
