PRO plot_neutrals,path=path,halo=halo,beam=beam

    if not keyword_set(path) then path=dialog_pickfile(path='~/FIDASIM/RESULTS/',/directory)
    print,path

	load_results,path,results

	grid=results.grid
	plasma=results.plasma
    neutrals=results.neutrals
	if grid.err eq 1 then begin
		print,'MISSING GRID FILE'
		goto,GET_OUT
	endif
	if neutrals.err eq 1 then begin
		print,'MISSING NEUTRALS FILE'
		goto,GET_OUT
	endif

	halodens=total(neutrals.halodens[*,*,*,*],4)
	fdens=total(neutrals.fdens[*,*,*,*],4)
	hdens=total(neutrals.hdens[*,*,*,*],4)
	tdens=total(neutrals.tdens[*,*,*,*],4)
	nbidens=fdens+hdens+tdens

	if keyword_set(beam) then halodens*=0
	if keyword_set(halo) then nbidens*=0
	dens=halodens+nbidens

    xmin=min(grid.u_grid) & ymin=min(grid.v_grid) & zmin=min(grid.w_grid)
    xmax=max(grid.u_grid) & ymax=max(grid.v_grid) & zmax=max(grid.w_grid)
    if xmin lt 0 then xmin1=1.2*xmin else xmin1=0.8*xmin
    if xmax lt 0 then xmax1=0.8*xmax else xmax1=1.2*xmax
    if ymin lt 0 then ymin1=1.2*ymin else ymin1=0.8*ymin
    if ymax lt 0 then ymax1=0.8*ymax else ymax1=1.2*ymax
    if zmin lt 0 then zmin1=1.2*zmin else zmin1=0.8*zmin
    if zmax lt 0 then zmax1=0.8*zmax else zmax1=1.2*zmax

    x_range=[xmin1,xmax1] & y_range=[ymin1,ymax1] & z_range=[zmin1,zmax1]
	r_range=[0.8*min(grid.r_grid),1.2*max(grid.r_grid)]

	loadct,39,/silent
	!p.multi=[0,2,1,0,1]
	window,0 & wset,0
	ind=long(grid.nz/2.)
	contour,dens[*,*,ind],grid.u_grid[*,*,ind],grid.v_grid[*,*,ind],/fill,nlevels=60 $
		   ,title='PLANE VIEW',xtitle='U [cm]',ytitle='V [cm]'$
		   ,color=0,background=255
    oplot,grid.u_grid,grid.v_grid,psym=3,color=0

	ind=long(grid.ny/2.)
	contour,reform(dens[*,ind,*]),abs(reform(grid.x_grid[*,ind,*])),reform(grid.z_grid[*,ind,*]),/fill,nlevels=60 $
           ,title='BEAM VIEW',xtitle='X [cm]',ytitle='Z [cm]'$
		   ,color=0,background=255
    oplot,abs(grid.x_grid[*,ind,*]),grid.z_grid[*,ind,*],psym=3,color=0

	GET_OUT:
END
