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
	!p.multi=[0,3,1,0,1]
	loadct,39,/silent
	window,0,xsize=1200,ysize=400 & wset,0
	ind=long(grid.nz/2.)
	contour,dens[*,*,ind],grid.u_grid[*,*,ind],grid.v_grid[*,*,ind],/fill,nlevels=60,$
			title='PLANE VIEW',xtitle='U [cm]',ytitle='V [cm]',$
			color=0,background=255
    oplot,grid.u_grid[*,*,ind],grid.v_grid[*,*,ind],psym=3,color=0
	
	ind=long(grid.ny/2.)
	contour,reform(dens[*,ind,*]),reform(grid.x_grid[*,ind,*]),reform(grid.z_grid[*,ind,*]),/fill,nlevels=60 $
           ,title='BEAM VIEW',xtitle='X [cm]',ytitle='Z [cm]'$
		   ,color=0,background=255
    oplot,grid.x_grid[*,ind,*],grid.z_grid[*,ind,*],psym=3,color=0
	
	ind=long(grid.nx/2.)
    contour,dens[ind,*,*],grid.y_grid[ind,*,*],grid.z_grid[ind,*,*],/fill,nlevels=60 $
           ,title='BEAM VIEW',xtitle='Y [cm]',ytitle='Z [cm]'$
           ,color=0,background=255
    oplot,grid.y_grid[ind,*,*],grid.z_grid[ind,*,*],psym=3,color=0


	GET_OUT:
END
