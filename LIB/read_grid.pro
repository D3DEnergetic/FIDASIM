PRO read_grid,file,grid,save=save

	if file_test(file) then begin
	    ncid=ncdf_open(file,/nowrite)
	    ncdf_varget,ncid,'Nx',nx
	    ncdf_varget,ncid,'Ny',ny
	    ncdf_varget,ncid,'Nz',nz
		ncdf_varget,ncid,'u_grid',u_grid
		ncdf_varget,ncid,'v_grid',v_grid
		ncdf_varget,ncid,'w_grid',w_grid
		ncdf_varget,ncid,'x_grid',x_grid
		ncdf_varget,ncid,'y_grid',y_grid
		ncdf_varget,ncid,'z_grid',z_grid
		ncdf_varget,ncid,'r_grid',r_grid
		ncdf_varget,ncid,'phi_grid',phi_grid

	    ncdf_close,ncid

		grid={nx:nx,ny:ny,nz:nz,u_grid:u_grid,v_grid:v_grid,w_grid:w_grid,$
			  r_grid:r_grid,phi_grid:phi_grid,x_grid:x_grid,y_grid:y_grid,z_grid:z_grid,err:0}
		if keyword_set(save) then save,grid,filename='grid.sav'
	endif else begin
		grid={err:1}
	endelse

END
