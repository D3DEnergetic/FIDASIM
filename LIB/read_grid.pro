PRO read_grid,file,grid

	if file_test(file) then begin
		;;READ ARRAY SIZES
		nx=0L & ny=0L & nz=0L
	    openr,lun,file,/get_lun
	    readu,lun, nx
	    readu,lun, ny
	    readu,lun, nz

		;;DEFINE ARRAYS
		u_grid=dblarr(nx,ny,nz)
		v_grid=dblarr(nx,ny,nz)
		w_grid=dblarr(nx,ny,nz)
		r_grid=dblarr(nx,ny,nz)
		phi_grid=dblarr(nx,ny,nz)
		x_grid=dblarr(nx,ny,nz)
		y_grid=dblarr(nx,ny,nz)
		z_grid=dblarr(nx,ny,nz)

		;;READ IN ARRAYS
	    readu,lun, u_grid
	    readu,lun, v_grid
	    readu,lun, w_grid
	    readu,lun, r_grid
	    readu,lun, phi_grid
	    readu,lun, x_grid
	    readu,lun, y_grid
	    readu,lun, z_grid
	    close,lun
	    free_lun,lun
		
		grid={nx:nx,ny:ny,nz:nz,u_grid:u_grid,v_grid:v_grid,w_grid:w_grid,$
			  r_grid:r_grid,phi_grid:phi_grid,$
			  x_grid:x_grid,y_grid:y_grid,z_grid:z_grid,err:0}
	endif else begin
		grid={err:1}
	endelse

END
