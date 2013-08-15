pro oplot_fidasim_grid,coords,thick=thick,pol=pol,tor=tor,cm=cm,color=color
  ;; -------------------------------------------------
  ;; -- ROUTINE TO PLOT THE FIDASIM SIMULATION GRID --
  ;; -------------------------------------------------
  if keyword_set(cm) then unit=1. else unit=0.01
  if not keyword_set(color) then tvlct,180,180,180,0
  if keyword_set(tor) then begin
     for ii=0,coords.nx-1 do $
        oplot,[coords.xx[ii],coords.xx[ii]]*unit $
               ,[coords.yy[0],coords.yy[coords.ny-1]+coords.dy]*unit $
               ,color=color,thick=thick
     oplot,[coords.xx[coords.nx-1]+coords.dx $
            ,coords.xx[coords.nx-1]+coords.dx]*unit $
           ,[coords.yy[0],coords.yy[coords.ny-1]+coords.dy ]*unit $
           ,color=color,thick=thick
     
     for ii=0,coords.ny-1 do $
        oplot,[coords.xx[0],coords.xx[coords.nx-1]+coords.dx]*unit $
              ,[coords.yy[ii],coords.yy[ii]]*unit $
              ,color=color,thick=thick
     oplot,[coords.xx[0],coords.xx[coords.nx-1]+coords.dx]*unit $
           ,[coords.yy[coords.ny-1]+coords.dy $
             ,coords.yy[coords.ny-1]+coords.dy]*unit $
           ,color=color,thick=thick
  endif

 if keyword_set(pol) then begin
     rrmin=min(coords.rrc_grid)
     rrmax=max(coords.rrc_grid)
     nr=round(sqrt(coords.nx^2+coords.ny^2))
     rrc=dindgen(nr)/(nr-1.)*(rrmax-rrmin)+rrmin
     dr=rrc[1]-rrc[0]
     rr=rrc-0.5*dr
     for ii=0,nr-1 do $
        oplot,[rr[ii],rr[ii]]*unit $
              ,[coords.zz[0],coords.zz[coords.nz-1]+coords.dz]*unit $
              ,thick=thick
     oplot,[rr[nr-1]+dr,rr[nr-1]+dr]*unit $
           ,[coords.zz[0],coords.zz[coords.nz-1]+coords.dz]*unit $
           ,thick=thick
     for ii=0,coords.nz-1 do $
        oplot,[rr[0],rr[nr-1]+dr]*unit $
              ,[coords.zz[ii],coords.zz[ii]]*unit $
              ,thick=thick
     oplot,[rr[0],rr[nr-1]+dr]*unit $
           ,[coords.zz[coords.nz-1]+coords.dz $
             ,coords.zz[coords.nz-1]+coords.dz]*unit $
           ,thick=thick  
  endif
 if not keyword_set(color) then tvlct,0,0,0,0

end
