PRO printc,str,f=f,b=b,s=s
    if not keyword_set(f) then f='w' ;;Foreground Color
    if not keyword_set(s) then s='b' ;;Style
    esc=string(27b)
    back=esc+"[0m"

    ;; Style Format Codes
    ;; b = Bright
    ;; d = Dim
    ;; u = Underline
    ;; r = Reverse (switch foreground and background colors)
    ;; h = Hidden

    style = {b:'1',d:'2',u:'4',r:'7',h:'8'}
    sTags = ["b","d","u","r","h"]

    ;; Color Codes
    ;; k = Black
    ;; r = Red
    ;; g = Green
    ;; y = Yellow
    ;; b = Blue
    ;; m = Magenta
    ;; c = Cyan
    ;; w = White

    fColors = { k:'30',r:'31',g:'32',y:'33',$
                b:'34',m:'35',c:'36',w:'37'}
    fTags = ["k","r","g","y","b","m","c","w"]

    sIndex = where(s eq sTags)
    fIndex = where(f eq fTags)
    if sIndex eq -1 then sIndex=1
    if fIndex eq -1 then fIndex=7

    print,esc+"["+style.(sIndex)+";"+fColors.(fIndex)+"m"+str+back

END

PRO rotate_uvw,uvw,Arot,Brot,Crot,updown,xyz    
    ;;rotate uvz by alpa alround z axis
    if updown lt 0 then qrz=MATRIX_MULTIPLY(Arot,uvw)
    if updown ge 0 then qrz=MATRIX_MULTIPLY(Brot,uvw)
    ;; rotate qrz_ray by phi_box on xyz_coordinates
    xyz=MATRIX_MULTIPLY(Crot,qrz)
END

PRO rotate_points,x0,y0,z0,Arot,Brot,Crot,xp,yp,zp
    if size(x0,/n_dimensions) ne 1 then begin
        x=reform(x0,n_elements(x0))
        y=reform(y0,n_elements(y0))
        z=reform(z0,n_elements(z0))
    endif else begin
        x=x0 & y=y0 & z=z0
    endelse
    xyz=transpose([[x],[y],[z]])
    rotate_uvw,xyz,Arot,Brot,Crot,1,xyz_p
    xp=transpose(xyz_p[0,*])
    yp=transpose(xyz_p[1,*])
    zp=transpose(xyz_p[2,*])
    if size(x0,/n_dimensions) ne 1 then begin
        xp=reform(xp,size(x0,/dimensions))
        yp=reform(yp,size(y0,/dimensions))
        zp=reform(zp,size(z0,/dimensions))
    endif
END

PRO make_rot_mat,ALPHA,BETA,Arot,Brot,Crot
    zero=0.d0
    one=1.d0
    ;;transformation matrix to rotate on NBI box axis by ALPHA
    Arot=dblarr(3,3)
    Arot[0,0]= cos(BETA)   & Arot[0,1]= zero    & Arot[0,2]= sin(BETA)
    Arot[1,0]= zero        & Arot[1,1]= one     & Arot[1,2]= zero
    Arot[2,0]=-sin(BETA)   & Arot[2,1]= zero    & Arot[2,2]= cos(BETA)
    ;; transformation matrix to rotate in vertical direction by beta
    Brot=dblarr(3,3)
    Brot[0,0]= cos(BETA)   & Brot[0,1]= zero    & Brot[0,2]= sin(BETA)
    Brot[1,0]= zero        & Brot[1,1]= one     & Brot[1,2]= zero
    Brot[2,0]=-sin(BETA)   & Brot[2,1]= zero    & Brot[2,2]= cos(BETA)
    ;; Arot and Brot are exactly the same. I dont know why Ben Gieger had it
    ;; but it will stay for now
    ;; transformation matrix to rotate towards the xy-coordinate system
    Crot=dblarr(3,3)
    Crot[0,0]= cos(ALPHA) & Crot[0,1]=-sin(ALPHA) & Crot[0,2]= zero
    Crot[1,0]= sin(ALPHA) & Crot[1,1]= cos(ALPHA) & Crot[1,2]= zero
    Crot[2,0]= zero       & Crot[2,1]= zero       & Crot[2,2]= one
END

PRO make_grid,inputs,grid,err

    err=1

    if inputs.alpha gt 2*!DPI or inputs.beta gt 2*!DPI then begin
        printc,'ERROR: Angles must be in radians',f='r'
        goto, GET_OUT
    endif
    make_rot_mat,inputs.alpha,inputs.beta,Arot,Brot,Crot

    nx=inputs.nx
    ny=inputs.ny
    nz=inputs.nz
    xmin=inputs.xmin
    xmax=inputs.xmax
    ymin=inputs.ymin
    ymax=inputs.ymax
    zmin=inputs.zmin
    zmax=inputs.zmax

    ;Basic grid points
    dx= (xmax-xmin) / double(nx)
    dy= (ymax-ymin) / double(ny)
    dz= (zmax-zmin) / double(nz)
    ng=long(nx)*long(ny)*long(nz)     ;; nr of cells

    ;;cell borders
    xx=dx*dindgen(nx)+xmin
    yy=dy*dindgen(ny)+ymin
    zz=dz*dindgen(nz)+zmin

    dr=[dx,dy,dz]
    drmin=min(dr)
    dv=dx*dy*dz

    ;;cell centers
    xxc=xx+0.5d0*dx & yyc=yy+0.5d0*dy & zzc=zz+0.5d0*dz
    ;; Put the basic grid into 1D array (useful for libkk routines)
    x=dblarr(ng) & y=dblarr(ng) & z=dblarr(ng)

    for i=0L,nx-1 do for j=0L,ny-1 do for k=0L,nz-1 do begin
        l=i+nx*j+nx*ny*k
        x[l]=xx[i] & y[l]=yy[j] & z[l]=zz[k]
    endfor

    ;; Make the corresponding grid center arrays
    xc=x+0.5d0*dx & yc=y+0.5d0*dy & zc=z+0.5d0*dz
    
    ;;Rotate all grid points to machine coordinates
    rotate_points,x,y,z,Arot,Brot,Crot,u,v,w
    rotate_points,xc,yc,zc,Arot,Brot,Crot,uc,vc,wc
    
    ;;Change origin for rotated points
    u+=inputs.origin[0] & uc+=inputs.origin[0]
    v+=inputs.origin[1] & vc+=inputs.origin[1] 
    w+=inputs.origin[2] & wc+=inputs.origin[2]
    
    rgrid=sqrt(uc^2+vc^2)
    phigrid=atan(vc,uc)
    
    r_grid=dblarr(nx,ny,nz) & phi_grid=r_grid
    w_grid=dblarr(nx,ny,nz) & u_grid=w_grid & v_grid=w_grid
    z_grid=dblarr(nx,ny,nz) & x_grid=z_grid & y_grid=z_grid
    for i=0L,nx-1 do for j=0L,ny-1 do for k=0L,nz-1 do begin
        l=i+nx*j+nx*ny*k
        r_grid[i,j,k]=rgrid[l] & phi_grid[i,j,k]=phigrid[l]
        u_grid[i,j,k]=uc[l] & v_grid[i,j,k]=vc[l] & w_grid[i,j,k]=wc[l]
        x_grid[i,j,k]=xc[l] & y_grid[i,j,k]=yc[l] & z_grid[i,j,k]=zc[l]
    endfor
    
    grid={nx:nx,ny:ny,nz:nz,x:x,y:y,z:z,xx:xx,yy:yy,zz:zz,xc:xc,yc:yc,zc:zc,xxc:xxc,yyc:yyc,zzc:zzc,$
          dx:dx,dy:dy,dz:dz,dr:dr,drmin:drmin,dv:dv,ng:ng,u:u,v:v,w:w,$
          uc:uc,vc:vc,wc:wc,r_grid:r_grid,phi_grid:phi_grid,$
          w_grid:w_grid,v_grid:v_grid,u_grid:u_grid,$
          x_grid:x_grid,y_grid:y_grid,z_grid:z_grid}
    err=0   
    GET_OUT:
END
    
PRO prepare_beam,inputs,nbi,nbgeom
    
    nbgeom={err:1}
    isource=inputs.isource

    if nbi.pinj le 0. then begin
        printc, 'ERROR: The selected source #'+string(isource)+' is not on',f='r'
        goto, GET_OUT
    endif
    uvw_src=nbi.xyz_src-inputs.origin
    uvw_pos=nbi.xyz_pos-inputs.origin
    uvw_origin=[0,0,0]-inputs.origin
    
    make_rot_mat,-inputs.alpha,inputs.beta,Drot,Erot,Frot
    rotate_uvw,uvw_src,Drot,Erot,Frot,1,xyz_src ;;rotate from machine to beam coordinates
    rotate_uvw,uvw_pos,Drot,Erot,Frot,1,xyz_pos
    rotate_uvw,uvw_origin,Drot,Erot,Frot,1,xyz_origin
    
    xs=xyz_src[0] & ys=xyz_src[1] & zs=xyz_src[2]
    xp=xyz_pos[0] & yp=xyz_pos[1] & zp=xyz_pos[2]

    dis=sqrt( (xs-xp)^2.0d +(ys-yp)^2.0d + (zs-zp)^2.0d)
    BETA=double(asin((zp-zs)/dis))
    ALPHA=double(atan((yp-ys),(xp-xs))-!DPI)
    print,'Beam injection start point in machine coordinates'
    print, nbi.xyz_src
    print,'Beam injection end point in machine coordinates'
    print, nbi.xyz_pos
    print,'Machine center in beam coordinates'
    print, xyz_origin
    print,'Beam injection start point in beam coordinates'
    print, xyz_src
    print,'Beam injection end point in beam coordinates'
    print, xyz_pos
    print,'Beam rotation angles as defined by fidasim.f90'
    print,'Alpha: '
    print,ALPHA/!pi*180,FORMAT='(F20.10)'
    print,'Beta:'
    print,BETA,FORMAT='(F20.10)'
    ;;MAKE ROTATION MATRICES 
    make_rot_mat,ALPHA,BETA,Arot,Brot,Crot

    nbgeom={isource:isource,alpha:ALPHA,beta:BETA,Arot:Arot,Brot:Brot,Crot:Crot,xyz_src:xyz_src,xyz_pos:xyz_pos,err:0}
    GET_OUT:
END

PRO los_track,coords,xyz_los_vec,xyspt,ri,tcell,cell,ncell
    ri=xyspt
    vn=xyz_los_vec
    ;; zeros are not good!
    index=where(vn eq 0,nind)
    if nind gt 0 then vn[index]=0.0001
    p=[0,0,0]
    dummy = min( abs( coords.xxc - ri[0] ), index )
    p[0]=index
    dummy = min( abs( coords.yyc - ri[1] ), index ) 
    p[1]=index
    dummy = min( abs( coords.zzc - ri[2] ), index ) 
    p[2]=index

    tcellh=fltarr(1000)
    cellh=fltarr(3,1000)
    m=0
    cellh[*,m]=p
    ;;loop along line of sight
    while(m lt 1000) do begin
        l=p
        index=where(vn gt 0.d0)
        if index[0] ne -1 then begin
            l(index)=p(index)+1
        endif 
        if l[0] gt coords.nx-1 or l[1] gt coords.ny-1 or $
           l[2] gt coords.nz-1 then goto, out
        ;;time needed to go into next cell
        dt_arr=fltarr(3)
        dt_arr[0] = ( coords.xx(l[0]) - ri[0] ) /vn[0]
        dt_arr[1] = ( coords.yy(l[1]) - ri[1] ) /vn[1]
        dt_arr[2] = ( coords.zz(l[2]) - ri[2] ) /vn[2] 
        dt=min(dt_arr,index)
        ri[*] = ri[*] + vn[*]*dt 
        if vn[index] gt 0.d0 then begin 
            p[index]=p[index]+1
        endif else begin
            p[index]=p[index]-1
        endelse
        if p[0] lt 0 or p[1] lt 0 or p[2] lt 0 then goto, out
        tcellh[m]=dt
        m=m+1
        cellh[*,m]=p[*] 
    endwhile
    out:
    ;; Store results into compressed arrays!
    if m gt 1 then begin
        tcell= tcellh[0:m-2]
        cell = cellh[*,0:m-2]
    endif else begin
        tcell = -1
        cell = -1
    endelse
    ncell=m-1
END

PRO fida_los_wght,grid,xlens,ylens,zlens,xlos,ylos,zlos,weight,err_arr

    nx=grid.nx
    ny=grid.ny
    nz=grid.nz
    nchan=n_elements(xlens)
    err_arr=dblarr(nchan)
    weight  = replicate(0.d0,nx,ny,nz,nchan)

    for chan=0L, n_elements(xlens)-1 do  begin
        xyzlens = [xlens[chan],ylens[chan],zlens[chan]]
        xyzlos  = [xlos[chan], ylos[chan], zlos[chan]]
        vi    = xyzlos-xyzlens
        dummy = max(abs(vi),ic)
        nstep = fix(700./grid.dr[ic])
        vi    = vi/sqrt(vi[0]^2+vi[1]^2+vi[2]^2) ;; unit vector
        xyz_pos = xyzlens

      ; find first grid cell
        for i=0L,nstep do begin
            xyz_pos[0] = xyz_pos[0] + grid.dr[ic] * vi[0]/abs(vi[ic])
            xyz_pos[1] = xyz_pos[1] + grid.dr[ic] * vi[1]/abs(vi[ic])
            xyz_pos[2] = xyz_pos[2] + grid.dr[ic] * vi[2]/abs(vi[ic])
            if xyz_pos[0] gt grid.xx[0] and xyz_pos[0] lt grid.xx[nx-1]+grid.dx and $ 
                xyz_pos[1] gt grid.yy[0] and xyz_pos[1] lt grid.yy[ny-1]+grid.dy and $
                xyz_pos[2] gt grid.zz[0] and xyz_pos[2] lt grid.zz[nz-1]+grid.dz then begin
                goto, out
            endif  
        endfor
        out:

      ; determine cells along the LOS
        if i lt nstep then begin
            los_track,grid,vi,xyz_pos,rout,tcell,cell,ncell
            if ncell gt 1 then begin
                for jj=0L,ncell-1 do begin
                    if finite(tcell[jj]) eq 0 then stop
                    ;;  tcell is the length of the track (cm) as v is 1cm/s
                    weight[cell[0,jj],cell[1,jj],cell[2,jj],chan]=tcell[jj]
                endfor
            endif else begin
                printc, 'WARNING: Channel #'+strtrim(string(chan),1)+' only crosses <= 1 cells',f='y'
                err_arr[chan]=1
            endelse
        endif else begin
            err_arr[chan]=1
        endelse
    endfor
    index=where(finite(weight) eq 0,nind)
    if nind gt 0 then begin
        printc,'ATTENTION: FIDA los weight at index '+strcompress(string(index),/remove_all)+$
               ' set to 0.0 as it was NAN or Infinite!',f='y'
        weight[index]=0.
    endif
END

PRO chord_coor,pi,pf,u,v,z,xp,yp,zp

    x0=pi[0] & xf=pf[0]
    y0=pi[1] & yf=pf[1]
    z0=pi[2] & zf=pf[2]

    phi=atan((yf-y0),(xf-x0))
    theta=-atan(SQRT((xf-x0)^2.0d + (yf-y0)^2.0d),(zf-z0))

    ;;Change coordinance system to chord view (where z in entering the plasma)
    xp=(u-x0)*cos(phi)+(v-y0)*sin(phi)
    yp=-(u-x0)*sin(phi)+(v-y0)*cos(phi)
    zp=z-z0
    xpp=xp*cos(theta)+zp*sin(theta)
    zpp=-xp*sin(theta)+zp*cos(theta)
    xp=xpp
    zp=zpp
END

FUNCTION npa_prob,x,y,xp,yp,zp,dx=dx,dy=dy
    if not keyword_set(dx) then dx=abs(x[1]-x[0])
    if not keyword_set(dy) then dy=abs(y[1]-y[0])

    r=(x-xp)^2.0 + (y-yp)^2.0+zp^2.0
    p=((4*!DPI)^(-1.0))*zp*r^(-1.5)
    return, p
END

PRO npa_los_wght,los,grid,weight,err_arr

    w=where(los.chan_id eq 1,nchan)
    rd=los.rd[w]
    ra=los.ra[w]
    h=los.h[w]
    xlens=los.xlens[w] & xlos=los.xlos[w]
    ylens=los.ylens[w] & ylos=los.ylos[w]
    zlens=los.zlens[w] & zlos=los.zlos[w]

    weight  = replicate(0.d0,grid.nx,grid.ny,grid.nz,nchan)
    err_arr=replicate(0.0,nchan)
    ny=200L
    nx=200L

    for chan=0,nchan-1 do begin
        xyzlens = [xlens[chan],ylens[chan],zlens[chan]]
        xyzlos  = [xlos[chan], ylos[chan], zlos[chan]]

        ymin=-1.1d0*rd[chan]
        xmin=-1.1d0*rd[chan]
        ymax= 1.1d0*rd[chan]
        xmax= 1.1d0*rd[chan]
        x = xmin + dindgen(nx)*(xmax-xmin)/nx
        y = ymin + dindgen(ny)*(ymax-ymin)/ny
        dx=abs(x[1]-x[0])
        dy=abs(y[1]-y[0])
        xd = reform(rebin(x,nx,ny,/sample),nx*ny)
        yd = reform(transpose(rebin(y,ny,nx,/sample)),nx*ny)    
        rrd = sqrt(xd^2 + yd^2)

        chord_coor,xyzlens,xyzlos,grid.u_grid,grid.v_grid,grid.w_grid,xp,yp,zp
        zp=zp+h[chan]
        ww=where(zp gt h[chan],nw)
        alpha=zp*0
        xcen=zp*0
        ycen=xcen
        rs=xcen+rd[chan]
        if nw ne 0 then begin
            alpha[ww]=zp[ww]/(h[chan]-zp[ww])
            rs[ww]=abs(ra[chan]*alpha[ww])
            xcen[ww]=-xp[ww]-xp[ww]*alpha[ww]
            ycen[ww]=-yp[ww]-yp[ww]*alpha[ww]
        endif
    
        for xi=0, grid.nx-1 do begin
            for yi=0, grid.ny-1 do begin
                for zi=0, grid.nz-1 do begin
                    if sqrt(xcen[xi,yi,zi]^2 + ycen[xi,yi,zi]^2) gt rd[chan]+rs[xi,yi,zi] then continue
                    xs=xd+xcen[xi,yi,zi]
                    ys=yd+ycen[xi,yi,zi]
                    rrs=sqrt(xs^2.0 + ys^2.0)
                    p=npa_prob(xd,yd,xp[xi,yi,zi],yp[xi,yi,zi],zp[xi,yi,zi],dx=dx,dy=dy)
                    www=where(rrs ge rs[xi,yi,zi] or rrd ge rd[chan],nw)
                    if nw ne 0 then p[www]=0
                    weight[xi,yi,zi,chan]=total(p*dx*dy)    
                endfor
            endfor
        endfor
        
        if total(weight[*,*,*,chan]) le 0 then err_arr[chan]=1
    endfor
            
END

PRO prepare_chords,inputs,grid,chords,fida

    nx=grid.nx
    ny=grid.ny
    nz=grid.nz

    ;;DECLARE ARRAYS
    err_arr=dblarr(chords.nchan) + 1
    weight  = replicate(0.d0,nx,ny,nz,chords.nchan)

    ;;CALCULATE RADIUS IN MACHINE COORDINATES
    rlos=sqrt(chords.xlos^2.0 + chords.ylos^2.0)

    ;;ROTATE CHORDS INTO BEAM COORDINATES
    make_rot_mat,-inputs.alpha,inputs.beta,Arot,Brot,Crot
    ulens=chords.xlens-inputs.origin[0] & ulos=chords.xlos-inputs.origin[0]
    vlens=chords.ylens-inputs.origin[1] & vlos=chords.ylos-inputs.origin[1]
    wlens=chords.zlens-inputs.origin[2] & wlos=chords.zlos-inputs.origin[2]
    rotate_points,ulens,vlens,wlens,Arot,Brot,Crot,xlens,ylens,zlens
    rotate_points,ulos,vlos,wlos,Arot,Brot,Crot,xlos,ylos,zlos

    ;;CALCULATE FIDA LOS WEIGHTS
    w=where(chords.chan_id eq 0,nw)
    if nw ne 0 and (inputs.calc_spec or inputs.calc_fida_wght) then begin
        print,'Calculating FIDA line of sight weights'
        fida_los_wght,grid,xlens[w],ylens[w],zlens[w],xlos[w],ylos[w],zlos[w],fida_wght,fida_err
        weight[*,*,*,w]=fida_wght
        err_arr[w]=fida_err
    endif

    ;;CALCULATE NPA LOS WEIGHTS
    w=where(chords.chan_id eq 1,nw)
    if nw ne 0 and (inputs.calc_npa or inputs.calc_npa_wght) then begin
        print,'Calculating NPA line of sight weights'
        npa_los_wght,chords,grid,npa_wght,npa_err
        weight[*,*,*,w]=npa_wght
        err_arr[w]=npa_err
    endif

    ;;GET RID OF LOS THAT DONT CROSS THE GRID
    los=where(err_arr eq 0,nw,complement=miss_los,ncomplement=nww)
    if nw eq 0 then begin
        printc,'ERROR: No chords crossed the simulation grid',f='r'
        err=1
    endif else begin
        print,strtrim(string(nw),1)+' out of '+strtrim(string(chords.nchan),1)+' chords crossed the simulation grid'
        if nww ne 0 then begin
            printc, 'ATTENTION: Missed indices:  '+strcompress(strjoin(string(miss_los))),f='y'
            printc, 'ATTENTION: Missed chan_ids: '+strcompress(strjoin(string(chords.chan_id[miss_los],f='(i2)'))),f='y'
        endif
        weight=weight[*,*,*,los]
        err=0
    endelse

    fida={nchan:n_elements(los),xlens:xlens[los],ylens:ylens[los],zlens:zlens[los],sigma_pi_ratio:chords.sigma_pi_ratio[los],$
          xlos:xlos[los],ylos:ylos[los],zlos:zlos[los],rlos:rlos[los],chan_id:chords.chan_id[los],$
          ra:chords.ra[los],rd:chords.rd[los],h:chords.h[los],los:los,weight:weight,err:err}
END

PRO transp_fbeam,inputs,grid,denf,fbm_struct,err

    !P.charsize=1.
    !P.background=255 & !P.color=0
    ;doplot=1
    ;print, 'reading fast ion distribution function from transp output'
    cdftest=findfile(inputs.cdf_file)
    ;print, '======================='
    if cdftest[0] eq '' then begin
        printc,'ERROR: '+inputs.cdf_file+' was not found',f='r'
        err=1
        goto,GET_OUT
    endif
    cdfid=NCDF_Open(inputs.cdf_file,/nowrite)
    ;; Retrieve signals
    ;; --------------------------------------------------------
    ncdf_varget, cdfid,'TRANSP_RUNID', runid
    ncdf_varget, cdfid,'TIME' , cdf_time    
    ncdf_varget, cdfid,'R2D'  , r2d     ; rposition of cells
    ncdf_varget, cdfid,'Z2D'  , z2d     ;zposition of cells
    ncdf_varget, cdfid,'BMVOL', bmvol   ; plasma volume
    ncdf_varget, cdfid,'E_D_NBI', energy ; central values
    ncdf_varget, cdfid,'A_D_NBI', pitch  ; central values
    ncdf_varget, cdfid,'F_D_NBI', FBM    ; fast-ion distribution function
    ncdf_varget, cdfid,'NTOT_D_NBI',ntot ; total number of fast ions
    ncdf_varget, cdfid,'RSURF', rsurf    ; flux surface
    ncdf_varget, cdfid,'ZSURF', zsurf    ; flux surface
    NCDF_Close,cdfid
    ngrid=n_elements(r2d)
    ;;================================
    ;; get tranpped -passing boundary
    ;;================================
    rmin=fltarr(ngrid)
    for i=0,ngrid -1 do begin
        dummy=min((rsurf-r2d[i])^2+(zsurf-z2d[i])^2,index)
        index = array_indices(rsurf, index)
        rmin[i]=min(rsurf[*,index[1]])
    endfor
    pitch_boundary=sqrt(1.-rmin[*]/r2d[*])
    
    ;; ----------- Check the time
    if abs(inputs.time-cdf_time) gt 0.02 then begin
        print, ' CDF file time:',cdf_time 
        printc, 'WARNING: Time of CDF file and simulation disagree!',f='y'
    endif     
    ;;-------------Convert eV-> keV
    energy=energy*1.0d-3          ;; fidasim needs energy in kev  
    fbm=fbm*1.0d3                 ;; now, this needs to be corrected
    ;; as we now calculate with fast-ions/omega/keV/cm^3 
    ;;------------Convert d_omega --> pitch
    ;; Fast-ion distribution is given as a function of cm^3, energy
    ;; and d_omega/4Pi. omega is the solild angle in 3D velocity space. In
    ;; order to transform this to a function depending on pitch instead
    ;; of d_omega/4PI, one has to multiply by 0.5!
    fbm*=0.5  
    ;; make sure that fbm is >=0:
    fbm>=0.
    ;;loading finished

    ;; TRANSP defines the pitch along the current direction. In
    ;; contrast, FIDASIM uses pitch along the B-field! Therefore,
    ;; reverse the pitch coordinate in fbm!
    if inputs.btipsign lt 0 then begin
        npitch=n_elements(pitch)
        index=npitch-(indgen(npitch)+1)
        fbm[*,*,*]=fbm[*,index,*]
    endif
    ;;----------- select energy range -------
    index=where(energy ge inputs.emin and energy le inputs.emax,nenergy)
    energy=energy[index]     
    fbm=fbm[index,*,*]
    dE      = energy[1] - energy[0]
    emin=(float(energy[0])         - float(0.5*dE))>0.
    emax=float(energy[nenergy-1]) + float(0.5*dE)
    print, 'Energy min/max:', emin,emax
    ;; --------- select Pitch range --------
    index=where(pitch ge inputs.pmin and pitch le inputs.pmax,npitch)
    pitch=pitch[index] 
    fbm=fbm[*,index,*]   
    dP  = abs(pitch[1]  - pitch[0])
    pmin=(float(pitch[0])       - float(0.5*dP))>(-1)
    pmax=(float(pitch[npitch-1])+ float(0.5*dP))<1
    print, 'Pitch  min/max:', pmin,pmax
    
    fbm_struct={cdf_file:inputs.cdf_file,cdf_time:cdf_time,ngrid:ngrid,r2d:r2d,z2d:z2d,bmvol:bmvol,$
                nenergy:nenergy,emin:emin,emax:emax,energy:energy,npitch:npitch,$
                pmin:pmin,pmax:pmax,pitch:pitch,fbm:FBM}
    
    ;; ------map fdens on FIDASIM grid and sort out
    ;; ------points outside the separatrix
    fdens=total(reform(total(fbm,1)),1)*dE*dP
    a=dblarr(grid.ng)
    if ngrid le 220 then width=6. else width=4.
    rout=reform(grid.r_grid,grid.ng)
    zout=reform(grid.w_grid,grid.ng)
    TRIANGULATE, r2d, z2d, tr
    fdens2=griddata(r2d,z2d,fdens,xout=rout,yout=zout,/linear,triangles=tr)

    ;;set negative values to zero
    fdens2=fdens2 >0.

    ;; only write fdens2 if it is close to r2d,z2d grid 
    ;;(this produced empty spots so i turned it off)
;   for i=0L,grid.ng-1 do a[i]=min(sqrt((z2d-zout[i])^2+(r2d-rout[i])^2))
;   ww=where(a gt width,nw)
;   if nw ne 0 then fdens2[ww]=0.
    denf=reform(fdens2,grid.nx,grid.ny,grid.nz)
    err=0
    GET_OUT:
END

PRO map_profiles,inputs,grid,equil,profiles,plasma,err
    ;;-------------------------------------------------
    ;; MAP kinetic profiles on FIDASIM grid
    ;;------------------------------------------------- 

    rhomax=max(profiles.rho)
    ww=where(equil.rho_grid gt rhomax,nww)
    if nww ne 0 then begin
        printc,'ATTENTION: ',f='y'
        printc,'    dene[rho > '+string(rhomax,f='(1f8.3)')+'] = 1e10 [cm^-3]',f='y'
        printc,'    te[  rho > '+string(rhomax,f='(1f8.3)')+'] = 1e-3 [keV]',f='y'
        printc,'    ti[  rho > '+string(rhomax,f='(1f8.3)')+'] = 1e-3 [keV]',f='y'
        printc,'    vtor[rho > '+string(rhomax,f='(1f8.3)')+'] = 0 [cm/s]',f='y'
        printc,'    zeff[rho > '+string(rhomax,f='(1f8.3)')+'] = 1 ',f='y'
    endif 

    ;;Electron density
    dene = 1.d-6 * interpol(profiles.dene,profiles.rho,equil.rho_grid) > 0. ;[1/cm^3]
    if nww ne 0 then dene[ww]=0.001*1d13

    ;;Zeff
    zeff = interpol(profiles.zeff,profiles.rho,equil.rho_grid) > 1.0 
    zeff = zeff < inputs.impurity_charge
    if nww ne 0 then zeff[ww]=1.0

    ;;Impurity density
    deni = (zeff-1.)/(inputs.impurity_charge*(inputs.impurity_charge-1))*dene
    
    ;;Proton density
    denp = dene-inputs.impurity_charge*deni
    print,'Percent impurity: '+string(total(deni)/total(denp)*100.,f='(1f8.3)')+' %'
    
    ;;Fast-ion density
    if inputs.load_fbm then begin
        transp_fbeam,inputs,grid,denf,fbm_struct,terr
        if terr eq 1 then begin
            printc,'ERROR: Failed to map fast ion density',f='r'
            err=1
            goto,GET_OUT
        endif
    endif else begin
        denf=dene*0.d0
        fbm_struct={err:1}
    endelse

    ;;Electron temperature
    te = 1.d-3 * interpol(profiles.te,profiles.rho,equil.rho_grid) > 0.001 ;keV
    if nww ne 0 then te[ww]=0.001
    
    ;;Ion temperature   
    ti = 1.d-3 * interpol(profiles.ti,profiles.rho,equil.rho_grid) > 0.001 ;keV
    if max(ti) gt 20. or max(te) gt 20. then begin
        printc, 'WARNING:',f='y'
        printc, '    Electron or Ion temperature greater than 10 keV',f='y'
        printc, '    Look at the tables, they might only consider',f='y'
        printc, '    temperatures less than 10keV!',f='y'
    endif
    if nww ne 0 then ti[ww]=0.001

    ;;Plasma rotation   
    vtor      =   interpol(profiles.vtor,profiles.rho,equil.rho_grid)*grid.r_grid ; [cm/s]  
    if nww ne 0 then vtor[ww]  =   replicate(0.0,nww)*grid.r_grid[ww]

    vrotu = - sin(grid.phi_grid)*vtor 
    vrotv =   cos(grid.phi_grid)*vtor
    vrotw =   0.d0*vrotu 

    ;;Rotate vector quantities to beam coordinates 
    make_rot_mat,-inputs.alpha,-inputs.beta,Arot,Brot,Crot
    rotate_points,vrotu,vrotv,vrotw,Arot,Brot,Crot,vrotx,vroty,vrotz ;;machine basis to beam basis
    rotate_points,equil.bx,equil.by,equil.bz,Arot,Brot,Crot,bx,by,bz
    rotate_points,equil.ex,equil.ey,equil.ez,Arot,Brot,Crot,ex,ey,ez

    ;; test if there are NANs or Infinites in the input profiels
    index=where(finite([ti,te,dene,denp,zeff,denp,deni]) eq 0,nind)
    if nind gt 0 then stop
    ;;-------SAVE-------
    plasma={fbm:fbm_struct,rho_grid:equil.rho_grid,$
            bx:bx,by:by,bz:bz,ex:ex,ey:ey,ez:ez,$
            bu:equil.bx,bv:equil.by,bw:equil.bz,eu:equil.ex,ev:equil.ey,ew:equil.ez,$
            ab:inputs.ab,ai:inputs.ai,$
            vrotx:vrotx,vroty:vroty,vrotz:vrotz,$
            vrotu:vrotu,vrotv:vrotv,vrotw:vrotw,$
            te:te,ti:ti,vtor:vtor,dene:dene,denp:denp,deni:deni,denf:denf,zeff:zeff}
    err=0
    GET_OUT: 
END

FUNCTION sinterpol,v,x,u,sortt=sortt,_extra=_extra
    if n_elements(sortt) lt 1 then sortt=0

    if sortt then begin
        ind=sort(X)
    endif else begin
        ind=lindgen(n_elements(x))
    endelse

    return,interpol(v[ind],x[ind],u,_extra=_extra)
END

PRO brems,inputs,det,profiles,equil,vbline
    ;; Calculates visible bremsstrahlung along FIDA sightlines
    ;; WWH 6/2013

    ;; INPUT
    ;; result_dir directory to write output
    ;; det       structure with detector lines of sight
    ;; profiles   structure with plasma profiles vs. rho

    ;; OUTPUT
    ;; file with the surface radiance (ph/s-m2-nm-sr) for
    ;; each sightline

    ;*******************************************************************
    ;*****************************************
    ; Plasma parameters
    rho=profiles.rho
    te=profiles.te              ; eV
    dene=profiles.dene*1.e-6    ; cm^-3
    zeff=profiles.zeff
    
    ; Require non-zero values for te and ne
    w=where(te le 0. or dene le 0.,nw)
    if nw gt 0 then begin
        rho=rho[0:w[0]-1]
        te=te[0:w[0]-1]
        dene=dene[0:w[0]-1]
        zeff=zeff[0:w[0]-1]
    endif
    w=where(zeff lt 1.,nw) & if nw gt 0 then zeff[w]=1.
    w=where(zeff gt 6.,nw) & if nw gt 0 then zeff[w]=6.
    rhomax=max(rho,nr) & nr+=1

    ;**********************************************
    ; Constants in calculation
    lambda=6561.    ; average wavelength (Angstroms)
    h_planck=4.135667e-15  ; [eV/s]
    c0=2.9979e8 ; [m/s]
    
    ; Visible bremsstrahlung emissivity versus rho
    gaunt=5.542-(3.108-alog(te/1000.))*(0.6905-0.1323/zeff)
    emisrho=10.*7.57d-9*gaunt*dene^2*zeff/(lambda*sqrt(te)) $
               *exp(-h_planck*c0/(lambda*te))

    ;***********************************************
    ;NOW do line integration to get surface radiance
    ;***********************************************
    nchan=det.nchan
    los=det.los
    vbline=replicate(0.,nchan)

    for i=0,nchan-1 do begin
        rhospath=equil.rho_chords.rhos[*,los[i]]
        rhospath=rhospath[where(finite(rhospath))]
        vbepath=sinterpol(emisrho,rho,rhospath,/sort)
        wgtr1=where(rhospath ge rhomax,nwgtr1)
        ;set emission at radii outside of max rho to zero
        if nwgtr1 ge 0 then vbepath[wgtr1]=0.0
        vbline[i]=total(vbepath)*equil.rho_chords.ds*0.001*inputs.dlambda*(4*!DPI)*1.d-4 > 0.001  ; (ph/s-m2-bin)
    endfor  ; channel loop
END

PRO write_namelist,inputs

    printc,'Writing namelist file...',f='c'
    spawn,'which git',git_command
    git_hash = ''
    if git_command ne '' then begin 
        spawn,git_command+' --git-dir='+inputs.install_dir+'.git rev-parse HEAD',git_hash
    endif
    filename = inputs.result_dir+inputs.runid+'/inputs.nml'
    openw,55,filename
    printf,55,'!! Created: ', systime()
    if git_hash ne '' then begin
        printf,55,'!! FIDASIM git commit: ',git_hash
    endif else begin
        printf,55,'!! FIDASIM version: '
    endelse

    printf,55,'&fidasim_inputs'
    printf,55,''
    printf,55,'!! Shot Info'
    printf,55,'!! Diagnostic: ',inputs.diag
    printf,55,f='("shot = ", i6 ,"    !! Shot Number")',inputs.shot    
    printf,55,f='("time = ", 1f8.5 ,"    !! Time")',inputs.time    
    printf,55,"runid = '" + inputs.runid + "'    !! runID"
    printf,55,''
    printf,55,'!! Simulation Switches'
    printf,55,f='("calc_spec = ",i2 , "    !! Calculate Spectra")',inputs.calc_spec
    printf,55,f='("calc_npa = ",i2 , "   !! Calculate NPA")',inputs.calc_npa
    printf,55,f='("calc_birth = ",i2 , "    !! Calculate Birth Profile")',inputs.calc_birth
    printf,55,f='("calc_fida_wght = ",i2 , "    !! Calculate FIDA weights")',inputs.calc_fida_wght
    printf,55,f='("calc_npa_wght = ",i2 , "    !! Calculate NPA weights")',inputs.calc_npa_wght
    printf,55,f='("calc_brems = ",i2,"    !! Calculate Bremsstrahlung else load from inputs")',inputs.calc_brems
    printf,55,f='("load_neutrals = ",i2,"    !! Load Neutrals")',inputs.load_neutrals
    printf,55,f='("load_fbm = ",i2,"    !! Load FBM")',inputs.load_fbm
    printf,55,''
    printf,55,'!! Wavelength Grid Settings'
    printf,55,f='("nlambda = ",1i5,"    !! Number of Wavelengths")',inputs.nlambda
    printf,55,f='("lambdamin = ",1f9.3,"    !! Minimum Wavelength")',inputs.lambdamin
    printf,55,f='("lambdamax = ",1f9.3,"    !! Maximum Wavelength")',inputs.lambdamax
    printf,55,''
    printf,55,'!! Monte Carlo Settings'
    printf,55,f='("nr_fast = ",i9,"    !! Number of FAST mc particles")',inputs.nr_fast
    printf,55,f='("nr_nbi = ",i9,"    !! Number of NBI mc particles")',inputs.nr_nbi
    printf,55,f='("nr_halo = ",i9,"    !! Number of HALO mc particles")',inputs.nr_halo
    printf,55,''
    printf,55,'!! Weight Function Settings'
    printf,55,f='("ne_wght = ",i9,"    !! Number of Energies")',inputs.ne_wght
    printf,55,f='("np_wght = ",i9,"    !! Number of Pitches")',inputs.np_wght
    printf,55,f='("nphi_wght = ",i9,"    !! Number of Gyro-angles")',inputs.nphi_wght
    printf,55,f='("ichan_wght = ",i3,"    !! Channel for weight function")',inputs.ichan_wght
    printf,55,f='("emax_wght = ",1f12.2,"    !! Maximum Energy for Weights")',inputs.emax_wght
    printf,55,f='("dwav_wght = ",1f12.5,"    !! Wavelength Seperation for Weights ")',inputs.dwav_wght
    printf,55,f='("wavel_start_wght = ",1f12.5,"    !! Wavelength Start for Weights ")',inputs.wavel_start_wght
    printf,55,f='("wavel_end_wght = ",1f12.5,"    !! Wavelength End for Weights ")',inputs.wavel_end_wght
    printf,55,''
    printf,55,'/'
    close,55
    printc, 'SUCCESS: Namelist file created: '+filename,f='g'
END

PRO check_inputs,inputs,err

    printc,'Checking input file...',f='c'
    vars=["shot","time","runid","device","install_dir","result_dir","cdf_file","profile_dir",$
          "emin","emax","pmin","pmax","isource","diag","einj","pinj","equil","btipsign","ab",$
          "ai","impurity_charge","lambdamin","lambdamax","nlambda","dlambda","nx","ny","nz",$
          "xmin","xmax","ymin","ymax","zmin","zmax","origin","alpha","beta","nr_fast","nr_nbi",$
          "nr_halo","ne_wght","np_wght","nphi_wght","emax_wght","ichan_wght","dwav_wght",$
          "wavel_start_wght","wavel_end_wght","calc_npa","calc_spec","calc_birth","calc_fida_wght",$
          "calc_npa_wght","calc_brems","load_neutrals","load_fbm"]

    inVars=strlowcase(TAG_NAMES(inputs))

    err=0
    for i=0,n_elements(vars)-1 do begin
       w=where(vars[i] eq inVars,nw)
       if nw eq 0 then begin
           printc,'ERROR: '+vars[i]+' missing from input file',f='r'
           err=1
       endif
    endfor
    
    if (inputs.calc_spec or inputs.calc_npa) and (not inputs.load_fbm) then begin 
       printc,'WARNING: load_fbm needs to be set',f='y'
       printc,'ATTENTION: setting load_fbm=1',f='y'
       inputs.load_fbm=1
    endif

    if inputs.calc_spec and inputs.calc_npa then begin
        printc,'WARNING: calc_spec and calc_npa cannot both be set',f='y'
        printc,'ATTENTION: setting calc_spec=1 & calc_npa = 0',f='y'
        inputs.calc_npa=0
    endif

    if err ne 0 then begin
        printc,'ERROR: Invalid input file. Exiting...',f='r'
    endif else begin
        printc,'SUCCESS: Input file is valid',f='g'
    endelse

END

PRO prefida,input_file,plot=plot,save=save

    COMPILE_OPT DEFINT32

    if n_elements(input_file) eq 0 then begin
        printc,'ERROR: Input file not specified',f='r'
        goto,GET_OUT
    endif

    if FILE_TEST(input_file) then begin
        ;;READ JSON INPUT FILE
        input_proc=0
        inputs=read_json(input_file)
        inputs=create_struct('install_dir',GETENV('FIDASIM_DIR'),inputs)
    endif else begin
        ;;CALL INPUT PROCEDURE/FILE
        input_proc=1
        printc,'WARNING: Input procedure is depreciated. Use JSON input file if possible',f='y'
        CALL_PROCEDURE,input_file,inputs
    endelse
    
    ;;CHECK INPUTS
    check_inputs,inputs,err
    if err ne 0 then goto,GET_OUT

    ;;CHECK FOR SLASH
    slash=strmid(inputs.result_dir,0,1,/reverse_offset)
    if slash ne '/' then inputs.result_dir+='/'
    slash=strmid(inputs.install_dir,0,1,/reverse_offset)
    if slash ne '/' then inputs.install_dir+='/'
    slash=strmid(inputs.profile_dir,0,1,/reverse_offset)
    if slash ne '/' then inputs.profile_dir+='/'
    
    ;;MAKE DIRECTORIES IF THEY DONT EXIST
    if file_test(inputs.result_dir,/directory) eq 0 then begin
        spawn,'mkdir '+inputs.result_dir
    endif
    if file_test(inputs.result_dir+inputs.runid,/directory) eq 0 then begin
        spawn,'mkdir '+inputs.result_dir+inputs.runid
    endif

    ;;ADD DEVICE DIRECTORY TO PATH
    !path = !path + ":" + expand_path("+"+inputs.install_dir+inputs.device)
    
    ;;ADD INSTALL DIRECTORY TO PATH
    !path = !path + ":" + expand_path(inputs.install_dir)

    ;;MAKE GRID
    printc,'Making grid...',f='c'
    make_grid,inputs,grid,err
    if err eq 1 then begin
        printc,'ERROR: Grid creation failed. Exiting...',f='r'
        goto,GET_OUT
    endif else begin
        printc,'SUCCESS: Grid creation completed',f='g'
        err=0
    endelse

    ;;CALL DEVICE ROUTINES THAT GET BEAM GEOMETRY, FIDA DIAGNOSTIC INFO, PROFILES, and the grid in flux coord.
    printc,'Calling device routines...',f='c'
    CALL_PROCEDURE, strlowcase(inputs.device)+'_routines',inputs,grid, nbi, chords, profiles, equil,err
    if err eq 1 then begin
        printc, 'ERROR: Device routines failed. Exiting...',f='r'
        goto,GET_OUT
    endif else begin
        printc,'SUCCESS: Device routines completed',f='g'
        err=0
    endelse

    ;;BEAM PRE PROCESSING
    printc,'Pre-processing beams...',f='c'
    prepare_beam,inputs,nbi,nbgeom
    if nbgeom.err eq 1 then begin
        printc,'ERROR: Beam pre-processing failed. Exiting...',f='r'
        goto, GET_OUT
    endif else begin
        printc,'SUCCESS: Beam pre-processing completed',f='g'
        err=0
    endelse

    ;;FIDA PRE PROCESSING 
    if inputs.calc_spec or inputs.calc_npa or inputs.calc_fida_wght or inputs.calc_npa_wght then begin
        printc,'Pre-processing chords...',f='c'
        prepare_chords,inputs,grid,chords,fida
        if fida.err eq 1 then begin
            printc,'ERROR: Chord pre-processing failed. Exiting...',f='r'
            goto, GET_OUT
        endif else begin
            printc,'SUCCESS: Chord pre-processing completed',f='g'
        endelse
    endif else begin
        err=0
        weight=replicate(0.d0,inputs.nx,inputs.ny,inputs.nz,1)
        fida={weight:weight,err:err}
    endelse

    ;;MAP PROFILES ONTO GRID
    printc,'Mapping profiles onto grid...',f='c'
    map_profiles,inputs,grid,equil,profiles,plasma,err
    if err eq 1 then begin
        printc,'ERROR: Profile mapping failed. Exiting...',f='r'
        goto,GET_OUT
    endif else begin
        printc,'SUCCESS: Profile mapping completed',f='g'
        err=0
    endelse

    ;; Calculate bremsstrahlung if desired
    if inputs.calc_brems eq 0 then begin
        printc,'Calculating bremsstrahlung...',f='c'
        brems,inputs,fida,profiles,equil,brems
        printc,'SUCCESS: Bremsstrahlung calculation completed',f='g'
    endif

    plot_file=inputs.install_dir+strupcase(inputs.device)+'/'+strlowcase(inputs.device)+'_plots.pro'
    ;; Plot grid, beam, sightlines, and equilibrium
    if keyword_set(plot) and FILE_TEST(plot_file) then begin
        CALL_PROCEDURE, strlowcase(inputs.device)+'_plots',inputs,grid,nbi,chords,fida,equil,nbgeom,plasma
    endif

    ;;SAVE STRUCTURES 
    if keyword_set(save) then begin
        file = inputs.result_dir+inputs.runid+'/'+inputs.runid+'.sav'
        save,inputs,grid,profiles,chords,nbi,equil,nbgeom,fida,plasma,filename=file,/compress
    endif

    ;;COPY INPUT PROCEDURE/FILE TO RESULT DIRECTORY
    if input_proc then begin
        file_info=ROUTINE_INFO(input_file,/source)
        file_path=file_info.path
        file_name=input_file+'.pro'
    endif else begin
        file_path=input_file
        file_name=FILE_BASENAME(input_file)
    endelse

    FILE_COPY,file_path,$
              inputs.result_dir+inputs.runid+'/'+file_name,$
              /overwrite,/allow_same

    ;;WRITE FIDASIM INPUT FILES
    write_namelist,inputs

    printc,'Writing input data file...',f='c'
    ;;WRITE TO FILE
    file =inputs.result_dir+inputs.runid+'/'+inputs.runid+'_inputs.cdf' 
    ncid = ncdf_create(file,/clobber)

    ;;DEFINE DIMENSIONS
    ncdf_control,ncid
    one_id = ncdf_dimdef(ncid,'dim001',1)
    three_id = ncdf_dimdef(ncid,'dim003',3)

    ndiag_id=ncdf_dimdef(ncid,'ndiag',n_elements(inputs.diag))
    strmax_id=ncdf_dimdef(ncid,'maxstr',max(strlen(inputs.diag)))

    if inputs.load_fbm then begin
        fbm_gdim= ncdf_dimdef(ncid,'fbm_grid',plasma.fbm.ngrid)
        fbm_edim=ncdf_dimdef(ncid,'fbm_energy',plasma.fbm.nenergy)  
        fbm_pdim=ncdf_dimdef(ncid,'fbm_pitch',plasma.fbm.npitch)    
    endif

    xid = ncdf_dimdef(ncid,'x',grid.nx)
    yid = ncdf_dimdef(ncid,'y',grid.ny)
    zid = ncdf_dimdef(ncid,'z',grid.nz)

    if inputs.calc_spec or inputs.calc_npa or inputs.calc_fida_wght or inputs.calc_npa_wght then begin
        chan_id=ncdf_dimdef(ncid,'chan',n_elements(fida.los))
        xyz_dim=[three_id,chan_id]
    endif

    griddim=[xid,yid,zid]
    dim3d=[three_id,three_id]
    diagstr_dim=[strmax_id,ndiag_id]
    ;;DEFINE VARIABLES
    shot_varid=ncdf_vardef(ncid,'shot',one_id,/long)
    time_varid=ncdf_vardef(ncid,'time',one_id,/float)
    diag_varid=ncdf_vardef(ncid,'diagnostic',diagstr_dim,/char)

    ;;SIZE VARIABLES
    nx_varid=ncdf_vardef(ncid,'Nx',one_id,/long)
    ny_varid=ncdf_vardef(ncid,'Ny',one_id,/long)
    nz_varid=ncdf_vardef(ncid,'Nz',one_id,/long)
    if inputs.load_fbm then begin
        gdim_varid=ncdf_vardef(ncid,'FBM_Ngrid',one_id,/long)
        edim_varid=ncdf_vardef(ncid,'FBM_Nenergy',one_id,/long)
        pdim_varid=ncdf_vardef(ncid,'FBM_Npitch',one_id,/long)
    endif
    if inputs.calc_spec or inputs.calc_npa or inputs.calc_fida_wght or inputs.calc_npa_wght then $
        nchan_varid=ncdf_vardef(ncid,'Nchan',one_id,/long)

    ;;DEFINE GRID VARIABLES
    ugrid_varid=ncdf_vardef(ncid,'u_grid',griddim,/double)
    vgrid_varid=ncdf_vardef(ncid,'v_grid',griddim,/double)
    wgrid_varid=ncdf_vardef(ncid,'w_grid',griddim,/double)
    rgrid_varid=ncdf_vardef(ncid,'r_grid',griddim,/double)
    phigrid_varid=ncdf_vardef(ncid,'phi_grid',griddim,/double)
    xgrid_varid=ncdf_vardef(ncid,'x_grid',griddim,/double)
    ygrid_varid=ncdf_vardef(ncid,'y_grid',griddim,/double)
    zgrid_varid=ncdf_vardef(ncid,'z_grid',griddim,/double)
    alpha_varid=ncdf_vardef(ncid,'alpha',one_id,/double)
    beta_varid=ncdf_vardef(ncid,'beta',one_id,/double)
    origin_varid=ncdf_vardef(ncid,'origin',three_id,/double)
    xx_varid=ncdf_vardef(ncid,'xx',xid,/double)
    yy_varid=ncdf_vardef(ncid,'yy',yid,/double)
    zz_varid=ncdf_vardef(ncid,'zz',zid,/double)

    ;;DEFINE BEAM VARIABLES
    bn_varid=ncdf_vardef(ncid,'beam',one_id,/long)
    ab_varid=ncdf_vardef(ncid,'ab',one_id,/double)
    divy_varid=ncdf_vardef(ncid,'divy',three_id,/double)    
    divz_varid=ncdf_vardef(ncid,'divz',three_id,/double)
    focy_varid=ncdf_vardef(ncid,'focy',one_id,/double)    
    focz_varid=ncdf_vardef(ncid,'focz',one_id,/double)
    bmwidra_varid=ncdf_vardef(ncid,'bmwidra',one_id,/double)
    bmwidza_varid=ncdf_vardef(ncid,'bmwidza',one_id,/double)
    einj_varid=ncdf_vardef(ncid,'einj',one_id,/double)
    pinj_varid=ncdf_vardef(ncid,'pinj',one_id,/double)
    sm_varid=ncdf_vardef(ncid,'species_mix',three_id,/double)
    xyzsrc_varid=ncdf_vardef(ncid,'xyz_src',three_id,/double)
    arot_varid=ncdf_vardef(ncid,'Arot',dim3d,/double)
    brot_varid=ncdf_vardef(ncid,'Brot',dim3d,/double)
    crot_varid=ncdf_vardef(ncid,'Crot',dim3d,/double)

    ;;DEFINE FBM VARIABLES
    if inputs.load_fbm then begin
        r2d_varid=ncdf_vardef(ncid,'FBM_r2d',fbm_gdim,/double)
        z2d_varid=ncdf_vardef(ncid,'FBM_z2d',fbm_gdim,/double)
        bmvol_varid=ncdf_vardef(ncid,'FBM_bmvol',fbm_gdim,/double)
        energy_varid=ncdf_vardef(ncid,'FBM_energy',fbm_edim,/double)
        pitch_varid=ncdf_vardef(ncid,'FBM_pitch',fbm_pdim,/double)
        emin_varid=ncdf_vardef(ncid,'FBM_emin',one_id,/double)
        emax_varid=ncdf_vardef(ncid,'FBM_emax',one_id,/double)
        pmin_varid=ncdf_vardef(ncid,'FBM_pmin',one_id,/double)
        pmax_varid=ncdf_vardef(ncid,'FBM_pmax',one_id,/double)
        cdftime_varid=ncdf_vardef(ncid,'FBM_time',one_id,/double)   
        fbm_varid=ncdf_vardef(ncid,'FBM',[fbm_edim,fbm_pdim,fbm_gdim],/double)
    endif

    ;;DEFINE PLASMA VARIABLES
    ai_varid=ncdf_vardef(ncid,'ai',one_id,/double)
    impc_varid=ncdf_vardef(ncid,'impurity_charge',one_id,/float)
    btip_varid=ncdf_vardef(ncid,'btipsign',one_id,/float)
    te_varid=ncdf_vardef(ncid,'te',griddim,/double)
    ti_varid=ncdf_vardef(ncid,'ti',griddim,/double)
    dene_varid=ncdf_vardef(ncid,'dene',griddim,/double)
    deni_varid=ncdf_vardef(ncid,'deni',griddim,/double)
    denp_varid=ncdf_vardef(ncid,'denp',griddim,/double)
    denf_varid=ncdf_vardef(ncid,'denf',griddim,/double)
    vx_varid=ncdf_vardef(ncid,'vrotx',griddim,/double)
    vy_varid=ncdf_vardef(ncid,'vroty',griddim,/double)
    vz_varid=ncdf_vardef(ncid,'vrotz',griddim,/double)
    zeff_varid=ncdf_vardef(ncid,'zeff',griddim,/double)
    bx_varid=ncdf_vardef(ncid,'bx',griddim,/double)
    by_varid=ncdf_vardef(ncid,'by',griddim,/double)
    bz_varid=ncdf_vardef(ncid,'bz',griddim,/double)
    ex_varid=ncdf_vardef(ncid,'ex',griddim,/double)
    ey_varid=ncdf_vardef(ncid,'ey',griddim,/double)
    ez_varid=ncdf_vardef(ncid,'ez',griddim,/double)
    rho_varid=ncdf_vardef(ncid,'rho_grid',griddim,/double)

    ;;DEFINE BREMSTRUHLUNG VARIABLES
    brems_varid=ncdf_vardef(ncid,'brems',chan_id,/double)

    ;;LOS VARIABLE DEFINITION
    if inputs.calc_spec or inputs.calc_npa or inputs.calc_fida_wght or inputs.calc_npa_wght then begin
        xlens_varid=ncdf_vardef(ncid,'xlens',chan_id,/double)
        ylens_varid=ncdf_vardef(ncid,'ylens',chan_id,/double)
        zlens_varid=ncdf_vardef(ncid,'zlens',chan_id,/double)
        xlos_varid=ncdf_vardef(ncid,'xlos',chan_id,/double)
        ylos_varid=ncdf_vardef(ncid,'ylos',chan_id,/double)
        zlos_varid=ncdf_vardef(ncid,'zlos',chan_id,/double)
        rlos_varid=ncdf_vardef(ncid,'rlos',chan_id,/double)
        ra_varid=ncdf_vardef(ncid,'ra',chan_id,/double)
        rd_varid=ncdf_vardef(ncid,'rd',chan_id,/double)
        h_varid=ncdf_vardef(ncid,'h',chan_id,/double)
        chan_id_varid=ncdf_vardef(ncid,'chan_id',chan_id,/double)
        sig_varid=ncdf_vardef(ncid,'sigma_pi',chan_id,/double)
        wght_varid=ncdf_vardef(ncid,'los_wght',[xid,yid,zid,chan_id],/double)
    endif
    ;;END DEFINITION
    ncdf_control,ncid,/ENDEF

    ;;WRITE VARIABLES TO FILE
    ;;WRITE ARRAY SIZES
    ncdf_varput,ncid,shot_varid,long(inputs.shot)
    ncdf_varput,ncid,time_varid,double(inputs.time)
    ncdf_varput,ncid,diag_varid,inputs.diag
    ncdf_varput,ncid,nx_varid,long(inputs.nx)
    ncdf_varput,ncid,ny_varid,long(inputs.ny)
    ncdf_varput,ncid,nz_varid,long(inputs.nz)
    if inputs.calc_spec or inputs.calc_npa or inputs.calc_fida_wght or inputs.calc_npa_wght then $
        ncdf_varput,ncid,nchan_varid,long(n_elements(fida.los))
    if inputs.load_fbm then begin
        ncdf_varput,ncid,gdim_varid,long(plasma.fbm.ngrid)
        ncdf_varput,ncid,edim_varid,long(plasma.fbm.nenergy)
        ncdf_varput,ncid,pdim_varid,long(plasma.fbm.npitch)
    endif

    ;;WRITE GRID VARIABLES
    ncdf_varput,ncid,ugrid_varid,double(grid.u_grid)    
    ncdf_varput,ncid,vgrid_varid,double(grid.v_grid)    
    ncdf_varput,ncid,wgrid_varid,double(grid.w_grid)    
    ncdf_varput,ncid,rgrid_varid,double(grid.r_grid)    
    ncdf_varput,ncid,phigrid_varid,double(grid.phi_grid)    
    ncdf_varput,ncid,xgrid_varid,double(grid.x_grid)    
    ncdf_varput,ncid,ygrid_varid,double(grid.y_grid)    
    ncdf_varput,ncid,zgrid_varid,double(grid.z_grid)    
    ncdf_varput,ncid,alpha_varid,double(inputs.alpha)
    ncdf_varput,ncid,beta_varid,double(inputs.beta)
    ncdf_varput,ncid,origin_varid,double(inputs.origin)
    ncdf_varput,ncid,xx_varid,double(grid.xx)
    ncdf_varput,ncid,yy_varid,double(grid.yy)
    ncdf_varput,ncid,zz_varid,double(grid.zz)

    ;;WRITE BEAM VARIABLES
    ncdf_varput,ncid,bn_varid,long(inputs.isource[0])
    ncdf_varput,ncid,ab_varid,double(inputs.ab)
    ncdf_varput,ncid,divy_varid,double(nbi.divy)
    ncdf_varput,ncid,divz_varid,double(nbi.divz)
    ncdf_varput,ncid,focy_varid,double(nbi.focy)
    ncdf_varput,ncid,focz_varid,double(nbi.focz)
    ncdf_varput,ncid,bmwidra_varid,double(nbi.BMWIDRA)
    ncdf_varput,ncid,bmwidza_varid,double(nbi.BMWIDZA)
    ncdf_varput,ncid,einj_varid,double(nbi.einj)
    ncdf_varput,ncid,pinj_varid,double(nbi.pinj)
    ncdf_varput,ncid,sm_varid,double([nbi.full,nbi.half,nbi.third])
    ncdf_varput,ncid,xyzsrc_varid,double(nbgeom.xyz_src)
    ncdf_varput,ncid,arot_varid,double(nbgeom.Arot)
    ncdf_varput,ncid,brot_varid,double(nbgeom.Brot)
    ncdf_varput,ncid,crot_varid,double(nbgeom.Crot)

    ;;WRITE FBM VARIABLES
    if inputs.load_fbm then begin
        ncdf_varput,ncid,r2d_varid,double(plasma.fbm.r2d)
        ncdf_varput,ncid,z2d_varid,double(plasma.fbm.z2d)
        ncdf_varput,ncid,bmvol_varid,double(plasma.fbm.bmvol)
        ncdf_varput,ncid,energy_varid,double(plasma.fbm.energy)
        ncdf_varput,ncid,pitch_varid,double(plasma.fbm.pitch)
        ncdf_varput,ncid,emin_varid,double(plasma.fbm.emin)
        ncdf_varput,ncid,emax_varid,double(plasma.fbm.emax)
        ncdf_varput,ncid,pmin_varid,double(plasma.fbm.pmin)
        ncdf_varput,ncid,pmax_varid,double(plasma.fbm.pmax)
        ncdf_varput,ncid,cdftime_varid,double(plasma.fbm.cdf_time)
        ncdf_varput,ncid,fbm_varid,double(plasma.fbm.fbm)
    endif

    ;;WRITE PLASMA VARIABLES
    ncdf_varput,ncid,ai_varid, double(inputs.ai)
    ncdf_varput,ncid,btip_varid,float(inputs.btipsign)
    ncdf_varput,ncid,impc_varid,float(inputs.impurity_charge)
    ncdf_varput,ncid,te_varid, double(plasma.te)
    ncdf_varput,ncid,ti_varid, double(plasma.ti)
    ncdf_varput,ncid,dene_varid, double(plasma.dene)
    ncdf_varput,ncid,denp_varid, double(plasma.denp)
    ncdf_varput,ncid,deni_varid, double(plasma.deni)
    ncdf_varput,ncid,denf_varid, double(plasma.denf)
    ncdf_varput,ncid,vx_varid, double(plasma.vrotx)
    ncdf_varput,ncid,vy_varid, double(plasma.vroty)
    ncdf_varput,ncid,vz_varid, double(plasma.vrotz)
    ncdf_varput,ncid,zeff_varid, double(plasma.zeff)
    ncdf_varput,ncid,bx_varid, double(plasma.bx)
    ncdf_varput,ncid,by_varid, double(plasma.by)
    ncdf_varput,ncid,bz_varid, double(plasma.bz)
    ncdf_varput,ncid,ex_varid, double(plasma.ex)
    ncdf_varput,ncid,ey_varid, double(plasma.ey)
    ncdf_varput,ncid,ez_varid, double(plasma.ez)
    ncdf_varput,ncid,rho_varid, double(plasma.rho_grid)

    ;;WRITE BREMS
    if n_elements(brems) ne 0 then ncdf_varput,ncid,brems_varid,double(brems)
    
    ;;WRITE LINE OF SIGHT (LOS)
    if inputs.calc_spec or inputs.calc_npa or inputs.calc_fida_wght or inputs.calc_npa_wght then begin
        los=fida.los
        ncdf_varput,ncid,xlens_varid,double(fida.xlens)
        ncdf_varput,ncid,ylens_varid,double(fida.ylens)
        ncdf_varput,ncid,zlens_varid,double(fida.zlens)
        ncdf_varput,ncid,xlos_varid,double(fida.xlos)
        ncdf_varput,ncid,ylos_varid,double(fida.ylos)
        ncdf_varput,ncid,zlos_varid,double(fida.zlos)
        ncdf_varput,ncid,rlos_varid,double(fida.rlos)
        ncdf_varput,ncid,ra_varid,double(fida.ra)
        ncdf_varput,ncid,rd_varid,double(fida.rd)
        ncdf_varput,ncid,h_varid,double(fida.h)
        ncdf_varput,ncid,chan_id_varid,double(fida.chan_id)
        ncdf_varput,ncid,sig_varid,double(fida.sigma_pi_ratio)
        ncdf_varput,ncid,wght_varid,double(fida.weight)
    endif
    ncdf_close,ncid
    printc,'SUCCESS: Input data file created: '+file,f='g'

    print,''
    print,''
    printc, 'SUCCESS: FIDASIM pre-processing completed',f='g'
    print, 'To run FIDASIM use the following command'
    print, inputs.install_dir+'fidasim '+inputs.result_dir+inputs.runid
    print,''
    print,''
    GET_OUT:
END

