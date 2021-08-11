function sunflower,n,alpha=alpha,doplot=doplot
    ;+#sunflower
    ;+Generates 2D array of `n` cartesian points in (x,y) that are approximately uniformly spaced
    ;+***
    ;+##Arguments
    ;+    **n**: Number of desired points
    ;+
    ;+##Keyword Arguments
    ;+    **alpha**: Parameter that determines boundary points
    ;+
    ;+    **doplot**: Plot samples
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> n = 10000
    ;+IDL> da = sunflower(n)
    ;+```

    if ~keyword_set(alpha) then alpha=1

    b=round(alpha*sqrt(n))  ; number of boundary pts
    phi=(sqrt(5)+1)/2       ; golden ratio
    r=fltarr(n) & theta=r
    for i=1,n do begin
      if i gt n-b then r[i-1]=1. else r[i-1]=sqrt(i-.5)/sqrt(n-0.5*(b+1))
      theta[i-1]=2*!pi*i/phi^2
    end

    if keyword_set(doplot) then begin
      plot,r*cos(theta),r*sin(theta),psym=2,xrange=1.1*[-1,1],yrange=1.1*[-1,1]
      phi=2*!pi*findgen(101)/100
      oplot,cos(phi),sin(phi),color=100
    end

    a=fltarr(2,n)
    a[0,*]=r*cos(theta) & a[1,*]=r*sin(theta)
    return,a

end


PRO orb_collimator,g,ri,vi,d,a,naperture,vsave,frac,norm,step=step,time_reverse=time_reverse,nsteps=nsteps,e0=e0,amu=amu,z=z,narea=narea,straight=straight
    ;+#orb_collimator
    ;+Calculates solid-angle weights of different velocity vectors that exit the collimator
    ;+***
    ;+##Arguments
    ;+    **g**: GEQDSK file
    ;+
    ;+    **ri**: Centered on detector in (R,phi,z) coordinates [m,radians,m]
    ;+
    ;+    **vi**: Axis orientation of collimator, i.e. velocity components in (R,phi,z) coordinates
    ;+
    ;+    **d**: Collimator length [m]
    ;+
    ;+    **a**: Collimator radius [m]
    ;+
    ;+    **naperture**: Number of velocities to launch
    ;+
    ;+##Keyword Arguments
    ;+    **step**: Step length [m]
    ;+
    ;+    **time_reverse**: Indicates reversed time orbits
    ;+
    ;+    **nsteps**: Number of steps
    ;+
    ;+    **e0**: Energy [keV]
    ;+
    ;+    **amu**: Atomic mass unit
    ;+
    ;+    **z**: Particle charge
    ;+
    ;+    **narea**: Number of positions to launch from
    ;+
    ;+    **straight**: Bypass orbit calculation and use straight rays for testing
    ;+
    ;+##Outputs
    ;+    **vsave**: Launch velocity vectors -- size(3,naperture)
    ;+
    ;+    **frac**:  Fraction of area that clears collimator -- size(naperture)
    ;+
    ;+    **norm**:  norm*total(frac) equals A*Omega
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> g = 'g000001.01000'
    ;+IDL> detector_aperture_geometry,0,rdist,zdist,v,d,rc
    ;+IDL> ri = [rdist[0],0,zdist[0]]
    ;+IDL> vi = -reform(v[*,0])
    ;+IDL> orb_collimator,g,ri,vi,d,rc[0],50,vsave,frac,norm
    ;+```

    common bcom,b0,r0,br,bphi,bz,gr0,gz0,dr,dz

    ; Ion orbit parameters
    if not keyword_set(step) then step=0.01	; step length in m
    if not keyword_set(nsteps) then nsteps=110
    if not keyword_set(e0) then e0=3030.			; keV
    if not keyword_set(amu) then amu=1. & mp=1.67e-27
    if not keyword_set(z) then z=1.
    if not keyword_set(narea) then narea=10000
    if not keyword_set(time_reverse) then time_reverse=1

    ;-------------
    ; Same for every calculation
    ; Use eqdsk to get wall location and magnetic field grid
    b0=abs(g.bcentr) & r0=g.rmaxis
    finewall,g,rwall,zwall
    calculate_bfield,bp,br,bphi,bz,g
    br=double(br) & bphi=double(bphi) & bz=double(bz)
    gr0=double(g.r(0)) & gz0=double(g.z(0))
    dr=double(g.r(1)-g.r(0)) & dz=double(g.z(1)-g.z(0))
    br=double(br) & bphi=double(bphi) & bz=double(bz)
    pphisgn=-g.cpasma/abs(g.cpasma)

    ; Normalization constants
    omega=z*1.6e-19*b0/(amu*mp)
    v0=sqrt(2*e0*1.e3*1.6e-19/(amu*mp))
    vconstant=v0/omega

    ;------
    ; Coordinate system along collimator tube
    ; unit vector along tube
    vr=vi[0] & vphi=vi[1] & vz=vi[2]

    ; transverse unit vector
    if vr eq 0 and vz eq 0 then begin
      ar=1. & az=0. & aphi=0.
    end else begin
      aphi=0.
      ar=vz/sqrt(vz^2+vr^2)
      az=-vr/sqrt(vz^2+vr^2)
    end

    ; toroidal unit vector
    cr=-az*vphi
    cphi=az*vr-ar*vz
    cz=ar*vphi

    ;------------
    ; delta positions on detector area
    ; use a modified sunflower arrangement based on the golden rule
    da=a*sunflower(narea)  ; (x,y) displacement array
    eps=1.02   ; fudge factor on whether ray clears

    ;---------
    ; Velocity angles to use
    vsave=fltarr(3,naperture)
    frac=replicate(0.,naperture)

    ; Find the maximum amount the curved orbit can expand the aperture
    rho=vconstant*sqrt(vr^2 + vz^2)
    if rho lt d then return
    expand=rho - sqrt(rho^2 - d^2)
    aperture=a + expand

    ; Points on effective aperture to aim velocities at
    daperture=aperture*sunflower(naperture) ; (x,y) target array

    ; Orbit initialization
    vi0=vi

    first=1
    ;-----------------
    ; Loop over initial velocities
    for iaperture=0,naperture-1 do begin
    ;  ivel=launch_vector(ivel0,tan(daperture[0,iaperture]/d),tan(daperture[1,iaperture]/d))
        xx=daperture[0,iaperture] & yy=daperture[1,iaperture]
        ivel=[ar*xx + cr*yy + vr*d, $
            aphi*xx + cphi*yy + vphi*d, $
            az*xx + cz*yy + vz*d]
        ivel/=sqrt(ivel[0]^2 + ivel[1]^2 + ivel[2]^2)
        vsave[*,iaperture]=ivel

        ; Preparing normal orb_mast calculation from center of detector
        Velocity=vconstant*ivel
        y=dblarr(6) & y(0:2)=Velocity & y(3:5)=ri
        h=double(step*omega/v0) & if time_reverse then h=-h
        yout=dblarr(6,nsteps)
        if keyword_set(straight) then begin
          s=sqrt(d^2+a^2)*findgen(nsteps)/(nsteps-1)
          uu=y[3]*cos(y[4]) + s*ivel[0]
          vv=y[3]*sin(y[4]) + s*ivel[1]
          zz=y[5] + s*ivel[2]
          yout[3,*]=sqrt(uu^2 + vv^2)
          yout[4,*]=atan(vv,uu)
          yout[5,*]=zz
        end else begin
            yout(*,0)=y(*)

            ; Orbit loop
            i=0
            lwall=1		; logical to stop if hits wall
            while i lt nsteps-1 do begin
                i=i + 1
              ; Van Zeeland's integrator
                ddeabm,'derivs',0.d,y,h,epsabs=1.e-8
              ; force energy conservation
              ;  y(0:2)=y(0:2)*vconstant/sqrt(y(0)^2+y(1)^2+y(2)^2)
                yout(*,i)=y(*)
            end
        end ; not straight keyword

        ;---------
        ; Find fraction of detector area that misses walls for this orbit
        ; Use an (xx,yy,zz) coordinate system with zz axis along collimator
        ;
        u=reform(yout[3,*]*cos(yout[4,*])) & u-=u[0]
        v=reform(yout[3,*]*sin(yout[4,*])) & v-=v[0]
        w=reform(yout[5,*]) & w-=w[0]

        xx=ar*u + aphi*v + az*w
        yy=cr*u + cphi*v + cz*w
        zz=vr*u + vphi*v + vz*w
        zz*=-1.

        ; Truncate orbit at end of collimator
        w=where(zz le d)
        xx=xx[w] & yy=yy[w] & zz=zz[w]

        ; Start this trajectory from various positions on detector area
        ; and store how many don't hit the wall

        plotarea=0 ; Optional plotting of which part of detector area works
        if plotarea then plot,a*cos(2*!pi*findgen(101)/100),a*sin(2*!pi*findgen(101)/100), $
                              title=string(1e3*daperture[0,iaperture])+' '+string(1e3*daperture[1,iaperture])
        for j=0,narea-1 do begin
            w=where((xx + da[0,j])^2 + (yy + da[1,j])^2 gt eps*a^2,nw)
            if nw lt 1 then frac[iaperture]+=1.
            if plotarea then begin
                if nw lt 1 then begin
                    col=250
                endif else begin
                    col=100
                endelse
                oplot,[da[0,j]],[da[1,j]],color=col,psym=1
            end
        end

        if plotarea then wait,0.5

        ; Optional plotting of trajectories
        doplot=0
        if doplot then begin
            if first then begin
                !p.multi=[0,0,1]
                plot,zz,sqrt(xx^2+yy^2),yrange=1.2*[min(sqrt(xx^2+yy^2)),max(sqrt(xx^2+yy^2))]
                first=0
            endif else begin
                oplot,zz,sqrt(xx^2+yy^2),color=100+fix(frac[iaperture])
            endelse
        end ; doplot

        ; cos(theta) factor
        frac[iaperture]*=d/sqrt(d^2 + daperture[0,iaperture]^2 + daperture[1,iaperture]^2)

    end ; iaperture loop

    frac/=narea

    norm=(!pi*a^2)*(!pi*aperture^2)/(4*!pi*d^2*naperture)

  end



FUNCTION orb_cfpd, g, rdist, zdist, v, d, rc, e0=e0, amu=amu, z=z, nrays=nrays, step=step, time_reverse=time_reverse, nsteps=nsteps, plot_show=plot_show
    ;+#orb_cfpd
    ;+Returns all of the orbit-related quantities for a given energy
    ;+***
    ;+##Arguments
    ;+    **g**: GEQDSK file
    ;+
    ;+    **rdist**: Radial coordinates of detector; (nchan), [m]
    ;+
    ;+    **zdist**: Vertical coordinates of detector; (nchan), [m]
    ;+
    ;+    **v**: Detector orientations; (3,nchan), [m/s, rad/s, m/s]
    ;+
    ;+    **d**: Collimator length; (nchan), [m]
    ;+
    ;+    **rc**: Outer collimator spacing; (nchan), [m]
    ;+
    ;+##Keyword Arguments
    ;+    **e0**: Energy (keV)
    ;+
    ;+    **amu**: Atomic mass unit
    ;+
    ;+    **z**: Particle charge
    ;+
    ;+    **nrays**: Number of orbits
    ;+
    ;+    **nsteps**: Number of steps
    ;+
    ;+    **step**: Step length [m]
    ;+
    ;+    **time_reverse**: Indicates reversed time orbits
    ;+
    ;+    **plot_show**: Plot orbit bundles -- (1) on (0) off
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> g = 'g000001.01000'
    ;+IDL> detector_aperture_geometry,0,rdist,zdist,v,d,rc
    ;+IDL> orb = orb_cfpd(g,rdist,zdist,v,d,rc,e0=3030,nrays=50,nsteps=110)
    ;+```

  common bcom,b0,r0,br,bphi,bz,gr0,gz0,dr,dz

  if ~keyword_set(e0) then e0=3030.; keV
  if ~keyword_set(amu) then amu=1
  if ~keyword_set(z) then z=1.
  if ~keyword_set(nrays) then nrays=50
  if ~keyword_set(nsteps) then nsteps=110
  if ~keyword_set(step) then step=0.01
  if ~keyword_set(time_reverse) then time_reverse=1
  if ~keyword_set(plot_show) then plot_show=0

  max_steps = 48 ;Skips inital steps along trajectory to avoid boundary check

  ; storage arrays
  nch=n_elements(rc)
  daomega=fltarr(nrays,nch)
  sightline=fltarr(6,nsteps,nrays,nch)
  initial_velocities=fltarr(3,nrays,nch)
  nactual=intarr(nrays,nch)

  ;TRANSMISSION FACTORS
  ; collimator
  for ich=0,(nch-1) do begin
    orb_collimator,g,[rdist[ich],0,zdist[ich]],-reform(v[*,ich]),d[ich],rc[ich],nrays,vsave,frac,norm, $
      e0=e0
    initial_velocities[*,*,ich]=vsave
    daomega[*,ich]=norm*frac
  end

  ;---------------
  ; Use orb to get sightlines
  ; As in orb_collimator, embed it to avoid unnecessary repetition

  ; Ion orbit parameters
  mp=1.67e-27

  ;-------------
  ; Same for every calculation
  ; Use eqdsk to get wall location and magnetic field grid
  b0=abs(g.bcentr) & r0=g.rmaxis
  finewall,g,rwall,zwall
  calculate_bfield,bp,br,bphi,bz,g
  br=double(br) & bphi=double(bphi) & bz=double(bz)
  gr0=double(g.r(0)) & gz0=double(g.z(0))
  dr=double(g.r(1)-g.r(0)) & dz=double(g.z(1)-g.z(0))
  br=double(br) & bphi=double(bphi) & bz=double(bz)
  pphisgn=-g.cpasma/abs(g.cpasma)

  ; Normalization constants
  omega=z*1.6e-19*b0/(amu*mp)
  v0=sqrt(2*e0*1.e3*1.6e-19/(amu*mp))
  vconstant=v0/omega

  ;-----------------
  ; Loop over sightlines
  for ich=0,nch-1 do begin
    ri=[rdist[ich],0,zdist[ich]]
    for iray=0,nrays-1 do if daomega[iray,ich] gt 0. then begin

  Velocity=vconstant*reform(initial_velocities[*,iray,ich])
  y=dblarr(6) & y(0:2)=Velocity & y(3:5)=ri
  h=double(step*omega/v0) & if time_reverse then h=-h
  yout=dblarr(6,nsteps)
  dydx=derivs(0.,y)
  yout(*,0)=y(*)

  ; Orbit loop
  i=0
  lwall=1		; logical to stop if hits wall
  inside = check_bdry(g,rwall,zwall,y(3),y(5))
  while i lt nsteps-1 and inside do begin
    i=i + 1
  ; Van Zeeland's integrator
    ddeabm,'derivs',0.d,y,h,epsabs=1.e-8
    yout(*,i)=y(*)
    dydx=derivs(0.,y)
    if i gt max_steps then inside = check_bdry(g,rwall,zwall,y(3),y(5))
  end
  nactual[iray,ich]=i+1
  for j=0,5 do sightline[j,0:i,iray,ich]=reform(yout[j,0:i])

  end ; iray loop
  end ; ich loop

  sightline[0:2,*,*,*]*=omega

  ;--------------
  ; Plot results
  if plot_show eq 1 then begin
    !p.multi=[0,2,0]
    !p.background=0 ;black background
    ; (R,z) elevation
      xmin=min(g.r) & xmax=max(g.r) & ymin=min(g.z) & ymax=max(g.z)

    contour,g.psirz,g.r,g.z,nlevels=10,color=255, $
      xrange=[xmin,xmax],yrange=[ymin,ymax], $
      xtitle='X [m]',ytitle='Y [m]',charsize=2
    oplot,g.lim(0,*),g.lim(1,*),color=40
    oplot,g.bdry(0,0:g.nbdry-1),g.bdry(1,0:g.nbdry-1),color=60
    colors=[255,100,200,300]
    for ich=0,nch-1 do for iray=0,nrays-1 do if daomega[iray,ich] gt 0. then begin
     nact=nactual[iray,ich]-1
     oplot,sightline[3,0:nact,iray,ich],sightline[5,0:nact,iray,ich],psym=3,color=colors[ich]
    end

    ; (R,phi) plan
    theta=2.*!pi*findgen(31)/30.
      xmin=-max(g.r) & xmax=max(g.r) & ymin=xmin & ymax=xmax
    plot,max(rwall)*cos(theta),max(rwall)*sin(theta), $
      xrange=[-2.5,2.5],yrange=[-4,4],xstyle=1,$
      xtitle='X [m]',ytitle='Y [m]',color=255,charsize=2
    ;  xrange=[xmin,xmax],yrange=[ymin,ymax]
    oplot,min(rwall)*cos(theta),min(rwall)*sin(theta),color=255
    for ich=0,nch-1 do for iray=0,nrays-1 do if daomega[iray,ich] gt 0. then begin
      nact=nactual[iray,ich]-1
      oplot,sightline[3,0:nact,iray,ich]*cos(sightline[4,0:nact,iray,ich]), $
            sightline[3,0:nact,iray,ich]*sin(sightline[4,0:nact,iray,ich]),color=colors[ich]
    end
  endif


  return,{sightline:sightline,nrays:nrays,nch:nch,nactual:nactual,daomega:daomega}

end
