FUNCTION nstx_beams,inputs,doplot=doplot
    ; Returns the nbi structure for one NBI source
    ;D. Liu 5/18/2014, for the existings on-axis NBs on NSTX
    ;;  IDL> help,nbi
    ;;  ** Structure <1d475af8>, 13 tags, length=168, data length=168, refs=1:
    ;;     EINJ            DOUBLE           80.775734
    ;;     PINJ            DOUBLE           2.4117758
    ;;     FULL            DOUBLE          0.54850105
    ;;     HALF            DOUBLE          0.28972649
    ;;     THIRD           DOUBLE          0.16177245
    ;;     XYZ_SRC         DOUBLE    Array[3]
    ;;     XYZ_POS         DOUBLE    Array[3]
    ;;     BMWIDRA         DOUBLE           6.0000000
    ;;     BMWIDZA         DOUBLE           24.000000
    ;;     DIVY            DOUBLE    Array[3]
    ;;     DIVZ            DOUBLE    Array[3]
    ;;     FOCY            DOUBLE           999999.90
    ;;     FOCZ            DOUBLE           1000.0000


    ; machine coordinate (u,v,z): typical N-S-E-W type definitions of machine hall
    
    ;Cross-over point of the on-axis NBI
    ;See Lane Roquemore's SSNPA coordinate graph
    ;    and Ron Bell's Neutral Beam Geometry graph
    
    rc_onaxis=158.58             ;major radius at crossover point of three NB lines
    phic_onaxis=79.298/180.*!pi  ;phi angle at crossover point of three NB lines
    uc_onaxis=rc_onaxis*cos(phic_onaxis)
    vc_onaxis=rc_onaxis*sin(phic_onaxis)
    print,'Machine coordinates of cross-over point of on-axis NBIs is ',uc_onaxis, vc_onaxis, 0.
    
    ;the angle difference of TRANSP coordinates and FIDA (u,v,z) coordinates
    angle_transp_uvz=142.2/180.*!pi ;It seems this value is not accurate, 139 instead of 142? D. Liu
    
    ;rc_offaxis=158.58             ;major radius at crossover point of three NB lines
    ;phic_offaxis=79.298/180.*!pi  ;phi angle at crossover point of three NB lines
    ;uc_offaxis=rc_offaxis*cos(phic_offaxis)
    ;vc_offaxis=rc_offaxis*sin(phic_offaxis)
    ;print,'Machine coordinates of cross-over point of on-axis NBIs is ',uc_offaxis, vc_offaxis, 0.
    
    ;Define the coordinates of beam apertures in (R, theta,z) coordinate
    ;angle (radians) in (u,v) of injected beam
    alpha_NB=fltarr(3)           
    ;(u,v) coordinates (cm) of neutral beam aperture
    rs_NB=fltarr(3)
    phis_NB=fltarr(3)
    us_NB=fltarr(3)
    vs_NB=fltarr(3)
    ;(u,v) coordinates (cm) of neutral beam ion souce
    phis_NB_source=fltarr(3)
    us_NB_source=fltarr(3)
    vs_NB_source=fltarr(3)
    
    ; Geometrical factors from TRANSP namelist
    ;-------------------------------
    ;  NB     1A       1B      1C
    XBZETA=[-85.5,    -82,    -78.5] ;toroidal angle of beam source in (R,zeta,Z)
    XLBAPA=[510.3,    509,    510.3] ;distance of ion source to beam aperture
    XLBTNA=[1134.,   1138.,   1142.] ;distance of ion source to beam tangency radius
    RTCENA=[69.4,    59.2,    49.7 ] ;beam tangency radius
    ysource=6.0     ; ion source half width in cm
    zsource=21.5    ; ion source half height in cm
    yedge=11.       ; beam aperture half width in cm
    zedge=30.25     ; beam aperture half width in cm
    xedge=XLBAPA    ; distance from ion source to beam aperture
    
    ;;for NSTX there is no focus on y and direction
    ;;Set focy and focz to be a very large value
    focy=0.988E3; (TRANSP)  ; horizontal focal length (infinity)
    focz=0.988E3; (TRANSP)  ; vertical focal length
    divy=4.94E-3    ; horizontal divergence in radians (TRANSP)
    divz=1.36E-2    ; vertical divergence in radians (TRANSP)
    ;divy=0.43/2./180.*!pi ;radians, Ron Bell's nominal value
    ;divz=1.08/2./180.*!pi ;radians, Ron Bell's nominal value
    
    ;;;Calculate the angle difference of TRANSP coordinates and FIDA (u,v,z) coordinates
    ;angle_transp_uvz=replicate(0.,n_elements(XBZETA))
    ;rc=rc_onaxis
    ;phic=phic_onaxis
    ;for i=0,n_elements(XBZETA)-1 do begin
    ;    angle_transp_uvz[i]=abs(XBZETA[i])/180.*!pi+phic-$
    ;                       ( asin(RTCENA[i]/rc)-$
    ;                         atan(RTCENA[i]/XLBTNA[i]) )
    ;endfor
    ;print,'Angle between TRANSP and machine coordiantes',angle_transp_uvz/!pi*180,mean(angle_transp_uvz/!pi*180)
    ;angle_transp_uvz=mean(angle_transp_uvz)
        
    for i=0,2 do begin
        ;;alpha_NB: angle of NB injection line in machine coordinates (FIDA uvz cordinates)
        alpha_NB[i]=!pi+XBZETA[i]/180.*!pi+angle_transp_uvz
        phis_NB[i]=XBZETA[i]/180.*!pi+angle_transp_uvz+$
                   atan(abs(RTCENA[i])/(XLBTNA[i]-XLBAPA[i]))          
        rs_NB[i]=sqrt(RTCENA[i]^2.+(XLBAPA[i]-XLBTNA[i])^2.)
        ;;(us_NB, vs_NB): beam aperture location in machine coordinates
        us_NB[i]=rs_NB[i]*cos(phis_NB[i])
        vs_NB[i]=rs_NB[i]*sin(phis_NB[i])
        ;;(us_NB, vs_NB): location of beam ion source in machine coordinates
        us_NB_source[i]=us_NB[i]+xedge[i]*cos(alpha_NB[i]-!pi)
        vs_NB_source[i]=vs_NB[i]+xedge[i]*sin(alpha_NB[i]-!pi)
        phis_NB_source[i]=XBZETA[i]/180.*!pi+angle_transp_uvz+$
                          atan(abs(RTCENA[i])/(XLBTNA[i]))
        us_NB_source[i]=(sqrt(RTCENA[i]^2.+(XLBTNA[i])^2.))*cos(phis_NB_source[i])
        vs_NB_source[i]=(sqrt(RTCENA[i]^2.+(XLBTNA[i])^2.))*sin(phis_NB_source[i])
        print,'Angle of on-axis NB injection line for NBI source ',strtrim(i,2), ' is ',alpha_NB[i]/!pi*180 
    endfor

    alpha_NB_onaxis=alpha_NB
    us_NB_onaxis=us_NB
    vs_NB_onaxis=vs_NB
    us_NB_source_onaxis=us_NB_source
    vs_NB_source_onaxis=vs_NB_source
    
    ;------------------------------
    
    
    ;------------------------------
    ;off axis NBI
    ;  NB     2A       2B      2C
    XBZETA=[-39.2,   -35.7,    -32.2] ;toroidal angle of beam source in (R,zeta,Z)
    XLBAPA=[510.3,    510.3,   510.3] ;distance of ion source to beam aperture
    XLBTNA=[1146.,    1146.,   1146.] ;distance of ion source to beam tangency radius
    RTCENA=[130.,     120.,    110. ] ;beam tangency radius
    ysource=6.0     ; ion source half width in cm
    zsource=21.5    ; ion source half height in cm
    yedge=11.       ; beam aperture half width in cm
    zedge=30.25     ; beam aperture half width in cm
    xedge=XLBAPA    ; distance from ion source to beam aperture
    
    ;;for NSTX there is no focus on y and direction
    ;;Set focy and focz to be a very large value
    focy=0.988E3; (TRANSP)  ; horizontal focal length (infinity)
    focz=0.988E3; (TRANSP)  ; vertical focal length
    divy=4.94E-3    ; horizontal divergence in radians (TRANSP)
    divz=1.36E-2    ; vertical divergence in radians (TRANSP)
    ;divy=0.43/2./180.*!pi ;radians, Ron Bell's nominal value
    ;divz=1.08/2./180.*!pi ;radians, Ron Bell's nominal value
    
    for i=0,2 do begin
        ;;alpha_NB: angle of NB injection line in machine coordinates (FIDA uvz cordinates)
        alpha_NB[i]=!pi+XBZETA[i]/180.*!pi+angle_transp_uvz
        phis_NB[i]=XBZETA[i]/180.*!pi+angle_transp_uvz+$
                   atan(abs(RTCENA[i])/(XLBTNA[i]-XLBAPA[i]))          
        rs_NB[i]=sqrt(RTCENA[i]^2.+(XLBAPA[i]-XLBTNA[i])^2.)
        ;;(us_NB, vs_NB): beam aperture location in machine coordinates
        us_NB[i]=rs_NB[i]*cos(phis_NB[i])
        vs_NB[i]=rs_NB[i]*sin(phis_NB[i])
        ;;(us_NB, vs_NB): location of beam ion source in machine coordinates
        us_NB_source[i]=us_NB[i]+xedge[i]*cos(alpha_NB[i]-!pi)
        vs_NB_source[i]=vs_NB[i]+xedge[i]*sin(alpha_NB[i]-!pi)
        phis_NB_source[i]=XBZETA[i]/180.*!pi+angle_transp_uvz+$
                          atan(abs(RTCENA[i])/(XLBTNA[i]))
        us_NB_source[i]=(sqrt(RTCENA[i]^2.+(XLBTNA[i])^2.))*cos(phis_NB_source[i])
        vs_NB_source[i]=(sqrt(RTCENA[i]^2.+(XLBTNA[i])^2.))*sin(phis_NB_source[i])
        print,'Angle of off-axis NB injection line for NBI source ',strtrim(i,2), ' is ',alpha_NB[i]/!pi*180 
    endfor
    
    alpha_NB_offaxis=alpha_NB
    us_NB_offaxis=us_NB
    vs_NB_offaxis=vs_NB
    us_NB_source_offaxis=us_NB_source
    vs_NB_source_offaxis=vs_NB_source
    
    alpha_NB=[alpha_NB_onaxis,alpha_NB_offaxis]
    us_NB=[us_NB_onaxis,us_NB_offaxis]
    vs_NB=[vs_NB_onaxis,vs_NB_offaxis]
    us_NB_source=[us_NB_source_onaxis,us_NB_source_offaxis]
    vs_NB_source=[vs_NB_source_onaxis,vs_NB_source_offaxis]
    
    ;;coordinates of beam apertures and cross over points in machine coordinate
    nsources=n_elements(us_NB)
    
    ;xyz_src=[[double(us_NB)],[double(vs_NB)],[replicate(0.0d,nsources)]]
    ;u_co=[replicate(uc_onaxis,3),replicate(uc_offaxis,3)]
    ;v_co=[replicate(vc_onaxis,3),replicate(vc_offaxis,3)]
    ;xyz_pos=[[double(u_co)],[double(v_co)],[replicate(0.0d,nsources)]]
    
    xyz_src=[[double(us_NB_source)],[double(vs_NB_source)],[replicate(0.0d,nsources)]]
    xyz_pos=[[double(us_NB)],[double(vs_NB)],[replicate(0.0d,nsources)]]
    
    bmwidra=replicate(ysource,nsources)   ;array of ion source half width in cm
    bmwidza=replicate(zsource,nsources)   ;array of ion source half height in cm
    divy=replicate(divy,nsources)
    divz=replicate(divz,nsources)
    focy=replicate(focy,nsources)
    focz=replicate(focz,nsources)
    
    einj=double(inputs.einj) & pinj=double(inputs.pinj)
    
    ;------------------------------------------
    ;;Get beam energy,power and fractions
    ;;D. Liu, need the subroutine get_nstx_beam 5/18/2014
    ;if inputs.einj eq 0. or inputs.pinj eq 0. then $
    ;   a=get_beam_power(inputs.shot,inputs.time*1000.,inputs.isource)
    ;if inputs.einj eq 0. then einj=a.einj else einj=inputs.einj
    ;if inputs.pinj eq 0. then pinj=a.pinj else pinj=inputs.pinj
    einj=double(einj) & pinj=double(pinj)
    
    ; Based on the subroutine /u/rbell/cfrac.pro
    ; Data for deuterium 
    denergy=[0., 80., 100.,  120., 200.]
    dmix1=[.467,.467,.440,.392,.392]
    dmix2=[.374,.374,.386,.395,.395]
    dmix3=[.159,.159,.174,.213,.213]
    curfrc=fltarr(3)
    ; interpolate beam current fraction for deuterium
    curfrc[0]=interpol(dmix1,denergy,einj)
    curfrc[1]=interpol(dmix2,denergy,einj)
    curfrc[2]=interpol(dmix3,denergy,einj)
    ;;ffull[i]=curfrc[0] & fhalf[i]=curfrc[1] & fthird[i]=curfrc[2]
    ;; FIDASIM uses current fractions   
    ffracs=curfrc[0] & hfracs=curfrc[1] & tfracs=1.0-ffracs-hfracs
        
    ;;SAVE IN NBI STRUCTURE
    nbi={einj:einj,pinj:pinj,full:ffracs,half:hfracs,third:tfracs,$
         xyz_src:reform(xyz_src[inputs.isource,*]),xyz_pos:reform(xyz_pos[inputs.isource,*]),$
         bmwidra:bmwidra[inputs.isource],bmwidza:bmwidza[inputs.isource],$
         divy:divy[inputs.isource],divz:divz[inputs.isource],$
         focy:focy[inputs.isource],focz:focz[inputs.isource] }
    
    ;plot geometry to check!
    if keyword_set(doplot) then begin
        window, 0, retain=2, xsize=800, ysize=800
        ;xrange=[-175,175] & yrange=[-175,175]
        xrange=[-180,180] & yrange=[-180,180]
        ;xrange=[-125,175] & yrange=[0,175]
        ;;aspect_ratio.pro from http://www.idlcoyote.com/programs/aspect.pro
        aspect_ratio=aspect((yrange[1]-yrange[0])/float(xrange[1]-xrange[0]) )
        npt=101
        phi = findgen(npt)/(npt-1.0)*2*!pi
        u = (findgen(npt)-npt*0.5)/(npt*0.5)*200.
        plot, 170.*cos(phi), 170.*sin(phi),xtitle='!6X[cm]',ytitle='!6Y[cm]',$
        xrange=xrange,yrange=yrange,position=aspect_ratio,/xstyle,/ystyle
        rc=rc_onaxis
        phic=phic_onaxis
        oplot, rc*cos(phi), rc*sin(phi),lineStyle=2
        oplot, [rc*cos(phic), rc*cos(phic)],[rc*sin(phic), rc*sin(phic)],psym=1
        oplot,[0.,0.],[0.,0.],psym=1
        for i=0,n_elements(us_NB)-1 do begin
            oplot, u, vs_NB[i] + (u-us_NB[i])*tan(alpha_NB[i]),thick=1.5
            ;for j=0L,npt-1 do begin
            ;print, 'x',u[j], ' v',vs_NB[i] + (u[j]-us_NB[i])*tan(alpha_NB[i]),' r',$
            ;       sqrt(u[j]^2.+ (vs_NB[i] + (u[j]-us_NB[i])*tan(alpha_NB[i]))^2.)
            ;endfor        
        endfor
        x=nbi.xyz_src[0]+(nbi.xyz_pos[0]-nbi.xyz_src[0])*findgen(npt)/(npt-1.0)
        y=nbi.xyz_src[1]+(nbi.xyz_pos[1]-nbi.xyz_src[1])*findgen(npt)/(npt-1.0)
        ;oplot,x,y,color=!blue,thick=2   
    endif
    
    return,nbi

end
