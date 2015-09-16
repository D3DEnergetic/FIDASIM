FUNCTION nstx_beams,inputs,doplot=doplot,ape=ape
    ; Returns the nbi structure for ONE neutral beam source
    ; D. Liu 6/26/2015, updated for NSTX-U NBI line#1 and line #2
    ;                   use conventional  N-S-E-W type definitions of machine hall
    ;IDL>nstx_input,inputs
    ;IDL>nbi=nstx_beams(inputs)
    ;IDL>help,nbi
    ;** Structure <1eb88468>, 15 tags, length=128, data length=124, refs=1:
    ;EINJ            DOUBLE           90.300003
    ;PINJ            DOUBLE           2.0000000
    ;FULL            FLOAT          0.453095
    ;HALF            FLOAT          0.380180
    ;THIRD           FLOAT          0.166725
    ;XYZ_SRC         FLOAT     Array[3]
    ;XYZ_POS         FLOAT     Array[3]
    ;BMWIDRA         FLOAT           6.00000
    ;BMWIDZA         FLOAT           21.5000
    ;EDGERA          FLOAT           11.0000
    ;EDGEZA          FLOAT           30.2500
    ;DIVY            DOUBLE    Array[3]
    ;DIVZ            DOUBLE    Array[3]
    ;FOCY            FLOAT           988.000
    ;FOCZ            FLOAT           988.000
     
    ; select one specific neutral beam source    
    isource=inputs.isource
    
    ;NB line #1
    rc_line1=158.58             ;major radius at crossover point of NB line #1
    phic_line1=79.298/180.*!pi  ;toroidal angle at crossover point
    uc_line1=rc_line1*cos(phic_line1)
    vc_line1=rc_line1*sin(phic_line1)
    
    ;NB line #2
    rc_line2=188.973             ;major radius at crossover point of NB line #2
    phic_line2=139.921/180.*!pi  ;toroidal angle at crossover point   
    uc_line2=rc_line2*cos(phic_line2)
    vc_line2=rc_line2*sin(phic_line2)
    print,'Machine coordinates of cross-over point of NBI line #1 is '
    print,uc_line1, vc_line1, 0.
    print,'Machine coordinates of cross-over point of NBI line #2 is '
    print,uc_line2, vc_line2, 0.
   
    ; Geometrical factors from TRANSP namelist
    ;NB line #1
    ;-------------------------------
    ;  NB        1A         1B      1C
    XBZETA_Line1=[ 56.91,   60.38,  63.85] ;toroidal angle of beam source in (R,zeta,Z)
    XLBAPA_Line1=[ 510.3,  510.3,  510.3 ] ;distance of ion source to beam aperture
    XLBTNA_Line1=[1130.7, 1135.2, 1138.9 ] ;distance of ion source to beam tangency radius
    RTCENA_Line1=[69.3 ,   59.1 ,   48.7 ] ;beam tangency radius [cm]
    XYBSCA_Line1=[0.0  ,    0.0 ,    0.0 ] ;Elevation above/below vacuum vessel midplane of center of beam source grid [cm]
    XYBAPA_Line1=[0.0  ,    0.0 ,    0.0 ] ;Elevation above/below vacuum vessel midplane of beam centerline at aperture [cm]
    NLCO_Line1  =[1.0  ,    1.0 ,    1.0 ] ;1 for counter-clockwise, -1 for clockwise
 
    ;NB line #2
    ;  NB          2A       2B      2C
    XBZETA_Line2=[103.09,  106.54 , 109.98  ] ;toroidal angle of beam source in (R,zeta,Z)
    XLBAPA_Line2=[510.3,   510.3,   510.3   ] ;distance of ion source to beam aperture
    XLBTNA_Line2=[1125.3, 1134.0,  1142.0   ] ;distance of ion source to beam tangency radius
    RTCENA_Line2=[129.9,   120.0,   109.5   ] ;beam tangency radius
    XYBSCA_Line2=[0.0  ,    0.0 ,    0.0 ] ;Elevation above/below vacuum vessel midplane of center of beam source grid [cm]
    XYBAPA_Line2=[0.0  ,    0.0 ,    0.0 ] ;Elevation above/below vacuum vessel midplane of beam centerline at aperture [cm]
    NLCO_Line2  =[1.0  ,    1.0 ,    1.0 ] ;1 for counter-clockwise, -1 for clockwise

    XBZETA=[XBZETA_Line1,XBZETA_Line2]
    XLBAPA=[XLBAPA_Line1,XLBAPA_Line2]
    XLBTNA=[XLBTNA_Line1,XLBTNA_Line2]
    RTCENA=[RTCENA_Line1,RTCENA_Line2]
    XYBSCA=[XYBSCA_Line1,XYBSCA_Line2]
    XYBAPA=[XYBAPA_Line1,XYBAPA_Line2]
    NLCO  =[NLCO_Line1,NLCO_Line2]

    nsources=n_elements(XBZETA)

    bmwidra=replicate(6.0,nsources)      ;array of ion source half width in cm
    bmwidza=replicate(21.5,nsources)     ;array of ion source half height in cm
    edgera=replicate(11.0,nsources)      ;array of beam aperture half width in cm
    edgeza=replicate(30.25,nsources)     ;array of beam aperture half height in cm
    
    divy=replicate(4.94d-3,3,nsources)   ;horizontal divergence in radians
    divz=replicate(1.36d-2,3,nsources)   ;vertical divergence in radians
    focy=replicate(988.0,nsources)       ;horizontal focal length in cm
    focz=replicate(988.0,nsources)       ;vertical focal length in cm

    angle_transp_uvz=0.0                   ;for NSTX-U TRANSP namelsit, in degrees
    print,'Angle between TRANSP and machine coordiantes',angle_transp_uvz
    transp_geometry,RTCENA,XLBAPA,XLBTNA,XBZETA,$
                    XYBAPA,XYBSCA,NLCO,xyz_src,xyz_ape,xyz_mid,angle=angle_transp_uvz
    print,'TRANSP NB ion source location:'
    for i=0,nsources-1 do begin
        print,i,xyz_src[0,i],xyz_src[1,i],xyz_src[2,i]
    endfor
    print,'TRANSP aperture location:'
    for i=0,nsources-1 do begin
        print,i,xyz_ape[0,i],xyz_ape[1,i],xyz_ape[2,i]
    endfor      

    ;ape=0;1: launch beam neutrals at beam aperture, 0: launch beam neutrals at ion source  
    if keyword_set(ape) then begin
        xyz_pos=xyz_mid
    xyz_src=xyz_ape    
    endif else begin
        xyz_pos=xyz_ape
        xyz_src=xyz_src
    endelse
            
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
    ;; FIDASIM uses current fractions   
    ffracs=curfrc[0] & hfracs=curfrc[1] & tfracs=1.0-ffracs-hfracs
    print,'current fractions:',ffracs,hfracs,tfracs
    print,'power fractions:',[ffracs,hfracs/2.,tfracs/3.]/$
                              total(ffracs+hfracs/2.+tfracs/3.)        
    
    ;;SAVE IN NBI STRUCTURE
    nbi={einj:einj,pinj:pinj,full:ffracs,half:hfracs,third:tfracs,$
         xyz_src:xyz_src[*,isource],xyz_pos:xyz_pos[*,isource],$
         bmwidra:bmwidra[isource],bmwidza:bmwidza[isource],$
     edgera:edgera[isource],edgeza:edgeza[isource],$
         divy:divy[*,isource],divz:divz[*,isource],$
         focy:focy[isource],focz:focz[isource] }     
           
     
   if keyword_set(doplot) then begin
      if !D.Name ne 'PS' then window, /free, retain=2, xsize=800, ysize=800
      xrange=[-210,210] & yrange=[-210,210]
      ;;aspect_ratio.pro from http://www.idlcoyote.com/programs/aspect.pro
      aspect_ratio=aspect((yrange[1]-yrange[0])/float(xrange[1]-xrange[0]) )
      npt=501
      phi = findgen(npt)/(npt-1.0)*2*!pi
      u = (findgen(npt)-npt*0.5)/(npt*0.5)*200.
      plot, 170.*cos(phi), 170.*sin(phi),xtitle='!6X[cm]',ytitle='!6Y[cm]',$
      xrange=xrange,yrange=yrange,position=aspect_ratio,/xstyle,/ystyle

      raxis=100.0
      rlcfs=145.
      
      rc=rc_line1      
      phic=phic_line1            
      oplot, rc*cos(phi), rc*sin(phi),lineStyle=2
      oplot, [rc*cos(phic), rc*cos(phic)],[rc*sin(phic), rc*sin(phic)],psym=1,color=!darkgreen
      
      rc=rc_line2
      phic=phic_line2            
      oplot, rc*cos(phi), rc*sin(phi),lineStyle=2
      oplot, [rc*cos(phic), rc*cos(phic)],[rc*sin(phic), rc*sin(phic)],psym=1,color=!darkgreen
      
      oplot,[0.,0.],[0.,0.],psym=1
      oplot, 30.*cos(phi), 30.*sin(phi)
      oplot, rlcfs*cos(phi), rlcfs*sin(phi),color=!ltblue
      oplot, raxis*cos(phi),raxis*sin(phi),color=!ltblue
   
      for i=0,nsources-1 do begin
          dir=reform(xyz_pos[*,i]-xyz_src[*,i])/$
              sqrt(total((xyz_pos[*,i]-xyz_src[*,i])^2.))
      k=(xyz_pos[1,i]-xyz_src[1,i])/(xyz_pos[0,i]-xyz_src[0,i])    
          angle=atan((xyz_src[1,i]-xyz_pos[1,i]),(xyz_src[0,i]-xyz_pos[0,i]))
      print,'Angle of NB line relative to East',angle/!pi*180.
      oplot,u,xyz_src[1,i]+(u-xyz_src[0,i])*k
      endfor
      
      ;plot norminal NB source 1B centerline based on R. Bell's notes
      k=(49.4706*2.54-36.9588*2.54)/(4*2.54-(-4)*2.54)
      angle_NB_1B=atan(k); 57.4053 degree /180. * !pi
      ds=0.001
      xmax=300.
      nx=round(xmax/ds)
      xbeam=-4.*2.54+(-nx/2.+lindgen(nx))*ds
      ybeam=36.9588*2.54+(-nx/2.+lindgen(nx))*ds*k
      w=where(ybeam gt 0.)
      dum=min((xbeam^2.+ybeam^2.),w)
      x_rtan=xbeam[w]
      y_rtan=ybeam[w]
   
      print,'Norminal Rtan for source 1B',sqrt(x_rtan^2.+y_rtan^2.)
      print,'Norminal angle for source 1B',angle_NB_1B/!pi*180.
      oplot,xbeam,ybeam,color=!blue

      ;plot norminal NB source 2B centerline
      angle_NB_2B=100.5/180.*!dpi
      k=tan(angle_NB_2B)
      ds=0.001
      xmax=300.
      nx=round(xmax/ds)
      xbeam=-144.594+(-nx/2.+lindgen(nx))*ds
      ybeam=121.669+(-nx/2.+lindgen(nx))*ds*k
      w=where(ybeam gt 0.)
      dum=min((xbeam^2.+ybeam^2.),w)
      x_rtan=xbeam[w]
      y_rtan=ybeam[w]
      print,'Norminal Rtan for source 2B',sqrt(x_rtan^2.+y_rtan^2.)
      print,'Norminal angle for source 2B',angle_NB_1B/!pi*180.
      oplot,xbeam,ybeam,color=!blue      
      
   endif        
   
   return,nbi

end
