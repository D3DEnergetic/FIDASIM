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

    isource=inputs.isource
    ; machine coordinate (u,v,z): typical N-S-E-W type definitions of machine hall
   
    ; Geometrical factors from TRANSP namelist
    ;-------------------------------
    ;  NB        1A       1B      1C
    XBZETA_on=[-85.5,    -82 ,   -78.5] ;toroidal angle of beam source in (R,zeta,Z)
    XLBAPA_on=[510.3,    509 ,   510.3] ;distance of ion source to beam aperture
    XLBTNA_on=[1134.,   1138.,   1142.] ;distance of ion source to beam tangency radius
    RTCENA_on=[69.4 ,   59.2 ,   49.7 ] ;beam tangency radius [cm]
    XYBSCA_on=[0.0  ,    0.0 ,    0.0 ] ;Elevation above/below vacuum vessel midplane of center of beam source grid [cm]
    XYBAPA_on=[0.0  ,    0.0 ,    0.0 ] ;Elevation above/below vacuum vessel midplane of beam centerline at aperture [cm]
    NLCO_on  =[1.0  ,    1.0 ,    1.0 ] ;1 for co-beams, -1 for counter-beams
 
    ;off axis NBI
    ;  NB          2A       2B      2C
    XBZETA_off=[-39.2,   -35.7 ,   -32.2] ;toroidal angle of beam source in (R,zeta,Z)
    XLBAPA_off=[510.3,    510.3,   510.3] ;distance of ion source to beam aperture
    XLBTNA_off=[1146.,    1146.,   1146.] ;distance of ion source to beam tangency radius
    RTCENA_off=[130. ,     120.,    110.] ;beam tangency radius
    XYBSCA_off=[0.0  ,     0.0 ,    0.0 ] ;Elevation above/below vacuum vessel midplane of center of beam source grid [cm]
    XYBAPA_off=[0.0  ,     0.0 ,    0.0 ] ;Elevation above/below vacuum vessel midplane of beam centerline at aperture [cm]
    NLCO_off  =[1.0  ,     1.0 ,    1.0 ] ;1 for co-beams, -1 for counter-beams

    XBZETA=[XBZETA_on,XBZETA_off]
    XLBAPA=[XLBAPA_on,XLBAPA_off]
    XLBTNA=[XLBTNA_on,XLBTNA_off]
    RTCENA=[RTCENA_on,RTCENA_off]
    XYBSCA=[XYBSCA_on,XYBSCA_off]
    XYBAPA=[XYBAPA_on,XYBAPA_off]
    NLCO  =[NLCO_on,NLCO_off]

    nsources=n_elements(XBZETA)

    bmwidra=replicate(6.0,nsources)
    bmwidza=replicate(21.5,nsources)   ;array of ion source half height in cm
    divy=replicate(4.94d-3,nsources)
    divz=replicate(1.36d-2,nsources)
    focy=replicate(988.0,nsources)
    focz=replicate(988.0,nsources)
    
    transp_geometry,RTCENA[isource],XLBAPA[isource],XLBTNA[isource],XBZETA[isource],$
                    XYBAPA[isource],XYBSCA[isource],NLCO[isource],xyz_src,xyz_pos,angle=145

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
         xyz_src:xyz_src,xyz_pos:xyz_pos,$
         bmwidra:bmwidra[isource],bmwidza:bmwidza[isource],$
         divy:divy[isource],divz:divz[isource],$
         focy:focy[isource],focz:focz[isource] }

    return,nbi

end
