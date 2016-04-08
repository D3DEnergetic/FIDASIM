PRO transp_geometry,RTCENA,XLBAPA,XLBTNA,XBZETA,XYBAPA,XYBSCA,NLCO,src,axis,angle=angle,plot=plot
    ;+#transp_geometry
    ;+Calculates src and axis from beam geometry variables in the TRANSP namelist
    ;+***
    ;+##Arguments
    ;+     **RTCENA**: Radius of tangency point [cm]
    ;+
    ;+     **XLBAPA**: Distance from center of beam source grid to aperture [cm]
    ;+
    ;+     **XLBTNA**: Distance from center of beam source grid to tangency point
    ;+
    ;+     **XBZETA**: Torodial angle [deg] Positive angles defined to be in the counter-clockwise direction
    ;+
    ;+     **XYBAPA**: Elevation above/below vacuum vessel midplane of beam centerline at aperture [cm]
    ;+
    ;+     **XYBSCA**: Elevation above/below vacuum vessel midplane of center of beam source grid [cm]
    ;+
    ;+     **NLCO**: 1 for Co-beam, -1 for Counter-beam
    ;+
    ;+##Keyword Arguments
    ;+     **angle**: Angle to add to XBZETA to rotate the beams into correct coordinates Ex. 90 for D3D,142.243 for NSTX
    ;+
    ;+     **plot**: Plot the beam
    ;+
    ;+##Outputs:
    ;+     **src**: Beam source position in machine coordinates [cm]
    ;+     **axis**: Direction vector
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> transp_geometry, RTCENA, XLBAPA, XLBTNA, XBZETA, XYBAPA, XYBSCA, NLCO, src, axis, angle=142.243
    ;+```

    if not keyword_set(angle) then angle=0
    if n_elements(RTCENA) gt 1 then begin
        print,'TRANSP_GEOMETRY does not take array arguements'
        goto,GET_OUT
    endif

    phi_s=(XBZETA+angle)*!DPI/180.0
    zs=XYBSCA
    za=XYBAPA
    alpha=asin((zs-za)/XLBAPA)
    pdst=XLBTNA*cos(alpha)
    rs=sqrt(RTCENA^2.0 + pdst^2.0)
    dat=XLBTNA-XLBAPA
    pdat=dat*cos(alpha)
    ra=sqrt(RTCENA^2.0 + pdat^2.0)
    beta_s=acos(RTCENA/rs)
    beta_a=acos(RTCENA/ra)
    phi_a=phi_s+NLCO*(beta_s-beta_a)

    src=[rs*cos(phi_s),rs*sin(phi_s),zs]
    aper_src=[ra*cos(phi_a),ra*sin(phi_a),za]
    axis=(aper_src-transp_src)
    axis=axis/sqrt(total(axis^2))
    pos=src+dir*(dat/2.0)

    if keyword_set(plot) then begin
        plot,[src[0],pos[0]],[src[1],pos[1]],$
             xrange=1.05*[-aper_src[0],aper_src[0]],$
             yrange=1.05*[-aper_src[1],aper_src[1]]
    endif
    GET_OUT:
END
