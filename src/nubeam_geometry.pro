FUNCTION nubeam_geometry, nubeam, angle=angle, verbose=verbose,plot=plot
    ;+#nubeam_geometry
    ;+Calculates the FIDASIM beam geometry from the beam geometry variables in the TRANSP/NUBEAM namelist
    ;+***
    ;+##Arguments
    ;+     **NUBEAM**: Structure containing the following
    ;+
    ;+     **NUBEAM.NAME**: Ion source name
    ;+
    ;+     **NUBEAM.NBSHAP**: Ion source shape 1=rectangular, 2=circular
    ;+
    ;+     **NUBEAM.FOCLZ**: Vertical focal length [cm]
    ;+
    ;+     **NUBEAM.FOCLR**: Horizontal focal length [cm]
    ;+
    ;+     **NUBEAM.DIVZ**: Vertical divergence [rad]
    ;+
    ;+     **NUBEAM.DIVR**: Horizontal divergence [rad]
    ;+
    ;+     **NUBEAM.BMWIDZ**: Ion source half height [cm]
    ;+
    ;+     **NUBEAM.BMWIDR**: Ion source half width [cm]
    ;+
    ;+     **NUBEAM.RTCENA**: Radius of tangency point [cm]
    ;+
    ;+     **NUBEAM.XLBTNA**: Distance from center of beam source grid to tangency point [cm]
    ;+
    ;+     **NUBEAM.XBZETA**: Torodial angle [deg] Positive angles defined to be in the counter-clockwise direction
    ;+
    ;+     **NUBEAM.XYBSCA**: Elevation above/below vacuum vessel midplane of center of beam source grid [cm]
    ;+
    ;+     **NUBEAM.NLCO**: 1 for Co-beam, 0 or -1 for Counter-beam
    ;+
    ;+     **NUBEAM.NBAPSHA**: Vector of aperture shapes 1=rectangular, 2=circular
    ;+
    ;+     **NUBEAM.XLBAPA**: Vector of distances from center of beam source grid to the aperture plane [cm]
    ;+
    ;+     **NUBEAM.XYBAPA**: Vector of elevation above/below vacuum vessel midplane of beam centerline at aperture [cm]
    ;+
    ;+     **NUBEAM.RAPEDGA**: Vector of aperture half-widths [cm]
    ;+
    ;+     **NUBEAM.XZPEDGA**: Vector of aperture half-heights [cm]
    ;+
    ;+     **NUBEAM.XRAPOFFA**: Vector of horizontal (y) offsets relative to the +x aligned beam centerline [cm]
    ;+
    ;+     **NUBEAM.XZAPOFFA**: Vector of vertical (z) offsets relative to the +x aligned beam centerline [cm]
    ;+
    ;+##Keyword Arguments
    ;+     **angle**: Angle to add to XBZETA to rotate the beams into correct coordinates [deg]
    ;+
    ;+     **verbose**: Print out positions
    ;+
    ;+     **plot**: Plot the beam
    ;+
    ;+##Return Value
    ;+ Neutral beam structure
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> nbi = nubeam_geometry(nubeam)
    ;+```

    if not keyword_set(angle) then angle=0
    if nubeam.nlco eq 0 then nubeam.nlco = -1

    phi_s=(nubeam.XBZETA+angle)*!DPI/180.0
    zs=nubeam.XYBSCA
    za=nubeam.XYBAPA[0]
    alpha=asin((zs-za)/nubeam.XLBAPA[0])
    pdst=nubeam.XLBTNA*cos(alpha)
    rs=sqrt(nubeam.RTCENA^2.0 + pdst^2.0)
    dat=nubeam.XLBTNA-nubeam.XLBAPA[0]
    pdat=dat*cos(alpha)
    ra=sqrt(nubeam.RTCENA^2.0 + pdat^2.0)
    beta_s=acos(nubeam.RTCENA/rs)
    beta_a=acos(nubeam.RTCENA/ra)
    phi_a=phi_s+nubeam.NLCO*(beta_s-beta_a)

    src=[rs*cos(phi_s),rs*sin(phi_s),zs]
    aper_src=[ra*cos(phi_a),ra*sin(phi_a),za]
    axis=(aper_src-src)
    axis=axis/sqrt(total(axis^2))
    pos=src+axis*nubeam.XLBTNA
  
    if keyword_set(verbose) then begin
        print,'Source position: ',src
        print,'1st Aperture position: ',aper_src
        print,'Tangency position: ', pos
    endif

    if keyword_set(plot) then begin
        plot,[src[0],pos[0]],[src[1],pos[1]]
    endif

    nbi = {data_source:"TRANSP/NUBEAM namelist",name:nubeam.name, $
           shape:nubeam.nbshap,src:src,axis:axis, $
           focy:double(nubeam.foclr),focz:double(nubeam.foclz), $
           divy:replicate(double(nubeam.divr),3), $
           divz:replicate(double(nubeam.divz),3), $
           widy:double(nubeam.bmwidr), widz:double(nubeam.bmwidz) $
           naperture:n_elements(nubeam.nbapsha),ashape:nubeam.nbapsha, $
           awidy:double(nubeam.rapedga),awidz:double(nubeam.xzpedga), $
           aoffy:double(nubeam.rapoffa),aoffz:double(nubeam.xzaoffa), $
           adist:double(nubeam.xlbapa) }
           } 
    
    return, nbi
END
