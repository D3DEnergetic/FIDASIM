PRO detector_aperture_geometry,g,wn,rdist,zdist,v,d,rc
    ;+#detector_aperture_geometry
    ;+ Defines 3 MeV proton diagnostic geometry
    ;+ Source: detector_aperture_geometry.dat and MAST_29908_221_g_ensemble.nml
    ;+***
    ;+##Arguments
    ;+    **g**: GEQDSK file
    ;+
    ;+    **wn**: indicates Werner (1) or Netepneko (0) diagnostic geometry
    ;+
    ;+##Outputs
    ;+    **rdist**: Detector-aperture geometry parameter
    ;+
    ;+    **zdist**: Detector-aperture geometry parameter
    ;+
    ;+    **v**: Detector-aperture geometry parameter
    ;+
    ;+    **d**: Detector-aperture geometry parameter
    ;+
    ;+    **rc**: Detector-aperture geometry parameter
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> g = 'g99999K26'
    ;+IDL> detector_aperture_geometry,g,0,rdist,zdist,v,d,rc
    ;+```

    if wn eq 1 then begin
        info, 'Using MAST geometry data 1 (Werner) and flipping the direction of B poloidal'
        ; Convert Alexander's input data file into IDL

        ;# detectors orientation
        theta_port=fltarr(4)
            theta_port(1-1)= 40.0
            theta_port(2-1)= 40.0
            theta_port(3-1)= 40.0
            theta_port(4-1)= 40.0
        theta_port*=!pi/180

        ;In MAST_29908_221_g_ensemble.nml, phi_port < 0, but it is positive here for
        ;comparison with Netepenko's data
        phi_port=fltarr(4)
            phi_port(1-1)=  30.0
            phi_port(2-1)=  35.0
            phi_port(3-1)=  40.0
            phi_port(4-1)=  45.0
        phi_port*=!pi/180

        v=fltarr(3,4)
        ;    # radial direction
            v(0,*) = - sin(theta_port)*cos(phi_port)

        ;    # toroidal direction
            v(1,*) = - sin(phi_port)

        ;    # vertical direction
            v(2,*) = cos(theta_port)*cos(phi_port)

        ;&detection (m)

        ;    #detectors position

        ;    #vertical coordinate
        ZDist=fltarr(4)
            ZDist(1-1)= 0.030014
            ZDist(2-1)= 0.038311
            ZDist(3-1)= 0.013981
            ZDist(4-1)= 0.024414

        ;    #radial coordinate
        RDist=fltarr(4)
            RDist(1-1)= 1.668328
            RDist(2-1)= 1.661375
            RDist(3-1)= 1.648754
            RDist(4-1)= 1.64

        ;    #toroidal angle
        PHDangle=fltarr(4)
            PHDangle(1-1)= 1.38277674641
            PHDangle(2-1)= 1.39294543453
            PHDangle(3-1)= 1.37658569418
            PHDangle(4-1)= 1.3855051501

            D = 0.0357      ;                       ! detector-collimator spacing (m)

        RC=fltarr(4)
            RC(1-1) = 2.5e-3  ;                     ! outer collimator radius (m)
            RC(2-1) = 2.5e-3
            RC(3-1) = 2.5e-3
            RC(4-1) = 2.5e-3

        RCD=fltarr(4)
            RCD(1-1) = 2.5e-3 ;             ! inner collimator radius (detector radius) (m)
            RCD(2-1) = 2.5e-3
            RCD(3-1) = 2.5e-3
            RCD(4-1) = 2.5e-3
    endif else begin
        info, 'Using MAST geometry data 0 (Netepenko)'
        ; Convert Alexander's input data file into IDL

        ;# detectors orientation
        theta_port=fltarr(4)
            theta_port(1-1)= 46.6875370324
            theta_port(2-1)= 47.6648339458
            theta_port(3-1)= 50.031360382
            theta_port(4-1)= 51.5837100275
        theta_port*=!pi/180

        phi_port=fltarr(4)
            phi_port(1-1)= 38.8895519328
            phi_port(2-1)= 44.2509710868
            phi_port(3-1)= 48.3160975078
            phi_port(4-1)= 53.6006103875
        phi_port*=!pi/180

        v=fltarr(3,4)
        ;    # radial direction
            v(0,*) = - sin(theta_port)*cos(phi_port)

        ;    # toroidal direction
            v(1,*) = - sin(phi_port)

        ;    # vertical direction
            v(2,*) = cos(theta_port)*cos(phi_port)

        ;&detection (m)

        ;    #detectors position

        ;    #vertical coordinate
        ZDist=fltarr(4)
            ZDist(1-1)= 0.036701870576
            ZDist(2-1)= 0.0409530530375
            ZDist(3-1)= 0.0232888146809
            ZDist(4-1)= 0.0301116448993

        ;    #radial coordinate
        RDist=fltarr(4)
            RDist(1-1)= 1.66830563971
            RDist(2-1)= 1.67508495819
            RDist(3-1)= 1.68813419386
            RDist(4-1)= 1.69658132076

        ;    #toroidal angle
        PHDangle=fltarr(4)
            PHDangle(1-1)= 79.8774705956
            PHDangle(2-1)= 79.2421358615
            PHDangle(3-1)= 80.3144072462
            PHDangle(4-1)= 79.7395308235
        PHDangle*=!pi/180

            D = 0.04      ;                       ! detector-collimator spacing (m)

        RC=fltarr(4)
            RC(1-1) = 1.288e-3  ;                     ! outer collimator radius (m)
            RC(2-1) = 1.294e-3
            RC(3-1) = 1.318e-3
            RC(4-1) = 1.343e-3

        RCD=fltarr(4)
            RCD(1-1) = 1.288e-3 ;                    ! inner collimator radius (detector radius) (m)
            RCD(2-1) = 1.294e-3
            RCD(3-1) = 1.318e-3
            RCD(4-1) = 1.343e-3
    endelse

    ;orb_mast,g,[rdist[0],0,zdist[0]],-reform(v[*,0])

    ; Looks qualitatively like Ramona's figures but not quantitatively
    ; Get midplane crossings of [1.015,1.055,1.134,1.188] for the four
    ; channels.  Seems reasonable.

end
