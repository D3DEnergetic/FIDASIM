PRO detector_aperture_geometry,g,wn,rdist,zdist,v,d,rc
    ;+#detector_aperture_geometry
    ;+ Defines charged fusion product diagnostic geometry
    ;+ Sources: detector_aperture_geometry.dat and MAST_29908_221_g_ensemble.nml
    ;+ NOTE: This file is hardcoded to be useful for developers. General users could modify parameters under
    ;+       wn = 0 or 1, and correspondingly set the wn flag when calling detector_aperture_geometry.pro
    ;+       Alternatively, the user could manually define the geometry, see run_tests.pro and test_cfpd.pro
    ;+***
    ;+##Arguments
    ;+    **g**: GEQDSK file
    ;+
    ;+    **wn**: indicates Werner (1) or Netepneko (0) diagnostic geometry
    ;+
    ;+##Outputs
    ;+    **rdist**: Radial coordinates of detector [m]
    ;+
    ;+    **zdist**: Vertical coordinates of detector [m]
    ;+
    ;+    **v**: Detector orientations [m/s, rad/s, m/s]
    ;+
    ;+    **d**: Collimator length [m]
    ;+
    ;+    **rc**: Outer collimator spacing [m]
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> g = 'g000001.01000'
    ;+IDL> detector_aperture_geometry,g,0,rdist,zdist,v,d,rc
    ;+```

    if wn eq 1 then begin ; Source: MAST_29908_221_g_ensemble.nml
        print, 'Using MAST geometry 1 (Werner)'
        ;; Detector orientation
        theta_port=fltarr(4)
            theta_port(1-1)= 40.0
            theta_port(2-1)= 40.0
            theta_port(3-1)= 40.0
            theta_port(4-1)= 40.0
        theta_port*=!pi/180

        ;In MAST_29908_221_g_ensemble.nml, phi_port < 0, but it is made positive here for
        ;comparison with detector_aperture_geometry.dat
        phi_port=fltarr(4)
            phi_port(1-1)=  30.0
            phi_port(2-1)=  35.0
            phi_port(3-1)=  40.0
            phi_port(4-1)=  45.0
        phi_port*=!pi/180

        v=fltarr(3,4)
            v(0,*) = - sin(theta_port)*cos(phi_port);    radial direction
            v(1,*) = - sin(phi_port)                ;    toroidal direction
            v(2,*) = cos(theta_port)*cos(phi_port)  ;    vertical direction

        ;; Detector position
        ZDist=fltarr(4)         ;    vertical coordinate
            ZDist(1-1)= 0.030014
            ZDist(2-1)= 0.038311
            ZDist(3-1)= 0.013981
            ZDist(4-1)= 0.024414

        RDist=fltarr(4)         ;    radial coordinate
            RDist(1-1)= 1.668328
            RDist(2-1)= 1.661375
            RDist(3-1)= 1.648754
            RDist(4-1)= 1.64

        PHDangle=fltarr(4)      ;    toroidal angle
            PHDangle(1-1)= 1.38277674641
            PHDangle(2-1)= 1.39294543453
            PHDangle(3-1)= 1.37658569418
            PHDangle(4-1)= 1.3855051501

        D = 0.0357      ;                        detector-collimator spacing (m)

        RC=fltarr(4)
            RC(1-1) = 2.5e-3  ;                      outer collimator radius (m)
            RC(2-1) = 2.5e-3
            RC(3-1) = 2.5e-3
            RC(4-1) = 2.5e-3

        RCD=fltarr(4)
            RCD(1-1) = 2.5e-3 ;              inner collimator radius (detector radius) (m)
            RCD(2-1) = 2.5e-3
            RCD(3-1) = 2.5e-3
            RCD(4-1) = 2.5e-3

    endif else begin ; Source: detector_aperture_geometry.dat
        print, 'Using MAST geometry 0 (Netepenko)'

        ;; Detector orientation
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
            v(0,*) = - sin(theta_port)*cos(phi_port);    radial direction
            v(1,*) = - sin(phi_port)                ;    toroidal direction
            v(2,*) = cos(theta_port)*cos(phi_port)  ;    vertical direction

        ;; Detector position
        ZDist=fltarr(4)     ;    vertical coordinate
            ZDist(1-1)= 0.036701870576
            ZDist(2-1)= 0.0409530530375
            ZDist(3-1)= 0.0232888146809
            ZDist(4-1)= 0.0301116448993

        RDist=fltarr(4)     ;    radial coordinate
            RDist(1-1)= 1.66830563971
            RDist(2-1)= 1.67508495819
            RDist(3-1)= 1.68813419386
            RDist(4-1)= 1.69658132076

        PHDangle=fltarr(4)  ;    toroidal angle
            PHDangle(1-1)= 79.8774705956
            PHDangle(2-1)= 79.2421358615
            PHDangle(3-1)= 80.3144072462
            PHDangle(4-1)= 79.7395308235
        PHDangle*=!pi/180

        D = 0.04      ;                        detector-collimator spacing (m)

        RC=fltarr(4)
            RC(1-1) = 1.288e-3  ;                      outer collimator radius (m)
            RC(2-1) = 1.294e-3
            RC(3-1) = 1.318e-3
            RC(4-1) = 1.343e-3

        RCD=fltarr(4)
            RCD(1-1) = 1.288e-3 ;                     inner collimator radius (detector radius) (m)
            RCD(2-1) = 1.294e-3
            RCD(3-1) = 1.318e-3
            RCD(4-1) = 1.343e-3

    endelse

end
