PRO detector_aperture_geometry,g,wn,rdist,zdist,v,d,rc
    ;+#detector_aperture_geometry
    ;+ Defines charged fusion product diagnostic geometry
    ;+ Sources: detector_aperture_geometry.dat and MAST_29908_221_g_ensemble.nml
    ;+ NOTE: This file is hardcoded to be useful for developers. General users could modify parameters under
    ;+       wn = 0 or 1, and correspondingly set the wn flag when calling detector_aperture_geometry.pro
    ;+       Alternatively, the user could manually define the geometry, see run_tests.pro and test_cfpd.pro
    ;+***
    ;+##Arguments
    ;+    **g**: GEQDSK structure
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
    ;+IDL> g = readg('g000001.01000')
    ;+IDL> detector_aperture_geometry,g,0,rdist,zdist,v,d,rc
    ;+```

    if wn eq 1 then begin ; Source: MAST_29908_221_g_ensemble.nml
        print, 'Using MAST geometry 1 (Werner)'
        ;; Detector orientation
        theta_port=replicate(40.d0, 4) * !dpi/180

        ;In MAST_29908_221_g_ensemble.nml, phi_port < 0, but it is made positive here for
        ;comparison with detector_aperture_geometry.dat
        phi_port = [30.d0, 35.d0, 40.d0, 45.d0] * !pi/180

        v=fltarr(3,4)
        v[0,*] = - sin(theta_port)*cos(phi_port);    radial direction
        v[1,*] = - sin(phi_port)                ;    toroidal direction
        v[2,*] = cos(theta_port)*cos(phi_port)  ;    vertical direction

        ;; Detector position
        ZDist = [0.030014, 0.038311, 0.013981, 0.024414]

        RDist = [1.668328, 1.661375, 1.648754, 1.64]

        PHDangle = [1.38277674641, 1.39294543453, 1.37658569418, 1.3855051501]

        D = 0.0357      ;                        detector-collimator spacing (m)

        RC = replicate(2.5d-3, 4) ; outer collimator radius (m)

        RCD = replicate(2.5d-3, 4) ; inner collimator radius (m)

    endif else begin ; Source: detector_aperture_geometry.dat
        print, 'Using MAST geometry 0 (Netepenko)'

        ;; Detector orientation
        theta_port = [46.6875370324, 47.6648339458, 50.031360382, 51.5837100275] * !pi/180

        phi_port = [38.8895519328, 44.2509710868, 48.3160975078, 53.6006103875] * !pi/180

        v=fltarr(3,4)
        v[0,*] = - sin(theta_port)*cos(phi_port);    radial direction
        v[1,*] = - sin(phi_port)                ;    toroidal direction
        v[2,*] = cos(theta_port)*cos(phi_port)  ;    vertical direction

        ;; Detector position
        ZDist = [0.036701870576, 0.0409530530375, 0.0232888146809, 0.0301116448993]

        RDist = [1.66830563971, 1.67508495819, 1.68813419386, 1.69658132076]

        PHDangle = [79.8774705956, 79.2421358615, 80.3144072462, 79.7395308235] * !pi/180

        D = 0.04      ;                        detector-collimator spacing (m)

        RC = [1.288d-3, 1.294d-3, 1.318d-3, 1.343d-3]; outer collimator radius (m)

        RCD = [1.288d-3, 1.294d-3, 1.318d-3, 1.343d-3]; inner collimator radius (detector radius) (m)

    endelse

end
