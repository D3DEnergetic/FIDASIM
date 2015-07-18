FUNCTION nstx_chords,shot,fida_diag,isource=isource,doplot=doplot,year=year
    
    ;; fida structure (15 == number of chords/channels)
    ;;** Structure <88d87f8>, 9 tags, length=800, data length=792, refs=1:
    ;;   SIGMA_PI_RATIO  DOUBLE          0.90000000 ;;COULD BE ARRAY
    ;;   NCHAN           LONG                15
    ;;   XLOS            DOUBLE    Array[15]
    ;;   YLOS            DOUBLE    Array[15]
    ;;   ZLOS            DOUBLE    Array[15]
    ;;   XHEAD           DOUBLE    Array[15]
    ;;   YHEAD           DOUBLE    Array[15]
    ;;   ZHEAD           DOUBLE    Array[15]
    ;;   ra              FLOAT     Array[15]
    ;;   rd              FLOAT     Array[15]
    ;;   h               FLOAT     Array[15]
    
    compile_opt defint32,strictarr,strictarrsubs
    
    ;;**********************************************************************
    ;; Vertial FIDA system, refer to fida_mario.pro in IDL FIDAsim4.0
    ;; Detector#1 at BAY A, 10 active tangential views
    ;; inner views (i.e. "outer" lens) at Bay A
    ;; (includes 10 fibers on a 200mm f/1.8 lens, diameter 600um)
    ;;**********************************************************************
    ;;lens and fiber parameter saved for refrence          
    ;;det_area = !pi*10.0^2
    ;;det_normal=-[0.0876, 0.1131, 0.9896]
    ;;fnumber=1.8
    ;;dfiber = replicate(0.6, 10)
    ;;dspot=[0.7082, 0.6596, 0.668, 0.6832, 0.6604,$
    ;;       0.6766, 0.6784, 0.6684, 0.6784, 0.6386]
    
    if n_elements(year) eq 0 then year='2008_new'
    case strlowcase(year) of
    '2008_old':begin
    ;;detector location
    detR_vfida1=127.794
    detZ_vfida1=191.757
    detphi_vfida1=1.4869  ; angle already in [rad] !!
    
    detp_vfida1=[detR_vfida1*cos(detphi_vfida1), $
                 detR_vfida1*sin(detphi_vfida1), $
             detZ_vfida1]
    
    ;intersection of fiber view with midplane
    rmid_vfida1=[  85.81,     89.99,     94.28,     98.72,  103.26,$
                  107.87,    112.51,    117.23,    122.04,  126.84]  
    xmid_vfida1=[ -20.34,    -17.17,    -14.12,    -10.93,   -7.80, $
                   -4.76,     -1.74,      1.33,      4.32,    7.31]
    ymid_vfida1=[  83.36,     88.34,     93.22,     98.11,  102.96, $
                  107.77,    112.50,    117.22,    121.96,  126.63]
    zmid_vfida1=replicate(0.,n_elements(xmid_vfida1))
    end
    '2008_new':begin; calculated by liu_compute_fida_midplane_position.pro
    ;;detector location
    detR_vfida1=127.130
    detZ_vfida1=188.075
    detphi_vfida1=1.48446  ; angle already in [rad] !!
    
    detp_vfida1=[detR_vfida1*cos(detphi_vfida1), $
                 detR_vfida1*sin(detphi_vfida1), $
                 detZ_vfida1]
    
    ;intersection of fiber view with midplane
    rmid_vfida1=[  86.298,    90.387,    94.697,    99.085, 103.570,$
                  108.134,   112.728,   117.441,   122.204, 126.995] 
    xmid_vfida1=[ -20.850,   -17.633,   -14.550,   -11.337,  -8.214,$
                   -5.123,    -2.118,     0.963,     3.981,   6.972]
    ymid_vfida1=[  83.741,    88.651,    93.572,    98.434, 103.244,$
                  108.012,   112.708,   117.437,   122.140, 126.803]
    zmid_vfida1=replicate(0.,n_elements(xmid_vfida1))
    end
    else:
    endcase
    
    
    ;; Detector#2, outer views (i.e. "inner" lens) at Bay A
    ;; (includes 6 fibers on a 135mm f/1.8 lens, diameter 800um)
    ;;lens and fiber parameter saved for refrence
    ;;det_area=!pi*6.75^2
    ;;det_normal=-[-0.1071, -0.1868, 0.9764]
    ;;fnumber=1.8
    ;;dfiber=replicate(0.8, 6)
    ;;dspot=[0.4058, 0.4270, 0.3722, 0.4122, 0.3952, 0.4066]
    case strlowcase(year) of
    '2008_old':begin
    ;;detector location
    detR_vfida2=104.350
    detZ_vfida2=196.716
    detphi_vfida2=1.5978  ; angle already in [rad] !!
    
    detp_vfida2=[detR_vfida2*cos(detphi_vfida2),$
                 detR_vfida2*sin(detphi_vfida2),$ 
                 detZ_vfida2]
    
    ;intersection of fiber view with midplane
    rmid_vfida2=[130.782, 135.567, 140.383, 145.357, 150.371, 155.429]
    xmid_vfida2=[ 10.82,   13.69,   16.58,   19.49,   22.45,   25.42 ]
    ymid_vfida2=[ 130.33, 134.87,  139.40,  144.05,  148.69,   153.34]
    zmid_vfida2=replicate(0.,n_elements(xmid_vfida2))
    end
    '2008_new':begin; calculated by liu_compute_fida_midplane_position.pro
    ;;detector location
    detR_vfida2=103.790
    detZ_vfida2=200.097
    detphi_vfida2=1.60961  ; angle already in [rad] !!
    
    detp_vfida2=[detR_vfida2*cos(detphi_vfida2),$
                 detR_vfida2*sin(detphi_vfida2),$ 
                 detZ_vfida2]
    
    ;intersection of fiber view with midplane
    rmid_vfida2=[130.848, 135.613, 140.445, 145.425, 150.422, 155.496]
    xmid_vfida2=[ 10.540,  13.406,  16.301,  19.243,  22.207,  25.181]
    ymid_vfida2=[130.423, 134.949, 139.496, 144.146, 148.774 ,153.443]
    zmid_vfida2=replicate(0.,n_elements(xmid_vfida2))
    end
    else:
    endcase
    
    
    detp_active_vfida=[detp_vfida1,detp_vfida2]
    rmid_active_vfida=[rmid_vfida1,rmid_vfida2]  
    xmid_active_vfida=[xmid_vfida1,xmid_vfida2]
    ymid_active_vfida=[ymid_vfida1,ymid_vfida2]
    zmid_active_vfida=[zmid_vfida1,zmid_vfida2]
    
    ;;**********************************************************************
    ;; Detector#3 at BAY B, 10 active tangential views
    ;; inner views (i.e. "outer" lens) at Bay B
    ;; (includes 10 fibers on a 200mm f/1.8 lens, diameter 600um)
    ;; Refer to fida_mario.pro in IDL FIDAsim4.0
    ;;**********************************************************************
    ;;lens and fiber parameter saved for refrence          
    ;;det_area = !pi*10.0^2
    ;;det_normal=-[0.0942, 0.0730, 0.9928]
    ;;fnumber=1.8
    ;;dfiber = replicate(0.6, 10)
    ;;dspot=[0.6420, 0.6458, 0.685, 0.6864, 0.6952,
    ;;       0.7448, 0.7664, 0.8298, 0.7838, 0.811]
    
    ;;detector location
    detR_vfida3=126.2458
    detZ_vfida3=183.587
    detphi_vfida3=0.7834  ; angle already in [rad] !!
    
    detp_vfida3=[detR_vfida3*cos(detphi_vfida3), $
                 detR_vfida3*sin(detphi_vfida3), $
                 detZ_vfida3]
    
    ;intersection of fiber view with midplane
    rmid_vfida3=[ 83.148,    87.895,  92.160,  97.410, 102.084,$
                 106.747,   111.391, 116.014, 120.607, 125.211]
    xmid_vfida3=[  57.09,    60.42,   63.73,   67.09,  70.41, $
                   73.65,    76.91,   80.15,   83.38,  86.62 ]
    ymid_vfida3=[  60.45,    63.84,   67.19,   70.63,  73.92,$
                   77.27,    80.57,   83.87,   87.14,  90.42]
    zmid_vfida3=replicate(0.,n_elements(xmid_vfida3))
    
    ;; Detector#4, outer views (i.e. "inner" lens) at Bay B
    ;; (includes 6 fibers on a 135mm f/1.8 lens, diameter 800um)
    ;;det_area = !pi*6.75^2
    ;;det_normal=[-0.4662, 0.0147, 0.8845]
    ;;fnumber=1.8
    ;;dfiber = replicate(0.8, 6)
    ;;dspot=[0.4502, 0.47, 0.4514, 0.435, 0.5622, 0.4744]
    ;;detector location
    detR_vfida4=103.8
    detZ_vfida4=223.209
    detphi_vfida4=0.7845  ; angle already in [rad] !!
    
    detp_vfida4=[detR_vfida4*cos(detphi_vfida4),$
                 detR_vfida4*sin(detphi_vfida4),$ 
                 detZ_vfida4]
    
    ;intersection of fiber view with midplane
    rmid_vfida4=[130.151, 135.032, 139.943, 144.841, 149.825, 154.827]
    xmid_vfida4=[ 91.69,   95.09,   98.54,  101.94,  105.40,  108.87 ]
    ymid_vfida4=[ 92.37,   95.87,   99.37,  102.90,  106.49,  110.09]
    zmid_vfida4=replicate(0.,n_elements(xmid_vfida4))
    
    detp_passive_vfida=[detp_vfida3,detp_vfida4]
    rmid_passive_vfida=[rmid_vfida3,rmid_vfida4]     
    xmid_passive_vfida=[xmid_vfida3,xmid_vfida4]
    ymid_passive_vfida=[ymid_vfida3,ymid_vfida4]
    zmid_passive_vfida=[zmid_vfida3,zmid_vfida4]
    
    ;;**********************************************************************
    ;; Tangential FIDA system, refer to fida_alessandro.pro in IDL FIDAsim4.0
    ;; Detectors at BAY L, 16 active tangential views
    ;; (includes 16 fibers on a 20mm f/1.7 lens, diameter 600um)
    ;;**********************************************************************
    ;;lens and fiber parameter saved for refrence          
    ;det_area = !pi*0.5882^2
    ;det_normal=[0.6125,   -0.5278,   -0.5884]
    ;fnumber=1.7
    ;dfiber = replicate(0.6, 16)
    ;dspot=2.5 * replicate(2,16)
    
    ;;detector location
    detR_active_tfida=168.91
    detZ_active_tfida=48.26
    detphi_active_tfida=1.8326  ; angle already in [rad] !!
    detp_active_tfida=[detR_active_tfida*cos(detphi_active_tfida),$
                       detR_active_tfida*sin(detphi_active_tfida),$
                       detZ_active_tfida]
    
    ;intersection of fiber view with midplane
    rmid_active_tfida=[ 0.8500,    0.8967,    0.9433,    0.9900, $
                        1.0367,    1.0833,    1.1300,    1.1767, $
                        1.2233,    1.2700,    1.3167,    1.3633, $
                        1.4100,    1.4567,    1.5033,    1.5500]*100     
    xmid_active_tfida=[-0.1698,   -0.1354,   -0.1024,   -0.0705, $
                       -0.0394,   -0.0089,    0.0209,    0.0504, $
                        0.0794,    0.1080,    0.1364,    0.1646, $
                        0.1925,    0.2202,    0.2477,    0.2751]*100
    ymid_active_tfida=[ 0.8329,    0.8864,    0.9378,    0.9875, $
                        1.0359,    1.0833,    1.1298,    1.1756, $
                        1.2208,    1.2654,    1.3096,    1.3534, $
                        1.3968,    1.4399,    1.4828,    1.5254]*100
    zmid_active_tfida=replicate(0.,n_elements(xmid_active_tfida))
    
    
    ;;**********************************************************************
    ;; Detectors at BAY F, 16 passive tangential views
    ;; (includes 16 fibers on a 20mm f/1.7 lens, diameter 600um)
    ;; Refer to fida_alessandro.pro in IDL FIDAsim4.0
    ;;**********************************************************************
    ;lens and fiber parameter saved for refrence           
    ;det_area = !pi*0.5882^2
    ;det_normal=[ -0.5722,    0.5712,   -0.5884]
    ;fnumber=1.7
    ;dfiber = replicate(0.6, 16)
    ;dspot=2.5 * replicate(2,16)
    
    ;;detector location
    detR_passive_tfida=168.91
    detZ_passive_tfida=48.26
    detphi_passive_tfida=-1.3823  ; angle already in [rad] !!
    detp_passive_tfida=[detR_passive_tfida*cos(detphi_passive_tfida),$
                        detR_passive_tfida*sin(detphi_passive_tfida),$
                        detZ_passive_tfida]
    
    ;intersection of fiber view with midplane
    xmid_passive_tfida=[ 0.8500,    0.8967,    0.9433,    0.9900, $
                         1.0367,    1.0833,    1.1300,    1.1767, $
                         1.2233,    1.2700,    1.3167,    1.3633, $
                         1.4100,    1.4567,    1.5033,    1.5500]*100
    xmid_passive_tfida=[ 0.1083,    0.0701,    0.0334,   -0.0021, $
                        -0.0367,   -0.0705,   -0.1037,   -0.1364, $
                        -0.1686,   -0.2005,   -0.2320,   -0.2633, $
                        -0.2943,   -0.3251,   -0.3557,   -0.3861]*100
    ymid_passive_tfida=[-0.8431,   -0.8939,   -0.9427,   -0.9900, $
                        -1.0360,   -1.0810,   -1.1252,   -1.1687, $
                        -1.2117,   -1.2541,   -1.2961,   -1.3377, $
                        -1.3789,   -1.4199,   -1.4607,   -1.5011]*100
    zmid_passive_tfida=replicate(0.,n_elements(xmid_passive_tfida)) 
    
    
    
    ;;from TRANSP namelist
    ;;NPA/SSNPA settings
    ;;CXRSTA=221.7, 221.7, 221.7, 221.7, 221.7, 183.6, 183.6, 183.6, 183.6
    ;;                              ; RADIUS FROM MACHINE CENTERLINE TO CX PIVOT
    ;;CXYSTA=9*0.0          ; ELEVATION OF PIVOT ABOVE MIDPLANE
    ;;CXZETA=0.0,   0.0,   0.0,   0.0,   0.0,   30.,   30.,    30.,   30.
    ;;                              ; TOROIDAL LOCATION OF PIVOT IN DEGREES
    ;;CXYTAN=9*0.0          ; HEIGHT ABOVE MIDPLANE OF PIVOT
    ;;CXRTAN=120.,  100.,  90.,   80.,    65., 120.73, 101.63, 88.52, 63.8
    ;;                              ; CX RTAN (+CO, -CTR)
    
    ;the angle difference of TRANSP coordinates and FIDA (u,v,z) coordinates
    angle_transp_uvz=0.0/180.*!pi ;check angle_transp_uvz within nstxu_beams.pro
    
    CXRSTA_nstx= [replicate(221.7,4),replicate(183.6,8)]
    nch_nstx_npa=n_elements(CXRSTA_nstx)
    CXYSTA_nstx=replicate(0.,nch_nstx_npa)
    CXZETA_nstx= [replicate(0.,4),   replicate(30.,8)]/180.*!pi+142.38/180.*!pi
    CXTHETA_nstx=[replicate(0.,4),   replicate(0.,8)]/180./!pi
    Rtan_nstx_npa= [70., 80., 100., 120., 63.8,  88.5,   101.6, 110.,$
                                          63.8,  88.5,   101.6, 110. ]  
                            
    npap_nstx=fltarr(nch_nstx_npa,3)         ;detector center position in (u,v,z) (cm)
    npa_mid_nstx=fltarr(nch_nstx_npa,3)      ;intersection of NPA and midplane in (u,v,z) (cm)
    det_radius_nstx=fltarr(nch_nstx_npa)     ;radius of NPA dector (cm)  
    
    npa_area_nstx=0.04d
    npa_l_nstx=26.d
    npa_aperture_nstx=sqrt(npa_area_nstx/!pi)
    
    ;;ssnpa flight tube geometry
    ;;ssnpa_l: flight tube length
    ;;ssnpa_hole: radius of at input side flight tube
    ;;ssnpa_det:  effective detector radius (which equals the aperture radius) 
    ssnpa_l_nstx=[     25.4,      25.4,       25.4,        25.4      ] 
    ssnpa_hole_nstx=[0.3/2.,    0.27/2.,    0.31/2.,     0.27/2.     ]
    ssnpa_det_nstx=[ 0.005/2.,  0.0035/2.,  0.0035/2.,   0.0025/2.   ] 
    
    det_radius_nstx=[replicate(sqrt(npa_area_nstx/!pi),4),   $   ;E||B NPA    
                     ssnpa_det_nstx,$                            ;actual SSNPA 
                     5.0d,  5.0d,    5.0d,   5.0d]               ;SSNPA in TRANSP
    
    aperture_radius_nstx=[replicate(sqrt(npa_area_nstx/!pi),4),   $   ;E||B NPA    
                          ssnpa_hole_nstx,$                           ;actual SSNPA 
                      5.d,  5.d,    5.0d,   5.0d]                 ;SSNPA in TRANSP
    
    dcollimator_nstx=[replicate(sqrt(npa_l_nstx),4),$
                      ssnpa_l_nstx,$
                      160.d, 160.d,200.d,200.d]
    
    npa_solid_angle_nstx=!pi*aperture_radius_nstx^2/dcollimator_nstx^2.
      
    for iDet=0,nch_nstx_npa-1 do begin
        ;detector position in uvz coordinate
        npap_nstx[iDet,*]=[CXRSTA_nstx[iDet]*cos(CXZETA_nstx[iDet]),$
                           CXRSTA_nstx[iDet]*sin(CXZETA_nstx[iDet]),$
                           CXYSTA_nstx[iDet] ]  
        ;intersection of npa detector and midplane in uvz coordinate
        angle=CXZETA_nstx[iDet]-sign(Rtan_nstx_npa[iDet])*$
                                acos(abs(Rtan_nstx_npa[iDet])/CXRSTA_nstx[iDet])     
        npa_mid_nstx[iDet,*]=[abs(Rtan_nstx_npa[iDet])*cos(angle),$
                              abs(Rtan_nstx_npa[iDet])*sin(angle),$
                              tan(CXTHETA_nstx[iDet])*$
                              sqrt(CXRSTA_nstx[iDet]^2.-Rtan_nstx_npa[iDet]^2.)+$
                              CXYSTA_nstx[iDet] ] 
    endfor
             
    
    ;New SSNPA at Bay I (16 cm below midplane)
    ;CXRSTA=5*188.0                   ! RADIUS FROM MACHINE CENTERLINE OF CX PIVOT
    ;CXYSTA=5*0.0             ! ELEVATION OF PIVOT ABOVE MIDPLANE
    ;CXZETA=5*61.7                    ! TOROIDAL LOCATION OF PIVOT IN DEGREES
    ;CXRTAN=140., 130.,120.,110.,100. ! CX RTAN (+CO, -CTR)
    ;CXTHEA=5*0.0                 ! ANGLE OF CHORD ABOVE MIDPLANE IN DEGREES
    ;CXYTAN=5*0.0                     ! HEIGHT ABOVE MIDPLANE OF TANGENCY RADIUS
    CXRTAN_tssnpa=[145., 142., 139., 136.,$
                   133., 130., 127., 124.,$
                   121., 118., 115., 112.,$
                   109., 106., 103., 100.,97.]  ;! CX RTAN (+CO, -CTR)
    nch_tssnpa=n_elements(CXRTAN_tssnpa)
    
    
    CXRSTA_tssnpa=replicate(188.0,nch_tssnpa)        ;! RADIUS FROM MACHINE CENTERLINE OF CX PIVOT
    CXYSTA_tssnpa=replicate(0.0,nch_tssnpa)          ;! ELEVATION OF PIVOT ABOVE MIDPLANE
    ;! TOROIDAL LOCATION OF PIVOT IN DEGREES
    CXZETA_tssnpa=replicate(61.7,nch_tssnpa)/180*!pi+142.38/180.*!pi
    CXTHEA_tssnpa=replicate(0.0,nch_tssnpa)      ;! ANGLE OF CHORD ABOVE MIDPLANE IN DEGREES
    CXYTAN_tssnpa=replicate(0.0,nch_tssnpa)          ;! HEIGHT ABOVE MIDPLANE OF TANGENCY RADIUS
    CXTHETA_tssnpa=replicate(0.0,nch_tssnpa)  
    
    tssnpap=fltarr(nch_tssnpa,3)             ;detector center position in (u,v,z) (cm)
    tssnpa_mid=fltarr(nch_tssnpa,3)          ;intersection of NPA and midplane in (u,v,z) (cm)
    det_radius_tssnpa=fltarr(nch_tssnpa)     ;radius of NPA dector (cm)  
    
    for iDet=0,nch_tssnpa-1 do begin
        ;detector position in uvz coordinate
        tssnpap[iDet,*]=[CXRSTA_tssnpa[iDet]*cos(CXZETA_tssnpa[iDet]),$
                         CXRSTA_tssnpa[iDet]*sin(CXZETA_tssnpa[iDet]),$
                         CXYSTA_tssnpa[iDet] ]  
        ;intersection of npa detector and midplane in uvz coordinate
        angle=CXZETA_tssnpa[iDet]-sign(CXRTAN_tssnpa[iDet])*$
                                  acos(abs(CXRTAN_tssnpa[iDet])/CXRSTA_tssnpa[iDet]) 
        tssnpa_mid[iDet,*]=[abs(CXRTAN_tssnpa[iDet])*cos(angle),$
                            abs(CXRTAN_tssnpa[iDet])*sin(angle),$
                            tan(CXTHETA_tssnpa[iDet])*$
                            sqrt(CXRSTA_tssnpa[iDet]^2.-CXRTAN_tssnpa[iDet]^2.)+$
                            CXYSTA_tssnpa[iDet] ] 
    endfor
    tssnpa_area=0.01d
    tssnpa_l=2.54d*2.75
    tssnpa_aperture=sqrt(0.03*0.085/!pi)
    
    det_radius_tssnpa=[replicate(sqrt(tssnpa_area/!pi),nch_tssnpa)] 
    aperture_radius_tssnpa=[replicate(sqrt(tssnpa_area/!pi),nch_tssnpa)]
    dcollimator_tssnpa=[replicate(sqrt(tssnpa_l),nch_tssnpa)]     
    
    tssnpa_solid_angle=!pi*aperture_radius_tssnpa^2/dcollimator_tssnpa^2.   
    
    
    ;New SSNPA at Bay L (+2.5cm above midplance)
    ;CXRSTA=5*203.                  ! RADIUS FROM MACHINE CENTERLINE OF CX PIVOT
    ;CXYSTA=5*0.0           ! ELEVATION OF PIVOT ABOVE MIDPLANE
    ;CXZETA=5*-39.4                 ! TOROIDAL LOCATION OF PIVOT IN DEGREES
    ;CXRTAN=90.,80.,70.,60.,50.     ! CX RTAN (+CO, -CTR)
    ;CXTHEA=5*0.0               ! ANGLE OF CHORD ABOVE MIDPLANE IN DEGREES
    ;CXYTAN=5*0.0                   ! HEIGHT ABOVE MIDPLANE OF TANGENCY RADIUS
    CXRTAN_rssnpa=[90., 85., 80., 75.,$
                   70., 65., 60., 55.,$
                   50., 45., 40., 35.,$
                   30., 25., 20., 15.,$
                   10., 5.]  
    nch_rssnpa=n_elements(CXRTAN_rssnpa)
    CXRSTA_rssnpa=replicate(203.0,nch_rssnpa)        ;! RADIUS FROM MACHINE CENTERLINE OF CX PIVOT
    CXYSTA_rssnpa=replicate(0.0,nch_rssnpa)          ;! ELEVATION OF PIVOT ABOVE MIDPLANE
    CXZETA_rssnpa=replicate(-39.4,nch_rssnpa)/180*!pi+142.38/180.*!pi
                                                     ;! TOROIDAL LOCATION OF PIVOT IN DEGREES
                  ;! CX RTAN (+CO, -CTR)
    CXTHEA_rssnpa=replicate(0.0,nch_rssnpa)      ;! ANGLE OF CHORD ABOVE MIDPLANE IN DEGREES
    CXYTAN_rssnpa=replicate(0.0,nch_rssnpa)          ;! HEIGHT ABOVE MIDPLANE OF TANGENCY RADIUS
    CXTHETA_rssnpa=replicate(0.0,nch_rssnpa)  
    
    rssnpap=fltarr(nch_rssnpa,3)             ;detector center position in (u,v,z) (cm)
    rssnpa_mid=fltarr(nch_rssnpa,3)          ;intersection of NPA and midplane in (u,v,z) (cm)
    det_radius_rssnpa=fltarr(nch_rssnpa)     ;radius of NPA dector (cm)  
    
    for iDet=0,nch_rssnpa-1 do begin
        ;detector position in uvz coordinate
        rssnpap[iDet,*]=[CXRSTA_rssnpa[iDet]*cos(CXZETA_rssnpa[iDet]),$
                         CXRSTA_rssnpa[iDet]*sin(CXZETA_rssnpa[iDet]),$
                         CXYSTA_rssnpa[iDet] ]  
        ;intersection of npa detector and midplane in uvz coordinate
        angle=CXZETA_rssnpa[iDet]-sign(CXRTAN_rssnpa[iDet])*$
                                  acos(abs(CXRTAN_rssnpa[iDet])/CXRSTA_rssnpa[iDet]) 
        rssnpa_mid[iDet,*]=[abs(CXRTAN_rssnpa[iDet])*cos(angle),$
                            abs(CXRTAN_rssnpa[iDet])*sin(angle),$
                            tan(CXTHETA_rssnpa[iDet])*$
                            sqrt(CXRSTA_rssnpa[iDet]^2.-CXRTAN_rssnpa[iDet]^2.)+$
                            CXYSTA_rssnpa[iDet] ] 
    endfor
    rssnpa_area=0.01d
    rssnpa_l=2.54d*2.35
    rssnpa_aperture=sqrt(0.03*0.085/!pi)
    
    det_radius_rssnpa=[replicate(sqrt(rssnpa_area/!pi),nch_rssnpa)] 
    aperture_radius_rssnpa=[replicate(sqrt(rssnpa_area/!pi),nch_rssnpa)]
    dcollimator_rssnpa=[replicate(sqrt(rssnpa_l),nch_rssnpa)]     
    
    rssnpa_solid_angle=!pi*aperture_radius_rssnpa^2/dcollimator_rssnpa^2. 
    
    
    ;New SSNPA at Bay B (+13.3cm above midplane)
    ;CXRSTA=5*170.                   ! RADIUS FROM MACHINE CENTERLINE OF CX PIVOT
    ;CXYSTA=5*0.0                ! ELEVATION OF PIVOT ABOVE MIDPLANE
    ;CXZETA=5*-88.6                  ! TOROIDAL LOCATION OF PIVOT IN DEGREES
    ;CXRTAN=40.,25.,10.,0.,-10.      ! CX RTAN (+CO, -CTR)
    ;CXTHEA=5*0.0                ! ANGLE OF CHORD ABOVE MIDPLANE IN DEGREES
    ;CXYTAN=5*0.0                    ! HEIGHT ABOVE MIDPLANE OF TANGENCY RADIUS
    CXRSTA_pssnpa=replicate(170.0,5)                   ;! RADIUS FROM MACHINE CENTERLINE OF CX PIVOT
    CXYSTA_pssnpa=replicate(0.0,5)             ;! ELEVATION OF PIVOT ABOVE MIDPLANE
    CXZETA_pssnpa=replicate(-88.6,5)/180*!pi+142.38/180.*!pi
                                                       ;! TOROIDAL LOCATION OF PIVOT IN DEGREES
    CXRTAN_pssnpa=[40.,25.,10.,0.,-10.]                ;! CX RTAN (+CO, -CTR)
    CXTHEA_pssnpa=replicate(0.0,5)                 ;! ANGLE OF CHORD ABOVE MIDPLANE IN DEGREES
    CXYTAN_pssnpa=replicate(0.0,5)                     ;! HEIGHT ABOVE MIDPLANE OF TANGENCY RADIUS
    CXTHETA_pssnpa=replicate(0.0,5)  
    nch_pssnpa=n_elements(CXRTAN_pssnpa)
    pssnpap=fltarr(nch_pssnpa,3)             ;detector center position in (u,v,z) (cm)
    pssnpa_mid=fltarr(nch_pssnpa,3)          ;intersection of NPA and midplane in (u,v,z) (cm)
    det_radius_pssnpa=fltarr(nch_pssnpa)     ;radius of NPA dector (cm)  
    
    for iDet=0,nch_pssnpa-1 do begin
        ;detector position in uvz coordinate
        pssnpap[iDet,*]=[CXRSTA_pssnpa[iDet]*cos(CXZETA_pssnpa[iDet]),$
                         CXRSTA_pssnpa[iDet]*sin(CXZETA_pssnpa[iDet]),$
                         CXYSTA_pssnpa[iDet] ]  
        ;intersection of npa detector and midplane in uvz coordinate
        angle=CXZETA_pssnpa[iDet]-sign(CXRTAN_pssnpa[iDet])*$
                                  acos(abs(CXRTAN_pssnpa[iDet])/CXRSTA_pssnpa[iDet]) 
        pssnpa_mid[iDet,*]=[abs(CXRTAN_pssnpa[iDet])*cos(angle),$
                            abs(CXRTAN_pssnpa[iDet])*sin(angle),$
                            tan(CXTHETA_pssnpa[iDet])*$
                            sqrt(CXRSTA_pssnpa[iDet]^2.-CXRTAN_pssnpa[iDet]^2.)+$
                            CXYSTA_pssnpa[iDet] ] 
    endfor
    pssnpa_area=0.01d
    pssnpa_l=2.54d*3.0
    pssnpa_aperture=sqrt(pssnpa_area/!pi)
    
    det_radius_pssnpa=[replicate(sqrt(pssnpa_area/!pi),5)] 
    aperture_radius_pssnpa=[replicate(sqrt(pssnpa_area/!pi),5)]
    dcollimator_pssnpa=[replicate(sqrt(pssnpa_l),5)]     
    
    pssnpa_solid_angle=!pi*aperture_radius_pssnpa^2/dcollimator_pssnpa^2.
    
    
    
    ;CXRSTA=5*170.                   ! RADIUS FROM MACHINE CENTERLINE OF CX PIVOT
    ;CXYSTA=5*0.0                ! ELEVATION OF PIVOT ABOVE MIDPLANE
    ;CXZETA=5*-88.6                  ! TOROIDAL LOCATION OF PIVOT IN DEGREES
    ;CXRTAN=40.,25.,10.,0.,-10.      ! CX RTAN (+CO, -CTR)
    ;CXTHEA=5*0.0                ! ANGLE OF CHORD ABOVE MIDPLANE IN DEGREES
    ;CXYTAN=5*0.0                    ! HEIGHT ABOVE MIDPLANE OF TANGENCY RADIUS
    CXRSTA_nstxu_npa=replicate(190.0,5)                   ;! RADIUS FROM MACHINE CENTERLINE OF CX PIVOT
    CXYSTA_nstxu_npa=replicate(0.0,5)              ;! ELEVATION OF PIVOT ABOVE MIDPLANE
    CXZETA_nstxu_npa=replicate(8.0,5)/180*!pi+142.38/180.*!pi 
                                                       ;! TOROIDAL LOCATION OF PIVOT IN DEGREES
    CXRTAN_nstxu_npa=[120.,100.,90.,80.,70.]                ;! CX RTAN (+CO, -CTR)
    CXTHEA_nstxu_npa=replicate(0.0,5)                  ;! ANGLE OF CHORD ABOVE MIDPLANE IN DEGREES
    CXYTAN_nstxu_npa=replicate(0.0,5)                     ;! HEIGHT ABOVE MIDPLANE OF TANGENCY RADIUS
    CXTHETA_nstxu_npa=replicate(0.0,5)  
    nch_nstxu_npa=n_elements(CXRTAN_nstxu_npa)
    nstxu_npap=fltarr(nch_nstxu_npa,3)             ;detector center position in (u,v,z) (cm)
    nstxu_npa_mid=fltarr(nch_nstxu_npa,3)          ;intersection of NPA and midplane in (u,v,z) (cm)
    det_radius_nstxu_npa=fltarr(nch_nstxu_npa)     ;radius of NPA dector (cm)  
    
    for iDet=0,nch_nstxu_npa-1 do begin
        ;detector position in uvz coordinate
        nstxu_npap[iDet,*]=[CXRSTA_nstxu_npa[iDet]*cos(CXZETA_nstxu_npa[iDet]),$
                            CXRSTA_nstxu_npa[iDet]*sin(CXZETA_nstxu_npa[iDet]),$
                            CXYSTA_nstxu_npa[iDet] ]  
        ;intersection of npa detector and midplane in uvz coordinate
        angle=CXZETA_nstxu_npa[iDet]-sign(CXRTAN_nstxu_npa[iDet])*$
                                  acos(abs(CXRTAN_nstxu_npa[iDet])/CXRSTA_nstxu_npa[iDet]) 
        nstxu_npa_mid[iDet,*]=[abs(CXRTAN_nstxu_npa[iDet])*cos(angle),$
                               abs(CXRTAN_nstxu_npa[iDet])*sin(angle),$
                               tan(CXTHETA_nstxu_npa[iDet])*$
                               sqrt(CXRSTA_nstxu_npa[iDet]^2.-CXRTAN_nstxu_npa[iDet]^2.)+$
                               CXYSTA_nstxu_npa[iDet] ] 
    endfor
    nstxu_npa_area=0.04d
    nstxu_npa_l=26.d
    nstxu_npa_aperture=sqrt(nstxu_npa_area/!pi)
    
    det_radius_nstxu_npa=[replicate(sqrt(nstxu_npa_area/!pi),5)] 
    aperture_radius_nstxu_npa=[replicate(sqrt(nstxu_npa_area/!pi),5)]
    dcollimator_nstxu_npa=[replicate(sqrt(nstxu_npa_l),5)]     
    
    nstxu_npa_solid_angle=!pi*aperture_radius_nstxu_npa^2/dcollimator_nstxu_npa^2.
    
                        
    ;; SELECT THE VIEWS
    for i=0,n_elements(fida_diag)-1 do begin
        CASE strupcase(fida_diag[i]) OF
            'ACTIVE_TFIDA': begin
             ;coordinate of intersection of this fiber view with midplane
             xlens_active_tfida=replicate(detp_active_tfida[0],n_elements(xmid_active_tfida))
             ylens_active_tfida=replicate(detp_active_tfida[1],n_elements(xmid_active_tfida)) 
             zlens_active_tfida=replicate(detp_active_tfida[2],n_elements(xmid_active_tfida))    
        
             xlos=xmid_active_tfida
             ylos=ymid_active_tfida
             zlos=zmid_active_tfida
            
             xhead=xlens_active_tfida
             yhead=ylens_active_tfida
             zhead=zlens_active_tfida
            
             nchan=n_elements(xlos)
             
             sigma_pi=replicate(1.d0,nchan)
             ra=replicate(0.d0,nchan)
             rd=replicate(0.d0,nchan)
             h=replicate(0.d0,nchan)
             chan_id=replicate(0.d0,nchan) ;;0 for fida chords
            end
            'PASSIVE_TFIDA': begin
             ;coordinate of intersection of this fiber view with midplane
             xlens_passive_tfida=replicate(detp_passive_tfida[0],n_elements(xmid_passive_tfida))
             ylens_passive_tfida=replicate(detp_passive_tfida[1],n_elements(xmid_passive_tfida)) 
             zlens_passive_tfida=replicate(detp_passive_tfida[2],n_elements(xmid_passive_tfida)) 
        
             xlos=xmid_passive_tfida
             ylos=ymid_passive_tfida
             zlos=zmid_passive_tfida
            
             xhead=xlens_passive_tfida
             yhead=ylens_passive_tfida
             zhead=zlens_passive_tfida
            
             nchan=n_elements(xlos)
             
             sigma_pi=replicate(1.d0,nchan)
             ra=replicate(0.d0,nchan)
             rd=replicate(0.d0,nchan)
             h=replicate(0.d0,nchan)
             chan_id=replicate(0.d0,nchan) ;;0 for fida chords
            end
            'ACTIVE_VFIDA': begin
             ;coordinate of intersection of this fiber view with midplane
             xlens_active_vfida=[replicate(detp_active_vfida[0],n_elements(xmid_vfida1)),$
                         replicate(detp_active_vfida[3],n_elements(xmid_vfida2))]                      
             ylens_active_vfida=[replicate(detp_active_vfida[1],n_elements(xmid_vfida1)),$
                         replicate(detp_active_vfida[4],n_elements(xmid_vfida2))]        
             zlens_active_vfida=[replicate(detp_active_vfida[2],n_elements(xmid_vfida1)),$
                         replicate(detp_active_vfida[5],n_elements(xmid_vfida2))]    
        
             xlos=xmid_active_vfida
             ylos=ymid_active_vfida
             zlos=zmid_active_vfida
         
             xhead=xlens_active_vfida
             yhead=ylens_active_vfida
             zhead=zlens_active_vfida
            
             nchan=n_elements(xlos)
            
             sigma_pi=replicate(1.d0,nchan)
             ra=replicate(0.d0,nchan)
             rd=replicate(0.d0,nchan)
             h=replicate(0.d0,nchan)
             chan_id=replicate(0.d0,nchan) ;;0 for fida chords
            end
            'PASSIVE_VFIDA': begin
             ;coordinate of intersection of this fiber view with midplane
             xlens_passive_vfida=[replicate(detp_passive_vfida[0],n_elements(xmid_vfida3)),$
                          replicate(detp_passive_vfida[3],n_elements(xmid_vfida4))]                    
             ylens_passive_vfida=[replicate(detp_passive_vfida[1],n_elements(xmid_vfida3)),$
                          replicate(detp_passive_vfida[4],n_elements(xmid_vfida4))]      
             zlens_passive_vfida=[replicate(detp_passive_vfida[2],n_elements(xmid_vfida3)),$
                          replicate(detp_passive_vfida[5],n_elements(xmid_vfida4))]
        
             xlos=xmid_passive_vfida
             ylos=ymid_passive_vfida
             zlos=zmid_passive_vfida
            
             xhead=xlens_passive_vfida
             yhead=ylens_passive_vfida
             zhead=zlens_passive_vfida
            
             nchan=n_elements(xlos)
            
             sigma_pi=replicate(1.d0,nchan)
             ra=replicate(0.d0,nchan)
             rd=replicate(0.d0,nchan)
             h=replicate(0.d0,nchan)
             chan_id=replicate(0.d0,nchan) ;;0 for fida chords
            end
            'NSTX_NPA': begin 
             ;;Intersection of NPA sightline and midplane
             xlos=npa_mid_nstx[*,0]
             ylos=npa_mid_nstx[*,1]
             zlos=npa_mid_nstx[*,2]
             ;;NPA detector location
             xhead=npap_nstx[*,0]
             yhead=npap_nstx[*,1]
             zhead=npap_nstx[*,2]
             
             nchan=n_elements(xlos)        
             sigma_pi=replicate(1.d0,nchan)
             ra=aperture_radius_nstx ;;aperture radius
             rd=det_radius_nstx      ;;detector radius
             h=dcollimator_nstx      ;;length collimator length
             chan_id=replicate(1.d0,nchan) ;;1 for npa chords
            end
            'NSTXU_TSSNPA': begin       
             ;;Intersection of NPA sightline and midplane
             xlos=tssnpa_mid[*,0]
             ylos=tssnpa_mid[*,1]
             zlos=tssnpa_mid[*,2]
             ;;NPA detector location
             xhead=tssnpap[*,0]
             yhead=tssnpap[*,1]
             zhead=tssnpap[*,2]
             
             nchan=n_elements(xlos)        
             sigma_pi=replicate(1.d0,nchan)
             ra=aperture_radius_tssnpa ;;aperture radius
             rd=det_radius_tssnpa      ;;detector radius
             h=dcollimator_tssnpa      ;;length collimator length
             chan_id=replicate(1.d0,nchan) ;;1 for npa chords
            end
            'NSTXU_RSSNPA': begin       
             ;;Intersection of NPA sightline and midplane
             xlos=rssnpa_mid[*,0]
             ylos=rssnpa_mid[*,1]
             zlos=rssnpa_mid[*,2]
             ;;NPA detector location
             xhead=rssnpap[*,0]
             yhead=rssnpap[*,1]
             zhead=rssnpap[*,2]
             
             nchan=n_elements(xlos)        
             sigma_pi=replicate(1.d0,nchan)
             ra=aperture_radius_rssnpa ;;aperture radius
             rd=det_radius_rssnpa      ;;detector radius
             h=dcollimator_rssnpa      ;;length collimator length
             chan_id=replicate(1.d0,nchan) ;;1 for npa chords
            end
            'NSTXU_PSSNPA': begin       
             ;;Intersection of NPA sightline and midplane
             xlos=pssnpa_mid[*,0]
             ylos=pssnpa_mid[*,1]
             zlos=pssnpa_mid[*,2]
             ;;NPA detector location
             xhead=pssnpap[*,0]
             yhead=pssnpap[*,1]
             zhead=pssnpap[*,2]
             
             nchan=n_elements(xlos)        
             sigma_pi=replicate(1.d0,nchan)
             ra=aperture_radius_pssnpa ;;aperture radius
             rd=det_radius_pssnpa      ;;detector radius
             h=dcollimator_pssnpa      ;;length collimator length
             chan_id=replicate(1.d0,nchan) ;;1 for npa chords
            end
            'NSTXU_NPA': begin      
             ;;Intersection of NPA sightline and midplane
             xlos=nstxu_npa_mid[*,0]
             ylos=nstxu_npa_mid[*,1]
             zlos=nstxu_npa_mid[*,2]
             ;;NPA detector location
             xhead=nstxu_npap[*,0]
             yhead=nstxu_npap[*,1]
             zhead=nstxu_npap[*,2]
             
             nchan=n_elements(xlos)        
             sigma_pi=replicate(1.d0,nchan)
             ra=aperture_radius_nstxu_npa ;;aperture radius
             rd=det_radius_nstxu_npa      ;;detector radius
             h=dcollimator_nstxu_npa      ;;length collimator length
             chan_id=replicate(1.d0,nchan) ;;1 for npa chords
            end     
            else: begin
              print, '% Diagnostic unknown'
              stop
            end
        endcase
    
        if i eq 0 then begin
           xloss=xlos & yloss=ylos & zloss=zlos
           xheads=xhead & yheads=yhead & zheads=zhead
           nchans=nchan & sigma_pis=sigma_pi & ras=ra & rds=rd
           hs=h & chan_ids=chan_id
        endif else begin
           xloss=[xloss,xlos] & yloss=[yloss,ylos] & zloss=[zloss,zlos]
           xheads=[xheads,xhead] & yheads=[yheads,yhead] & zheads=[zheads,zhead]
           nchans+=nchan & sigma_pis=[sigma_pis,sigma_pi] & ras=[ras,ra] & rds=[rds,rd] 
           hs=[hs,h] & chan_ids=[chan_ids,chan_id]
        endelse
    endfor
        
    ;;SAVE IN FIDA STRUCTURE,xloss
    fida={nchan:nchans,diag:fida_diag,$
          xlos:double(xloss),ylos:double(yloss),zlos:double(zloss),$
          xlens:double(xheads),ylens:double(yheads),zlens:double(zheads),$
          sigma_pi_ratio:sigma_pis,ra:ras,rd:rds,h:hs,chan_id:chan_ids}
    
    if keyword_set(doplot) then begin
       nstx_input,inputs
       print,'Isource',inputs.isource
       nbi=nstx_beams(inputs,/doplot)
       npt=1001
       nbi_x=nbi.xyz_src[0]+(nbi.xyz_pos[0]-nbi.xyz_src[0])*findgen(npt)/(npt-1.0)*2.5
       nbi_y=nbi.xyz_src[1]+(nbi.xyz_pos[1]-nbi.xyz_src[1])*findgen(npt)/(npt-1.0)*2.5
       min_i=-1
       min_j=-1
       oplot,[0.,200*cos(angle_transp_uvz)],[0,200.*sin(angle_transp_uvz)]
       for ich=0,fida.nchan-1 do begin
           x=fida.xlens[ich]+(fida.xlos[ich]-fida.xlens[ich])*findgen(npt)/(npt-1.)*1.1
           y=fida.ylens[ich]+(fida.ylos[ich]-fida.ylens[ich])*findgen(npt)/(npt-1.)*1.1
           ;oplot,x,y,color=!red,thick=2
           oplot,fida.xlens,fida.ylens,color=!darkgreen,psym=1,symsize=2
           oplot,fida.xlos,fida.ylos,color=!red,psym=4
           dmin=999999999.
           for i=0L,npt-1 do begin
              for j=0L,npt-1 do begin
                  dd=sqrt((x[i]-nbi_x[j])^2.+ (y[i]-nbi_y[j])^2)
              if dd lt dmin then begin
                 dmin=dd
             min_i=i
             min_j=j
              endif
              endfor
           endfor
           ;oplot,[x[min_i],x[min_i]],[y[min_j],y[min_j]],color=!red,psym=4
           print,'Detector head (u,v,z)',fida.xlens[ich],fida.ylens[ich],fida.zlens[ich]
           print,'Intersectional point at midplane (u,v,z)',$
                  fida.xlos[ich],fida.ylos[ich],fida.zlos[ich],sqrt(min((nbi_x-fida.xlos[ich])^2.0+(nbi_y-fida.ylos[ich])^2.0))
           print,'Rmajor at intersectional point'
           print,sqrt(x[min_i]^2.+y[min_i]^2.),sqrt(nbi_x[min_j]^2.+nbi_y[min_j]^2.),dmin
           print,'Distance from det head',sqrt( (x[min_i]-fida.xlens[ich])^2.+$
                                                (y[min_i]-fida.ylens[ich])^2.)
           rmajor=sqrt(x^2.+y^2.)
           raxis=105.684
           dum=min(abs(rmajor-raxis),w)
           print,'Distance from det head to magnetic axis',sqrt( (x[w]-fida.xlens[ich])^2.+$
                                                             (y[w]-fida.ylens[ich])^2.),rmajor[w]
           
           ;print,'Rtan'
           ;print,CXRTAN_rssnpa[ich]
           ;print,CXRTAN_tssnpa[i]
       
       endfor
    endif
    
    return,fida

END
