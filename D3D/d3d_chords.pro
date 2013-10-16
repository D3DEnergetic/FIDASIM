FUNCTION d3d_chords,shot,fida_diag

	fida_diag=strupcase(fida_diag)
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
    ;;   HEADSIZE        FLOAT     Array[15]

	;;FROM fida_vanzeeland DIII-D 2-D camera at 90 degrees
	xlens1=[276.326]
	ylens1=[-5.45047]
	zlens1=[7.62]
	xmid1=[19.6787]
	ymid1=[179.223]
	zmid1=[0.]	
	
	;;From fida_mckee for 105 port for linear BES array

	xmid2=[ 133.617,  133.490,  133.350,  133.222,  133.094, $
			132.965,  132.837,  132.708,  132.579,  132.449, $
			132.320,  132.190,  132.073,  131.942,  131.811, $
			131.693,  131.575,  131.443,  131.324,  131.205, $
			131.086,  130.966,  130.846,  130.725,  130.605, $
			130.498,  130.376,  130.255,  130.146,  130.024, $
			129.915,  129.806]

	ymid2=[-168.692, -167.515, -166.219, -165.040, -163.856, -162.669, -161.482, -160.293, -159.100, -157.904 $
		 , -156.706, -155.506, -154.424, -153.217, -152.009, -150.919, -149.828, -148.610, -147.511, -146.411 $
		 , -145.309, -144.204, -143.095, -141.983, -140.870, -139.879, -138.758, -137.635, -136.635, -135.507 $
		 , -134.502, -133.493]
	
	zmid2=replicate(0.,32)
	
	xlens2=replicate(261.7,32)
	ylens2=replicate(-70.1,32)
	zlens2=replicate(15.0 ,32)
	
	;;from fida_chris for tangential viewing FIDA system
	xmid3=[-135.21,-134.667,-133.97,-133.239,-132.486,-131.677,-130.913,-130.249,-129.474,-128.621,-127.879]
	ymid3=[-177.627,-172.757,-166.495,-159.936,-153.178,-145.923,-139.065,-133.102,-126.145,-118.493,-111.834]
	zmid3=replicate(0.0,n_elements(xmid3))
	xlens3=replicate(-46.02,n_elements(xmid3))
	ylens3=replicate(-198.5,n_elements(xmid3))
	zlens3=replicate(122.0,n_elements(xmid3))

    ;;from fida_grierson_MU1_MU2 (30LT/RT from 315T0)
    xlens4=replicate(-182.609,8)
    ylens4=replicate(203.626,8)
    zlens4=replicate(3.4302,8)
    xmid4=[34.803, 42.116, 48.071, 54.495, 67.172, 68.501, 75.616, 82.214] ;; cm
    ymid4=[167.035, 172.633, 177.193, 182.110, 187.223, 192.833, 198.280, 203.331] ;; cm
    zmid4=[-3.510, -3.710, -3.870, -4.270, -4.700, -4.830, -4.940,-5.310] ;; cm

	;;fida_yadong
	if shot le 143721 then begin
		xmid5=[-128.4,-129.3,-130.0,-130.8,-131.3,-131.9,-132.7,-128.9,-131.0,-133.2]
		ymid5=[120.9,129.1,135.1,142.2,147.2,152.4,160.2,125.1,144.9,165.2]
		nchan5=n_elements(xmid5)
		zmid5=replicate(0.,nchan5)
		xlens5=replicate(-142.55,nchan5)
		ylens5=replicate( 142.55,nchan5)
		zlens5=replicate(-152.00,nchan5)

		xmid6=replicate(0.0,nchan5) & ymid6=xmid5
		zmid6=replicate(0.0,nchan5)
		xmid6[0]=-128.9
		ymid6[0]=125.1
		xlens6=replicate(-163.0,nchan5)
		ylens6=replicate(163.0,nchan5)
		zlens6=replicate(0.0,nchan5)
	endif else begin
		name_cv1=['v17','v18','v19','v20']
		radii_cv1=[176.3,182.48,187.75,192.77]
		phi_cv1=[313.18,314.93,316.05,317.27]
		umid_cv1=radii_cv1*sin(!pi*phi_cv1/180.)
		vmid_cv1=radii_cv1*abs(cos(!pi*phi_cv1/180.))

		; Note that v08 and v24 are on u3 but they are collected by 330 R+1 lens
		name_u3=['v01','v02','v06','v21','v22','v23']
		radii_u3=[185.3,190.03,216.42,196.49,201.09,207.63]
		phi_u3=[315.5,316.64,321.77,318.16,319.06,320.33]
		umid_u3=radii_u3*sin(!pi*phi_u3/180.)
		vmid_u3=radii_u3*abs(cos(!pi*phi_u3/180.))

		; Using old values for yadong's system
		name_1G=['v1','v2','v3']
		umid_1G=[-128.9,-131.0,-133.2]
		vmid_1G=[125.1,144.9,165.2]
		radii_1G=[179.625,195.338,212.210]

		xmid5=[umid_cv1,umid_u3,umid_1G]
		ymid5=[vmid_cv1,vmid_u3,vmid_1G]
		nchan5=n_elements(xmid5)
		zmid5=replicate(0.0,nchan5)
		xlens5=replicate(-142.55,nchan5)
        ylens5=replicate( 142.55,nchan5)
        zlens5=replicate(-152.00,nchan5)

        xmid6=replicate(0.0,nchan5) & ymid6=xmid5
        zmid6=replicate(0.0,nchan5)
		radii=[224.67,223.]
		phi=[325.62,325.5]
		xmid6[0:1]=radii*sin(!pi*phi/180.)
		ymid6[0:1]=radii*abs(cos(!pi*phi/180.))

        xlens6=replicate(235*sin(!pi*330./180.),nchan5)
        ylens6=replicate(235*abs(cos(!pi*330./180.)),nchan5)
        zlens6=replicate(100.0,nchan5)
	endelse

    ;; from fida_grierson_MU3_MU4 210LT/RT from 165R0
    xlens7=replicate(83.51,8)
    ylens7=replicate(-239.2,8)
    zlens7=replicate(0.,8)
    xmid7=[-126.145, -128.594, -129.522, -130.514, -131.476, -132.440, -133.568, -134.655]
    ymid7=[-99.732,  -121.144, -129.260, -137.932, -146.338, -154.773, -164.629, -174.137]
    zmid7=[-4.950,   -4.620,   -4.450,   -4.450,   -4.320,   -3.94 , -4.010  , -4.270]

	;;NPA CHORDS
    detxyz=fltarr(3,3)
    detxyz[0,*]=[2.52,.08,.81]
    detxyz[1,*]=[2.5,.05,.86]
    detxyz[2,*]=[2.49,.09,.89]
	npap=detxyz*0.
    for i=0,2 do begin
        npap[i,0]=-(detxyz[i,0] + detxyz[i,1])/sqrt(2.)
        npap[i,1]=(detxyz[i,1] - detxyz[i,0])/sqrt(2.)
        npap[i,2]=-detxyz[i,2]
    endfor
    npap*=100   ; convert from meters to centimeters
    
    iwall=fltarr(3,3)
    iwall[0,*]=[0.95,1.15,32.4]
    iwall[1,*]=[0.95,0.70,4.19]
    iwall[2,*]=[0.95,0.49,-10.96]
	npa_mid=iwall*0.
    for i=0,2 do begin
        npa_mid[i,0]=iwall[i,0]*cos(!pi*(225+iwall[i,2])/180)
        npa_mid[i,1]=iwall[i,0]*sin(!pi*(225+iwall[i,2])/180)
        npa_mid[i,2]=iwall[i,1]
    endfor
    npa_mid*=100    ; cm

 ;; SELECT THE VIEWS
	for i=0,n_elements(fida_diag)-1 do begin
    	CASE (fida_diag[i]) OF
      	   'VERTICAL': begin
        	 xlos=xmid5
        	 ylos=ymid5
        	 zlos=zmid5
        	 xhead=xlens5
        	 yhead=ylens5
        	 zhead=zlens5
        	 nchan=n_elements(xlos)
        	 sigma_pi=replicate(1.d0,nchan)
        	 headsize=replicate(0.d0,nchan)
        	 opening_angle=replicate(0.d0,nchan)
       		end
      		'OBLIQUE': begin
        	 xlos=xmid3
        	 ylos=ymid3
        	 zlos=zmid3
        	 xhead=xlens3
        	 yhead=ylens3
        	 zhead=zlens3
        	 nchan=n_elements(xlos)
        	 sigma_pi=replicate(1.d0,nchan)
        	 headsize=replicate(0.d0,nchan)
         	 opening_angle=replicate(0.d0,nchan)
       		end
       		'TANGENTIAL30': begin
        	 xlos=xmid4
        	 ylos=ymid4
        	 zlos=zmid4
        	 xhead=xlens4
        	 yhead=ylens4
        	 zhead=zlens4
        	 nchan=n_elements(xlos)
        	 sigma_pi=replicate(1.d0,nchan)
        	 headsize=replicate(0.d0,nchan)
        	 opening_angle=replicate(0.d0,nchan)
      		end
      		'TANGENTIAL210': begin
        	 xlos=xmid7
        	 ylos=ymid7
        	 zlos=zmid7
        	 xhead=xlens7
        	 yhead=ylens7
        	 zhead=zlens7
        	 nchan=n_elements(xlos)
        	 sigma_pi=replicate(1.d0,nchan)
        	 headsize=replicate(0.d0,nchan)
        	 opening_angle=replicate(0.d0,nchan)
      		end
      	 	'ALL': begin
			 xlos=[xmid1,xmid2,xmid3,xmid4,xmid5,xmid6,xmid7]
			 ylos=[ymid1,ymid2,ymid3,ymid4,ymid5,ymid6,ymid7]
			 zlos=[zmid1,zmid2,zmid3,zmid4,zmid5,zmid6,zmid7]
			 xhead=[xlens1,xlens2,xlens3,xlens4,xlens5,xlens6,xlens7]
			 yhead=[ylens1,ylens2,ylens3,ylens4,ylens5,ylens6,ylens7]
			 zhead=[zlens1,zlens2,zlens3,zlens4,zlens5,zlens6,zlens7]
			 nchan=n_elements(xlos)
        	 sigma_pi=replicate(1.d0,nchan)
        	 headsize=replicate(0.d0,nchan)
        	 opening_angle=replicate(0.d0,nchan)
	   	 	end
      		'NPA': begin
        	 xlos=npa_mid[*,0]
        	 ylos=npa_mid[*,1]
        	 zlos=npa_mid[*,2]
        	 xhead=npap[*,0]
        	 yhead=npap[*,1]
        	 zhead=npap[*,2]
        	 nchan=n_elements(xlos)        
        	 sigma_pi=replicate(1.d0,nchan)
        	 headsize=replicate(1.5d0,nchan)
        	 opening_angle=replicate(0.0261799388d0,nchan)
       		end
        	ELSE: begin
         	 PRINT, '% Diagnostic unknown'
         	 STOP
            end
    	ENDCASE

		if i eq 0 then begin
			xloss=xlos & yloss=ylos & zloss=zlos
			xheads=xhead & yheads=yhead & zheads=zhead
			nchans=nchan & sigma_pis=sigma_pi & headsizes=headsize 
			opening_angles=opening_angle
		endif else begin
            xloss=[xloss,xlos] & yloss=[yloss,ylos] & zloss=[zloss,zlos]
            xheads=[xheads,xhead] & yheads=[yheads,yhead] & zheads=[zheads,zhead]
            nchans+=nchan & sigma_pis=[sigma_pis,sigma_pi] & headsizes=[headsizes,headsize] 
    	    opening_angles=[opening_angles,opening_angle]
		endelse
	endfor
	
	;;SAVE IN FIDA STRUCTURE
	fida={nchan:nchans,diag:fida_diag,xlos:double(xloss),ylos:double(yloss),zlos:double(zloss),$
		  xlens:double(xheads),ylens:double(yheads),zlens:double(zheads),$
		  sigma_pi_ratio:sigma_pis,headsize:headsizes,opening_angle:opening_angles}
	return,fida
END
