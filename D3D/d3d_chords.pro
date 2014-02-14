FUNCTION get_cer_geom,shot,isource

	a=GET_CERGEOM(shot)
	b=GET_CER_BEAM_ORDER(shot)
	beams=['30LT','30RT','150LT','150RT','210LT','210RT','330LT','330RT']
	nchan=n_elements(a.labels)
	xlens=dblarr(nchan) & xlos=xlens
	ylens=dblarr(nchan)	& ylos=ylens
	zlens=dblarr(nchan)	& zlos=zlens
	err_arr=intarr(nchan)

	whb=where(b eq beams[isource])
	for i=0,nchan-1 do begin
      	rl=a.rcers[i]*100.
      	phil=a.phicers[i]
      	rb=a.rcere[i,whb]*100.
      	phib=a.phicere[i,whb]

		if phib le 360.0 then begin
      		zlens[i]=a.zcers[i]*100.
    	  	xlens[i]=a.rcers[i]*COS((90. - a.phicers[i])*!DTOR)*100.
	      	ylens[i]=a.rcers[i]*SIN((90. - a.phicers[i])*!DTOR)*100.
      		zlos[i]=a.zcere[i,whb]*100.
      		xlos[i]=a.rcere[i,whb]*COS((90 - a.phicere[i,whb])*!DTOR)*100.
      		ylos[i]=a.rcere[i,whb]*SIN((90 - a.phicere[i,whb])*!DTOR)*100.
		endif else err_arr[i]=1
  	endfor
	w=where(err_arr eq 0)
	output={chords:a.labels[w],xlens:xlens[w],ylens:ylens[w],zlens:zlens[w],xlos:xlos[w],ylos:ylos[w],zlos:zlos[w]}
	return,output
END

FUNCTION d3d_chords,shot,fida_diag,isource=isource

	if n_elements(isource) eq 0 then isource=6
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
    ;;   ra              FLOAT     Array[15]
    ;;   rd              FLOAT     Array[15]
    ;;   h               FLOAT     Array[15]

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

    cer_chords=get_cer_geom(shot,isource)
	w=where(strmid(cer_chords.chords,0,1) eq 'V',nw)
	if nw ne 0 then begin
		xmid5=cer_chords.xlos[w]
		ymid5=cer_chords.ylos[w]
		zmid5=cer_chords.zlos[w]
		xlens5=cer_chords.xlens[w]
		ylens5=cer_chords.ylens[w]
		zlens5=cer_chords.zlens[w]
	endif else begin
        xmid5=0. & ymid5=0. & zmid5=100.
        xlens5=0. & ylens5=0. & zlens5=0.
    endelse

    ;; from fida_grierson_MU3_MU4 210LT/RT from 165R0
    xlens7=replicate(83.51,8)
    ylens7=replicate(-239.2,8)
    zlens7=replicate(0.,8)
    xmid7=[-126.145, -128.594, -129.522, -130.514, -131.476, -132.440, -133.568, -134.655]
    ymid7=[-99.732,  -121.144, -129.260, -137.932, -146.338, -154.773, -164.629, -174.137]
    zmid7=[-4.950,   -4.620,   -4.450,   -4.450,   -4.320,   -3.94 , -4.010  , -4.270]

	;;CER TANGENTIAL CHORDS
	w=where(strmid(cer_chords.chords,0,1) ne 'V',nw)
	if nw ne 0 then begin
		xmid8=cer_chords.xlos[w]
		ymid8=cer_chords.ylos[w]
		zmid8=cer_chords.zlos[w]
		xlens8=cer_chords.xlens[w]
		ylens8=cer_chords.ylens[w]
		zlens8=cer_chords.zlens[w]
	endif else begin
		xmid8=0. & ymid8=0. & zmid8=100.
		xlens8=0. & ylens8=0. & zlens8=0.
	endelse

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
        	 ra=replicate(0.d0,nchan)
        	 rd=replicate(0.d0,nchan)
        	 h=replicate(0.d0,nchan)
        	 chan_id=replicate(0.d0,nchan)
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
        	 ra=replicate(0.d0,nchan)
           	 rd=replicate(0.d0,nchan)
		  	 h=replicate(0.d0,nchan)
        	 chan_id=replicate(0.d0,nchan)
       		end
            'TANGENTIAL': begin
             xlos=xmid8
             ylos=ymid8
             zlos=zmid8
             xhead=xlens8
             yhead=ylens8
             zhead=zlens8
             nchan=n_elements(xlos)
             sigma_pi=replicate(1.d0,nchan)
             ra=replicate(0.d0,nchan)
        	 rd=replicate(0.d0,nchan)
             h=replicate(0.d0,nchan)
        	 chan_id=replicate(0.d0,nchan)
            end
       		'MAIN_ION30': begin
        	 xlos=xmid4
        	 ylos=ymid4
        	 zlos=zmid4
        	 xhead=xlens4
        	 yhead=ylens4
        	 zhead=zlens4
        	 nchan=n_elements(xlos)
        	 sigma_pi=replicate(1.d0,nchan)
        	 ra=replicate(0.d0,nchan)
        	 rd=replicate(0.d0,nchan)
        	 h=replicate(0.d0,nchan)
        	 chan_id=replicate(0.d0,nchan)
      		end
      		'MAIN_ION210': begin
        	 xlos=xmid7
        	 ylos=ymid7
        	 zlos=zmid7
        	 xhead=xlens7
        	 yhead=ylens7
        	 zhead=zlens7
        	 nchan=n_elements(xlos)
        	 sigma_pi=replicate(1.d0,nchan)
        	 ra=replicate(0.d0,nchan)
        	 rd=replicate(0.d0,nchan)
        	 h=replicate(0.d0,nchan)
        	 chan_id=replicate(0.d0,nchan)
      		end
      	 	'ALL': begin
			 xlos=[xmid1,xmid2,xmid3,xmid4,xmid5,xmid7,xmid8]
			 ylos=[ymid1,ymid2,ymid3,ymid4,ymid5,ymid7,ymid8]
			 zlos=[zmid1,zmid2,zmid3,zmid4,zmid5,zmid7,zmid8]
			 xhead=[xlens1,xlens2,xlens3,xlens4,xlens5,xlens7,xlens8]
			 yhead=[ylens1,ylens2,ylens3,ylens4,ylens5,ylens7,ylens8]
			 zhead=[zlens1,zlens2,zlens3,zlens4,zlens5,zlens7,zlens8]
			 nchan=n_elements(xlos)
        	 sigma_pi=replicate(1.d0,nchan)
        	 ra=replicate(0.d0,nchan)
        	 rd=replicate(0.d0,nchan)
        	 chan_id=replicate(0.d0,nchan)
        	 h=replicate(0.d0,nchan)
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
        	 ra=replicate(0.5d0,nchan)
        	 rd=replicate(0.5d0,nchan)
        	 h=replicate(25.4d0,nchan)
             ;;Rev. Sci. Instrum. 83, 10D304 (2012) <-NPA Geometry
        	 chan_id=replicate(1.d0,nchan)
       		end
        	ELSE: begin
         	 PRINT, '% Diagnostic unknown'
         	 STOP
            end
    	ENDCASE

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
	
	;;SAVE IN FIDA STRUCTURE
	fida={nchan:nchans,diag:fida_diag,xlos:double(xloss),ylos:double(yloss),zlos:double(zloss),$
		  xlens:double(xheads),ylens:double(yheads),zlens:double(zheads),$
		  sigma_pi_ratio:sigma_pis,ra:ras,rd:rds,h:hs,chan_id:chan_ids}
	return,fida
END
