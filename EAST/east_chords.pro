FUNCTION east_chords,shot,fida_diag

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

	
	;;From Tangential O-port

	xmid1=[ 234.4987,  230.1616,  225.7936,  221.3921, $
			    216.954,  212.478,  207.959,  203.3939, $
			    198.7784,  194.1077,  189.376,  184.5769, $
			    179.703,  174.7453,  169.6934,  164.5347]

	ymid1=[ 34.3418,  36.2577,  38.1873,  40.1316, $
          42.0919,  44.0694,  46.0657,  48.0823, $
          50.1212,  52.1845,  54.2747,  56.3946, $
          58.5477,  60.7377,  62.9694,  65.2483]
	
	zmid1=replicate(0.,16)
	
	xlens1=replicate(186.4,16)
	ylens1=replicate(-187.0,16)
	zlens1=replicate(-12.63 ,16)
	
	;;From Vertical B-port

  xmid2=[ 198.7784,  194.1077,  189.376,  184.5769, $
          179.703,  174.7453,  169.6934,  164.5347]

  ymid2=[ 50.1212,  52.1845,  54.2747,  56.3946, $
          58.5477,  60.7377,  62.9694,  65.2483]
  
  zmid2=replicate(0.,8)
  
  xlens2=replicate(183.1,8)
  ylens2=replicate(80.3,8)
  zlens2=replicate(130.4 ,8)

  ;;From Vertical N-port

  xmid3=[ 198.7784,  194.1077,  189.376,  184.5769, $
          179.703,  174.7453,  169.6934,  164.5347]

  ymid3=[ 50.1212,  52.1845,  54.2747,  56.3946, $
          58.5477,  60.7377,  62.9694,  65.2483]
  
  zmid3=replicate(0.,8)
  
  xlens3=replicate(80.3,8)
  ylens3=replicate(-183.1,8)
  zlens3=replicate(130.4 ,8)
  
  ;;From Vertical A-port
  
  
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
    CASE (fida_diag) OF
      'TANGENTIAL-O': begin
        xlos=xmid1
        ylos=ymid1
        zlos=zmid1
        xhead=xlens1
        yhead=ylens1
        zhead=zlens1
        nchan=n_elements(xlos)
        sigma_pi=replicate(1.d0,nchan)
        headsize=replicate(0.d0,nchan)
        opening_angle=replicate(0.d0,nchan)
       end
      'VERTICAL-B': begin
        xlos=xmid2
        ylos=ymid2
        zlos=zmid2
        xhead=xlens2
        yhead=ylens2
        zhead=zlens2
        nchan=n_elements(xlos)
        sigma_pi=replicate(1.d0,nchan)
        headsize=replicate(0.d0,nchan)
        opening_angle=replicate(0.d0,nchan)
       end
      'VERTICAL-N': begin
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
      'VERTICAL-A': begin
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
      'ALL': begin
		xlos=[xmid1,xmid2,xmid3]
		ylos=[ymid1,ymid2,ymid3]
		zlos=[zmid1,zmid2,zmid3]
		xhead=[xlens1,xlens2,xlens3]
		yhead=[ylens1,ylens2,ylens3]
		zhead=[zlens1,zlens2,zlens3]
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
       END
    ENDCASE

	;;SAVE IN FIDA STRUCTURE
	fida={nchan:nchan,diag:fida_diag,xlos:double(xlos),ylos:double(ylos),zlos:double(zlos),$
		  xlens:double(xhead),ylens:double(yhead),zlens:double(zhead),$
		  sigma_pi_ratio:sigma_pi,headsize:headsize,opening_angle:opening_angle}
	return,fida
END
