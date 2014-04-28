FUNCTION nstx_chords,shot,fida_diag,isource=isource

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
        	 chan_id=replicate(0.d0,nchan) ;;0 for fida chords
           end
      	   'TANGENTIAL': begin
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
        	 chan_id=replicate(0.d0,nchan) ;;0 for fida chords
           end
      		'ssNPA': begin 
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
        	 chan_id=replicate(1.d0,nchan) ;;1 for npa chords
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
