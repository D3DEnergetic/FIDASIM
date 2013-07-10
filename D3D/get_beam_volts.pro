FUNCTION get_beam_volts,shot,time,beam_num,all=all
	
	if keyword_set(all) then beam_num=indgen(8)
	voltage=fltarr(n_elements(beam_num))
	for i=0,n_elements(beam_num)-1 do begin
		if beam_num[i] eq 0 then begin
			gadat,t,tmp,'nbvolt_30l',shot		
			w=where(abs(t-time) lt 5.0)
			volts=max(tmp[w])/1000.0
		endif

		if beam_num[i] eq 1 then begin
			gadat,t,tmp,'nbvolt_30r',shot		
			w=where(abs(t-time) lt 5.0)
			volts=max(tmp[w])/1000.0
		endif

		if beam_num[i] eq 2 then begin
			gadat,t,tmp,'nbvolt_15l',shot		
			w=where(abs(t-time) lt 5.0)
			volts=max(tmp[w])/1000.0
		endif

		if beam_num[i] eq 3 then begin
			gadat,t,tmp,'nbvolt_15r',shot		
			w=where(abs(t-time) lt 5.0)
			volts=max(tmp[w])/1000.0
		endif

		if beam_num[i] eq 4 then begin
			gadat,t,tmp,'nbvolt_21l',shot		
			w=where(abs(t-time) lt 5.0)
			volts=max(tmp[w])/1000.0
		endif

		if beam_num[i] eq 5 then begin
			gadat,t,tmp,'nbvolt_21r',shot		
			w=where(abs(t-time) lt 5.0)
			volts=max(tmp[w])/1000.0
		endif

		if beam_num[i] eq 6 then begin
			gadat,t,tmp,'nbvolt_33l',shot		
			w=where(abs(t-time) lt 5.0)
			volts=max(tmp[w])/1000.0
		endif

		if beam_num[i] eq 7 then begin
			gadat,t,tmp,'nbvolt_33r',shot		
			w=where(abs(t-time) lt 5.0)
			volts=max(tmp[w])/1000.0
		endif
			
		voltage[i]=volts
		t=0 & tmp=0 & volts=0
	endfor
	return,voltage
END
		
