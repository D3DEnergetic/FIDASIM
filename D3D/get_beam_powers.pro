FUNCTION get_beam_powers,shot,time,beam_num,all=all
	
	if keyword_set(all) then beam_num=indgen(8)
	powers=fltarr(n_elements(beam_num))
	for i=0,n_elements(beam_num)-1 do begin
		if beam_num[i] eq 0 then begin
			gadat,t,tmp,'pinj_30l',shot		
			w=where(abs(t-time) lt 5.0)
			pow=max(tmp[w])/1000.0
		endif

		if beam_num[i] eq 1 then begin
			gadat,t,tmp,'pinj_30r',shot		
			w=where(abs(t-time) lt 5.0)
			pow=max(tmp[w])/1000.0
		endif

		if beam_num[i] eq 2 then begin
			gadat,t,tmp,'pinj_15l',shot		
			w=where(abs(t-time) lt 5.0)
			pow=max(tmp[w])/1000.0
		endif

		if beam_num[i] eq 3 then begin
			gadat,t,tmp,'pinj_15r',shot		
			w=where(abs(t-time) lt 5.0)
			pow=max(tmp[w])/1000.0
		endif

		if beam_num[i] eq 4 then begin
			gadat,t,tmp,'pinj_21l',shot		
			w=where(abs(t-time) lt 5.0)
			pow=max(tmp[w])/1000.0
		endif

		if beam_num[i] eq 5 then begin
			gadat,t,tmp,'pinj_21r',shot		
			w=where(abs(t-time) lt 5.0)
			pow=max(tmp[w])/1000.0
		endif

		if beam_num[i] eq 6 then begin
			gadat,t,tmp,'pinj_33l',shot		
			w=where(abs(t-time) lt 5.0)
			pow=max(tmp[w])/1000.0
		endif

		if beam_num[i] eq 7 then begin
			gadat,t,tmp,'pinj_33r',shot		
			w=where(abs(t-time) lt 5.0)
			pow=max(tmp[w])/1000.0
		endif
			
		powers[i]=pow
		t=0 & tmp=0 & pow=0
	endfor
	return,powers
END
		
