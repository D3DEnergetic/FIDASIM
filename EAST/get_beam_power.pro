FUNCTION get_beam_power,shot,time,beam_num,avg=avg
; WWH 7/11/13
; INPUTS
; shot
; time (ms)
; beam_num   integer
; OUTPUT
; pinj (for source beam_num) (MW)
; einj (kV)
; KEYWORD
; avg   average power from time through time+avg (ms)
;       default is 5 ms

if not keyword_set(avg) then avg=5.
n=get_nbi(shot,/fast)
w=where(1000.*n.time ge time and 1000.*n.time le time+avg,nw)
if nw gt 0 then pinj=1.e-6*total(n.pbeam[w,beam_num])/nw else begin
  print,'Failed to retrieve beam power',time,avg,beam_num
  pinj=0.
end

einj=1.e-3*n.volts[beam_num]
return,{pinj:pinj,einj:einj}

end
		
