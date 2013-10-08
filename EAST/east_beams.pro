FUNCTION east_beams,inputs
; Returns the nbi structure

	;; Here from design division 
	;; This places the aperture as actually located, and the source is
	;; then "xedge" cm further away from the vessel.
	;;     A port, F port
	us_NB=[871.4, 899.5, $
	       -449.1, -561.6]
	vs_NB=[-247.0, -167.7,$
	       766.8, 710.5]


	; Cross-over points from drawing
	u_co=[362.2, 362.2, -159.0, -159.0]
	v_co=[-22.1, -22.1, 326.2, 326.2]

	nsources=n_elements(vs_NB)
	xyz_src=[[double(us_NB)],[double(vs_NB)],[replicate(0.0d,nsources)]]
	xyz_pos=[[double(u_co)],[double(v_co)],[replicate(0.0d,nsources)]]

;	focy=replicate(1d33,nsources)      ; horizontal focal length
;	(infinity)
  focy=replicate(999999.9d0,nsources) ; so f90 can read input 
	focz=replicate(1000d0,nsources)      ; vertical focal length is 10 m
	divy=replicate(1.0472d-2,3, nsources);(8.73d-3,3,nsources)    ; horizontal divergence in radians
	divz=replicate(2.0944d-2,3, nsources);(2.27d-2,3,nsources)

	bmwidra=6d0     ; ion source half width in cm
	bmwidza=24d0    ; ion source half height in cm
 

	;------------------------------------------
	;;Get beam energy,power and fractions

	if inputs.einj eq 0. or inputs.pinj eq 0. then $
		a=get_beam_power(inputs.shot,inputs.time*1000.,inputs.isource)
	if inputs.einj eq 0. then einj=a.einj else einj=inputs.einj
	if inputs.pinj eq 0. then pinj=a.pinj else pinj=inputs.pinj
	einj=double(einj) & pinj=double(pinj)

    ;;SPECIES_MIX
    ;; Power fractions from Hu's conceptual design paper
    pfractions=[0.8,0.14,0.06]
   ;; FIDASIM uses current fractions
     const=1./(pfractions[0]+2*pfractions[1]+3*pfractions[2])
     ffracs=const*pfractions[0]
     hfracs=2*const*pfractions[1]
     tfracs=3*const*pfractions[2]
	

	;;SAVE IN NBI STRUCTURE
	nbi={einj:einj,pinj:pinj,full:ffracs,half:hfracs,third:tfracs,$
		 xyz_src:reform(xyz_src[inputs.isource,*]),xyz_pos:reform(xyz_pos[inputs.isource,*]),bmwidra:bmwidra,bmwidza:bmwidza,$
		 divy:divy[*,inputs.isource],divz:divz[*,inputs.isource],focy:focy[inputs.isource],focz:focz[inputs.isource]}
	return,nbi
END
