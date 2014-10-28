;;RENAME TO "DEVICE"_ROUTINES I.E. D3D_ROUTINES AND RENAME FILE ACCORDINGLY
PRO d3d_routines,inputs,grid,$ 			;;INPUT: INPUTS AND GRID
					   nbi,$ 			;;OUTPUT: NEUTRAL BEAM INJECTION INFO STRUCTURE
					   chords,$ 		;;OUTPUT: CHORDS INFO STRUCTURE
					   profiles,$		;;OUTPUT: PROFILES STRUCTURE
					   equil,$			;;OUTPUT: MAGNETIC GRID STRUCTURE
					   err				;;OUTPUT: ERROR STATUS ERR=1 == SOMETHING WENT WRONG

	err=0
 
    ;; CHECK FOR D3D SPECIFIC INPUTS
    inVars=strlowcase(TAG_NAMES(inputs))
    if where('use_transp' eq inVars) eq -1 then begin
        inputs = create_struct(inputs,'use_transp',0) ;;DEFAULT TO GAPROFILES
    endif

	;;GET BEAM GEOMETRY
	nbi=d3d_beams(inputs)
	
	;;GET CHORD GEOMETRY
	chords=d3d_chords(inputs.shot,inputs.diag,isource=inputs.isource)

	;;GET PROFILES
    if inputs.use_transp eq 1 then begin
        profiles=create_transp_profiles(inputs)
    endif else begin
	    profiles=d3d_profiles(inputs)
    endelse

	if profiles.err eq 1 then begin
		print,'FAILED TO GET PROFILES'
		err=1
		goto,GET_OUT
	endif

	;;GET E&M FIELDS AT GRID POINTS
	equil=d3d_equil(inputs,grid,chords)
	if equil.err eq 1 then begin
		print,'FAILED TO GET EQUILIBRIUM'
		err=1
		goto,GET_OUT
	endif

	GET_OUT:
END 
