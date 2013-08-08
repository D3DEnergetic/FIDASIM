;;RENAME TO "DEVICE"_ROUTINES I.E. D3D_ROUTINES AND RENAME FILE ACCORDINGLY
PRO d3d_routines,inputs,grid,$ 			;;INPUT: INPUTS AND GRID
					   nbi,$ 			;;OUTPUT: NEUTRAL BEAM INJECTION INFO STRUCTURE
					   fida,$ 			;;OUTPUT: FIDA DIAGNOSTIC INFO STRUCTURE
					   profiles,$		;;OUTPUT: PROFILES STRUCTURE
					   equil,$			;;OUTPUT: MAGNETIC GRID STRUCTURE
					   err				;;OUTPUT: ERROR STATUS ERR=1 == SOMETHING WENT WRONG

	err=0

	;;GET BEAM GEOMETRY
	nbi=d3d_beams(inputs)
	
	;;GET CHORD GEOMETRY
	fida=d3d_chords(inputs.shot,inputs.fida_diag)

	;;GET PROFILES
	profiles=d3d_profiles(inputs)
	if profiles.err eq 1 then begin
		print,'FAILED TO GET PROFILES'
		err=1
		goto,GET_OUT
	endif

	;;GET E&M FIELDS AT GRID POINTS
	equil=d3d_equil(inputs,grid,fida)
	if equil.err eq 1 then begin
		print,'FAILED TO GET EQUILIBRIUM'
		err=1
		goto,GET_OUT
	endif

	GET_OUT:
END 
