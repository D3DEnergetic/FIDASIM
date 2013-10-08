;;RENAME TO "DEVICE"_ROUTINES I.E. EAST_ROUTINES AND RENAME FILE ACCORDINGLY
PRO east_routines,inputs,grid,$ 			;;INPUT: INPUTS AND GRID
					   nbi,$ 			;;OUTPUT: NEUTRAL BEAM INJECTION INFO STRUCTURE
					   chords,$ 		;;OUTPUT: CHORDS INFO STRUCTURE
					   profiles,$		;;OUTPUT: PROFILES STRUCTURE
					   equil,$			;;OUTPUT: MAGNETIC GRID STRUCTURE
					   err				;;OUTPUT: ERROR STATUS ERR=1 == SOMETHING WENT WRONG

	err=0

	;;GET BEAM GEOMETRY
	nbi=east_beams(inputs)
	
	;;GET CHORD GEOMETRY
	chords=east_chords(inputs.shot,inputs.diag)

	;;GET PROFILES
	profiles=create_transp_profiles(inputs,/no_omega)
	if profiles.err eq 1 then begin
		print,'FAILED TO GET PROFILES'
		err=1
		goto,GET_OUT
	endif

	;;GET E&M FIELDS AT GRID POINTS
	equil=east_equil(inputs,grid,chords)
	if equil.err eq 1 then begin
		print,'FAILED TO GET EQUILIBRIUM'
		err=1
		goto,GET_OUT
	endif

	GET_OUT:
END 
