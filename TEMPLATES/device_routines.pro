;;RENAME TO "DEVICE"_ROUTINES I.E. D3D_ROUTINES AND RENAME FILE ACCORDINGLY
PRO templete_routines,inputs,grid,$     ;;INPUT: INPUTS AND GRID POINTS DO NOT CHANGE
					   nbi,$ 			;;OUTPUT: NEUTRAL BEAM INJECTION INFO STRUCTURE
					   fida,$ 			;;OUTPUT: FIDA DIAGNOSTIC INFO STRUCTURE
					   profiles,$		;;OUTPUT: PROFILES STRUCTURE
					   equil,$			;;OUTPUT: MAGNETIC GRID STRUCTURE
					   err				;;OUTPUT: ERROR STATUS ERR=1 == SOMETHING WENT WRONG

	
	;;IN THIS SECTION YOU CAN USE WHATEVER ROUTINES 
	;;YOU WANT SO LONG AS YOU DEFINE THE OUTPUT STRUCTURES
	;;CONTAIN AT LEAST THE FOLLOWING TAGS
	
	;;** Structure <10d7c568>, 13 tags, length=168, data length=164, refs=1:
	;;   EINJ            DOUBLE           74.369202
	;;   PINJ            DOUBLE           2.0655587
	;;   FULL            DOUBLE          0.53365599
	;;   HALF            DOUBLE          0.27345340
	;;   THIRD           DOUBLE          0.19289061
	;;   XYZ_SRC         DOUBLE    Array[3]
	;;   XYZ_POS         DOUBLE    Array[3]
	;;   BMWIDRA         DOUBLE           6.0000000
	;;   BMWIDZA         DOUBLE           24.000000
	;;   DIVY            DOUBLE    Array[3]
	;;   DIVZ            DOUBLE    Array[3]
	;;   FOCY            FLOAT            1.e+06
	;;   FOCZ            DOUBLE           1000.0000

	;;** Structure <10d7aeb8>, 9 tags, length=288, data length=288, refs=1:
	;;   NCHAN           LONG                10
	;;   XMID            FLOAT     Array[10]
	;;   YMID            FLOAT     Array[10]
	;;   ZMID            FLOAT     Array[10]
	;;   XLENS           FLOAT     Array[10]
	;;   YLENS           FLOAT     Array[10]
	;;   ZLENS           FLOAT     Array[10]
	;;   SIGMA_PI_RATIO  FLOAT           1.00000
	;;   HEADSIZE        FLOAT     Array[10]

	;;** Structure <10dc9e28>, 5 tags, length=3380784, data length=3380772, refs=1:
	;;   RHO_GRID        FLOAT     Array[120000]
	;;	 RHO_CHORDS		 STRUC	   STRUCTURE
	;;   B               FLOAT     Array[3, 120000]
	;;   E               FLOAT     Array[3, 120000]

	;;** Structure <10dc9e28>, 5 tags, length=3380784, data length=3380772, refs=1:
	;;   RHOS			 DOUBLE    Array[4000,10] ;;10 = number of channels
	;;   DS				 FLOAT 	   0.30000		  ;;Step size in cm

	;;** Structure <10d955b8>, 8 tags, length=3648, data length=3642, refs=1:
	;;   RHO             FLOAT     Array[101]
	;;   TE              DOUBLE    Array[101]
	;;   TI              DOUBLE    Array[101]
	;;   VTOR            FLOAT     Array[101]
	;;   DENE            DOUBLE    Array[101]
	;;   ZEFF            FLOAT     Array[101]

	;;FOR CONVINIENCE HERE ARE THE MINIMUM STRUCTURE DEFINITIONS
	equil={rho_grid:rho_grid,$	   			;;FIDA GRID IN MAGNETIC FLUX COORDINATES (RHO)
		   rho_chords:rho_chords,$			;;STRUCTURE CONTAINING AN ARRAY OF RHO VALUES AND STEP SIZE IN [cm]
		   b:b,$					   		;;MAGNETIC FIELD COMPONENTS AT GRID POINTS
		   e:e }							;;ELECTRIC FIELD COMPONENTS AT GRID POINTS

	nbi={einj:einj,$				   		;;BEAM INJECTION ENERGY [keV]
		 pinj:pinj,$				   		;;BEAM INJECTION POWER  [MW]
		 full:full,$				   		;;FULL BEAM FRACTION
		 half:half,$				   		;;HALF BEAM FRACTION
		 third:third,$				   		;;THIRD BEAM FRACTION
		 xyz_src:xyz_src,$			   		;;POSITION OF BEAM SOURCE IN MACHINE COORDINATES [cm]
		 xyz_pos:xyz_pos,$			   		;;BEAM CROSSOVER POINT IN MACHINE COORDINATES [cm]
		 bmwidra:bmwidra,$			   		;;HORIZONTAL BEAM WIDTH [cm]
		 bmwidza:mbwidza,$			   		;;VERTICAL BEAM WIDTH   [cm]
		 focy:focy,$				   		;;HORIZONTAL FOCAL LENGTH [cm]
	     focz:focz,$						;;VERTICAL FOCAL LENGTH [cm]
		 divy:divy,$				   		;;HORIZONTAL BEAM DIVERGENCE [rad]
		 divz:divz }				   		;;VERTICAL BEAM DIVERGENCE [rad]

  	fida={sigma_pi_ratio:sigma_pi_ratio,$	;;RATIO OF SIGMA LINES TO PI LINES
		 nchan:nchan,$				  		;;NUMBER OF CHANNELS
		 xmid:xmid,$						;;X POS. OF WHERE CHORD CROSSES MIDPLANE [cm]
		 ymid:ymid,$						;;Y POS. OF WHERE CHORD CROSSES MIDPLANE [cm]
         zmid:zmid,$						;;Z POS. OF WHERE CHORD CROSSES MIDPLANE [cm]
		 xlens:xlens,$						;;X POS. OF LENS [cm]
		 ylens:ylens,$						;;Y POS. OF LENS [cm]
 		 zlens:zlens,$						;;Z POS. OF LENS [cm]
  		 headsize:headsize}					;;SIZE OF HEAD

	profiles={time:time,$					;;SHOT TIME
			  rho:rho,$						;;RHO VALUES
			  ti:ti,$						;;ION TEMPERATURE [eV]
			  vtor:vtor,$					;;TORODIAL ANGULAR VELOCITY [rad/s]
			  te:te,$						;;ELECTRON TEMPERATURE [eV]
			  dene:dene,$					;;ELECTRON DENSITY [m^-3]
			  zeff:zeff}					;;ZEFF
END 
