;;RENAME TO "DEVICE"_ROUTINES I.E. D3D_ROUTINES AND RENAME FILE ACCORDINGLY
PRO templete_routines,inputs,grid,$     ;;INPUT: INPUTS AND GRID POINTS DO NOT CHANGE
					   nbi,$ 			;;OUTPUT: NEUTRAL BEAM INJECTION INFO STRUCTURE
					   chords,$ 		;;OUTPUT: CHORDS INFO STRUCTURE
					   profiles,$		;;OUTPUT: PROFILES STRUCTURE
					   equil,$			;;OUTPUT: MAGNETIC GRID STRUCTURE
					   err				;;OUTPUT: ERROR STATUS ERR=1 == SOMETHING WENT WRONG

	
	;;IN THIS SECTION YOU CAN USE WHATEVER ROUTINES 
	;;YOU WANT SO LONG AS YOU DEFINE THE OUTPUT STRUCTURES
	;;CONTAIN AT LEAST THE FOLLOWING TAGS

	;;	IDL> help,chords 
	;;	** Structure <1d447c48>, 11 tags, length=728, data length=724, refs=1:
	;;	   NCHAN           LONG                11
	;;	   DIAG            STRING    'OBLIQUE'
	;;	   XLOS            DOUBLE    Array[11]
	;;	   YLOS            DOUBLE    Array[11]
	;;	   ZLOS            DOUBLE    Array[11]
	;;	   XLENS           DOUBLE    Array[11]
	;;	   YLENS           DOUBLE    Array[11]
	;;	   ZLENS           DOUBLE    Array[11]
	;;	   SIGMA_PI_RATIO  DOUBLE    Array[11]
	;;	   HEADSIZE        FLOAT     Array[11]
	;;	   OPENING_ANGLE   FLOAT     Array[11]

	;;	IDL> help,equil
	;;	** Structure <1d474638>, 10 tags, length=6636160, data length=6636138, refs=1:
	;;	   RHO_GRID        FLOAT     Array[40, 60, 50]
	;;	   RHO_CHORDS      STRUCT    -> <Anonymous> Array[1]
	;;	   BX              DOUBLE    Array[40, 60, 50]
	;;	   BY              DOUBLE    Array[40, 60, 50]
	;;	   BZ              DOUBLE    Array[40, 60, 50]
	;;	   EX              DOUBLE    Array[40, 60, 50]
	;;	   EY              DOUBLE    Array[40, 60, 50]
	;;	   EZ              DOUBLE    Array[40, 60, 50]
	;;	   ERR             INT              0

	;;	IDL> help,equil.rho_chords
	;;	** Structure <1d48bf08>, 2 tags, length=352008, data length=352004, refs=2:
	;;	   RHOS            DOUBLE    Array[4000, 11]
	;;	   DS              FLOAT          0.300000
	;;

	;;	IDL> help,profiles
	;;	** Structure <1d475698>, 7 tags, length=5816, data length=5810, refs=1:
	;;	   RHO             DOUBLE    Array[121]
	;;	   TE              DOUBLE    Array[121]
	;;	   TI              DOUBLE    Array[121]
	;;	   VTOR            DOUBLE    Array[121]
	;;	   DENE            DOUBLE    Array[121]
	;;	   ZEFF            DOUBLE    Array[121]
	;;	   ERR             INT              0

	;;	IDL> help,nbi
	;;	** Structure <1d475af8>, 13 tags, length=168, data length=168, refs=1:
	;;	   EINJ            DOUBLE           80.775734
	;;	   PINJ            DOUBLE           2.4117758
	;;	   FULL            DOUBLE          0.54850105
	;;	   HALF            DOUBLE          0.28972649
	;;	   THIRD           DOUBLE          0.16177245
	;;	   XYZ_SRC         DOUBLE    Array[3]
	;;	   XYZ_POS         DOUBLE    Array[3]
	;;	   BMWIDRA         DOUBLE           6.0000000
	;;	   BMWIDZA         DOUBLE           24.000000
	;;	   DIVY            DOUBLE    Array[3]
	;;	   DIVZ            DOUBLE    Array[3]
	;;	   FOCY            DOUBLE           999999.90
	;;	   FOCZ            DOUBLE           1000.0000

	;;FOR CONVINIENCE HERE ARE THE MINIMUM STRUCTURE DEFINITIONS
	equil={rho_grid:rho_grid,$	   			;;FIDA GRID IN MAGNETIC FLUX COORDINATES (RHO)
		   rho_chords:rho_chords,$			;;STRUCTURE CONTAINING AN ARRAY OF RHO VALUES AND STEP SIZE IN [cm]
		   bx:bx,$					   		;;X MAGNETIC FIELD COMPONENT AT GRID POINTS
		   by:by,$					   		;;Y MAGNETIC FIELD COMPONENT AT GRID POINTS
		   bz:bz,$					   		;;Z MAGNETIC FIELD COMPONENT AT GRID POINTS
		   ex:ex,$							;;X ELECTRIC FIELD COMPONENT AT GRID POINTS
		   ey:ey,$							;;Y ELECTRIC FIELD COMPONENT AT GRID POINTS
		   ez:ez }							;;Z ELECTRIC FIELD COMPONENT AT GRID POINTS

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

  	chords={sigma_pi_ratio:sigma_pi_ratio,$	;;RATIO OF SIGMA LINES TO PI LINES
		 nchan:nchan,$				  		;;NUMBER OF CHANNELS
		 xmid:xmid,$						;;X POS. OF WHERE CHORD CROSSES MIDPLANE [cm]
		 ymid:ymid,$						;;Y POS. OF WHERE CHORD CROSSES MIDPLANE [cm]
         zmid:zmid,$						;;Z POS. OF WHERE CHORD CROSSES MIDPLANE [cm]
		 xlens:xlens,$						;;X POS. OF LENS [cm]
		 ylens:ylens,$						;;Y POS. OF LENS [cm]
 		 zlens:zlens,$						;;Z POS. OF LENS [cm]
  		 headsize:headsize,$				;;SIZE OF HEAD
		 opening_angle:opening_angle}		;;OPENING ANGLE

	profiles={time:time,$					;;SHOT TIME
			  rho:rho,$						;;RHO VALUES
			  ti:ti,$						;;ION TEMPERATURE [eV]
			  vtor:vtor,$					;;TORODIAL ANGULAR VELOCITY [rad/s]
			  te:te,$						;;ELECTRON TEMPERATURE [eV]
			  dene:dene,$					;;ELECTRON DENSITY [m^-3]
			  zeff:zeff}					;;ZEFF
END 
