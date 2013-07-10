;;RENAME TO "DEVICE"_ROUTINES I.E. D3D_ROUTINES AND RENAME FILE ACCORDINGLY
PRO templete_routines,inputs,grid,$     ;;INPUT: INPUTS AND GRID POINTS DO NOT CHANGE
					   nbi,$ 			;;OUTPUT: NEUTRAL BEAM INJECTION INFO STRUCTURE
					   fida,$ 			;;OUTPUT: FIDA DIAGNOSTIC INFO STRUCTURE
					   profiles,$		;;OUTPUT: PROFILES STRUCTURE
					   equil,$			;;OUTPUT: MAGNETIC GRID STRUCTURE
					   err				;;OUTPUT: ERROR STATUS ERR=1 == SOMETHING WENT WRONG

	
	;;IN THIS SECTION YOU CAN USE WHATEVER ROUTINES YOU WANT SO LONG AS YOU DEFINE THE OUTPUTS AS FOLLOWS
	
	;; nbi structure (8 == number of beams)
	;;** Structure <8394cd8>, 13 tags, length=1104, data length=1104, refs=1:
	;;   EINJ            FLOAT     Array[8]
	;;   PINJ            DOUBLE    Array[8]
	;;   FULL            FLOAT     Array[8]
	;;   HALF            FLOAT     Array[8]
	;;   THIRD           FLOAT     Array[8]
	;;   XYZ_SRC         DOUBLE    Array[8, 3]
	;;   XYZ_POS         DOUBLE    Array[8, 3]
	;;   BMWIDRA         DOUBLE           11.000000
	;;   BMWIDZA         DOUBLE           25.000000
	;;   DIVY            DOUBLE    Array[3, 8]
	;;   DIVZ            DOUBLE    Array[3, 8]
	;;   FOCY            DOUBLE    Array[8]
	;;   FOCZ            DOUBLE    Array[8]

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
	;;   HEADSIZE        FLOAT     Array[15]

	;; profiles structure
	;;** Structure <83e8518>, 7 tags, length=5768, data length=5764, refs=1:
	;;   TIME            FLOAT           4.42000
	;;   RHO             DOUBLE    Array[120]
	;;   TI              DOUBLE    Array[120] [eV]
	;;   VTOR            DOUBLE    Array[120] [m/s]
	;;   TE              DOUBLE    Array[120] [eV]
	;;   DENE            DOUBLE    Array[120] [m^-3]
	;;   ZEFF            DOUBLE    Array[120]
	
	;; equil structure (3 == x,y,z components
	;;** Structure <83e79e8>, 2 tags, length=3360000, data length=3360000, refs=1:
	;;   RHO_GRID        FLOAT     Array[120000]
	;;   B               DOUBLE    Array[3, 120000]
	;;   E               DOUBLE    Array[3, 120000]

	equil={rho_grid:rho_grid,$	   			;;FIDA GRID IN MAGNETIC FLUX COORDINATES (RHO)
			  b:b,$					   		;;MAGNETIC FIELD COMPONENTS AT GRID POINTS
			  e:e}							;;ELECTRIC FIELD COMPONENTS AT GRID POINTS

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
		 xlos:xlos,$						;;X POS. OF LINE OF SIGHT [cm]
		 ylos:ylos,$						;;Y POS. OF LINE OF SIGHT [cm]
         zlos:zlos,$						;;Z POS. OF LINE OF SIGHT [cm]
		 xhead:xhead,$						;;X POS. OF HEAD [cm]
		 yhead:yhead,$						;;Y POS. OF HEAD [cm]
 		 zhead:zhead,$						;;Z POS. OF HEAD [cm]
  		 headsize:headsize}					;;SIZE OF HEAD

	profiles={time:time,$					;;SHOT TIME
			  rho:rho,$						;;RHO VALUES
			  ti:ti,$						;;ION TEMPERATURE [eV]
			  vtor:vtor,$					;;TORODIAL VELOCITY [m/s]
			  te:te,$						;;ELECTRON TEMPERATURE [eV]
			  dene:dene,$					;;ELECTRON DENSITY [m^-3]
			  zeff:zeff}					;;ZEFF
END 
