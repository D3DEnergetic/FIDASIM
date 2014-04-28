FUNCTION nstx_beams,inputs
; Returns the nbi structure
;;  IDL> help,nbi
    ;;  ** Structure <1d475af8>, 13 tags, length=168, data length=168, refs=1:
    ;;     EINJ            DOUBLE           80.775734
    ;;     PINJ            DOUBLE           2.4117758
    ;;     FULL            DOUBLE          0.54850105
    ;;     HALF            DOUBLE          0.28972649
    ;;     THIRD           DOUBLE          0.16177245
    ;;     XYZ_SRC         DOUBLE    Array[3]
    ;;     XYZ_POS         DOUBLE    Array[3]
    ;;     BMWIDRA         DOUBLE           6.0000000
    ;;     BMWIDZA         DOUBLE           24.000000
    ;;     DIVY            DOUBLE    Array[3]
    ;;     DIVZ            DOUBLE    Array[3]
    ;;     FOCY            DOUBLE           999999.90
    ;;     FOCZ            DOUBLE           1000.0000

	einj=double(inputs.einj) & pinj=double(inputs.pinj)

    ;;GET SPECIES_MIX
    cgfitf=[-0.109171,0.0144685,-7.83224e-5]
    cgfith=[0.0841037,0.00255160,-7.42683e-8]
    ;; Current fractions
    ffracs=cgfitf[0]+cgfitf[1]*einj+cgfitf[2]*einj^2
    hfracs=cgfith[0]+cgfith[1]*einj+cgfith[2]*einj^2
    tfracs=1.0-ffracs-hfracs
	
    ;;BEAM GEOMETRY STUFF GOES HERE

	;;SAVE IN NBI STRUCTURE
	nbi={einj:einj,pinj:pinj,full:ffracs,half:hfracs,third:tfracs,$
		 xyz_src:reform(xyz_src[inputs.isource,*]),xyz_pos:reform(xyz_pos[inputs.isource,*]),bmwidra:bmwidra,bmwidza:bmwidza,$
		 divy:divy[*,inputs.isource],divz:divz[*,inputs.isource],focy:focy[inputs.isource],focz:focz[inputs.isource]}
	return,nbi
END
