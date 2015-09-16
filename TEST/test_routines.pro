PRO test_routines,inputs,grid,$     ;;INPUT: INPUTS AND GRID POINTS DO NOT CHANGE
					   nbi,$ 			;;OUTPUT: NEUTRAL BEAM INJECTION INFO STRUCTURE
					   chords,$ 		;;OUTPUT: CHORDS INFO STRUCTURE
					   profiles,$		;;OUTPUT: PROFILES STRUCTURE
					   equil,$			;;OUTPUT: MAGNETIC GRID STRUCTURE
					   err				;;OUTPUT: ERROR STATUS ERR=1 == SOMETHING WENT WRONG
           
    profiles = read_ncdf(inputs.profile_dir)
    equil = test_equil(inputs,grid)

    ;; Chords
    ulens = dblarr(3)
    wlens = replicate(100.0,3)
    vlens = [-170.0,-170.0,-170.0]
    ulos = dblarr(3)
    wlos = dblarr(3)
    vlos = [-200.0,-170.0,-140.0]

    ulens = [ulens,ulens]
    vlens = [vlens,vlens]
    wlens = [wlens,wlens]
    ulos = [ulos,ulos]
    vlos = [vlos,vlos]
    wlos = [wlos,wlos]

    ra = [replicate(0.0,3),replicate(3.0,3)]
    rd = [replicate(0.0,3),replicate(3.0,3)]
    h  = [replicate(0.0,3),replicate(50.0,3)]
    sigma_pi = [replicate(1.0,3),replicate(0.0,3)]
    chan_id = [replicate(0,3),replicate(1,3)]
    
    chords = {nchan:6,diag:["SPECTRAL","NPA"],$
              ulos:ulos,vlos:vlos,wlos:wlos,xlos:ulos,ylos:vlos,zlos:wlos,$
              ulens:ulens,vlens:vlens,wlens:wlens,xlens:ulens,ylens:vlens,zlens:wlens,$
              sigma_pi_ratio:sigma_pi,ra:ra,rd:rd,h:h,chan_id:chan_id,err:0}

    uvw_src = [0.0,-230.0 - 2.0*cos(inputs.beta),2.0*sin(inputs.beta)]
    uvw_pos = [0.0,-230.0,0.0]
    einj = inputs.einj
    pinj = inputs.pinj
    focy = 999999.9d0
    focz = 1000d0
    divy = replicate(8.73d-3,3)
    divz = replicate(2.27d-2,3)
    bmwidra=6d0
    bmwidza=24d0

    cgfitf=[-0.109171,0.0144685,-7.83224e-5]
    cgfith=[0.0841037,0.00255160,-7.42683e-8]

    ffracs=cgfitf[0]+cgfitf[1]*einj+cgfitf[2]*einj^2
    hfracs=cgfith[0]+cgfith[1]*einj+cgfith[2]*einj^2
    tfracs=1.0-ffracs-hfracs

    nbi={einj:einj,pinj:pinj,full:ffracs,half:hfracs,third:tfracs,$
         uvw_src:uvw_src,uvw_pos:uvw_pos,bmwidra:bmwidra,bmwidza:bmwidza,$
         divy:divy,divz:divz,focy:focy,focz:focz,err:0}

    ;err=0
END
