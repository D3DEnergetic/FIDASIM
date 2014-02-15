PRO load_results,result_dir,results,save=save

	;;CHECK FOR SLASH
	slash=strmid(result_dir,0,1,/reverse_offset)
	if slash ne '/' then result_dir+='/'

	runid=file_basename(result_dir)
	input_file=runid+'_inputs.cdf'

    ;;READ INPUTS
    inputs=read_ncdf(result_dir+input_file,vars=['shot','time'])

	;;READ GRID
    gridvars=['Nx','Ny','Nz','alpha','beta','origin',$
              'x_grid','y_grid','z_grid',$
              'u_grid','v_grid','w_grid',$
              'phi_grid','r_grid','rho_grid','xx','yy','zz']
    grid=read_ncdf(result_dir+input_file,vars=gridvars)

	;;READ LOS
    losvars=['diagnostic','chan_id','Nchan','xlos','ylos','zlos','rlos',$
             'xlens','ylens','zlens','sigma_pi','h','rd','ra','los_wght']
	los=read_ncdf(result_dir+input_file,vars=losvars)

	;;READ PLASMA 
    pvars=['ti','te','dene','deni','denp','denf','vrotx','vroty','vrotz',$
           'zeff','bx','by','bz','ex','ey','ez','ai','impurity_charge','btipsign']
	plasma=read_ncdf(result_dir+input_file,vars=pvars)

    ;;READ BEAM
    bvars=['beam','ab','einj','pinj','divy','divz','focy','focz','bmwidra','bmwidza',$
           'xyz_src','species_mix','Arot','Brot','Crot']
    beam=read_ncdf(result_dir+input_file,vars=bvars)

	;;READ SPECTRA
	spectra=read_ncdf(result_dir+runid+'_spectra.cdf')

	;;READ NEUTRALS
	neutrals=read_ncdf(result_dir+runid+'_neutrals.cdf')

	;;READ NPA
	npa=read_ncdf(result_dir+runid+'_npa.cdf')

	;;READ FBM
    fbmvars=['FBM_time','FBM_Nenergy','FBM_Npitch','FBM_Ngrid',$
             'FBM_r2d','FBM_z2d','FBM_bmvol','FBM','FBM_emin',$
             'FBM_emax','FBM_energy','FBM_pmin','FBM_pmax','FBM_pitch']
	fbm=read_ncdf(result_dir+input_file,vars=fbmvars)
    fbm.FBM_pitch=fbm.FBM_pitch*plasma.btipsign

	;;READ FIDA WEIGHT FUNCTIONS
	fida_weights=read_ncdf(result_dir+runid+'_fida_weights.cdf')
    fida_weights.pitch=fida_weights.pitch*plasma.btipsign

	;;READ NPA WEIGHT FUNCTIONS
	npa_weights=read_ncdf(result_dir+runid+'_npa_weights.cdf')
    npa_weights.pitch=npa_weights.pitch*plasma.btipsign

	;;READ BIRTH PROFILE
	birth=read_ncdf(result_dir+runid+'_birth.cdf')

	results={inputs:inputs,grid:grid,los:los,plasma:plasma,$
			spectra:spectra,neutrals:neutrals,beam:beam,$
			npa:npa,fbm:fbm,fida_weights:fida_weights,npa_weights:npa_weights,birth:birth}

	if keyword_set(save) then begin
		save,inputs,grid,los,plasma,spectra,neutrals,$
		npa,fbm,fida_weights,npa_weights,birth,filename=inputs.fidasim_runid+'_results.sav',/compress
	endif
END
