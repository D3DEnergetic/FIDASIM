PRO load_results,result_dir,results,save=save

	;;CHECK FOR SLASH
	slash=strmid(result_dir,0,1,/reverse_offset)
	if slash ne '/' then result_dir+='/'

	;;READ INPUTS
	read_inputs,result_dir+'inputs.dat',inputs

	;;READ GRID
	read_grid,result_dir+'parameters.cdf',grid

	;;READ LOS
	read_los,result_dir+'parameters.cdf',los

	;;READ PLASMA 
	read_plasma,result_dir+'parameters.cdf',plasma

	;;READ FIDA SPECTRA
	read_fida,result_dir+'fida_spectra.cdf',fida

	;;READ HALO AND NBI SPECTRA
	read_nbi_halo,result_dir+'nbi_halo_spectra.cdf',nbi_halo

	;;READ NEUTRALS
	read_neutrals,result_dir+'neutrals.cdf',neutrals

	;;READ NPA
	read_npa,result_dir+'npa.cdf',npa

	;;READ FBM
	read_fbm,result_dir+'parameters.cdf',fbm

	;;READ WEIGHT FUNCTIONS
	read_weights,result_dir+'weight_function.cdf',weights

	;;READ BIRTH PROFILE
	read_birth,result_dir+'birth.cdf',birth

	results={inputs:inputs,grid:grid,los:los,plasma:plasma,$
			fida:fida,nbi_halo:nbi_halo,neutrals:neutrals,$
			npa:npa,fbm:fbm,weights:weights,birth:birth}
	if keyword_set(save) then begin
		save,inputs,grid,los,plasma,fida,nbi_halo,neutrals,$
		npa,fbm,weights,birth,filename=inputs.fidasim_runid+'_results.sav',/compress
	endif
END
