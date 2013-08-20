PRO load_results,result_dir,output

	;;READ INPUTS
	read_inputs,result_dir+'inputs.dat',inputs

	;;READ GRID
	read_grid,result_dir+'grid.bin',grid

	;;READ LOS
	read_los,result_dir+'los.bin',inputs,los

	;;READ PLASMA 
	read_plasma,result_dir+'plasma.bin',plasma

	;;READ FIDA SPECTRA
	read_fida,result_dir+'fida_spectra.bin',fida

	;;READ HALO AND NBI SPECTRA
	read_nbi_halo,result_dir+'nbi_halo_spectra.bin',nbi_halo

	;;READ NEUTRALS
	read_neutrals,result_dir+'neutrals.bin',neutrals

	;;READ NPA
	read_npa,result_dir+'npa.bin',npa

	;;READ FBM
	read_fbm,result_dir+'transp_fbm.bin',fbm

	;;READ WEIGHT FUNCTIONS
	read_weights,result_dir+'weight_function.bin',weights

	;;READ BIRTH PROFILE
	read_birth,result_dir+'birth.bin',birth

	output={inputs:inputs,grid:grid,los:los,plasma:plasma,$
			fida:fida,nbi_halo:nbi_halo,neutrals:neutrals,$
			npa:npa,fbm:fbm,weights:weights,birth:birth}
END
