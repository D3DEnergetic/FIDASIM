PRO calc_spectra,grid,fbm,weights,los,neutrals,lambda,radiance,mean_fbm

	nchan=los.nchan
	nwav=weights.nwav
	lambda=weights.central_wavel

    nen_transp=n_elements(fbm.energy)
    npitch_transp=n_elements(fbm.pitch) 
    dE_transp=fbm.energy[1]-fbm.energy[0]
    dpitch_transp=abs(fbm.pitch[1]-fbm.pitch[0])


    radiance= fltarr(nchan,nwav)

	calc_mean_fbm,grid,fbm,weights,los,neutrals,mean_fbm
        ;;------------------------------------------------------------
        ;; ------ CALCULATE SYNTHETIC SPECTRA AND PROFIELES ----------
        ;;------------------------------------------------------------
	for ichan=0,los.nchan-1 do begin
        for ii=0,weights.nwav-1 do begin
           ;; PRODUCT with fast-ion distribution funciton
           prod=replicate(0.,weights.nen,weights.npitch)
           for ie=0,weights.nen-1 do begin
              for ip=0,weights.npitch-1 do begin
                 prod[ie,ip]=mean_fbm[ie,ip,ichan]*weights.weight_tot[ii,ie,ip,ichan]
                 ;;--> [ph/(s cm^2 dP keV)]
              endfor
           endfor
           radiance[ichan,ii]=total(prod[*,*])*weights.dE * $
                              weights.dpitch/(4.d0*!pi)*1.d4 / $
                              weights.dwav ;;--> [ph/(s m^2 sr nm)]
        endfor
  endfor ;; loop over channels	
  radiance=transpose(radiance)
END
