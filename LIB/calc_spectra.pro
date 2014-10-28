PRO calc_spectra,grid,fbm,weights,los,neutrals,lambda,radiance,mean_fbm

	nchan=weights.nchan
	lambda=weights.lambda
    nwav=n_elements(lambda)
    nen=n_elements(weights.energy)
    npitch=n_elements(weights.pitch)
    dE=abs(weights.energy[1]-weights.energy[0])
    dpitch=abs(weights.pitch[1]-weights.pitch[0])
    dwav=abs(weights.lambda[1]-weights.lambda[0])

    radiance= fltarr(nchan,nwav)

	calc_mean_fbm,grid,fbm,weights,los,neutrals,mean_fbm,elevel=2
        ;;------------------------------------------------------------
        ;; ------ CALCULATE SYNTHETIC SPECTRA AND PROFIELES ----------
        ;;------------------------------------------------------------
	for ichan=0,nchan-1 do begin
        for ii=0,nwav-1 do begin
           ;; PRODUCT with fast-ion distribution funciton
           prod=replicate(0.,nen,npitch)
           for ie=0,nen-1 do begin
              for ip=0,npitch-1 do begin
                 prod[ie,ip]=mean_fbm[ie,ip,ichan]*weights.wfunct[ii,ie,ip,ichan]
                 ;;--> [ph/(s cm^2 dP keV)]
              endfor
           endfor
           radiance[ichan,ii]=total(prod[*,*])$
                    *dE*dpitch/(4.d0*!pi)*1.d4 / dwav ;;--> [ph/(s m^2 sr nm)]
        endfor
  endfor ;; loop over channels	
  radiance=transpose(radiance)
END
