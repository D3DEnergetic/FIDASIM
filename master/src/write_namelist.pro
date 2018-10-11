PRO write_namelist, filename, inputs
    ;+##`write_namelist, filename, inputs`
    ;+Writes namelist file
    ;+
    ;+###Input Arguments
    ;+     **filename**: Name of the namelist file
    ;+
    ;+     **inputs**: Input structure
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> write_namelist, filename, inputs
    ;+```
    info,'Writing namelist file...'

    fidasim_version = get_version(get_fidasim_dir())

    openw,55,filename
    printf,55,'!! Created: ', systime()
    printf,55,'!! FIDASIM version: '+fidasim_version
    printf,55,'!! Comment: '+inputs.comment
    printf,55,'&fidasim_inputs'
    printf,55,''
    printf,55,'!! Shot Info'
    printf,55,f='("shot = ", i6 ,"    !! Shot Number")',inputs.shot
    printf,55,f='("time = ", 1f8.5 ,"    !! Time [s]")',inputs.time
    printf,55,"runid = '" + inputs.runid + "'    !! runID"
    printf,55,"result_dir = '" + inputs.result_dir+"'    !! Result Directory"
    printf,55,''
    printf,55,'!! Input Files'
    printf,55,"tables_file = '" + inputs.tables_file+"'    !! Atomic Tables File"
    printf,55,"equilibrium_file = '" + inputs.equilibrium_file +"'    !! File containing plasma parameters and fields"
    printf,55,"geometry_file = '" + inputs.geometry_file +"'    !! File containing NBI and diagnostic geometry"
    printf,55,"distribution_file = '" + inputs.distribution_file +"'    !! File containing fast-ion distribution"
    printf,55,''
    printf,55,'!! Simulation Switches'
    printf,55,f='("calc_bes = ",i2 , "    !! Calculate NBI Spectra")',inputs.calc_bes
    printf,55,f='("calc_dcx = ",i2 , "    !! Calculate Direct CX Spectra")',inputs.calc_dcx
    printf,55,f='("calc_halo = ",i2 , "    !! Calculate Halo Spectra")',inputs.calc_halo
    printf,55,f='("calc_cold = ",i2 , "    !! Calculate Cold D-alpha Spectra")',inputs.calc_cold
    printf,55,f='("calc_brems = ",i2 , "    !! Calculate Bremsstrahlung")',inputs.calc_brems
    printf,55,f='("calc_fida = ",i2 , "    !! Calculate Active FIDA Spectra")',inputs.calc_fida
    printf,55,f='("calc_npa = ",i2 , "   !! Calculate Active NPA")',inputs.calc_npa
    printf,55,f='("calc_pfida = ",i2 , "    !! Calculate Passive FIDA Spectra")',inputs.calc_pfida
    printf,55,f='("calc_pnpa = ",i2 , "   !! Calculate Passive NPA")',inputs.calc_pnpa
    printf,55,f='("calc_neutron = ",i2 , "   !! Calculate B-T Neutron Rate")',inputs.calc_neutron
    printf,55,f='("calc_birth = ",i2 , "    !! Calculate Birth Profile")',inputs.calc_birth
    printf,55,f='("calc_fida_wght = ",i2 , "    !! Calculate FIDA weights")',inputs.calc_fida_wght
    printf,55,f='("calc_npa_wght = ",i2 , "    !! Calculate NPA weights")',inputs.calc_npa_wght
    printf,55,''
    printf,55,'!! Debugging Switches'
    printf,55,f='("no_flr = ",i2,"    !! Turn off Finite Larmor Radius effects")',inputs.no_flr
    printf,55,f='("load_neutrals = ",i2,"    !! Load neutrals from neutrals file")',inputs.load_neutrals
    printf,55,"neutrals_file = '" + inputs.neutrals_file +"'    !! File containing the neutral density"
    printf,55,f='("verbose = ",i2,"    !! Verbose")',inputs.verbose
    printf,55,''
    printf,55,'!! Monte Carlo Settings'
    printf,55,f='("n_fida = ",i9,"    !! Number of Active FIDA mc particles")',inputs.n_fida
    printf,55,f='("n_pfida = ",i9,"    !! Number of Passive FIDA mc particles")',inputs.n_pfida
    printf,55,f='("n_npa = ",i9,"    !! Number of Active NPA mc particles")',inputs.n_npa
    printf,55,f='("n_pnpa = ",i9,"    !! Number of Passive NPA mc particles")',inputs.n_pnpa
    printf,55,f='("n_nbi = ",i9,"    !! Number of NBI mc particles")',inputs.n_nbi
    printf,55,f='("n_halo = ",i9,"    !! Number of HALO mc particles")',inputs.n_halo
    printf,55,f='("n_dcx = ",i9,"     !! Number of DCX mc particles")',inputs.n_dcx
    printf,55,f='("n_birth = ",i9,"    !! Number of BIRTH mc particles")',inputs.n_birth
    printf,55,''
    printf,55,'!! Neutral Beam Settings'
    printf,55,f='("ab = ",1f9.5,"     !! Beam Species mass [amu]")',inputs.ab
    printf,55,f='("pinj = ",1f9.3,"     !! Beam Power [MW]")',inputs.pinj
    printf,55,f='("einj = ",1f9.3,"     !! Beam Energy [keV]")',inputs.einj
    printf,55,f='("current_fractions(1) = ",1f9.5," !! Current Fractions (Full component)")',inputs.current_fractions[0]
    printf,55,f='("current_fractions(2) = ",1f9.5," !! Current Fractions (Half component)")',inputs.current_fractions[1]
    printf,55,f='("current_fractions(3) = ",1f9.5," !! Current Fractions (Third component)")',inputs.current_fractions[2]
    printf,55,''
    printf,55,'!! Plasma Settings'
    printf,55,f='("ai = ",1f9.5,"     !! Ion Species mass [amu]")',inputs.ai
    printf,55,f='("impurity_charge = ",i3,"     !! Impurity Charge")',inputs.impurity_charge
    printf,55,''
    printf,55,'!! Beam Grid Settings'
    printf,55,f='("nx = ",i4,"    !! Number of cells in X direction (Into Plasma)")',inputs.nx
    printf,55,f='("ny = ",i4,"    !! Number of cells in Y direction")',inputs.ny
    printf,55,f='("nz = ",i4,"    !! Number of cells in Z direction")',inputs.nz
    printf,55,f='("xmin = ",1f9.3,"     !! Minimum X value [cm]")',inputs.xmin
    printf,55,f='("xmax = ",1f9.3,"     !! Maximum X value [cm]")',inputs.xmax
    printf,55,f='("ymin = ",1f9.3,"     !! Minimum Y value [cm]")',inputs.ymin
    printf,55,f='("ymax = ",1f9.3,"     !! Maximum Y value [cm]")',inputs.ymax
    printf,55,f='("zmin = ",1f9.3,"     !! Minimum Z value [cm]")',inputs.zmin
    printf,55,f='("zmax = ",1f9.3,"     !! Maximum Z value [cm]")',inputs.zmax
    printf,55,'!! Tait-Bryan Angles for z-y`-x`` rotation'
    printf,55,f='("alpha = ",1f9.5,"     !! Rotation about z-axis [rad]")',inputs.alpha
    printf,55,f='("beta  = ",1f9.5,"     !! Rotation about y`-axis [rad]")',inputs.beta
    printf,55,f='("gamma = ",1f9.5,"     !! Rotation about x``-axis [rad]")',inputs.gamma
    printf,55,'!! Beam Grid origin in machine coordinates (cartesian)'
    printf,55,f='("origin(1) = ",1f9.3,"     !! U value [cm]")',inputs.origin[0]
    printf,55,f='("origin(2) = ",1f9.3,"     !! V value [cm]")',inputs.origin[1]
    printf,55,f='("origin(3) = ",1f9.3,"     !! W value [cm]")',inputs.origin[2]
    printf,55,''
    printf,55,'!! Wavelength Grid Settings'
    printf,55,f='("nlambda = ",1i5,"    !! Number of Wavelengths")',inputs.nlambda
    printf,55,f='("lambdamin = ",1f9.3,"    !! Minimum Wavelength [nm]")',inputs.lambdamin
    printf,55,f='("lambdamax = ",1f9.3,"    !! Maximum Wavelength [nm]")',inputs.lambdamax
    printf,55,''
    printf,55,'!! Weight Function Settings'
    printf,55,f='("ne_wght = ",i9,"    !! Number of Energies for Weights")',inputs.ne_wght
    printf,55,f='("np_wght = ",i9,"    !! Number of Pitches for Weights")',inputs.np_wght
    printf,55,f='("nphi_wght = ",i9,"    !! Number of Gyro-angles for Weights")',inputs.nphi_wght
    printf,55,f='("emax_wght = ",1f9.2,"    !! Maximum Energy for Weights [keV]")',inputs.emax_wght
    printf,55,f='("nlambda_wght = ",1i5,"    !! Number of Wavelengths for Weights ")',$
              inputs.nlambda_wght
    printf,55,f='("lambdamin_wght = ",1f9.3,"    !! Minimum Wavelength for Weights [nm]")',$
              inputs.lambdamin_wght
    printf,55,f='("lambdamax_wght = ",1f9.3,"    !! Maximum Wavelength for Weights [nm]")',$
              inputs.lambdamax_wght
    printf,55,''
    printf,55,'/'
    printf,55,''
    close,55
    success,'Namelist file created: '+filename
END

