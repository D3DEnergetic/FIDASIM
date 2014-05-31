PRO nstx_input,inputs
;;=============================================================================
;;PREFIDA INPUT FILE FOR NSTX or NSTX-U
;;=============================================================================
shot=122631L                           ;; Shot Number
time=0.1                               ;; Time in second
runid='122631L21_fi1_tssnpa'  ;; runid of FIDASIM
device='NSTX'                          ;; Device name: D3D,NSTX,AUGD,MAST
install_dir='/u/dliu/FIDASIM_140304/'  ;; Location of FIDAsim code & executable
result_dir='/p/fida/NSTX/'             ;; Location where results will be stored,
                                       ;; /RESULTS/runid dir will be created
profile_dir='/p/nstxusr/nstx-users/dliu/FIDAsim_prep/NSTX/122631/'
                                       ;; Location of profile save file

;;----------------------------------------------------
;; Fast-ion distribution function from transp
;;----------------------------------------------------
cdf_file='/p/nstxusr/nstx-users/dliu/FIDAsim_prep/NSTX/122631/122631L21_fi_1.cdf'
                        ;; CDF file from transp with the distribution funciton
emin=10.                ;; Minimum energy used from the distribution function
emax=100.               ;; Maximum energy used from the distribution function
pmin=-1.                ;; Minimum pitch used from the distribution function
pmax=1.                 ;; Maximum pitch used from the distribution function

;;-----------------------------------------------------
;; Beam/FIDA/EQUILIBRIUM Selection
;;-----------------------------------------------------
isource=4                  ;; Beam source index, the index starts from 0 to 5
                           ;; FIDASIM only simulates one beam source per run
einj=90.5                  ;; Beam energy in [keV]. If 0, get data from MDS+
pinj=2.12                  ;; Beam power in [MW}.If 0, get data from MDS+

diag=['NSTXU_TSSNPA']   ;; Name of the FIDA diag
equil=''                   ;; Name of equilibrium. Ex. for NSTX EFIT02.
                           ;; If empty, look for gfile in profile_dir

;;-----------------------------------------------------
;; Discharge Parameters
;;-----------------------------------------------------
btipsign=-1.d0             ;; Bt and Ip are in the opposite direction   
ab=2.01410178d0            ;; Atomic mass of beam [u]
ai=2.01410178d0            ;; Atomic mass of hydrogenic plasma ions [u]
impurity_charge=6          ;; 5: BORON, 6: carbon, 7: Nitrogen

;;-----------------------------------------------------
;; Wavelength Grid
;;-----------------------------------------------------
lambdamin=648.8d0           ;; Minimum wavelength of wavelength grid[nm]
lambdamax=663.8d0           ;; Maximum wavelength of wavelength grid[nm]
nlambda=1500L               ;; Number of wavelengths
dlambda= abs(lambdamax-lambdamin)/double(nlambda)    ;; Wavelength seperation

;;---------------------------------------------------
;; Define FIDASIM origin and grids in machine coordinates(x,y,z)
;;---------------------------------------------------
;rc_onaxis=158.58             ;major radius
;phic_onaxis=79.298/180.*!pi  ;angle
;origin=[rc_onaxis*cos(phic_onaxis),rc_onaxis*sin(phic_onaxis),0.d]
rc_offaxis=158.58             ;major radius
phic_offaxis=286.2/180.*!pi-!pi/2-atan(120./rc_offaxis)
origin=[rc_offaxis*cos(phic_offaxis),rc_offaxis*sin(phic_offaxis),0.d]

;NB source       0        1         2         3         4        5
;angle of NBI 236.700,  240.200,  243.700,  283.000,  286.500,  290.000
;origin=[29.4485, 155.822,0.d] ;; Cross-over point was chosen as (x,y,z) origin 
alpha=286.2/180.*!pi         ;; Rotation angle in radians from +x about z axis
                              ;; that transforms machine coordinates to (x,y,z)
                              ;; The centerline of NB was chosen as x axis
beta=0.0                      ;; Rotation about z axis
nx=71                         ;; Number of cells in x direction
ny=41                         ;; Number of cells in y direction
nz=71                         ;; Number of cells in z direction
xmin=-20.                    ;; Minimum x value
xmax=120.                    ;; Maximum x value
ymin=-40.                    ;; Minimum y value
ymax=45.                     ;; Maximum y value
zmin=-70.                    ;; Minimum z value
zmax=70.                     ;; Maximum z value

;;--------------------------------------------------
;; Define number of Monte Carlo particles
;;--------------------------------------------------
nr_fast=5000000L                ;; FIDA
nr_nbi=500000L                   ;; Beam emission
nr_halo=500000L                 ;; Halo contribution

;;--------------------------------------------------
;; Calculation of the weight function
;;--------------------------------------------------
ne_wght=50                    ;; Number of Pitches, energyies and gyro angles 
np_wght=50                    ;; Number of Pitches, energyies and gyro angles 
nphi_wght=50                  ;; Number of Pitches, energyies and gyro angles
emax_wght=100.                ;; Maximum energy (keV)
ichan_wght=-1                 ;; -1 for all channels,
                              ;; otherwise a given channel index
dwav_wght=.02                 ;; Wavelength interval
wavel_start_wght=648.         ;; Minimum wavelength
wavel_end_wght=663.           ;; Maximum wavelength

;;-------------------------------------------------
;; Simulation switches
;;-------------------------------------------------
calc_npa=[1]            ;; (0 or 1) If 1 do a simulation for NPA
calc_spec=[1]           ;; (0 or 1) If 1 then spectra is calculated
calc_birth=[1]          ;; (0 or 1) If then the birth profile is calculated
calc_brems=[0]            ;; (0 or 1) If 0 use the IDL bremstrahlung calculation
calc_fida_wght=[1]      ;; (0 or 1) If 1 then weight functions are calculated
calc_npa_wght=[1]       ;; (0 or 1) If 1 then weight functions are calculated
load_neutrals=[0]       ;; (0 or 1) If 1 then the neutral density is loaded 
                        ;;                    from an existing FIDAsim run  
load_fbm=[1]            ;; (0 or 1) If 1 then the fbm is loaded 
                        ;; (calc_spec/npa overwrites)

;;------------------------------------------------
;; DO NOT MODIFY THIS PART
;;------------------------------------------------
install_dir=+GETENV('FIDASIM_DIR')
;;device=strupcase(device)

inputs={shot:shot,time:time,runid:runid,device:device,install_dir:install_dir,$
        result_dir:result_dir,cdf_file:cdf_file,profile_dir:profile_dir,$   
        emin:emin,emax:emax,pmin:pmin,pmax:pmax,isource:isource,diag:diag,$
        einj:einj,pinj:pinj,equil:equil,btipsign:btipsign,$
        ab:ab,ai:ai,impurity_charge:impurity_charge,$
        lambdamin:lambdamin,lambdamax:lambdamax,nlambda:nlambda,dlambda:dlambda,$
        origin:origin,alpha:alpha,beta:beta,nx:nx,ny:ny,nz:nz,$
        xmin:xmin,xmax:xmax,ymin:ymin,ymax:ymax,zmin:zmin,zmax:zmax,$
        nr_fast:nr_fast,nr_nbi:nr_nbi,nr_halo:nr_halo,$
        ne_wght:ne_wght,np_wght:np_wght,nphi_wght:nphi_wght,$
        emax_wght:emax_wght,ichan_wght:ichan_wght,$
        dwav_wght:dwav_wght,wavel_start_wght:wavel_start_wght,$
        wavel_end_wght:wavel_end_wght,$
        calc_npa:calc_npa,calc_spec:calc_spec,$
        calc_birth:calc_birth,calc_brems:calc_brems,calc_fida_wght:calc_fida_wght,$
        calc_npa_wght:calc_npa_wght,$
        load_neutrals:load_neutrals,load_fbm:load_fbm }

END
