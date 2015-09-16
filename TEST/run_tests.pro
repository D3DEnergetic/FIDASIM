PRO run_tests,result_dir,test_case=test_case


   if not keyword_set(test_case) then begin
     test_case = 'TEST_1A'
   endif else begin
     test_case = strupcase(test_case)
   endelse

   fida_dir = GETENV('FIDASIM_DIR')
   inputs = {device:"TEST",shot:1,time:1.0,$
             emin:0.0,emax:100.0,pmin:-1.0,pmax:1.0,$
             einj:72.5,pinj:1.7,equil:'g000001.01000',btipsign:-1.0,$
             ab:2.01410178,ai:2.0141078,impurity_charge:6,$
             lambdamin:647.0,lambdamax:667.0,nlambda:2000,dlambda:0.01,$
             nr_fast:5000000,nr_nbi:50000,nr_halo:500000,$
             ne_wght:50,np_wght:50,nphi_wght:100,emax_wght:100,ichan_wght:-1,$
             dwav_wght:0.02,wavel_start_wght:647.0,wavel_end_wght:667.0,$
             calc_npa:0,calc_spec:1,calc_birth:1,calc_brems:1,calc_fida_wght:1,$
             calc_npa_wght:1,load_neutrals:0,load_fbm:1,interactive:0,$
             install_dir:fida_dir,diag:["ALL"],isource:0,$
             result_dir:result_dir,cdf_file:fida_dir+'TEST/test_fi_1.cdf'}
 
   inter_grid = {nr:70,nw:100,$
                 rmin:100.0,rmax:240,$
                 wmin:-100.0,wmax:100.0}

   ;; Non-rotated, Non-tilted Grid
   grid01 = {nx:50,ny:60,nz:70,$
             xmin:-50.0,xmax:50.0,$
             ymin:-230.0,ymax:-110.0,$
             zmin:-70.0,zmax:70.0,$
             alpha:0.0,beta:0.0,$
             origin:[0.0,0.0,0.0]}

   ;; Rotated, Non-tilted
   grid02 = {nx:60,ny:50,nz:70,$
             xmin:0.0,xmax:120.0,$
             ymin:-50.0,ymax:50.0,$
             zmin:-70.0,zmax:70.0,$
             alpha:!DPI/2.0,beta:0.0,$
             origin:[0.0,-230.0,0.0]}
      
   ;; Rotated, Tilted down
   grid03 = {nx:60,ny:50,nz:70,$
             xmin:0.0,xmax:120.0,$
             ymin:-50.0,ymax:50.0,$
             zmin:-70.0,zmax:70.0,$
             alpha:!DPI/2.0,beta:0.25,$
             origin:[0.0,-230.0,0.0]}
   
   ;; Rotated, Tilted up
   grid04 = {nx:60,ny:50,nz:70,$
             xmin:0.0,xmax:120.0,$
             ymin:-50.0,ymax:50.0,$
             zmin:-70.0,zmax:70.0,$
             alpha:!DPI/2.0,beta:-0.25,$
             origin:[0.0,-230.0,0.0]}
   
   inputs01a = create_struct("profile_dir",fida_dir+'/TEST/test_profiles_a.cdf',$
                             "runid",'test_1a',$
                             "comment",'Non-rotated, Non-tilted grid; flat profiles',$
                             inter_grid,grid01,inputs)
                          
   inputs01b = create_struct("profile_dir",fida_dir+'/TEST/test_profiles_b.cdf',$
                             "runid",'test_1b',$
                             "comment",'Non-rotated, Non-tilted grid; realistic profiles',$
                             inter_grid,grid01,inputs)

   inputs02a = create_struct("profile_dir",fida_dir+'/TEST/test_profiles_a.cdf',$
                             "runid",'test_2a',$
                             "comment",'Rotated, Non-tilted grid; flat profiles',$
                             inter_grid,grid02,inputs)
                          
   inputs02b = create_struct("profile_dir",fida_dir+'/TEST/test_profiles_b.cdf',$
                             "runid",'test_2b',$
                             "comment",'Rotated, Non-tilted grid; realistic profiles',$
                             inter_grid,grid02,inputs)

   inputs03a = create_struct("profile_dir",fida_dir+'/TEST/test_profiles_a.cdf',$
                             "runid",'test_3a',$
                             "comment",'Rotated, Tilted down grid; flat profiles',$
                             inter_grid,grid03,inputs)
                          
   inputs03b = create_struct("profile_dir",fida_dir+'/TEST/test_profiles_b.cdf',$
                             "runid",'test_3b',$
                             "comment",'Rotated, Tilted down grid; realistic profiles',$
                             inter_grid,grid03,inputs)

   inputs04a = create_struct("profile_dir",fida_dir+'/TEST/test_profiles_a.cdf',$
                             "runid",'test_4a',$
                             "comment",'Rotated, Tilted up grid; flat profiles',$
                             inter_grid,grid04,inputs)
                          
   inputs04b = create_struct("profile_dir",fida_dir+'/TEST/test_profiles_b.cdf',$
                             "runid",'test_4b',$
                             "comment",'Rotated, Tilted down grid; realistic profiles',$
                             inter_grid,grid04,inputs)

   CASE test_case OF
       'TEST_1A': input_str = inputs01a
       'TEST_1B': input_str = inputs01b
       'TEST_2A': input_str = inputs02a
       'TEST_2B': input_str = inputs02b
       'TEST_3A': input_str = inputs03a
       'TEST_3B': input_str = inputs03b
       'TEST_4A': input_str = inputs04a    
       'TEST_4B': input_str = inputs04b
       ELSE: PRINT, 'Unknown test case: ',test_case
   ENDCASE

   PRINT, 'Preparing test case: ',test_case
   PRINT, input_str.comment
   prefida,input_str=input_str

end    
