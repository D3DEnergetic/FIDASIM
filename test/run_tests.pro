PRO run_tests,result_dir,test_case=test_case


   if not keyword_set(test_case) then begin
     test_case = 'TEST_1A'
   endif else begin
     test_case = strupcase(test_case)
   endelse

   fida_dir = GETENV('FIDASIM_DIR')
   test_dir = fida_dir+'/test/'
   einj = double(72.5)
   pinj = double(1.7)

   cgfitf=[-0.109171,0.0144685,-7.83224e-5]
   cgfith=[0.0841037,0.00255160,-7.42683e-8]
   ffracs=cgfitf[0]+cgfitf[1]*einj+cgfitf[2]*einj^2
   hfracs=cgfith[0]+cgfith[1]*einj+cgfith[2]*einj^2
   tfracs=1.0-ffracs-hfracs
   species_mix = double([ffracs,hfracs,tfracs])

   basic_inputs = {max_threads:16,device:"TEST",shot:1,time:1.d0,$
                   einj:einj,pinj:pinj,species_mix:species_mix,$
                   ab:2.01410178d0,ai:2.0141078d0,impurity_charge:6,$
                   lambdamin:647.0d0,lambdamax:667.0d0,nlambda:2000,$
                   n_fida:5000000L,n_npa:500000L,n_nbi:50000L, $
                   n_halo:500000L,n_dcx:500000L,n_birth:10000L,$
                   ne_wght:50,np_wght:50,nphi_wght:100,emax_wght:100.0d0,$
                   nlambda_wght:1000,lambdamin_wght:647.d0,lambdamax_wght:667.d0,$
                   calc_npa:1,calc_brems:1,calc_bes:1,calc_fida:1,$
                   calc_birth:1,calc_fida_wght:1,calc_npa_wght:1,$
                   load_neutrals:0,dump_dcx:1,verbose:1,$
                   install_dir:fida_dir,result_dir:result_dir,tables_file:fida_dir+'/tables/atomic_tables.h5'}

   ;; Non-rotated, Non-tilted Grid
   grid01 = {nx:50,ny:60,nz:70,$
             xmin:-50.d0,xmax:50.d0,$
             ymin:-230.d0,ymax:-110.d0,$
             zmin:-70.d0,zmax:70.d0,$
             alpha:0.d0,beta:0.d0,gamma:0.d0,$
             origin:[0.d0,0.d0,0.d0]}

   ;; Rotated, Non-tilted
   grid02 = {nx:60,ny:50,nz:70,$
             xmin:0.d0,xmax:120.d0,$
             ymin:-50.d0,ymax:50.d0,$
             zmin:-70.d0,zmax:70.d0,$
             alpha:!DPI/2.0,beta:0.d0,gamma:0.d0,$
             origin:[0.d0,-230.d0,0.d0]}
      
   ;; Rotated, Tilted down
   grid03 = {nx:60,ny:50,nz:70,$
             xmin:0.d0,xmax:120.d0,$
             ymin:-50.d0,ymax:50.d0,$
             zmin:-70.d0,zmax:70.d0,$
             alpha:!DPI/2.0,beta:0.25d0,gamma:0.d0,$
             origin:[0.d0,-230.d0,0.d0]}
   
   ;; Rotated, Tilted up
   grid04 = {nx:60,ny:50,nz:70,$
             xmin:0.d0,xmax:120.d0,$
             ymin:-50.d0,ymax:50.d0,$
             zmin:-70.d0,zmax:70.d0,$
             alpha:!DPI/2.0,beta:-0.25d,gamma:0.d0,$
             origin:[0.d0,-230.d0,0.d0]}
   
   inputs01a = create_struct("runid","test_1a",$
                             "comment",'Non-rotated, Non-tilted grid; flat profiles',$
                             grid01,basic_inputs)
                          
   inputs01b = create_struct("runid",'test_1b',$
                             "comment",'Non-rotated, Non-tilted grid; realistic profiles',$
                             grid01,basic_inputs)

   inputs02a = create_struct("runid",'test_2a',$
                             "comment",'Rotated, Non-tilted grid; flat profiles',$
                             grid02,basic_inputs)
                          
   inputs02b = create_struct("runid",'test_2b',$
                             "comment",'Rotated, Non-tilted grid; realistic profiles',$
                             grid02,basic_inputs)

   inputs03a = create_struct("runid",'test_3a',$
                             "comment",'Rotated, Tilted down grid; flat profiles',$
                             grid03,basic_inputs)
                          
   inputs03b = create_struct("runid",'test_3b',$
                             "comment",'Rotated, Tilted down grid; realistic profiles',$
                             grid03,basic_inputs)

   inputs04a = create_struct("runid",'test_4a',$
                             "comment",'Rotated, Tilted up grid; flat profiles',$
                             grid04,basic_inputs)
                          
   inputs04b = create_struct("runid",'test_4b',$
                             "comment",'Rotated, Tilted down grid; realistic profiles',$
                             grid04,basic_inputs)

   test_cases = 'TEST_' + ['1A','1B','2A','2B','3A','3B','4A','4B']
   input_str = [inputs01a,inputs01b, $
                inputs02a,inputs02b, $
                inputs03a,inputs03b, $
                inputs04a,inputs04b]

   run = intarr(8)
   CASE test_case OF
       'TEST_1A': run[0] = 1
       'TEST_1B': run[1] = 1
       'TEST_2A': run[2] = 1
       'TEST_2B': run[3] = 1
       'TEST_3A': run[4] = 1
       'TEST_3B': run[5] = 1
       'TEST_4A': run[6] = 1
       'TEST_4B': run[7] = 1
       'ALL': run[*] = 1
       ELSE: BEGIN
           PRINT, 'Unknown test case: ',test_case
           return
       END
   ENDCASE

   for i=0,7 do begin
       if run[i] eq 0 then continue
       PRINT, 'Preparing test case '+test_cases[i]
       PRINT, input_str[i].comment
       grid = rz_grid(100.d0,240.d0, 70, -100.d0,100.d0, 100) 
       
       fbm = read_nubeam(test_dir+'test_fi_1.cdf',grid,$
                         btipsign=-1.0,$
                         e_range=[67.0,77.0], $
                         p_range=[-0.1,0.1])

       spec = test_chords()
       npa = test_npa()

       equil = read_geqdsk(test_dir+'g000001.01000',grid,flux=flux)
       nbi = test_beam(input_str[i].beta)

       tcb = byte(strlowcase(test_cases[i]))
       pfile = test_dir+'test_profiles_'+string(tcb[-1])+'.cdf'
       plasma = test_profiles(pfile,grid,flux)
       
       prefida,input_str[i], grid, nbi, plasma, equil, fbm, spec=spec, npa=npa
   endfor

end    
