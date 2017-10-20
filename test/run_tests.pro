PRO run_tests,result_dir,test_case=test_case


   if not keyword_set(test_case) then begin
     test_case = 'TEST_1A'
   endif else begin
     test_case = strupcase(test_case)
   endelse

   fida_dir = get_fidasim_dir()
   test_dir = fida_dir+'/test/'
   einj = double(72.5)
   pinj = double(1.7)

   cgfitf=[-0.109171,0.0144685,-7.83224e-5]
   cgfith=[0.0841037,0.00255160,-7.42683e-8]
   ffracs=cgfitf[0]+cgfitf[1]*einj+cgfitf[2]*einj^2
   hfracs=cgfith[0]+cgfith[1]*einj+cgfith[2]*einj^2
   tfracs=1.0-ffracs-hfracs
   current_fractions = double([ffracs,hfracs,tfracs])

   basic_inputs = {device:"TEST",shot:1L,time:1.d0,$
                   einj:einj,pinj:pinj,current_fractions:current_fractions,$
                   ab:2.01410178d0,ai:2.0141078d0,impurity_charge:6,$
                   lambdamin:647.0d0,lambdamax:667.0d0,nlambda:2000,$
                   n_fida:5000000L,n_npa:5000000L,n_nbi:50000L, $
                   n_halo:500000L,n_dcx:500000L,n_birth:10000L,$
                   ne_wght:50,np_wght:50,nphi_wght:100,emax_wght:100.0d0,$
                   nlambda_wght:1000,lambdamin_wght:647.d0,lambdamax_wght:667.d0,$
                   calc_npa:1,calc_brems:1,calc_bes:1,calc_fida:1,calc_neutron:1,$
                   calc_birth:1,calc_fida_wght:1,calc_npa_wght:1,dump_dcx:1,$
                   result_dir:result_dir,tables_file:fida_dir+'/tables/atomic_tables.h5'}

   basic_bgrid = {nx:50,ny:60,nz:70,$
                  xmin:-50.d0,xmax:50.d0,$
                  ymin:-230.d0,ymax:-110.d0,$
                  zmin:-70.d0,zmax:70.d0,$
                  alpha:0.d0,beta:0.d0,gamma:0.d0,$
                  origin:[0.d0,0.d0,0.d0]}

   test_ids = ['1A','1B','2A','2B','3A','3B','4A','4B']
   test_cases = 'TEST_' + test_ids 
   test_cases = [test_cases,'ALL']
   run = strmatch(test_cases,test_case,/fold_case)
   if total(run) eq 0 then begin
       error, "Unknown test case",/halt
   endif
   if run[-1] eq 1 then run[*]=1

   for i=0,n_elements(test_ids)-1 do begin
       if run[i] eq 0 then continue
       PRINT, 'Preparing test case '+test_cases[i]
       CASE test_cases[i] OF
           "TEST_1A": BEGIN
               inputs = create_struct("runid","test_1a",$
                        "comment",'Non-rotated, Non-tilted grid; flat profiles',$
                        basic_inputs, basic_bgrid)
               nbi = test_beam(0.d0)
           END
           "TEST_1B": BEGIN
               inputs = create_struct("runid",'test_1b',$
                        "comment",'Non-rotated, Non-tilted grid; realistic profiles',$
                        basic_inputs, basic_bgrid)
               nbi = test_beam(0.d0)
           END
           "TEST_2A": BEGIN
               nbi = test_beam(0.d0)
               bgrid = beam_grid(nbi,230.d0, nx=60, ny=50, nz=70, $
                       length=120.d0,width=100.d0,height=140.d0)
               inputs= create_struct("runid",'test_2a',$
                       "comment",'Rotated, Non-tilted grid; flat profiles',$
                       basic_inputs, bgrid)
           END
           "TEST_2B": BEGIN
               nbi = test_beam(0.d0)
               bgrid = beam_grid(nbi,230.d0, nx=60, ny=50, nz=70, $
                       length=120.d0,width=100.d0,height=140.d0)
               inputs = create_struct("runid",'test_2b',$
                        "comment",'Rotated, Non-tilted grid; realistic profiles',$
                        basic_inputs, bgrid)
           END
           "TEST_3A": BEGIN
               nbi = test_beam(0.25)
               bgrid = beam_grid(nbi,230.d0, nx=60, ny=50, nz=70, $
                       length=120.d0,width=100.d0,height=140.d0)
               inputs = create_struct("runid",'test_3a',$
                        "comment",'Rotated, Tilted down grid; flat profiles',$
                        basic_inputs, bgrid)
           END
           "TEST_3B": BEGIN
               nbi = test_beam(0.25)
               bgrid = beam_grid(nbi,230.d0, nx=60, ny=50, nz=70, $
                       length=120.d0,width=100.d0,height=140.d0)
               inputs = create_struct("runid",'test_3b',$
                        "comment",'Rotated, Tilted down grid; realistic profiles',$
                        basic_inputs, bgrid)

           END
           "TEST_4A": BEGIN
               nbi = test_beam(-0.25)
               bgrid = beam_grid(nbi,230.d0, nx=60, ny=50, nz=70, $
                       length=120.d0,width=100.d0,height=140.d0)
               inputs = create_struct("runid",'test_4a',$
                        "comment",'Rotated, Tilted up grid; flat profiles',$
                        basic_inputs, bgrid)
                          
           END
           "TEST_4B": BEGIN
               nbi = test_beam(-0.25)
               bgrid = beam_grid(nbi,230.d0, nx=60, ny=50, nz=70, $
                       length=120.d0,width=100.d0,height=140.d0)
               inputs = create_struct("runid",'test_4b',$
                        "comment",'Rotated, Tilted down grid; realistic profiles',$
                        basic_inputs, bgrid)
           END
       ENDCASE 
       PRINT, inputs.comment

       grid = rz_grid(100.d0,240.d0, 70, -100.d0,100.d0, 100) 
       equil = read_geqdsk(test_dir+'g000001.01000',grid,flux=flux,g=g)
       equil = create_struct(equil,"geqdsk",g)
       fbm = read_nubeam(test_dir+'test_fi_1.cdf',grid,$
                         btipsign=-1.0,$
                         e_range=[67.0,77.0], $
                         p_range=[-0.1,0.1])
       spec = test_chords()
       npa = test_npa()

       tcb = byte(strlowcase(test_cases[i]))
       pfile = test_dir+'test_profiles_'+string(tcb[-1])+'.cdf'
       plasma = test_profiles(pfile,grid,flux)
       
       prefida,inputs, grid, nbi, plasma, equil, fbm, spec=spec, npa=npa
   endfor

end    
