PRO write_geometry, filename, nbi, spec=spec, npa=npa
    ;+##`write_geometry, filename, nbi, spec=spec, npa=npa`
    ;+Write geometry values to a HDF5 file
    ;+
    ;+###Input Arguments
    ;+     **filename**: Name of the geometry file
    ;+
    ;+     **nbi**: NBI geometry structure
    ;+
    ;+###Keyword Arguments
    ;+     **spec**: Optional, Spectral geometry structure
    ;+
    ;+     **npa**: Optional, NPA geometry structure
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> write_geometry, filename, nbi, spec=spec, npa=npa
    ;+```
    info,'Writing geometry file...'

    ;; Create attributes
    root_atts = {attribute,obj:'/', $
                name:'description',$
                data:'Geometric quantities for FIDASIM'}

    ;; NBI attributes
    nbi_desc = {attribute,obj:'/nbi', $
                name:'description',$
                data:'Neutral Beam Geometry'}

    nbi_cs = {attribute,obj:'/nbi', $
              name:'coordinate_system',$
              data:'Right-handed cartesian'}

    nbi_ds_desc = {attribute,obj:'/nbi/data_source', $
                   name:'description', $
                   data:'Source of the NBI geometry'}

    nbi_name_desc = {attribute,obj:'/nbi/name', $
                     name:'description', $
                     data:'Beam name'}

    nbi_src_desc = {attribute,obj:'/nbi/src', $
                    name:'description', $
                    data:'Position of the center of the beam source grid'}
    nbi_src_unit = {attribute,obj:'/nbi/src', $
                    name:'units', $
                    data:'cm'}

    nbi_axis_desc = {attribute,obj:'/nbi/axis', $
                     name:'description', $
                     data:'Axis of the beam centerline: Centerline(t) = src + axis*t '}
    nbi_axis_unit = {attribute,obj:'/nbi/axis', $
                     name:'units', $
                     data:'cm'}

    nbi_focy_desc = {attribute,obj:'/nbi/focy', $
                     name:'description', $
                     data:'Horizonal focal length of the beam'}
    nbi_focy_unit = {attribute,obj:'/nbi/focy', $
                     name:'units', $
                     data:'cm'}

    nbi_focz_desc = {attribute,obj:'/nbi/focz', $
                     name:'description', $
                     data:'Vertical focal length of the beam'}
    nbi_focz_unit = {attribute,obj:'/nbi/focz', $
                     name:'units', $
                     data:'cm'}
     
    nbi_divy_desc = {attribute,obj:'/nbi/divy', $
                     name:'description', $
                     data:'Horizonal divergences of the beam. One for each energy component'}
    nbi_divy_unit = {attribute,obj:'/nbi/divy', $
                     name:'units', $
                     data:'radians'}

    nbi_divz_desc = {attribute,obj:'/nbi/divz', $
                     name:'description', $
                     data:'Vertical divergences of the beam. One for each energy component'}
    nbi_divz_unit = {attribute,obj:'/nbi/divz', $
                     name:'units', $
                     data:'radians'}

    nbi_widy_desc = {attribute,obj:'/nbi/widy', $
                     name:'description', $
                     data:'Half width of the beam source grid'}
    nbi_widy_unit = {attribute,obj:'/nbi/widy', $
                     name:'units', $
                     data:'cm'}

    nbi_widz_desc = {attribute,obj:'/nbi/widz', $
                     name:'description', $
                     data:'Half height of the beam source grid'}
    nbi_widz_unit = {attribute,obj:'/nbi/widz', $
                     name:'units', $
                     data:'cm'}

    nbi_shape_desc = {attribute,obj:'/nbi/shape', $
                      name:'description', $
                      data:'Shape of the beam source grid: 1="rectangular", 2="circular"'}

    nbi_atts = [nbi_desc, nbi_cs, nbi_ds_desc,nbi_name_desc, $
                nbi_src_desc, nbi_src_unit, $
                nbi_axis_desc, nbi_axis_unit, $
                nbi_focy_desc, nbi_focy_unit, $
                nbi_focz_desc, nbi_focz_unit, $
                nbi_divy_desc, nbi_divy_unit, $
                nbi_divz_desc, nbi_divz_unit, $
                nbi_widy_desc, nbi_widy_unit, $
                nbi_widz_desc, nbi_widz_unit, $
                nbi_shape_desc]

    ;; Spectroscopic attributes
    spec_desc = {attribute,obj:'/spec', $
                 name:'description', $
                 data:'FIDA/BES Chord Geometry'}
    spec_cs = {attribute,obj:'/spec', $
               name:'coordinate_system', $
               data:'Right-handed cartesian'}

    spec_ds_desc = {attribute,obj:'/spec/data_source', $
                    name:'description', $
                    data:'Source of the chord geometry'}

    spec_nchan_desc = {attribute,obj:'/spec/nchan', $
                       name:'description', $
                       data:'Number of channels'}

    spec_system_desc = {attribute,obj:'/spec/system', $
                        name:'description', $
                        data:'Names of the different spectrocopic systems'}

    spec_id_desc = {attribute,obj:'/spec/id', $
                    name:'description', $
                    data:'Line of sight ID'}

    spec_lens_desc = {attribute,obj:'/spec/lens', $
                      name:'description', $
                      data:'Positions of the lenses'}
    spec_lens_unit = {attribute,obj:'/spec/lens', $
                      name:'units', $
                      data:'cm'}

    spec_axis_desc = {attribute,obj:'/spec/axis', $
                      name:'description', $
                      data:'Optical axis of the lines of sight: LOS(t) = lens + axis*t '}
    spec_axis_unit = {attribute,obj:'/spec/axis', $
                      name:'units', $
                      data:'cm'}

    spec_radius_desc = {attribute,obj:'/spec/radius', $
                        name:'description', $
                        data:'Line of sight radius at midplane or tangency point' }
    spec_radius_unit = {attribute,obj:'/spec/radius', $
                        name:'units', $
                        data:'cm'}

    spec_sigma_desc = {attribute,obj:'/spec/sigma_pi', $
                       name:'description', $
                       data:'Ratio of the intensities of the sigma and pi stark lines. Measured quantity'}

    spec_spot_desc = {attribute,obj:'/spec/spot_size', $
                      name:'description', $
                      data:'Radius of spot size'}
    spec_spot_unit = {attribute,obj:'/spec/spot_size', $
                      name:'units', $
                      data:'cm'}
 
    spec_atts = [spec_desc,spec_cs, spec_ds_desc, $
                 spec_nchan_desc, spec_system_desc, $
                 spec_id_desc, $
                 spec_lens_desc,spec_lens_unit, $
                 spec_axis_desc,spec_axis_unit, $
                 spec_radius_desc, spec_radius_unit, $
                 spec_sigma_desc, $
                 spec_spot_desc, spec_spot_unit]

    ;;NPA attributes
    npa_desc = {attribute,obj:'/npa', $
                name:'description', $
                data:'NPA Geometry'}

    npa_cs = {attribute,obj:'/npa', $
              name:'coordinate_system', $
              data:'Right-handed cartesian'}

    npa_ds_desc = {attribute,obj:'/npa/data_source', $
                   name:'description', $
                   data:'Source of the NPA geometry'}
    
    npa_nchan_desc = {attribute,obj:'/npa/nchan', $
                      name:'description', $
                      data:'Number of channels'}

    npa_system_desc = {attribute,obj:'/npa/system', $
                       name:'description', $
                       data:'Names of the different NPA systems'}

    npa_id_desc = {attribute,obj:'/npa/id', $
                   name:'description', $
                   data:'Line of sight ID'}

    npa_dshape_desc = {attribute,obj:'/npa/d_shape', $
                       name:'description', $
                       data:'Shape of the detector: 1="rectangular", 2="circular"'}

    npa_dcent_desc = {attribute,obj:'/npa/d_cent', $
                      name:'description', $
                      data:'Center of the detector'}
    npa_dcent_unit = {attribute,obj:'/npa/d_cent', $
                      name:'units', $
                      data:'cm'}

    npa_dtedge_desc = {attribute,obj:'/npa/d_tedge', $
                       name:'description', $
                       data:'Center of the detectors top edge'}
    npa_dtedge_unit = {attribute,obj:'/npa/d_tedge', $
                       name:'units', $
                       data:'cm'}

    npa_dredge_desc = {attribute,obj:'/npa/d_redge', $
                       name:'description', $
                       data:'Center of the detectors right edge'}
    npa_dredge_unit = {attribute,obj:'/npa/d_redge', $
                       name:'units', $
                       data:'cm'}

    npa_ashape_desc = {attribute,obj:'/npa/a_shape', $
                       name:'description', $
                       data:'Shape of the aperture: 1="rectangular", 2="circular"'}

    npa_acent_desc = {attribute,obj:'/npa/a_cent', $
                      name:'description', $
                      data:'Center of the aperture'}
    npa_acent_unit = {attribute,obj:'/npa/a_cent', $
                      name:'units', $
                      data:'cm'}

    npa_atedge_desc = {attribute,obj:'/npa/a_tedge', $
                       name:'description', $
                       data:'Center of the apertures top edge'}
    npa_atedge_unit = {attribute,obj:'/npa/a_tedge', $
                       name:'units', $
                       data:'cm'}

    npa_aredge_desc = {attribute,obj:'/npa/a_redge', $
                       name:'description', $
                       data:'Center of the apertures right edge'}
    npa_aredge_unit = {attribute,obj:'/npa/a_redge', $
                       name:'units', $
                       data:'cm'}

    npa_radius_desc = {attribute,obj:'/npa/radius', $
                       name:'description', $
                       data:'Line of sight radius at midplane or tangency point' }
    npa_radius_unit = {attribute,obj:'/npa/radius', $
                       name:'units', $
                       data:'cm'}

    npa_atts = [npa_desc, npa_cs, npa_ds_desc, $ 
                npa_nchan_desc, npa_system_desc, $
                npa_id_desc, $
                npa_dshape_desc, npa_ashape_desc, $
                npa_dcent_desc, npa_dcent_unit, $
                npa_acent_desc, npa_acent_unit, $
                npa_dtedge_desc, npa_dtedge_unit, $
                npa_atedge_desc, npa_atedge_unit, $
                npa_dredge_desc, npa_dredge_unit, $
                npa_aredge_desc, npa_aredge_unit, $
                npa_radius_desc, npa_radius_unit ]

    atts = [root_atts, nbi_atts]
    geom = {nbi:nbi}

    if keyword_set(spec) then begin
        geom = create_struct(geom, "spec", spec)
        atts = [atts, spec_atts]
    endif

    if keyword_set(npa) then begin
        geom = create_struct(geom, "npa", npa)
        atts = [atts,npa_atts]
    endif

    write_hdf5, geom, filename=filename, atts=atts, /clobber

    if file_test(filename) then begin
        success, 'Geometry file created: '+filename
    endif else begin
        error, 'Geometry file creation failed.'
    endelse

END

