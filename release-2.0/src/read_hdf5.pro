FUNCTION valid_name, name, bad_names=bad_names, post = post

    reserved=['AND','BEGIN','BREAK','CASE','COMMON','COMPILE_OPT',$
              'CONTINUE','DO','ELSE','END','ENDCASE','ENDELSE',$
              'ENDFOR','ENDFOREACH','ENDIF','ENDREP','ENDSWITCH','ENDWHILE',$
              'EQ','FOR','FOREACH','FORWARD_FUNCTION','FUNCTION','GE',$
              'GOTO','GT','IF','INHERITS','LE','LT','MOD','NE','NOT','OF',$
              'ON_IOERROR','OR','PRO','REPEAT','SWITCH','THEN','UNTIL',$
              'WHILE','X0R']

    if not keyword_set(bad_names) then bad_names = ''
    if not keyword_set(post) then post = ''

    bad_names = [reserved, bad_names]
    if total(strmatch(bad_names, name,/fold_case)) ne 0 then begin
        valid_name = name + '_' + post
    endif else begin
        valid_name = name
    endelse

    ; added 2016-09-15 by NGB to fix invalid structure tags
    valid_name = idl_validname(valid_name, /convert_all)
    
    return, valid_name

END

FUNCTION create_nested_struct,path,data

     rpath = reverse(strsplit(path,'/',/extract))
     d = data
     for i=0,n_elements(rpath) -1 do begin
         varname = valid_name(rpath[i])
         d = create_struct(varname,d)
     endfor

     return, d
END

FUNCTION hdf5_read_attributes,id,bad_names=bad_names

    ;; Get any attributes
    natts = h5a_get_num_attrs(id)
    atts = {}
    for i=0L, natts-1 do begin
        ;; Open attribute id
        attribute_id = h5a_open_idx(id,i)

        ;; Get attribute name and make sure its valid
        attribute_name = h5a_get_name(attribute_id)
        attribute_name = valid_name(attribute_name, $
                         bad_names=bad_names, $
                         post=strcompress(string(i),/remove_all))

        ;; Get the attribute data
        attribute_data = h5a_read(attribute_id)
        atts = create_struct(atts,attribute_name,attribute_data)
       
        ;; Close attribute id
        h5a_close,attribute_id
    endfor

    return, atts
END

FUNCTION hdf5_read_dataset,id,name,shallow=shallow

    ;; Get data
    dataset_id = h5d_open(id,name)
    data = h5d_read(dataset_id)

    ;; Get any attributes
    atts = hdf5_read_attributes(dataset_id,bad_names="data")

    ;; Close the dataset
    h5d_close, dataset_id

    if keyword_set(shallow) then begin
        return, data
    endif else begin
        return, create_struct(atts, "data", data)
    endelse
END

FUNCTION hdf5_read_group,id,shallow=shallow

    FORWARD_FUNCTION hdf5_read_group

    nobjs = h5g_get_num_objs(id)
    d = {}
    for i = 0, nobjs-1 do begin
        obj_name = h5g_get_obj_name_by_idx(id, i)
        var_name = valid_name(obj_name)

        obj_info = h5g_get_objinfo(id, obj_name)
        obj_type = obj_info.type
 
        CASE obj_type OF
            'GROUP': BEGIN
                gid = h5g_open(id, obj_name)
                var = hdf5_read_group(gid,shallow=shallow)
                h5g_close, gid
                if n_elements(var) ne 0 then begin
                    d = create_struct(d, var_name, var)
                endif
            END
            'DATASET': BEGIN
                var = hdf5_read_dataset(id,obj_name,shallow=shallow)
                if n_elements(var) ne 0 then begin
                    d = create_struct(d, var_name, var)
                endif
            END
            'TYPE': BEGIN
                tid = h5t_open(id, obj_name)
                var = hdf5_read_attributes(tid)
                h5t_close, tid
                if n_elements(var) ne 0 then begin
                    d = create_struct(d, var_name, var)
                endif
            END
            ELSE:
        ENDCASE
    endfor

    if not keyword_set(shallow) then begin
        atts = hdf5_read_attributes(id)
    endif

    if n_elements(atts) ne 0 then begin 
        return, create_struct(d, atts)
    endif else begin
        return, d
    endelse
    
END

FUNCTION hdf5_read_from_list, id, var_paths, flatten=flatten,shallow=shallow

    d = {}
    used_names = []
    i = 0L
    while i lt n_elements(var_paths) do begin
        catch,err_status
        if err_status ne 0 then begin
            print,'Error reading '+var_paths[i]
            print,!ERROR_STATE.MSG
            catch,/cancel
            i = i + 1
            continue
        endif
        path = var_paths[i]
        obj_info = h5g_get_objinfo(id, path)
        obj_type = obj_info.type

        CASE obj_type OF
            'LINK': BEGIN
                var_name = h5g_get_linkval(id,path)
                var_paths = [var_paths,var_name]
             END
            'GROUP': BEGIN
                gid = h5g_open(id,path)
                var = hdf5_read_group(gid,shallow=shallow)
                h5g_close, gid
                if n_elements(var) ne 0 then begin
                    if keyword_set(flatten) then begin
                        var_names = strsplit(path,'/',/extract)
                        var_name = valid_name(var_names[-1], $
                                   bad_names=used_names, $
                                   post=strcompress(string(i),/remove_all))
                        used_names = [used_names,var_name]
                        d = create_struct(d, var_name, var)
                    endif else begin
                        d = create_struct(d, create_nested_struct(path,var))
                    endelse
                endif
            END
            'DATASET': BEGIN
                var = hdf5_read_dataset(id,path,shallow=shallow)
                if n_elements(var) ne 0 then begin
                    if keyword_set(flatten) then begin
                        var_names = strsplit(path,'/',/extract)
                        var_name = valid_name(var_names[-1], $
                                   bad_names=used_names, $
                                   post=strcompress(string(i),/remove_all))
                        used_names = [used_names, var_name]
                        d = create_struct(d, var_name, var)
                    endif else begin
                        d = create_struct(d, create_nested_struct(path,var))
                    endelse
                endif
            END
            'TYPE': BEGIN
                tid = h5t_open(id, path)
                var = hdf5_read_attributes(tid)
                h5t_close, tid
                if n_elements(var) ne 0 then begin
                    if keyword_set(flatten) then begin
                        var_names = strsplit(path,'/',/extract)
                        var_name = valid_name(var_names[-1], $
                                   bad_names=used_names, $
                                   post=strcompress(string(i),/remove_all))
                        used_names = [used_names,var_name]
                        d = create_struct(d, var_name, var)
                    endif else begin
                        d = create_struct(d, create_nested_struct(path,var))
                    endelse
                endif
            END
            ELSE:
        ENDCASE
        i = i + 1
    endwhile

    return, d

END

FUNCTION read_hdf5,filename,paths=paths,flatten=flatten,shallow=shallow
    ;+#read_hdf5
    ;+Reads HDF5 file variables and attributes
    ;+***
    ;+##Arguments
    ;+    **filename**: HDF5 file
    ;+
    ;+##Keyword Arguments
    ;+    **paths**: Paths to variables to be read
    ;+
    ;+    **flatten**: Flatten tree structure
    ;+
    ;+    **shallow**: Performs a shallow read i.e. no dataset/group attributes
    ;+
    ;+##Return Value
    ;+Structure containing variables and attributes
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> a = read_hdf5("./test_1a_geometry.h5")
    ;+IDL> b = read_hdf5("./test_1a_geometry.h5",paths="/spec/lens",/flatten,/shallow)
    ;+```
    if file_test(filename) then begin
        ;; Open file
        fid = h5f_open(filename)

        if not keyword_set(paths) then begin
            ;; Read group and sub-groups
            d = hdf5_read_group(fid,shallow=shallow)
        endif else begin
            ;; Read datasets from list
            d = hdf5_read_from_list(fid,paths,flatten=flatten,shallow=shallow)
        endelse

        ;; Close file
        h5f_close,fid
    endif else begin
        print,'File does not exist'
        return, 0
    endelse

    return, d

END
