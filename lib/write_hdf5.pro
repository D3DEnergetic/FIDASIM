FUNCTION vars_to_struct,vars=vars,level=level
    ;; struct_from_list
    ;; creates structure from list of variables in specified scope
    ;; if vars not set then use all valid variables in specified scope
    if not keyword_set(level) then level=-1
    if not keyword_set(vars) then vars = scope_varname(level=level)

    nvars = n_elements(vars)

    for i=0,nvars-1 do begin
        catch, err_status
        if err_status ne 0 then begin
            catch,/cancel
            continue
        endif
        var = scope_varfetch(vars[i],level=level)
        if n_elements(var) ne 0 then begin
           if i eq 0 then begin
               s = create_struct(vars[i],var)
           endif else begin
               s = create_struct(s,vars[i],var)
           endelse
        endif
    endfor
    
    if nvars eq 1 and size(s.(0),/tname) eq 'STRUCT' then begin
        s = s.(0)
    endif
    return, s
END

PRO hdf5_write_struct, id, struct

    ntags = n_tags(struct)
    if ntags eq 0 then goto, GET_OUT

    tags = strlowcase(TAG_NAMES(struct))
    for i=0,ntags-1 do begin
        var = struct.(i)
        varInfo = size(var,/structure)
        typeName = varInfo.type_name

        if typeName eq 'STRUCT' then begin
            gid = h5g_create(id,tags[i])
            hdf5_write_struct,gid,var
            h5g_close, gid
        endif else begin
            data = var
            ndim = size(data,/n_dim)
            dims = size(data,/dim)

            datatype_id = h5t_idl_create(data)
            if ndim eq 0L then begin
                dataspace_id = h5s_create_scalar()
                dataset_id = h5d_create(id, tags[i], datatype_id, dataspace_id)
            endif else begin
                dataspace_id = h5s_create_simple(dims)
                dataset_id = h5d_create(id, tags[i], datatype_id, dataspace_id, $
                                        chunk_dimensions=dims, gzip=9,/shuffle)
            endelse

            h5d_write, dataset_id, data

            h5d_close, dataset_id
            h5s_close, dataspace_id
            h5t_close, datatype_id
        endelse
    endfor
    GET_OUT:
END

FUNCTION valid_attribute, att
    tags = ["obj","name","data"]
    att_tags = strlowcase(TAG_NAMES(att))

    is_valid = 1
    for i=0,n_elements(tags)-1 do begin
        w = where(tags[i] eq att_tags,nw)
        if nw eq 0 then begin
            print, 'ERROR: Structure tag "'+tags[i]+'" missing from attribute definition'
            is_valid = 0
        endif
    endfor
    if is_valid eq 0 then begin
        help, att
        goto, GET_OUT
    endif
    value_info = size(att.data,/structure)
    if value_info.type_name eq 'STRUCT' then begin
        print,'ERROR: attribute value cannot be a structure'
        is_valid = 0
    endif

    GET_OUT:    
    return, is_valid
END

PRO hdf5_write_att_data, id, name, data

    data_info = size(data,/structure)
    type_name = data_info.type_name

    if type_name eq 'STRING' then begin
        value = strjoin(data,", ",/single)
    endif else begin
        value = data
    endelse

    dims = size(value,/dim)
    ndims = size(value,/n_dim)
    datatype_id = h5t_idl_create(value)

    if ndims eq 0L then begin
        dataspace_id = h5s_create_scalar()
    endif else begin
        dataspace_id = h5s_create_simple(dims)
    endelse

    att_id = h5a_create(id, name, datatype_id, dataspace_id)
    h5a_write, att_id, value
    h5a_close, att_id
END

PRO hdf5_write_attributes,id,atts
    
    natts = n_elements(atts)

    for i=0, natts-1 do begin
        if not valid_attribute(atts[i]) then continue
        
        object_info = h5g_get_objinfo(id,atts[i].obj)
    
        CASE object_info.type OF
            'LINK': print,'ERROR: Can not handle an attribute of a reference'
            'GROUP': BEGIN
                gid = h5g_open(id, atts[i].obj)
                hdf5_write_att_data, gid, atts[i].name, atts[i].data
                h5g_close,gid
            END
            'DATASET': BEGIN
                did = h5d_open(id,atts[i].obj)
                hdf5_write_att_data, did, atts[i].name, atts[i].data
                h5d_close,did
            END
            'TYPE': BEGIN
                tid = h5t_open(id,atts[i].obj)
                hdf5_write_att_data, tid, atts[i].name, atts[i].data
                h5t_close,tid
            END
            ELSE: print,'ERROR: Unknown object'
        ENDCASE
    endfor

END

PRO write_hdf5,vars,atts=atts,filename=filename,clobber=clobber
    ;+#write_hdf5
    ;+Writes HDF5 files from variables in the local scope or a structure
    ;+***
    ;+##Arguments
    ;+    **vars**: List of variables or a structure
    ;+
    ;+##Keyword Arguments
    ;+    **atts**: Attributes to write
    ;+
    ;+    **filename**: Filename of output HDF5 file
    ;+
    ;+    **clobber**: Overwrite exisiting HDF5 file
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> a = [1,2,3]
    ;+IDL> write_hdf5,"a",atts={attribute,obj:"/a",name:"description",data:"example variable"}
    ;+IDL> s = {a:1,b:{a:2}}
    ;+IDL> write_hdf5,s,/clobber
    ;+```

    if not keyword_set(filename) then filename = 'idlsave.h5'

    if file_test(filename) and not keyword_set(clobber) then begin
        print,"File already exists. Use clobber keyword to overwrite"
        goto, GET_OUT
    endif

    nvars = n_elements(vars)
    if nvars eq 0 then goto, GET_OUT

    type = size(vars,/tname)
    if type ne 'STRUCT' and type ne 'STRING' then begin
       print, "Invalid argument type. Expected STRING or STRUCT"
       print, type
       goto, GET_OUT
    endif

    if nvars ne 1 and type eq 'STRUCT' then begin
        print, "Invalid argument type. Arrays of structs not permitted"
        goto, GET_OUT
    endif
 
    if nvars eq 1 and type eq 'STRUCT' then begin
        var_struct = vars
    endif else begin
        var_struct = vars_to_struct(vars=vars,level=-2)
    end

    file_id = h5f_create(filename)
    hdf5_write_struct, file_id, var_struct

    if keyword_set(atts) then begin
        hdf5_write_attributes, file_id, atts
    endif

    h5f_close, file_id

    GET_OUT:
END
