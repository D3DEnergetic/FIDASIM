PRO write_ncdf,vars,filename=filename,merge_struct=merge_struct,clobber = clobber

    if not keyword_set(filename) then filename = 'idlsave.cdf'

    if file_test(filename) and not keyword_set(clobber) then begin
        print,'File already exists. Use clobber keyword to overwrite'
        goto, GET_OUT
    endif

    if size(vars,/tname) ne 'STRING' then begin
        print,'Invalid argument type. Expected STRING'
        goto, GET_OUT
    endif

    nVars = n_elements(vars)
    ncid = ncdf_create(filename,/clobber)
    varids = []
    v = {}
    while n_elements(vars) ne 0 do begin
        ;;check to see if variable is set in local scope
        void = execute("isSet = n_elements(" + vars[0] + ")")
        ;;if not set pull variable from above scope
        if isSet eq 0 then begin
            var = scope_varfetch(vars[0],level=-1)
        endif else begin
            void = execute("var="+vars[0])
        endelse

        varInfo = size(var,/structure)
        typeName = varInfo.type_name
        CASE typeName OF
            'STRUCT': begin
                tagNames = strlowcase(tag_names(var))
                if keyword_set(merge_struct) then begin
                    varNames = tagNames
                endif else begin
                    varNames = vars[0]+'_'+tagNames
                endelse
                for i=0,n_elements(tagNames)-1 do begin
                    void = execute(varNames[i]+"= var."+tagNames[i])
                endfor
                if n_elements(vars) gt 1 then vars = vars[1:-1] else vars = []
                vars = [varNames,vars]
                goto,SKIP
            end
            'STRING': begin
                str_id = ncdf_dimdef(ncid,vars[0]+'_dim',n_elements(var))
                strmax_id = ncdf_dimdef(ncid,vars[0]+'_strlen_dim',max(strlen(var)))
                dimids = [strmax_id,str_id]
                void = execute("varids = [varids,ncdf_vardef(ncid,'"+vars[0]+"',dimids,/byte)]")
                var=byte(var)
            end
            ELSE: begin
                if varInfo.n_dimensions gt 0 then begin
                    dimids = []
                    for i=0,varInfo.n_dimensions-1 do begin
                        tmp = ncdf_dimdef(ncid,vars[0]+'_dim'+strcompress(string(i+1),/remove_all),varInfo.dimensions[i])
                        dimids = [dimids,tmp]
                    endfor
                endif else begin
                    dimids = ncdf_dimdef(ncid,vars[0]+'_dim',1)
                endelse
                if typeName eq "INT" then typeName = "LONG"
                void = execute("varids = [varids,ncdf_vardef(ncid,'" + vars[0] +"',dimids,/" +typeName+ ")]")
            end
        ENDCASE

        v = create_struct(v,vars[0],var)
        if n_elements(vars) gt 1 then vars = vars[1:-1] else vars = []
        SKIP:
    endwhile

    ncdf_control,ncid,/ENDEF

    for i=1,N_TAGS(v)-1 do begin
        ncdf_varput,ncid,varids[i-1],v.(i)
    endfor
    
    ncdf_close,ncid
    GET_OUT:
END
