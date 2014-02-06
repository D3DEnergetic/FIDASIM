FUNCTION read_ncdf,file,vars=vars
    d={err:1}
    if file_test(file) then begin
        d.err=0
        ncid=ncdf_open(file,/nowrite)
        info=ncdf_inquire(ncid)
        if keyword_set(vars) then nvars=n_elements(vars) else nvars=info.nvars
        for i=0,nvars-1 do begin
            if keyword_set(vars) then begin
                name=vars[i]
                ncdf_varget,ncid,name,tmp
            endif else begin
                var_info=ncdf_varinq(ncid,i)
                name=var_info.name
                ncdf_varget,ncid,i,tmp
            endelse
            d=create_struct(d,name,tmp)
        endfor
	ncdf_close,ncid
    endif else print,'FILE DOES NOT EXIST: '+file
    return,d
END

