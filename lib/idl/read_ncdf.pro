FUNCTION read_ncdf,file,vars=vars
    ;+#read_ncdf
    ;+Reads a flat NetCDF file
    ;+***
    ;+##Arguments
    ;+    **file**: NetCDF file
    ;+
    ;+##Keyword Arguments
    ;+    **vars**: List of variables to read
    ;+
    ;+##Return Value
    ;+Structure containing NetCDF variables
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> a = read_ncdf("./123324H01_fi_1.cdf")
    ;+```

    ;;List of reserved names
    reserved=['AND','BEGIN','BREAK','CASE','COMMON','COMPILE_OPT',$
              'CONTINUE','DO','ELSE','END','ENDCASE','ENDELSE',$
              'ENDFOR','ENDFOREACH','ENDIF','ENDREP','ENDSWITCH','ENDWHILE',$
              'EQ','FOR','FOREACH','FORWARD_FUNCTION','FUNCTION','GE',$
              'GOTO','GT','IF','INHERITS','LE','LT','MOD','NE','NOT','OF',$
              'ON_IOERROR','OR','PRO','REPEAT','SWITCH','THEN','UNTIL',$
              'WHILE','X0R']
    d={err:1}
    if file_test(file) then begin
        d.err=0
        ncid=ncdf_open(file,/nowrite)
        info=ncdf_inquire(ncid)
        if keyword_set(vars) then nvars=n_elements(vars) else nvars=info.nvars
        for i=0,nvars-1 do begin
            if keyword_set(vars) then begin
                name=vars[i]
                CATCH, err_status
                if err_status ne 0 then begin
                    CATCH,/CANCEL
                    continue
                endif
                ncdf_varget,ncid,name,tmp
            endif else begin
                var_info=ncdf_varinq(ncid,i)
                name=strjoin(strsplit(var_info.name,".",/extract),"_")
                ncdf_varget,ncid,i,tmp
            endelse
            if total(strmatch(reserved,name,/FOLD_CASE)) gt 0 then name+='_'
            d=create_struct(d,name,tmp)
        endfor
	ncdf_close,ncid
    endif else message,'FILE DOES NOT EXIST: '+file
    return,d
END

