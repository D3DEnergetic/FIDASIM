PRO check_struct_schema, schema, s, err_status, desc=desc
    ;+#check_struct_schema
    ;+ Check structure `s` is formatted according to `schema`
    ;+***
    ;+##Input Arguments
    ;+     **schema**: structure schema
    ;+
    ;+     **s**: structure to check
    ;+
    ;+##Output Arguments
    ;+     **err**: error code
    ;+
    ;+##Keyword Arguments
    ;+     **desc**: description of structure `s`
    ;+
    ;+##Example usage
    ;+```idl
    ;+IDL> s = {a:0, b:[1.d0,2.d0], c:"example"}
    ;+IDL> schema = {a:{dims:0,type:"INT"}, b:{dims:[2],type:"DOUBLE"}, c:{dims:0,type:"STRING"}  }
    ;+
    ;+IDL> check_struct_schema, schema, s, err, desc="Example structure"
    ;+IDL> print, err
    ;+    0
    ;+```
    if not keyword_set(struct_name) then desc = 'structure'

    err_status = 0
    schema_tags = strlowcase(TAG_NAMES(schema)) 
    tags = strlowcase(TAG_NAMES(s)) 

    for i=0,n_elements(tags)-1 do begin
        w=where(tags[i] eq schema_tags,nw)
        if nw eq 0 then begin
            info,'Extra variable "'+tags[i]+'" found in '+desc
        endif
    endfor

    for i=0,n_elements(schema_tags)-1 do begin
        w=where(schema_tags[i] eq tags,nw)
        if nw eq 0 then begin
            error,'"'+schema_tags[i]+'" is missing from the '+desc
            err_status = 1
        endif else begin
            ;; Check dimensions
            ww = where((size(s.(w),/dim) eq schema.(i).dims) ne 1,nww) 
            if nww ne 0 then begin
                error,'"'+schema_tags[i]+'" has the wrong dimensions. Expected ('+ $
                      strjoin(strcompress(string(schema.(i).dims),/remove_all),',')+')'
                print,'size('+schema_tags[i]+') = ',size(s.(w),/dim)
                err_status = 1
            endif
            ;; Check type
            tname = size(s.(w),/tname)
            if tname ne schema.(i).type then begin
                error,'"'+schema_tags[i]+'" has the wrong type. Expected '+schema.(i).type
                print,'type('+schema_tags[i]+') = '+tname
                err_status = 1
            endif
            ;; Check for NaNs or Inf
            if tname ne 'STRING' and tname ne 'STRUCT' then begin
                ww = where(finite(s.(w)) eq 0,nww) 
            endif else nww = 0
            if nww ne 0 then begin
                error,'NaN or Infinity detected in "'+schema_tags[i]+'"'
                err_status = 1
            endif
        endelse
    endfor

END

