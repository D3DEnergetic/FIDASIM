FUNCTION source_file,name
    ;+ Returns the source file of the routine `name`
    ;+ If name is not given then it returns the source file of caller routine
    if N_PARAMS() eq 0 then begin
        s = scope_traceback(/structure)
        nlevels = n_elements(s)
        sfile = s[nlevels-2].filename
        return, file_expand_path(sfile)
    endif else begin
        help,/source_files,output=csf ;all compiled source files
        nc = n_elements(csf)
        for i=2,nc-1 do begin
            has_name = stregex(csf[i],name,/fold_case) ne -1
            if has_name then begin
                sfile = stregex(csf[i],"\/[\/a-z0-9_\-]*.[a-z0-9_\-]*",/extract,/fold_case)
                return, file_expand_path(sfile)
            endif
        endfor
    endelse

    return,''

END
