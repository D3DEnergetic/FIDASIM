FUNCTION get_version, fidasim_dir
    ;+#get_version
    ;+ Gets FIDASIM version number from git. 
    ;+ Falls back to reading VERSION file when git is not available
    ;+***
    ;+##Input Arguments
    ;+    **fidasim_dir**: FIDASIM install directory
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> version = get_version(getenv("FIDASIM_DIR"))
    ;+```

    version = ''
    git_dir = fidasim_dir+'/.git'
    if strcmp(!VERSION.OS_FAMILY ,'windows', /fold_case) then begin
        spawn,'command -v git ',git_command
    endif else begin
        spawn,'command -v git ',git_command,/sh
    endelse    
    if file_test(git_command) and file_test(git_dir,/dir) then begin
        spawn,git_command+' --git-dir='+git_dir+' describe --tags --always --dirty',version,err_status
    endif else begin
        version_file = fidasim_dir+'/VERSION'
        version = ''
        if file_test(version_file) then begin
            openr, lun, version_file, /get_lun
            readf, lun, version
            free_lun, lun
        endif
    endelse

    return, version

END

