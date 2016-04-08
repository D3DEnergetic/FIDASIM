FUNCTION get_version, fidasim_dir
    ;+##`get_version(fidasim_dir)`
    ;+ Gets FIDASIM version number from git. 
    ;+ Falls back to reading VERSION file when git is not available
    ;+
    ;+###Input Arguments
    ;+    **fidasim_dir**: FIDASIM install directory
    ;+
    ;+###Example Usage
    ;+```idl
    ;+IDL> version = get_version("~/FIDASIM")
    ;+```

    version = ''
    git_dir = fidasim_dir+'/.git'
    spawn,'command -v git ',git_command,/sh
    if file_test(git_command) and file_test(git_dir,/dir) then begin
        spawn,git_command+' --git-dir='+git_dir+' describe --tags --always',version,err_status
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

