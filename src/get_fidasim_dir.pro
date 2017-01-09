FUNCTION get_fidasim_dir
    ;+#get_fidasim_dir
    ;+ Gets FIDASIM install directory
    ;+***
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> fida_dir = get_fidasim_dir()
    ;+``` 
    return, file_dirname(file_dirname(file_dirname(source_file())))
END
