PRO sav_to_hdf5,savefile,filename=filename
    ;+#save_to_hdf5
    ;+Converts an IDL save file to a HDF5 file
    ;+***
    ;+##Arguments
    ;+    **savefile**: IDL save file
    ;+
    ;+##Keyword Arguments
    ;+   **filename**: HDF5 file name defaults to IDL save file name
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> sav_to_hdf5,'./dne159245.00500'
    ;+```
    if not keyword_set(filename) then begin
        p = strsplit(savefile,'.',/extract)
        if strlowcase(p[-1]) eq 'sav' then c = 2 else c = 1
        filename = strjoin(p[0:n_elements(p)-c],'.')+'.h5'
    endif

    restore, savefile

    sObj = OBJ_NEW('IDL_Savefile',savefile)
    vars = strlowcase(sObj->Names())
    write_hdf5,vars,filename=filename

END
