; READ_NETCDF_DATA            C.M. Greenfield                      1998
;                             Last modification           July 12, 2001
;
; NAME:
;       READ_NETCDF_DATA
;
; PURPOSE:
;       READ_NETCDF_DATA reads data from a NetCDF file which has
;       already been opened with the NCDF_Open procedure. Note that
;       the structure returned contains pointers to the axis data
;       structures which should be freed with ptr_free when
;       finished. Otherwise, this procedure will be a memory leak.
;
;       If the variable turns out to be a multigraph, a structure
;       containing a list of the included variables is returned.
;
; CATEGORY:
;       TRANSP data handling
;
; CALLING SEQUENCE:
;       read_netcdf_data,cdfid,var[,axrd=axrd][,plist=plist][,/GINFO_ONLY]
;
; INPUTS:
;       CDFID   NetCDF ID returned by NCDF_Open
;       VAR     TRANSP output variable name
;
; KEYWORD PARAMETERS:
;       AXRD    Axis read flag used during recursive calls.
;               AXRD SHOULD NOT BE SET BY THE CALLER!
;               Eventually, this will be removed in favor of a 
;               call to HELP,/CALLS
;       PLIST   Data for return of axis data during recursive calls.
;               PLIST SHOULD NOT BE SET BY THE CALLER!
;       GINFO_ONLY      Flag indicating that read_netcdf_data should
;                       just read the globals and return.
; OTHER REQUIRED ROUTINES:
;       NONE
;
;-----------------------------------------------------------------------------
pro cleanup_pointers
common READ_NETCFD_DATA_PTR, pointerlist
for i=0,n_elements(pointerlist)-1 do Ptr_free,pointerlist[i]
end

function read_netcdf_data,cdfid,var,axrd=axrd,plist=plist,globals=globals, $
                          ginfo_only=ginfo_only
common READ_NETCFD_DATA_PTR, pointerlist
forward_function read_netcdf_data
name=strupcase(var)
;Get global stuff
info=ncdf_inquire(cdfid)
for i=0,info.ngatts-1 do begin
    gname=ncdf_attname(cdfid,i,/global)
    ncdf_attget,cdfid,gname,/global,gvalue
    if (size(gvalue,/type) eq 1) then gvalue=strtrim(string(gvalue),2)
    if (n_elements(globals) eq 0) then globals=create_struct(gname,gvalue) $
    else globals=create_struct(globals,gname,gvalue)
endfor
for i=0,info.ndims-1 do begin
    ncdf_diminq,cdfid,i,dname,dsize
    if (n_elements(dims) eq 0) then dims=create_struct(dname,dsize) $
    else dims=create_struct(dims,dname,dsize)
endfor
ginfo={globals:globals,dims:dims}
if keyword_set(ginfo_only) then return,ginfo
;Get the variable id, return if not found
varid=NCDF_Varid(cdfid,name) & if (varid lt 0) then return,{error:-1}
varinfo=NCDF_Varinq(cdfid,varid)
;Get the long name and units of the variable
NCDF_Attget,cdfid,varid,'units',units & units=strtrim(string(units),2)
NCDF_Attget,cdfid,varid,'long_name',lname & lname=strtrim(string(lname),2)
title=lname+' ('+units+')'
;Get the data
NCDF_Varget,cdfid,varid,data
;Check for multigraph definition
if (varinfo.ndims eq 0) then begin
    message,strupcase(var)+' is a multigraph definition',/information
    NCDF_Attget,cdfid,varid,'Fct_Ids',Fct_Ids
    n=n_elements(Fct_Ids)
    case data of
        1:offset=1+ginfo.globals.nfxt ;Function of t only
        0:offset=1 ;Function of x and t
        else:BEGIN
            help,data
            message,'Not understood in multigraph handling.',/continue
            return,{error:-99}
        END
    endcase
    multi=strarr(n)
    for i=0,n-1 do multi[i]=(NCDF_Varinq(cdfid,abs(fct_ids[i])+offset)).name
    return,{name:name,ndims:0, $
            lname:lname,units:units,title:title, $
            multi:multi,sign:fct_ids/abs(fct_ids),ginfo:ginfo,error:-100}
endif
s=size(data) & if (s[0] eq 0) then dim=1 else dim=s(1:s[0])
;read_netcdf_data calls itself recursively here to get the axis
;data. Then it builds a results structure and returns to the caller.
;axrd (axis read flag) is only set when read_netcdf_data is called to
;get axis data. If set, read_netcdf_data will not call itself again -
;this prevents an infinite loop
if (n_elements(axrd) eq 0) then begin 
    axdat=ptrarr(varinfo.ndims)
    for i=0,varinfo.ndims-1 do begin
        NCDF_Diminq,cdfid,varinfo.dim(i),aname,asize
        axis_data=read_netcdf_data(cdfid,aname,/axrd,plist=plist)
        axdat(i)=Ptr_new(axis_data(0))
        if (n_elements(plist) eq 0) then plist=[axdat(i)] else  $
          plist=[plist,axdat(i)]
    endfor        
    if (n_elements(pointerlist) eq 0) then pointerlist=axdat $
    else pointerlist=[pointerlist,axdat]
    return,{name:name,ndims:varinfo.ndims,dim:dim, $
            lname:lname,units:units,title:title, $
            data:data,axis:axdat,ginfo:ginfo,error:0}
endif else return,{name:name,ndims:varinfo.ndims,dim:dim, $
                   lname:lname,units:units,title:title, $
                   data:data,ginfo:ginfo,error:0}
end

