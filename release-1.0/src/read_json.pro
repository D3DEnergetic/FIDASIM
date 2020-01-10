FUNCTION tokenize,str,regex
	;;; Splits string into tokens according to regex ;;;
    json=str
    group=[]
    start=[]
    finish=[]
    cnt=0
    while stregex(json,regex,/BOOLEAN) do begin
        pos=stregex(json,regex,length=len)
        tmp=strmid(json,pos,len)
        group=[group,tmp]
        if cnt eq 0 then begin
            start=[start,pos]
            finish=[finish,pos+len]
        endif else begin
            start=[start,finish[cnt-1]+pos]
            finish=[finish,finish[cnt-1]+pos+len]
        endelse
        json=strmid(json,pos+len,strlen(json))
        cnt+=1
    endwhile
    return,{group:group,start:start,finish:finish}   
END

FUNCTION json_minify,str
    ;;; Strips out C like comments out of JSON string ;;;
    ;;; Based of JSON.minify https://github.com/getify/JSON.minify ;;;
    regex='"|(#=)|(=#)|(#)|'+string(10b)+'|'+string(13b)
    slashes='(\\)*$'

    in_string = 0
    in_multi = 0 
    in_single = 0
    new_str=''
    index=0
    tokenizer=tokenize(str,regex)
 
    for i=0,n_elements(tokenizer.group)-1 do begin
        
        if not (in_multi or in_single) then begin
             tmp=strmid(str,index,tokenizer.start[i]-index)
             if not in_string then begin
                 tmp=strcompress(tmp,/remove_all)
                 tmp=strjoin(strsplit(tmp,string(10b),/extract,/regex),'')
                 tmp=strjoin(strsplit(tmp,string(13b),/extract,/regex),'')                
             endif
             new_str+=tmp
        endif

        index=tokenizer.finish[i]
        val=tokenizer.group[i]

        if val eq '"' and not (in_multi or in_single) then begin
            escaped = stregex(strmid(str,0,tokenizer.start[i]),slashes,/EXTRACT)
            if not in_string or (escaped eq '' or strlen(escaped) mod 2 eq 0) then in_string = not in_string
            index-=1

        endif else if not (in_string or in_multi or in_single) then begin
            if val eq '#=' then in_multi = 1 else if val eq '#' then in_single=1

        endif else if val eq '=#' and in_multi and not (in_string or in_single) then begin
            in_multi=0
        endif else if (stregex(val,string(10b),/boolean) or stregex(val,string(13b),/boolean)) $
                      and not (in_multi or in_string) and in_single then begin
            in_single=0
        endif else if not ((in_multi or in_single) or $
                      (stregex(val,string(10b),/boolean) or $ 
                       stregex(val,string(13b),/boolean) or $
                       stregex(val,string(9b),/boolean)  or $
                       stregex(val,string(32b),/boolean) )) then begin
            new_str+=val
        endif
    endfor            
    new_str+=strmid(str,index,strlen(str)-index)

    return,new_str
END

FUNCTION read_json, file
    ;+#read_json
    ;+ Reads a JSON file that can have YAML like comments
    ;+***
    ;+##Arguments
    ;+    **file**: JSON file
    ;+
    ;+##Return Value 
    ;+Structure containg JSON values
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> json_struct = read_json("./file.json")
    ;+```

    openr,lun,file,/GET_LUN
    json=''
    tmp=''
    while not EOF(lun) do begin
        readf,lun,tmp
        json+=tmp
        json+=string(10b)
    endwhile
    free_lun,lun
    stripped_json=json_minify(json)

    return,json_parse(stripped_json,/toarray,/tostruct)
END
