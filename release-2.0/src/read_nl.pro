;9/12/2016
;************************************************************************
;+
; NAME:
;	READ_NL
;
; PURPOSE:
;    	Reads Fortran style namelists, returns data as a structure.
;
; CATEGORY:
;	I/O
;
; CALLING SEQUENCE:
;	result = read_nl(filename)
;
; INPUTS:
;	filename - the name of a Fortran namelist file.
;
; KEYWORD PARAMETERS:
;	DBL	 - returns double precision floating point. Default is
;	float.
;       MODEL - if a
;   structure is passed in with the keyword "model" that is the
;   structure used as the starting point for the structure which will
;   be returned by this function.  Data in the namelist file are added to
;   this structure or used to overwrite data already in the
;   structure.   This model can be used to define an array.  Then,
;   entries in the namelist that only modify certain array elements
;   can exist without the array having been previously defined at its
;   full size somewhere else in the namelist.
;        ERROR - returned equal to 0 unless the file cannot be opened
;                in which case 1 is returned.
;
; OUTPUTS:
;	An IDL structure.  The HEADER tagname includes lines (usually 
;	description) in the file before the first namelist.  Each namelist 
;	is a sub-structure.  The namelist name and variable names appear 
;	as the tagnames.  If the number of varilabes of one namelist 
;	exceeds maximum tags (251) allowed in an IDL structure, the namelist
;	will be split into two sub-structures (the second has _ext attached
;	to the tagname).
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;
; EXAMPLE:
;
; HISTORY:
;   02-12-96 B.Rice - created
;   ?? 		- change the line b=strarr(8000) to handle larger files
;    		- each namelist structure is limited to 251 tag names, 
;		  additional names will be put in a separate _ext namelist 
;    		- quotes are left on strings to distinguish from boolean 
;		  variables
;   08-25-98    - IDL 5.1 allows more than 128 tag names, so code has been 
;                 modified to take advantage of this.  Tag_name limit has
;                 been increased to 251
;   10-09-97    - changed reading of strings to preserve white space
;               - also, embedded commas and double quotes '' are preserved  
;   05-31-99 QP - fixed a bug that caused infinite looping when parameters
;	          are separated by ',' instead of ' '.
;   09-12-00  JRF expand the capability to specify values for certain
;   elements of an array.  First, the format y( 1) = 0.0 is allowed;
;   i.e. there can be spaces within the parentheses.  Second, if a
;   structure is passed in with the keyword "model" that is the
;   structure used initially.  Data in the namelist file are added to
;   this structure or used to overwrite data already in the
;   structure.   This model can be used to define an array.  Then,
;   entries in the namelist that only modify certain array elements
;   can exist without the array having been previously defined at its
;   full size somewhere else in the namelist.
;   Finally, the array can be multidimensional.  This is
;   only possible if the array is previously defined multidimensional
;   in the model.  Up to 3 dimensions are allowed.  The format, for
;   example, y(1, 1) = 0.0 can be used.
;
;   Also, a bug was fixed: the array element number should match
;                           Fortran style with element numbers
;                           starting at 1.  Previously, if element "1"
;                           was specified, the second element of the
;                           array was modified.
;
;   Also, an empty namelist is handled.
;   Also, lines beginning with an exclamation point (!)  are
;    ignored.  These are comment lines.
;   08-01-02 QP - F90 on HP generates namelist using (&name,/) instead of
;                 the old ($name,$end), also strings are not in quotes; 
;		  Modified to read both the old and new format.
;   08-23-02 QP - For F90 namelist, get rid of lines after the last /,
;		  the old style did it naturally with splitting around $.
;   08-26-02 QP - convert array[1] to scalar after str_sep(sd,...). In IDL5.5,
;		  struct.array[1] is not treated as scalar as in 5.4 thus
;		  needs be subscripted to be used as a scalar.
;   03-18-03 QP - fixed a bug that removed one char from the last var instead
;		  of the intended last '/'.
;   06-05-03 QP - made it compatiable with Linux F90 namelist where all vars
;	          in one namelist are concatenated into one line.
;   06-09-03 QP - fixed bugs introduced last time
;   10-01-04 QP - fixed a bug that read 1e-8 as 1 (checked '.' only for float)
;                 now check both '.' and 'e'/'E'
;   06-17-05 QP - for tag name with parentheses, if it is in form of name(1),
;		  strip out (1), otherwise, encode using
;		  create_struct_h
;   02-10-06 JRF - In new format namelists which end with a /,
;                 allow for the possibility that there are characters
;                 following the /.  Apparently, some compilers allow
;                 the namelist name to be there.
;                 Comment out the debug printout: 
;                 print,'converting ',origparam,' to ',param
;   02-14-06 JRF -  Handle the case where the file isn't empty but
;                  there is no namelist in the file. Also, fix a bug
;                  that occurred when a namelist entry like "name = 1" 
;                  was encountered when the output structure already had a
;                  tag "name" and that tag was an array.  The namelist 
;                  entry should assign the value to the first element
;                  of the array.  However, because of the way IDL
;                  treats assignment to an array in a structure, all
;                  elements of the array were being altered.
;   08-4-06 JRF - Remove the conversion of the file text to lower
;                  case that was being done at the initial file read.
;                  This conversion was changing strings.  Then, when
;                  determining if a value is boolean, do the
;                  comparison in lower case so that, for instance,
;                  both t and T are allowed.
;   02-9-07 JRF - Merge the versions of this file from the PCS and
;                 hydra. This was done by adding the 08-4-06 change to 
;                 the PCS version.
;                 Added the ERROR keyword to return a flag if the file 
;                 cannot be opened.
;   10-8-07 JRF - 1. Add the comment_character keyword which works as follows:
;                 By default, any lines beginning with an
;                 exclamation point are ignored (this default behavior
;                 is retained for backward compatibility). In
;                 addition, the comment_character keyword can be set
;                 equal to a string with a single character.  This
;                 character is treated as a comment character. The
;                 portions of a line following this
;                 character are ignored. However, if the comment
;                 character is between quotes or apostrophes then it
;                 is ignored.
;                 2. If the keyword print_tokens is set, some
;                 debugging information is printed.  As the code
;                 processes the characters in the namelist file it
;                 splits the characters into strings separated by
;                 equals signs.  This keyword causes the individual
;                 strings to be printed as each one is processed.
;                 When the code chokes on some unexpected new format,
;                 this keyword allows the problem string to be identified.
;
;                 Fix a bunch of bugs:
;                 1.  The "header" was always being returned as a null 
;                 string.  Now, the header will be a string array
;                 containing the lines in the namelist file that
;                 precede the first line that begins a namelist.
;                 2.  The length of each namelist was being limited to 
;                 100000 characters.  There wasn't any fundamental
;                 reason for this and this limitation was removed.
;                 3.  The code that determined whether a namelist
;                 contained more than 250 variables was confused by
;                 multiple references to the same variable in the
;                 namelist file.  This was corrected.
;
;                 Add some comments to the code.
;   4-15-08 JRF - In the function quoted_reverse_find, replaced the
;                 function rstrpos with strpos(...,/reverse_search).
;                 The function rstrpos seemed to be choking on long
;                 strings (e.g. the problem occurred on a string with
;                 more than 203000 characters).  The IDL documentation 
;                 says that rstrpos is obsolete anyway.
;   4-21-09 JRF - Define n_line as a long so that the number of lines
;                 in the file can be larger than the largest positive 
;                 value an integer variable can hold.
;   9-02-09 JRF - Added the keyword line_count_estimate. This value can
;                 be an estimate of the number of lines that are in
;                 the namelist file that needs to be read. The value
;                 is used to initialize an array of strings that will
;                 hold each of the lines in the file. If this keyword
;                 is not specified, the default value is 8000. If the
;                 file contains significantly more than 8000 lines the
;                 process of reading the file goes very slowly. The
;                 file will be read much more quickly if the initial
;                 array of strings is  close to or larger than is
;                 necessary to hold
;                 the file. So, if it is expected that the input
;                 namelist file will have significantly more than 8000
;                 lines, use this keyword.
;   2-21-12 JRF - In a case where the input is a string variable,
;                 allow for a repeat count (e.g. 6*'test').
;   6-16-15 JRF - Fix a couple of issues that were revealed by quirks
;                 in a particular namelist file.
;                 1. Handle the case where the last character in a
;                    line of values is not a comma or a space and the
;                    next line in the file begins in the first column.
;                 2. Handle the case where there are tab characters
;                    following the name of the first namelist. (Now
;                    tab characters will be removed from all lines of
;                    the file before processing (except for the
;                    header), replaced by a space if necessary.
;   9-12-16 JRF - IDL gets very unhappy if an attempt is made to
;                 assign a tag name in a structure to one of IDL's
;                 reserved words (e.g.'ne'). Similarly to what is done
;                 in read_nc.pro, if a tag name needs to be equal to
;                 one of the reserved words, add one an
;                 underscore to the end of the tag in order to avoid
;                 conflict with the reserved word and with any other
;                 tags in the structure.
;                
;
;-
;************************************************************************
;If the input string variable is equal to one of IDL's reserved
;words, replace it with a similar name with an underscore
;at the end.
;
function replace_reserved,name
   new_name = name
   switch strupcase(name) of
      'AND':
      'BEGIN':
      'BREAK':
      'CASE':
      'COMMON':
      'COMPILE_OPT':
      'CONTINUE':
      'DO':
      'ELSE':
      'END':
      'ENDCASE':
      'ENDELSE':
      'ENDFOR':
      'ENDIF':
      'ENDREP':
      'ENDSWITCH':
      'ENDWHILE':
      'EQ':
      'FOR':
      'FORWARD_FUNCTION':
      'FUNCTION':
      'GE':
      'GOTO':
      'GT':
      'IF':
      'INHERITS':
      'LE':
      'LT':
      'MOD':
      'NE':
      'NOT':
      'OF':
      'ON_IOERROR':
      'OR':
      'PRO':
      'REPEAT':
      'SWITCH':
      'THEN':
      'UNTIL':
      'WHILE':
      'XOR': begin
         new_name = name + '_'
         break
      end
   endswitch
   return,new_name
end

PRO nl_remove_comments,input,comment_character,output
;
; From the string "input", remove all characters that follow the first 
; instance of the character specified by comment_character.  However, any
; instances of comment_character between quotes or apostrophes are ignored.

   start =  0
   length =  strlen(input)
   quote_started =  0
   apostrophe_started =  0

   WHILE(1) DO BEGIN
      IF(start EQ length) THEN BEGIN
         output =  input
         return
      ENDIF

      a =  strmid(input,start,1)
      IF( (a EQ comment_character) AND (quote_started EQ 0) AND $
          (apostrophe_started EQ 0)) THEN BEGIN
         output =  strmid(input,0,start)
         return
      ENDIF ELSE IF(a EQ '"') THEN BEGIN
         IF(quote_started EQ 0) THEN BEGIN
            quote_started =  1
         ENDIF ELSE BEGIN
            quote_started =  0
         ENDELSE
      ENDIF ELSE IF(a EQ "'") THEN BEGIN
         IF(apostrophe_started EQ 0) THEN BEGIN
            apostrophe_started =  1
         ENDIF ELSE BEGIN
            apostrophe_started =  0
         ENDELSE
      ENDIF
      start =  start + 1
   ENDWHILE

   return
END

FUNCTION quoted_reverse_find,string,target
; This function calls rstrpos, but, if the string searched for
; is located between a pair of parentheses, it is ignored and the
; next target string is located.  We assume here that there is only
; one pair of parentheses at most.
;
   length =  strlen(string)
   i =  strpos(string,target,/reverse_search)
   k1 =  strpos(string,")",/reverse_search)
   k2 =  strpos(string,"(",/reverse_search)
   IF( (k1 NE -1) AND (k2 NE -1) ) THEN BEGIN
      IF( (i GT k2) AND (i LT k1) ) THEN BEGIN
         substring =  strmid(string,0,k2)
         i =  strpos(substring,target,/reverse_search)
      ENDIF
   ENDIF
   return,i
END


function read_nl,fname,dbl=dbl,debug=debug,model=model,error = error_out,$
                 comment_character= comment_character_in,$
                 print_tokens= print_tokens,$
                 line_count_estimate = line_count_estimate

  dbl=keyword_set(dbl)
  x=''
  error_out = 0

  if (keyword_set(debug)) then error=0 else catch,error
  if (error ne 0) then begin
    message,/info,'Error reading namelist file: '+fname
    message,/info,!ERROR_STATE.msg
    message,/info,!ERROR_STATE.sys_msg
    if (n_elements(lun) gt 0) then free_lun,lun
    return,x
  endif

   k=file_search(fname,count=count)
   if count eq 0 then begin
     if (keyword_set(debug)) then print,'File not found: '+fname
     error_out = 1
     return,x
   endif

   openr,lun,fname,/get_lun 

;
; Read the file and create a string array with one element per line in 
; the file.  The initial array size of nmax is just an initial
; allocation.  If file is longer than this number of lines, the array
; will be extended.
;
   if(n_elements(line_count_estimate) ne 0) then begin
      nmax = line_count_estimate
   endif else begin
      nmax = 8000               ;assume < 8000 lines
   endelse
   b=strarr(nmax)
   n_line = long(0)
   oneline=''
;
; By default, ignore any lines beginning with an exclamation point
; (this default behavior is retained for backward compatibility).
; In addition, if the comment_character 
; keyword is specified, ignore the portions
; of a line following the specified comment character.
;
   while (not(eof(lun))) do begin
     readf,lun,oneline

     IF(n_elements(comment_character_in) NE 0) THEN BEGIN
        FOR i= 0,n_elements(comment_character_in) - 1 DO begin
           nl_remove_comments,oneline,comment_character_in[i],oneline
        ENDFOR
     ENDIF

     IF(strpos(oneline,"!") ne 0)  THEN begin
        if (n_line lt nmax) then $
         b[n_line]=oneline else b = [b, oneline]
        n_line = n_line+1
     ENDIF
   endwhile
   free_lun,lun
   if (n_line eq 0) then message,'empty namelist file'
;   b = strlowcase(b[0:n_line-1])

; Find the first line that starts a namelist.  Any lines prior to this 
; line are returned as the header.

   istart = min(where(strpos(b,'$') ne -1))
   if istart eq -1 then begin	;f90 namelist
      istart = min(where(strpos(b,'&') ne -1))
      if(istart ne -1) then begin
         index = where(strcompress(b,/remove) eq '/',count)
         if count ne 0 then n_line = index[count-1]+1 ;will rm lines after last /
      endif
   endif
   if(istart eq -1) then message,'File does not contain a namelist.'

   IF(istart GT 0) THEN header =  b[0:istart - 1] ELSE header =  ''
;
; Concatenate all of the lines into a single string beginning with the
; first line that starts a namelist. For economy of storage, compress
; blocks of white space except in lines that contain a string
; variable.
;
; Treat the "new line" character at the end of every line as white
; space by adding a space at the end of every line. The compression
; operation will take out this space if it is not needed. However,
; the space will be useful to make certain that lines that begin in the
; first column are separated from the previous line. (Without this, the
; code that handles an array of strings can get into an infinite loop
; if the line previous to the string array doesn't end in a comma and
; the next line begins in the first column.)
;
   s = strcompress(b(istart))
   for i=istart+1,n_line-1 DO BEGIN ;put all lines into one string
      IF (strpos(b[i],"'") EQ -1) THEN s=s+strcompress(b(i) + ' ')  $
      ELSE s = s+b(i) + ' '
   endfor
;
; Separate the single string into an array of strings, one for each
; namelist.  In the process, the dollar sign or ampersand that begins
; the namelist is removed and the $end or / that ends the namelist is
; also removed. Determine the count of namelists.
;
   s=strtrim(s,2)               ;remove leading and trailing blanks
   s=str_sep(s,'$')             ;separate into namelists
   if n_elements(s) eq 1 then begin
      s=str_sep(s,'&') 		;f90 namelist on HP
      newnl = 1
   endif else newnl = 0

   if not newnl then begin 
      nl=(n_elements(s)-1)/2    ;number of lists
      w=indgen(nl)*2+1
      s=s(w)                    ;remove $end from string array
      s=strtrim(s,2)
   endif else begin 		;f90 namelist on HP
      nl=n_elements(s)-1	;number of lists
      s=s[1:nl]			;remove initial string which should be null
      s=strtrim(s,2)
      for i=0,nl-1 do begin	;remove / and any following characters 
                                ;from the end of each namelist 
                                ;string. Apparently some compilers
                                ;accept
                                ;the namelist name after the /.
        temp = strpos(s[i],'/',/reverse_search)
        if(temp ne -1) then s[i] = strmid(s[i],0,temp) 
     endfor
  end
;
; At this point, the variable s is a string array with one element per 
; namelist.  Each string begins with the namelist name. Determine the
; names by finding the portion of each string that precedes the first space.
;
   k=strpos(s,' ')
   lst_nm=strarr(nl) 
   for j=0,nl-1 do BEGIN
      IF(k[j] EQ -1) THEN BEGIN
         lst_nm[j] =  s[j]
         s[j] =  ""
      ENDIF ELSE BEGIN
         lst_nm[j]=strmid(s[j],0,k[j]) ;extract list name
         s[j]=strmid(s[j],k[j]) ;eliminate list name from string
      ENDELSE
   endfor
;
; Loop through each of the namelists.
;
   final_structure_created =  0
   for j=0,nl-1 do begin        ;namelist loop
;
; If a model structure was provided, look for a model for the namelist 
; currently being processed.  If a model exists, use it as the
; structure that will be filled in here.  Otherwise, create a
; placeholder structure and set a flag indicating that a structure
; needs to be initialized when the first tag name is evaluated.
;
      IF(n_elements(model) GT 0) THEN BEGIN
         a =  where(strupcase(tag_names(model)) EQ strupcase(lst_nm[j]),count)
         IF(count EQ 1) THEN BEGIN
            new_structure =  0
            x1 =  model.(a(0))
         ENDIF ELSE BEGIN
            new_structure =  1
            x1 =  {empty_namelist:long(0)} ; default if namelist is empty
         ENDELSE
      ENDIF else BEGIN
         new_structure =  1
         x1 =  {empty_namelist:long(0)} ; default if namelist is empty
      ENDELSE
;
; Divide the single namelist string into an array of strings, breaking 
; the string at each equals sign.
;
      s1=s[j]
      s_sep=str_sep(s1,'=')
      s_sep=strtrim(s_sep,2)
      ns=n_elements(s_sep)-1
      k=0
;
; Each element in the new array should contain a variable name and a
; variable value.  Process each of the elements in the string array.
;
      for i=1,ns do BEGIN       ;loop to extract each variable in nl
         IF(keyword_set(print_tokens)) THEN print,'' ;for debugging
         IF(keyword_set(print_tokens)) THEN print,'***previous token: ',s_sep[i-1] ;for debugging
         IF(keyword_set(print_tokens)) THEN print,'***current token: ',s_sep[i] ;for debugging
         param=strtrim(strmid(s_sep[i-1],k,1000),2) ;get parameter name
         IF (i eq ns) then begin
	    k=quoted_reverse_find(s_sep[i],'/')
	    if k ne -1 then begin
	    if quoted_reverse_find(s_sep[i],' ') gt k then begin ; Linux style
	     ; f90 linux namelist has one namelist on one line, ending / needs
	     ; be trimmed for the last namelistg
	       s_sep[i] = strmid(s_sep[i],0,k)
	       ;print,param+" = '",s_sep[i]+"'"
	    endif
	    endif
	    s_sep[i]=s_sep[i]+' ' ;bug fix, add end space	    
         ENDIF
;    print,j,i,k,'   ',param
       ; 20030605 - previously checked here if vars are separated by ' ' first, then ','.
       ; Now check if vars are separated by ',' first then ' ', which fixes
       ; the problem with Linux f90 nl where the var name would be extracted 
       ; with leading ',' if the previous var is a string with ending spaces.
       ; 20030609 - either order fails for certain files. Now check both and 
       ; take the last one as the separator.
         k1=quoted_reverse_find(s_sep[i],',') ; vars separated by ','
         k2=quoted_reverse_find(s_sep[i],' ') ; vars separated by ' '
         k=max([k1,k2]) ; whichever comes last should be the separator
         if k ne -1 and k eq k1 then k=k+1
         if k eq -1 then k=strlen(s_sep[i]) ;last value in string
         sd=strmid(s_sep[i],0,k) ;extract data
         IF(keyword_set(print_tokens)) THEN print,'***data: ',sd

         kk = 0
         q0 = -1
         q1=strpos(sd,"'")
         q2=strpos(sd,"'",q1+1)
         q2k = q2
         q3 = strpos(sd,"'",q2+1) 
         if (q1 ge 0) and (q2 ge 0) then begin 
;This is a string
;
;            print,'1: ','q1,q2,q3=',q1,q2,q3
;            jcount = 0   ;used to debug an infinite loop problem
;
            REPEAT BEGIN
               WHILE q3 EQ q2k+1 DO BEGIN ;look for double quotes (JRF,huh?)
                  q2k = q3            
                  q3 = strpos(sd,"'",q2k+1)
                  q2 = q3
               ENDWHILE
;               print,'2: ','jcount,q0,q1,q2,q3=',jcount,q0,q1,q2,q3
               if(q1 - q0 - 1 gt 0) then begin
;Allow for input of the form: 6*'test','hello',3*'goodbye'
;
                  extra = strmid(sd,q0 + 1,q1 - q0 - 1)
                  if(strmid(extra,0,1) eq ',') then strput,extra,' ',0
                  extra = strtrim(strcompress(extra,/remove_all),2)
                  temp = strlen(extra)
                  if( (temp gt 0) and $
                      (strmid(extra,temp - 1,1) eq '*') ) then begin
                     temp = long(strmid(extra,0,temp - 1))
                     IF(kk EQ 0) then begin
                        sd1 = replicate(strmid(sd,q1,q2-q1+1),temp)
                     endif else begin
                        sd1 = [sd1,replicate(strmid(sd,q1,q2-q1+1),temp)]
                     endelse
                  endif else begin
                     IF(kk EQ 0) then begin
                        sd1 = strmid(sd,q1,q2-q1+1) 
                     endif else begin
                        sd1 = [sd1,strmid(sd,q1,q2-q1+1)]
                     endelse
                  endelse
               endif else begin
                  IF kk EQ 0 THEN sd1 = strmid(sd,q1,q2-q1+1)  $
                  ELSE sd1 = [sd1,strmid(sd,q1,q2-q1+1)]
               endelse
               q0 = q2
               kk = 1
               q1=strpos(sd,"'",q2+1)
               q2=strpos(sd,"'",q1+1)
               q2k = q2
               q3 = strpos(sd,"'",q2+1) ;look for double quotes
;               print,'3: ','jcount,q0,q1,q2,q3=',jcount,q0,q1,q2,q3
;               jcount = jcount + 1
;               if(jcount gt 50) then stop               
            ENDREP UNTIL q1 EQ -1
            sd = sd1
;            print,'sd: ',sd
;            stop
         endif else BEGIN       ;must be float, integer, or boolean
            repeat begin        ;remove commas if present
               kk=0
               kk=strpos(sd,',',kk+1)
               if kk ne -1 then strput,sd,' ',kk
            endrep until kk eq -1

            sd=str_sep(strtrim(strcompress(sd),2),' ')
            if n_elements(sd) eq 1 then sd=sd[0] ; convert array[1] to scalar

            c1=strmid(sd[0],0,1)+strmid(sd[0],strlen(sd[0])-1,1) ;check 1st/last char
            if c1 ne ".." and strlowcase(c1) ne "ff" and strlowcase(c1) ne "tt" then begin
               if max(strpos(sd,'*')) ne -1 then begin ;split * into array
                  for kk=0,n_elements(sd)-1 do begin  
                     y=str_sep(sd[kk],'*')                      
                     if n_elements(y) eq 1 then begin
                        if dbl then s1=double(y) else s1=float(y)
                        if strpos(y[0],'.') eq -1 then s1=long(s1)
                     endif else begin
                        if dbl then y1=double(y[1]) else y1=float(y[1])
                        if strpos(y[1],'.') eq -1 then y1=long(y1)
                        s1=replicate(y1,y[0])
                     endelse
                     if kk eq 0 then s2=s1 else s2=[s2,s1]
                  endfor
                  sd=s2
               endif else begin
		  on_ioerror,STRING	; handle type conversion error
                  if max(strpos(sd,'.')) eq -1 and $
		     max(strpos(strlowcase(sd),'e')) eq -1 then begin ;check for float
                     sd=long(sd) 	; convert to long
		     goto,NULL
                  endif else begin
                     if dbl then sd=double(sd) else sd=float(sd) ;convert to float
		     goto,NULL
                  endelse
                  STRING:sd="'"+sd+"'" ; convert to string
                  NULL:on_ioerror,NULL
               endelse
            endif ;; will be .T or .F or T or F etc. 
         endelse

	 ; if param is in form of name(1), strip out (1)  QP 17JUN05
	 if (pos = strpos(strcompress(param,/re),'(1)')) ge 0 then begin
	    origparam = param
	    param = strmid(param,0,pos)
;	    print,'converting ',origparam,' to ',param
	 endif
         IF(keyword_set(print_tokens)) THEN print,'***param: ',param

	 skip = 0
         ;;; NEW CODE to handle subscripting into existing array for param="foo(x)"
         if (strpos(param,'(') ge 0) then BEGIN
           subpieces = str_sep(param,'(')
           temp = replace_reserved(strtrim(subpieces[0],2))
           subindex = (where(tag_names(x1) eq $
                             strupcase(temp),subnindex))[0]
           if (subnindex eq 1) then BEGIN
	     skip = 1		; process here, skip in later part
;             print,"new code"
             subpieces =  str_sep(subpieces[1],')')
             subpieces =  str_sep(subpieces[0],',')
             subpieces =  strtrim(subpieces,2)
             subelements = long(subpieces) - 1
             subdata = x1.(subindex)
                                ; Need a test here to make certain
                                ; that the number of dimensions in the 
                                ; array subdata matches the number of
                                ; elements in subelements.
             IF(n_elements(subelements) EQ 1) THEN $
              for subi=0,n_elements(sd)-1 do $
              subdata[subi+subelements[0]] = sd[subi]
             IF(n_elements(subelements) EQ 2) THEN $
              for subi=0,n_elements(sd)-1 do $
              subdata[subi+subelements[0],subelements[1]] = sd[subi]
             IF(n_elements(subelements) EQ 3) THEN $
              for subi=0,n_elements(sd)-1 do $
              subdata[subi+subelements[0],subelements[1],subelements[2]] = $
              sd[subi]
             x1.(subindex) = subdata
           endif ;else message,/info,'Attempt to subscript into non-existent tag: '+subpieces[0]
         endif
         ;;; END NEW CODE

	 if skip eq 0 then begin
           if (new_structure EQ 1) then begin ;make structure
;              print,"old code 1"
             if (strpos(param,'(') ge 0) then begin
                x1=create_struct_h(param,sd)
             endif else begin
                param = replace_reserved(param)
                x1=create_struct(param,sd)
             endelse
             new_structure =  0
           endif else begin
                                ;If param contains an "(" it won't
                                ;match any of the reserved words.
             param = replace_reserved(param)
                                ;If param contains an "(" it won't
                                ;match any of the existing tags.
             w=where(strupcase(param) eq tag_names(x1))

             if w[0] eq -1 then begin 
;                print,"old code 2"
               if (strpos(param,'(') ge 0) then begin
                                ;If param contains an "(" and a
                                ;corresponding tag does not already
                                ;exist (this would have been detected
                                ;above in the "new code"), then create
                                ;the  new structure element with a
                                ;unique hex encoded name because
                                ;there is no way to fill in an array
                                ;element in an array that doesn't
                                ;exist yet.
                  x1=create_struct_h(x1, param,sd)
               endif else begin
                                ;param was already compared to the
                                ;reserved words above.
                  x1=create_struct(x1, param,sd)
               endelse 
             endif else begin   ;overwrite duplicates
;                print,"old code 3"
                                ;We know at this point that there is
                                ;no "(" within param because that
                                ;would have been handled in the "new
                                ;code" above.
               if n_elements(sd) eq 1 then x1.(w[0])[0]=sd[0] else x1.(w[0])=sd
            endelse
           endelse
           
           
           if n_elements(tag_names(x1)) eq 250 then begin ;250 tag name limit, need to start another name
             if final_structure_created eq 0 then begin ;make namelist struct
               x=create_struct('header',header, lst_nm[j],x1)
               final_structure_created =  1
             endif else begin
               x=create_struct(x, lst_nm[j],x1)
             endelse
             lst_nm[j]=lst_nm[j]+'_ext'
           endif
         endif ; skip==0


      endfor                    ;i loop

     if (final_structure_created eq 0) then begin ;make namelist struct
         x=create_struct('header',header, lst_nm[j],x1)
         final_structure_created =  1
      endif else begin
         x=create_struct(x, lst_nm[j],x1)
      endelse

   endfor                       ;j loop

jump1:

;   print,'systime = ',systime(1)-stime
   return,x
end


