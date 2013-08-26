PRO libdefine
  IF (!version.memory_bits eq 64) then begin  
     DEFSYSV, '!libkk', exist=exist
     if exist eq 0 then $
        defsysv,'!libkk'  ,'/afs/ipp/aug/ads/lib64/@sys/libkk.so'
     DEFSYSV, '!libddww', exist=exist
     if exist eq 0 then  $
        defsysv,'!libddww','/afs/ipp/aug/ads/lib64/@sys/libddww.so' 
  endif else begin
     DEFSYSV, '!LIBKK', exists=exists
     IF exists EQ 0L THEN $
        DEFSYSV, '!LIBKK', '/usr/ads/lib/libkk.so', 1
     DEFSYSV, '!LIBDDWW', exists=exists
     IF exists EQ 0L THEN $
        DEFSYSV, '!LIBDDWW', '/usr/ads/lib/libddww.so', 1 
  endelse
END
