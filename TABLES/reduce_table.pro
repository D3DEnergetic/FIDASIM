pro reduce_table,table,neb,nt,nmax_red,table_red
  ;; programm to reduce the size of the tables to nmax_red!
  sz=size(table)
  nmax=sz[2]
  table_red=dblarr(nmax_red+1,nmax_red,neb,nt)
  ;; reduce the array size to nmax_red:
  index=indgen(nmax_red)
  table_red[index,*,*,*]=table[index,index,*,*]
  ;; write the charge exchange and ionization rates into the
  ;; nmax_red+1 level:
  table_red[nmax_red,*,*,*]=table[nmax,index,*,*]
  ;; assume the excitation into higher states as a loss mechanism
  for n=0,nmax_red-1 do begin 
     for m=nmax_red,nmax-1 do begin
        table_red[nmax_red,n,*,*]+=table[m,n,*,*]
     endfor
  endfor
end
