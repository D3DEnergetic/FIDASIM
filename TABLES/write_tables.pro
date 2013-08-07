@all_tables.pro
pro write_tables
  low_res=1
  if low_res eq 1 then begin
     print, '!! ATTENTION, tables are calculated with a low resoulution!'
     print, '!! This is fast but not recommended!!'
     print, 'press ENTER to continue'
     dummy=''
     read,dummy
  endif

  einstein_table
  proton_table
  electron_table
  impurity_table,qimp=6 ;; for carbon
  neut_table
end
