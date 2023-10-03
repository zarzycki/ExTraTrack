undef("list_allocated")
procedure list_allocated(
  debug      : logical
  )

begin

if (debug) then
  list_vars()
  list_hlus()
end if

end

;-----------------------------------------------

undef("print_memory_info")
procedure print_memory_info()

begin

  print(systemfunc("ps u"))

end