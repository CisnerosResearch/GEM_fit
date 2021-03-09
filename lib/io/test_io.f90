  implicit none
  integer inf1,inf2,inf3,ios
  character*80 line
  include "io.fh"

  inf1 = File_open("check1","r")
  inf2 = File_open("check2","r")
  write(stdout,*)'inf1,inf2 = ',inf1,inf2
  write(stdout,*)'======================================='
  read(inf1,'(A)',iostat=ios)line
  write(stdout,*)line
  read(inf2,'(A)',iostat=ios)line
  write(stdout,*)line
  call File_close(inf1)
  call File_close(inf2)
  end

subroutine error_handler
end
