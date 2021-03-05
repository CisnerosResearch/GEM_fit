!---------------------------------------------------------------
integer function TRDPRM_get_next_line(inunit,line)
  implicit none
  integer inunit
  character(len=*) line
  integer j,k,length,ios

  read(inunit,'(A)',iostat=ios)line
  TRDPRM_get_next_line = ios
  if ( ios /= 0 )return
! remove comments (after exclamation points)
! using this in place of a direct read is important whenever num_fields
! remaining is to be used
  length = len(line)
  do j = 1,length
    if ( line(j:j) == '!')then
      do k = j,length
        line(k:k) = ' '
      enddo
      return
    endif
  enddo
  return
end function TRDPRM_get_next_line
!------------------------------------------------------
