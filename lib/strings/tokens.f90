!----------------------------------------------------------------------
! Return next "token", a white-space delimited string.

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine str_next_token here]
subroutine str_next_token(str,ptr,token,length)
  character(len=*) str    !input
  integer ptr             !in/out
  character(len=*) token  !output
  integer length          !output
  include "strings.fh"

  integer start

  length=0
  ! ptr is zero after all tokens are parsed; return now if already done.
  if (ptr == 0) return

  ! If the last token goes to the end of the string, ptr = len(str)+1
  if (ptr > len(str)) return

  !print *,'ptr token = ',ptr
  call str_next_nonblank(str,ptr)
  if (ptr == 0) return  ! Non-blank character not found.

  start=ptr
  call str_next_blank(str,ptr)
  if (ptr==0) then !No more tokens after this, but accept this one!
    ptr=len(str)+1
  end if
  length=ptr-start
  token=str(start:ptr-1)
end subroutine str_next_token
!----------------------------------------------------------------------
! Return next String; a white-space delimited string that accepts
! quoted strings containing white space.

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine str_next_string here]
subroutine str_next_string(str,ptr,token,length)
  character(len=*) str    !input
  integer ptr             !in/out
  character(len=*) token  !output
  integer length          !output
  include "strings.fh"
  character(len=4) start_chars,end_chars
  integer start,grouping

  data start_chars /'''"(['/
  data end_chars   /'''")]'/

  length=0
  ! ptr is zero after all tokens are parsed; return now if already done.
  if (ptr == 0) return

  call str_next_nonblank(str,ptr)
  if (ptr == 0) return  ! Non-blank character not found.

  grouping=index(start_chars,str(ptr:ptr))
  ! Handle grouped text: quotes, parenthesis, brackets
  if (grouping>0) then
    ptr=ptr+1
    length=index(str(ptr:),end_chars(grouping:grouping))
    token=str(ptr:ptr+length-2)
    ptr=ptr+length
  else ! (grouping>0)
    start=ptr
    call str_next_blank(str,ptr)
    length=ptr-start
    token=str(start:ptr-1)
  end if ! (grouping>0)
  if (ptr > Len_Trim(str)) then !No more tokens after this
    ptr=0
  end if ! (ptr > Len_Trim(str))
end subroutine str_next_string
!----------------------------------------------------------------------------


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine str_next_name here]
subroutine str_next_name(str,ptr,token,length)
  implicit none
  include "strings.fh"
  character(len=*) str, token
  integer ptr, length
  integer str_len, start

  str_len=Len_Trim(str)
  if (ptr >= str_len) then
    ptr=0
    return
  end if ! (ptr >= str_len)

  length=0
  call str_next_nonblank(str,ptr)
  if (ptr == 0 .or. .not. char_is_alpha(str(ptr:ptr)) ) then
    return
  end if ! (ptr == 0 .or. .not. char_is_alpha(str(ptr:ptr)) )
  start=ptr
  ptr=ptr+1
  do while ( ptr <= len(str) &
      .and. ( char_is_alpha(str(ptr:ptr)) &
      .or. char_is_digit(str(ptr:ptr)) &
      .or. str(ptr:ptr) == '$' &
      .or. str(ptr:ptr) == '_' ) )
    ptr=ptr+1
  end do !  while ( ptr <= len(str)
  token=str(start:ptr-1)
  length=ptr-start
end subroutine str_next_name

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ parse through a string returning PARAM=VALUE pairs.
function str_next_param(str,ptr,param,value) !logical
  implicit none
  include "strings.fh"
  include "io.fh"
  character(len=*) str ! Input: string to parse
  integer ptr       ! In/Out: current position in string
  character(len=*) param ! Output: next parameter
  character(len=*) value ! Output: value for the parameter

  integer param_len, value_len
  character(len=128) pair
  integer eq ! position of '='

  str_next_param = .false.
  param=" "; value=" "

  call str_next_name(str,ptr,param,param_len)
  if (param_len==0) return ! No more param name found; all done.
  call str_next_nonblank(str,ptr)
  if (str(ptr:ptr) /= '=') then
    write(stderr,*)'Param/Value Error: no "=" for param: ',param
    return
  end if

  ptr=ptr+1 !Skip over '='
  call str_next_nonblank(str,ptr)
  call str_next_string(str,ptr,value,value_len)
  str_next_param = .true.

end function str_next_param

!------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function str_lookup_token here]
function str_lookup_token(tok,list,ntok)
  implicit none
  include "strings.fh"
  character(len=*) tok,list(*)
  integer ntok

  integer keylen
  character(len=128) key

  if (ntok == 0) then
    str_lookup_token=0
    return
  end if ! (ntok == 0)
  key=tok
  keylen=len(list(1))
  do str_lookup_token=1,ntok
    if (key(1:keylen) == list(str_lookup_token)) then
      return
    end if ! (key(1:keylen) == list(str_lookup_token))
  end do !  str_lookup_token=1,ntok
  str_lookup_token=0
end function str_lookup_token 
!------------------------------------------------------
