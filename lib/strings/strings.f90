
!------------------------------------------------------------
! NOTE: F90 (and most F77) compilers should have the Len_Trim()
! intrinsic that does this. G77 and some others also have LNBLNK().

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of integer function len_trim here]
!integer function Len_Trim(str)
function Len_Trim(str)
  implicit none
  character(len=*) str
  integer Len_Trim
  integer LastNonBlank
  ! purpose is to get "true" length of string; i.e. w/o blanks
  ! INPUT a string str of unknown length with possibly blanks at end
  ! OUTPUT length up to last non-blank
  LastNonBlank = len(str)
  do while (str(lastnonblank:lastnonblank) == ' ' .and. &
      LastNonBlank > 0 )
    LastNonBlank = LastNonBlank - 1
  end do !  while (str(lastnonblank:lastnonblank) == ' ' .and.
!end integer function len_trim 
end function len_trim 

!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine str_next_nonblank here]
subroutine str_next_nonblank(str,ptr)
  implicit none
  character(len=*) str
  integer ptr
  include "strings.fh"
  character(len=2), parameter :: delimiters = " "//char(9)
  do while (ptr <= len(str) .and. &
      index(delimiters,str(ptr:ptr)) > 0 )
      !print *,'ptr nonblank = ',ptr
    ptr = ptr + 1
  end do !  while (ptr <= len(str) .and.
  if (ptr > len(str)) then
    ptr = 0
  end if ! (ptr > len(str))
end subroutine str_next_nonblank

!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine str_next_blank here]
subroutine str_next_blank(str,ptr)
  implicit none
  character(len=*) str
  integer ptr
  include "strings.fh"
  character(len=2), parameter :: delimiters = " "//char(9)
  do while (ptr <= len(str) .and. &
      index(delimiters,str(ptr:ptr)) == 0 )
    ptr = ptr + 1
  end do !  while (ptr <= len(str) .and.
  if (ptr > len(str)) then
    ptr = 0
  end if ! (ptr > len(str))
end subroutine str_next_blank

!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function char_is_alpha here]
function char_is_alpha(chr) !logical
  implicit none
  character chr
  include "strings.fh"
  integer i
  i=IAChar(chr)
  if (  (i >= 65 .and. i <= 90)   &! Upper case ASCII characters
      .or. (i >= 97 .and. i <= 122) ) then
    char_is_alpha = .true.
  else !  (  (i >= 65 .and. i <= 90)  
    char_is_alpha = .false.
  end if !  (  (i >= 65 .and. i <= 90)  
end function char_is_alpha

!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function char_is_digit here]
function char_is_digit(chr) !logical
  implicit none
  character chr
  include "strings.fh"
  integer i
  i=IAChar(chr)
  if (i >= 48 .and. i <= 57) then ! ASCII digits
    char_is_digit = .true.
  else ! (i >= 48 .and. i <= 57)
    char_is_digit = .false.
  end if ! (i >= 48 .and. i <= 57)
end function char_is_digit

!------------------------------------------------------------------------
! Same as the C strspn() function.
! Returns the number of characters in the  initial
! segment of str which are in the string accept.

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function str_spn here]
function str_spn(str,accept)
  implicit none
  include "strings.fh"
  character(len=*) str, accept
  str_spn=0
  do while ( str_cspn < len(str) .and. &
      index(accept,str(str_cspn+1:str_cspn+1)) > 0)
    str_spn = str_spn + 1
  end do !  while ( str_cspn < len(str) .and.
end function str_spn

!------------------------------------------------------------------------
! Same as the C strcspn() function.
! Returns the number of characters in the  initial
! segment of str which are not in the string reject.

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function str_cspn here]
function str_cspn(str,reject)
  implicit none
  include "strings.fh"
  character(len=*) str, reject
  str_cspn=0
  do while ( str_cspn < len(str) .and. &
      index(reject,str(str_cspn+1:str_cspn+1)) == 0)
    str_cspn = str_cspn + 1
  end do !  while ( str_cspn < len(str) .and.
end function str_cspn

!------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine str_toupper here]
subroutine str_toupper(str)
  implicit none
  include "strings.fh"
  character(len=*) str
  integer i

  do i=1,Len_Trim(str)
    str(i:i)=char_toupper(str(i:i))
  end do !  i=1,Len_Trim(str)
end subroutine str_toupper 
!------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function char_toupper here]
function char_toupper(char)
  implicit none
  include "strings.fh"
  character char
  character(len=1) achar !intrinsic
  integer iachar !intrinsic

  if (char >= 'a' .and. char <= 'z') then
    char_toupper=achar(iachar(char)-32)
  else ! (char >= 'a' .and. char <= 'z')
    char_toupper=char
  end if ! (char >= 'a' .and. char <= 'z')
end function char_toupper 
!------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine str_tolower here]
subroutine str_tolower(str)
  implicit none
  include "strings.fh"
  character(len=*) str
  integer i

  do i=1,Len_Trim(str)
    str(i:i)=char_tolower(str(i:i))
  end do !  i=1,Len_Trim(str)
end subroutine str_tolower 
!------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function char_tolower here]
function char_tolower(char)
  implicit none
  include "strings.fh"
  character char
  character(len=1) achar !intrinsic
  integer iachar !intrinsic

  if (char >= 'A' .and. char <= 'Z') then
    char_tolower=achar(iachar(char)+32)
  else ! (char >= 'A' .and. char <= 'Z')
    char_tolower=char
  end if ! (char >= 'A' .and. char <= 'Z')
end function char_tolower 
!------------------------------------------------------
! Some possibly useful C <string.h> functions to add:
! char *strchr (char *s, int c)
! char *strrchr (char *s, int c)
! char *strchrnul (char *s, int c)
! size_t strcspn (char *s, char *reject)
! size_t strspn (char *s, char *accept)
! char *strpbrk (char *s, char *accept)
! char *strstr (char *haystack, char *needle)
! char *strtok (char *restrict s, char *restrict delim)
! char *strcasestr (char *haystack, char *needle)
! char *rindex (char *s, int c)
! int strcasecmp (char *s1, char *s2)
! int strncasecmp (char *s1, char *s2, size_t n)
! char *strsep (char **restrict stringp,)
! char *strsignal (int sig)
! char *basename (char *filename)
