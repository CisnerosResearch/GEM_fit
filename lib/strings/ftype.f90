! Fortan version of standard ctype functions (no locale support)
!
! Character class functions are all
!   logical function ch_<test>(character(len=1))
! where test is one of:  isalnum, isalpha, iscntrl, isdigit,
! islower, isgraph, isprint, ispunct, isspace, isupper, isxdigit
!
! Also case conversion functions:
!   character(len=1) function ch_tolower(character(len=1))
!   character(len=1) function ch_toupper(character(len=1))
! Return upper or lower case for given any-case character.
!
! Also, string conversion utilities:
! subroutine str_toupper(character(len=*) in, character(len=*) out)
! subroutine str_tolower(character(len=*) in, character(len=*) out)


block data ftype_blockdata
  implicit none
!#define IN_BLOCKDATA
 include "ftype.inc"
!#undef IN_BLOCKDATA

  data ctype_class_table / &
    z'0200', z'0200', z'0200', z'0200', z'0200', z'0200', z'0200', z'0200', &
    z'0200', z'0314', z'0214', z'0214', z'0214', z'0214', z'0200', z'0200', &
    z'0200', z'0200', z'0200', z'0200', z'0200', z'0200', z'0200', z'0200', &
    z'0200', z'0200', z'0200', z'0200', z'0200', z'0200', z'0200', z'0200', &
    z'013C', z'0400', z'0400', z'0400', z'0400', z'0400', z'0400', z'0400', &
    z'0400', z'0400', z'0400', z'0400', z'0400', z'0400', z'0400', z'0400', &
    z'0800', z'0800', z'0800', z'0800', z'0800', z'0800', z'0800', z'0800', &
    z'0800', z'0800', z'0400', z'0400', z'0400', z'0400', z'0400', z'0400', &
    z'0400', z'0800', z'0800', z'0800', z'0800', z'0800', z'0800', z'0800', &
    z'0800', z'0800', z'0800', z'0800', z'0800', z'0800', z'0800', z'0800', &
    z'0800', z'0800', z'0800', z'0800', z'0800', z'0800', z'0800', z'0800', &
    z'0800', z'0800', z'0800', z'0400', z'0400', z'0400', z'0400', z'0400', &
    z'0400', z'0800', z'0800', z'0800', z'0800', z'0800', z'0800', z'0800', &
    z'0800', z'0800', z'0800', z'0800', z'0800', z'0800', z'0800', z'0800', &
    z'0800', z'0800', z'0800', z'0800', z'0800', z'0800', z'0800', z'0800', &
    z'0800', z'0800', z'0800', z'0400', z'0400', z'0400', z'0400', z'0200', &
    z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', &
    z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', &
    z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', &
    z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', &
    z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', &
    z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', &
    z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', &
    z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', &
    z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', &
    z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', &
    z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', &
    z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', &
    z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', &
    z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', &
    z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', &
    z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000', z'0000'  /
    
    
  data ch_toupper_table / &
    z'00', z'01', z'02', z'03', z'04', z'05', z'06', z'07', &
    z'08', z'09', z'0a', z'0b', z'0c', z'0d', z'0e', z'0f', &
    z'10', z'11', z'12', z'13', z'14', z'15', z'16', z'17', &
    z'18', z'19', z'1a', z'1b', z'1c', z'1d', z'1e', z'1f', &
    z'20', z'21', z'22', z'23', z'24', z'25', z'26', z'27', &
    z'28', z'29', z'2a', z'2b', z'2c', z'2d', z'2e', z'2f', &
    z'30', z'31', z'32', z'33', z'34', z'35', z'36', z'37', &
    z'38', z'39', z'3a', z'3b', z'3c', z'3d', z'3e', z'3f', &
    z'40', z'41', z'42', z'43', z'44', z'45', z'46', z'47', &
    z'48', z'49', z'4a', z'4b', z'4c', z'4d', z'4e', z'4f', &
    z'50', z'51', z'52', z'53', z'54', z'55', z'56', z'57', &
    z'58', z'59', z'5a', z'5b', z'5c', z'5d', z'5e', z'5f', &
    z'60', z'41', z'42', z'43', z'44', z'45', z'46', z'47', &
    z'48', z'49', z'4a', z'4b', z'4c', z'4d', z'4e', z'4f', &
    z'50', z'51', z'52', z'53', z'54', z'55', z'56', z'57', &
    z'58', z'59', z'5a', z'7b', z'7c', z'7d', z'7e', z'7f', &
    z'80', z'81', z'82', z'83', z'84', z'85', z'86', z'87', &
    z'88', z'89', z'8a', z'8b', z'8c', z'8d', z'8e', z'8f', &
    z'90', z'91', z'92', z'93', z'94', z'95', z'96', z'97', &
    z'98', z'99', z'9a', z'9b', z'9c', z'9d', z'9e', z'9f', &
    z'a0', z'a1', z'a2', z'a3', z'a4', z'a5', z'a6', z'a7', &
    z'a8', z'a9', z'aa', z'ab', z'ac', z'ad', z'ae', z'af', &
    z'b0', z'b1', z'b2', z'b3', z'b4', z'b5', z'b6', z'b7', &
    z'b8', z'b9', z'ba', z'bb', z'bc', z'bd', z'be', z'bf', &
    z'c0', z'c1', z'c2', z'c3', z'c4', z'c5', z'c6', z'c7', &
    z'c8', z'c9', z'ca', z'cb', z'cc', z'cd', z'ce', z'cf', &
    z'd0', z'd1', z'd2', z'd3', z'd4', z'd5', z'd6', z'd7', &
    z'd8', z'd9', z'da', z'db', z'dc', z'dd', z'de', z'df', &
    z'e0', z'e1', z'e2', z'e3', z'e4', z'e5', z'e6', z'e7', &
    z'e8', z'e9', z'ea', z'eb', z'ec', z'ed', z'ee', z'ef', &
    z'f0', z'f1', z'f2', z'f3', z'f4', z'f5', z'f6', z'f7', &
    z'f8', z'f9', z'fa', z'fb', z'fc', z'fd', z'fe', z'ff'  /
    
  data ch_tolower_table / &
    z'00', z'01', z'02', z'03', z'04', z'05', z'06', z'07', &
    z'08', z'09', z'0a', z'0b', z'0c', z'0d', z'0e', z'0f', &
    z'10', z'11', z'12', z'13', z'14', z'15', z'16', z'17', &
    z'18', z'19', z'1a', z'1b', z'1c', z'1d', z'1e', z'1f', &
    z'20', z'21', z'22', z'23', z'24', z'25', z'26', z'27', &
    z'28', z'29', z'2a', z'2b', z'2c', z'2d', z'2e', z'2f', &
    z'30', z'31', z'32', z'33', z'34', z'35', z'36', z'37', &
    z'38', z'39', z'3a', z'3b', z'3c', z'3d', z'3e', z'3f', &
    z'40', z'61', z'62', z'63', z'64', z'65', z'66', z'67', &
    z'68', z'69', z'6a', z'6b', z'6c', z'6d', z'6e', z'6f', &
    z'70', z'71', z'72', z'73', z'74', z'75', z'76', z'77', &
    z'78', z'79', z'7a', z'5b', z'5c', z'5d', z'5e', z'5f', &
    z'60', z'61', z'62', z'63', z'64', z'65', z'66', z'67', &
    z'68', z'69', z'6a', z'6b', z'6c', z'6d', z'6e', z'6f', &
    z'70', z'71', z'72', z'73', z'74', z'75', z'76', z'77', &
    z'78', z'79', z'7a', z'7b', z'7c', z'7d', z'7e', z'7f', &
    z'80', z'81', z'82', z'83', z'84', z'85', z'86', z'87', &
    z'88', z'89', z'8a', z'8b', z'8c', z'8d', z'8e', z'8f', &
    z'90', z'91', z'92', z'93', z'94', z'95', z'96', z'97', &
    z'98', z'99', z'9a', z'9b', z'9c', z'9d', z'9e', z'9f', &
    z'a0', z'a1', z'a2', z'a3', z'a4', z'a5', z'a6', z'a7', &
    z'a8', z'a9', z'aa', z'ab', z'ac', z'ad', z'ae', z'af', &
    z'b0', z'b1', z'b2', z'b3', z'b4', z'b5', z'b6', z'b7', &
    z'b8', z'b9', z'ba', z'bb', z'bc', z'bd', z'be', z'bf', &
    z'c0', z'c1', z'c2', z'c3', z'c4', z'c5', z'c6', z'c7', &
    z'c8', z'c9', z'ca', z'cb', z'cc', z'cd', z'ce', z'cf', &
    z'd0', z'd1', z'd2', z'd3', z'd4', z'd5', z'd6', z'd7', &
    z'd8', z'd9', z'da', z'db', z'dc', z'dd', z'de', z'df', &
    z'e0', z'e1', z'e2', z'e3', z'e4', z'e5', z'e6', z'e7', &
    z'e8', z'e9', z'ea', z'eb', z'ec', z'ed', z'ee', z'ef', &
    z'f0', z'f1', z'f2', z'f3', z'f4', z'f5', z'f6', z'f7', &
    z'f8', z'f9', z'fa', z'fb', z'fc', z'fd', z'fe', z'ff'  /

end block data ftype_blockdata

function ch_isalnum(ch)
  implicit none
  include "ftype.inc"
  ch_isalnum = (IAnd(ctype_class_table(IAchar(ch)),ISalnum_bit)>0)
end function ch_isalnum

function ch_isalpha(ch)
  implicit none
  include "ftype.inc"
  ch_isalpha = (IAnd(ctype_class_table(IAchar(ch)),ISalpha_bit)>0)
end function ch_isalpha

function ch_iscntrl(ch)
  implicit none
  include "ftype.inc"
  ch_iscntrl = (IAnd(ctype_class_table(IAchar(ch)),IScntrl_bit)>0)
end function ch_iscntrl

function ch_isdigit(ch)
  implicit none
  include "ftype.inc"
  ch_isdigit = (IAnd(ctype_class_table(IAchar(ch)),ISdigit_bit)>0)
end function ch_isdigit

function ch_islower(ch)
  implicit none
  include "ftype.inc"
  ch_islower = (IAnd(ctype_class_table(IAchar(ch)),ISlower_bit)>0)
end function ch_islower

function ch_isgraph(ch)
  implicit none
  include "ftype.inc"
  ch_isgraph = (IAnd(ctype_class_table(IAchar(ch)),ISgraph_bit)>0)
end function ch_isgraph

function ch_isprint(ch)
  implicit none
  include "ftype.inc"
  ch_isprint = (IAnd(ctype_class_table(IAchar(ch)),ISprint_bit)>0)
end function ch_isprint

function ch_ispunct(ch)
  implicit none
  include "ftype.inc"
  ch_ispunct = (IAnd(ctype_class_table(IAchar(ch)),ISpunct_bit)>0)
end function ch_ispunct

function ch_isspace(ch)
  implicit none
  include "ftype.inc"
  ch_isspace = (IAnd(ctype_class_table(IAchar(ch)),ISspace_bit)>0)
end function ch_isspace

function ch_isupper(ch)
  implicit none
  include "ftype.inc"
  ch_isupper = (IAnd(ctype_class_table(IAchar(ch)),ISupper_bit)>0)
end function ch_isupper

function ch_isxdigit(ch)
  implicit none
  include "ftype.inc"
  ch_isxdigit = (IAnd(ctype_class_table(IAchar(ch)),ISxdigit_bit)>0)
end function ch_isxdigit

!----------------------------------------------------------------
function ch_toupper(ch)
  implicit none
  include "ftype.inc"
  ch_toupper = ch_toupper_table(IAchar(ch)+1)
end function ch_toupper

function ch_tolower(ch)
  implicit none
  include "ftype.inc"
  ch_tolower = ch_tolower_table(IAchar(ch)+1)
end function ch_tolower

subroutine str_toupper(str_in,str_out)
  implicit none
  character(len=*) str_in, str_out
  include "ftype.inc"
  integer i,n

  if (Len(str_out)>Len(str_in)) then
    str_out = ' '
    n=Len(str_in)
  else
    n=Len(str_out)
  endif
  do i=1,n
    str_out(i:i)=ch_toupper(str_in(i:i))
  end do
end subroutine str_toupper

subroutine str_tolower(str_in,str_out)
  implicit none
  character(len=*) str_in, str_out
  include "ftype.inc"
  integer i,n

  if (Len(str_out)>Len(str_in)) then
    str_out = ' '
    n=Len(str_in)
  else
    n=Len(str_out)
  endif
  do i=1,n
    str_out(i:i)=ch_tolower(str_in(i:i))
  end do
end subroutine str_tolower

