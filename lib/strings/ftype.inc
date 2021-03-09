
 include "ftype.fh"

!#ifndef IN_BLOCKDATA
external ftype_blockdata
!#endif
character ch

integer*2 ctype_class_table(256)
character(len=1) ch_toupper_table(256), ch_tolower_table(256)

common /ftype_private/ ctype_class_table, ch_toupper_table, ch_tolower_table
save /ftype_private/

integer*2 ISupper_bit;  parameter (ISupper_bit = z'0001') ! UPPERCASE
integer*2 ISlower_bit;  parameter (ISlower_bit = z'0002') ! lowercase
integer*2 ISalpha_bit;  parameter (ISalpha_bit = z'0004') ! Alphabetic
integer*2 ISdigit_bit;  parameter (ISdigit_bit = z'0008') ! Numeric
integer*2 ISxdigit_bit; parameter (ISxdigit_bit= z'0010') ! Hexadecimal numeric
integer*2 ISspace_bit;  parameter (ISspace_bit = z'0020') ! Whitespace
integer*2 ISprint_bit;  parameter (ISprint_bit = z'0040') ! Printing
integer*2 ISgraph_bit;  parameter (ISgraph_bit = z'0080') ! Graphical
integer*2 ISblank_bit;  parameter (ISblank_bit = z'0100') ! Blank (usually SPC and TAB)
integer*2 IScntrl_bit;  parameter (IScntrl_bit = z'0200') ! Control character
integer*2 ISpunct_bit;  parameter (ISpunct_bit = z'0400') ! Punctuation
integer*2 ISalnum_bit;  parameter (ISalnum_bit = z'0800') ! Alphanumeric

