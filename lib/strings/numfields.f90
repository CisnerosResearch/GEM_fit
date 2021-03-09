integer function Num_fields(str)
  implicit none
  character(len=*) str

  character(len=120) word
  integer ptr,lenword

  Num_fields = 0
  ptr = 1
  do while ( ptr > 0 )
    call str_next_token(str,ptr,word,lenword)
    if ( lenword > 0 )Num_fields = Num_fields + 1
  enddo ! while ( ptr > 0 )

  return
end function Num_fields

