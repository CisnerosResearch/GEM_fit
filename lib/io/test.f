program test
implicit none
character(len=128) buffer_long_name

buffer_long_name="Hello World"
write(stdout,*)'First line of file:                                                                   ',buffer_long_name
end
