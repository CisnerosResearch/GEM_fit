
  integer MAXBASIS
  parameter(MAXBASIS=20)
  character(len=60) auxname(MAXBASIS),tmp_name 
  character(len=5) NAMEFMT
  parameter (NAMEFMT='(A60)')
  common/bnames/auxname
