
  integer num_AO_basis
  integer MAXBASIS
  parameter(MAXBASIS=20)
  character(len=60) auxname(MAXBASIS),tmp_name,basname(MAXBASIS)
  character(len=5) NAMEFMT
  parameter (NAMEFMT='(A60)')
  common/numbas/num_AO_basis
  common/bnames/auxname,basname
