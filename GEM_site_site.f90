!GEM_fit: Gaussian Electrostatic Moment fitting
!Copyright (C) 2012  G. Andres Cisneros
!
!This program is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with this program.  If not, see <http://www.gnu.org/licenses/>.
program GEM_site_site

use definition

  implicit none

  include "io.fh"
  type(auxiliary_basis)::atom_auxis
  type(sites_info)::site_info
  double precision ener
  integer argc,arg,l,allochk, i,iargc
  character(len=PATH_MAX) argv
  character(len=80)parmfile,auxfile,crdfile,loc_hermcoeffs
  logical readdens, debug, chgmtp

  readdens = .false.
  debug = .false.
  chgmtp = .false.
  parmfile = ''
  auxfile = ''
  crdfile = ''
  loc_hermcoeffs = ''
  arg=1
  argc=IArgC()
  do while (arg <= argc)
    call GetArg(arg,argv)
    if(argv == '-crd') then
      arg=arg+1
      call GetArg(arg,crdfile)
    elseif ( argv == '-parm')then
      arg=arg+1
      call GetArg(arg,parmfile)
    elseif ( argv == '-aux')then
      arg=arg+1
      call GetArg(arg,auxfile)
    elseif ( argv == '-lherm')then
      arg=arg+1
      call GetArg(arg,loc_hermcoeffs)
    elseif ( argv == '-debug')then
      arg=arg+1
      debug = .true.
      if (debug) open(12,file='GEMoutput.debug',status="unknown")
    endif
    arg = arg + 1
  enddo
  if ( (crdfile == '') .or. (parmfile == '') .or.  &
       (auxfile == '') .or. (loc_hermcoeffs == '') )then
    write(6,*)'usage: ', &
      'GEM_site_site -crd crdfile -parm parmfile -aux auxfile -lherm &
       &lhermcoeffs '
    stop
  endif
 
  call AHCRD_read_file(site_info,crdfile)
  call AHBASE_load_auxnames(atom_auxis,parmfile)
  call AHBASE_load_auxbasis(atom_auxis,auxfile)
  call AHTYPE_load_site_info_2(site_info,atom_auxis,parmfile)
  call AHCOEFF_load_local_hermite(site_info,loc_hermcoeffs)
  call AHFRAME_load_deflist(site_info,parmfile,readdens,debug)
  call AHFRAME_build_frames(site_info)
  call AHFRAME_local_to_global(site_info)

  inquire(file = 'CHARGE_IN_MTP', exist = chgmtp)
  if (.not. chgmtp) then
     call AHSITESITE_interact(site_info,ener)
  else
     call AHSITESITE_interact_CHG(site_info,ener)
  endif

  call H_exchange(site_info,1.00d0,ener) ! for overlap of density (regression)

  if(debug)close(12)
end program GEM_site_site
