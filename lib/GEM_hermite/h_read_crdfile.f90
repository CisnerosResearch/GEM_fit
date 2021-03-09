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
!--------------------------------------------------------
subroutine AHCRD_read_file(site_info,filename)
use definition
  implicit none
  type(sites_info)::site_info
  character(len=*) filename
  include "io.fh"
  include "site_names.inc"
  
  integer inunit, ios, allochk, ptr
  integer nsites,lenword,i,j,k,numf,res,save_res,num_res, natoms, mol_chg
  character(len=120)line,word
  character*4 nameat
  double precision crd(3),bohr,bohr_per_angstrom

  ! convert coord to siteic units
  bohr=0.529177249d0
  bohr_per_angstrom = 1.d0 / bohr

  inunit = File_Open(filename,"r")
  if ( inunit < 0 )then
    write(6,*)'failed to open sitefile ',filename(1:Len_Trim(filename))
    stop
  endif
  read(inunit,*,iostat=ios)nsites, natoms, mol_chg
  if ( ios /= 0 )then
    write(6,*)'read_site: bad file: ',filename
    !call error_handler; stop
    close(inunit)
    stop
  endif ! ( ios /= 0 )
! allocate 
  site_info%nsites = nsites
  site_info%natoms = natoms
  site_info%mol_chg = mol_chg

  allocate(site_info%residue_num(nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHCRD_read_file:could not allocate residue_num, exiting'
     stop
  endif
  allocate(site_info%site_crds(3*nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHCRD_read_file:could not allocate site_crds, exiting'
     stop
  endif
  allocate(site_info%local_crds(3*nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHCRD_read_file:could not allocate local_crds, exiting'
     stop
  endif

  ios = 0
  i = 0
  save_res = -1
  num_res = 0
  do while (ios==0)
    line = ' '  ! to avoid junk in line after some reads
    read(inunit,'(A)',iostat=ios)line
    if ( ios == 0 )then
      i = i + 1
      ptr = 1
      call str_next_token(line,ptr,word,lenword)
      sitename(i)(1:lenword) = word(1:lenword)
      do j = lenword+1,4
        sitename(i)(j:j) = ' '
      enddo
      call str_next_token(line,ptr,word,lenword)
      read(word,*)res
      site_info%residue_num(i) = res
      if ( res /= save_res )then
         num_res = num_res + 1
         save_res = res
      endif
      read(line(ptr:),*)(crd(k),k=1,3)
      ! get crds in siteic units
      do k = 1,3
        site_info%site_crds(3*(i-1)+k) = bohr_per_angstrom*crd(k)
      enddo 
    endif !(ios==0) 
  enddo !while (ios==0)
  !write(6,*)'num res = ',num_res
  site_info%num_residues = num_res

  allocate(site_info%residue_start(num_res), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHCRD_read_file:could not allocate residue_start, exiting'
     stop
  endif
  allocate(site_info%residue_end(num_res), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHCRD_read_file:could not allocate residue_end, exiting'
     stop
  endif

  ! read again to assign residue_start
  rewind(inunit)
  read(inunit,*,iostat=ios)nsites
  ios = 0
  i = 0
  save_res = -1
  num_res = 0
  do while (ios==0)
    line = ' '  ! to avoid junk in line after some reads
    read(inunit,'(A)',iostat=ios)line
    if ( ios == 0 )then
      i = i + 1
      ptr = 1
      call str_next_token(line,ptr,word,lenword)
      call str_next_token(line,ptr,word,lenword)
      read(word,*)res
      if ( res /= save_res )then
         num_res = num_res + 1
         save_res = res
         site_info%residue_start(num_res) = i
         if ( num_res > 1 )then
           site_info%residue_end(num_res-1) = i-1
         endif
      endif
    endif !(ios==0) 
  enddo !while (ios==0)
  !IH(p_reshi+num_res-1) = nsites
  site_info%residue_end(num_res) = nsites
 
  !do i = 1,nsites
    !write(6,*)'site name res,x,y,z = ',sitename(i),site_info%residue_num(i), &
              !site_info%site_crd(3*(i-1)+k),k=1,3)
  !enddo
  !do i = 1,num_res
    !write(6,*)'res_start,end = ',site_info%residue_start(i), &
    !                             site_info%residue_end(i)
  !enddo
  call File_close(inunit)
  return
end subroutine AHCRD_read_file
!-------------------------------------------------------------
