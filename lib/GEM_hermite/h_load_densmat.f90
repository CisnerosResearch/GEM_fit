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
!-------------------------------------------------------------------------
subroutine HLOAD_dens(site_info, filename, dtype, promolfit)
use definition
  implicit none
  type(sites_info)::site_info
  character(len=*) filename
  character(len=4) tmp_name
  character(len=3) dtype
  
  integer elements, scan_file
  integer i,natoms,num_stos
  logical file_there, promolfit

  file_there = .false.
  scan_file = 0

  inquire(file = filename, exist = file_there)
  if(.not. file_there) then
     write(6,*)'AHLOAD_dens: could not find density file',filename,'exiting'
     stop
  endif

  write (6,1001)filename
1001 format (T2,'Reading density matrix from file: ',a50)
  tmp_name = 'fchk'
  scan_file = index(filename,tmp_name)

  if (scan_file == 0) then !Read density matrix from cadpac or hondo
     !print *,'found denstity, calling load_p'
     call H_read_p(site_info, filename, promolfit)
  else !Read density matrix from gaussian formatted checkpoint
     call H_read_fchk(site_info, filename, dtype, promolfit)
  endif !filename .ne. ''

!! debug print to see if I read it right
!  print *,'------------------- DENS_MAT -----------------------'
!  num_stos = site_info%num_stos
!  elements = (num_stos*(num_stos+1))/2
!  do i = 1, elements
!     print *,site_info%densmat(i)
!  enddo
!  print *,'------------------- DENS_MAT -----------------------'
!! end debug print to see if I read it right

  end subroutine HLOAD_dens
!-------------------------------------------------------------------------
subroutine H_read_p(site_info, filename, promolfit)

use definition
  implicit none
  type(sites_info)::site_info
  !double precision, allocatable:: densmat(:)
  !double precision, densmat(:)
  character(len=*) filename
  include "io.fh"

  integer inunit, ios, allochk
  integer num_stos, elements, i, j, k 
  logical promolfit
  
  inunit = File_Open(filename,"r")
  if ( inunit < 0 )then
    write(6,*)'failed to open dens_mat ',filename(1:Len_Trim(filename))
    stop
  endif

  read(inunit,*,iostat=ios)num_stos
! allocate num_stos
  site_info%num_stos = num_stos 

! allocate density matrix
  elements = (num_stos*(num_stos+1))/2
  if (.not. promolfit) then
     allocate(site_info%densmat(elements), stat = allochk)
  else
     allocate(site_info%promol_dens(elements), stat = allochk)
  endif
  if (allochk .gt. 0 ) then
     write(6,*)'H_read_p:could not allocate densmat, exiting'
     stop
  endif

! read density matrix, need to check if it's from cadpac of hondo because
! it's different format.

    do i=1,elements
       if (.not. promolfit) then
          read(inunit,*,iostat=ios)site_info%densmat(i)
       else
          read(inunit,*,iostat=ios)site_info%promol_dens(i)
       endif
       if (ios<0) goto 50
50  enddo

  call File_close(inunit)

  end subroutine H_read_p
!-------------------------------------------------------------------------
subroutine H_read_fchk(site_info, filename, dtype, promolfit)

use definition
  implicit none
  type(sites_info)::site_info
  character(len=*) filename
  character(len=17) dens
  character(len=17) dens1
  character(len=17) dens2
  character(len=17) dens3
  character(len=16) dens4
  character(len=16) dens5
  character(len=3) dtype
  character(len=15) natm
  character(len=25) stos
  character(len=80) line
  include "io.fh"

  integer natoms, num_stos
  integer inunit, scan_res_dens, scan_res_atm, scan_res_stos, scan_res_dens1
  integer scan_res_dens2, scan_res_dens3, scan_res_dens4, scan_res_dens5
  integer ios, counter, ptr, endline, testend, allochk, posthf

  integer dim, elements, temp, j, k, l, stat, ngau, k1, temp2, k2
  real*8 tmpvec(5)
  logical end, close_file, end_of_file, promolfit

  inunit = File_Open(filename,"r")
  if ( inunit < 0 )then
    write(6,*)'failed to open fchk ',filename
    stop
  endif

  dens = 'Total SCF Density'
  dens1 = 'Total MP2 Density'
  dens2 = 'Total MP3 Density'
  dens3 = 'Total MP4 Density'
  dens4 = 'Total CC Density'
  dens5 = 'Total CI Density'
  natm = 'Number of atoms'
  stos = 'Number of basis functions'
  ios = 0
  scan_res_dens = 0
  scan_res_dens1 = 0
  scan_res_dens2 = 0
  scan_res_dens3 = 0
  scan_res_dens4 = 0
  scan_res_dens5 = 0
  scan_res_atm = 0
  scan_res_stos = 0
  counter = 0
  end_of_file = .false.
  dtype='SCF'
  posthf = 0
  do while (ios==0 .and. .not. end_of_file)
    line = ' '  ! to avoid junk in line after some reads
    read(inunit,'(A)',iostat=ios)line
    scan_res_dens = index(line(1:80),dens)
    scan_res_dens1 = index(line(1:80),dens1)
    scan_res_dens2 = index(line(1:80),dens2)
    scan_res_dens3 = index(line(1:80),dens3)
    scan_res_dens4 = index(line(1:80),dens4)
    scan_res_dens5 = index(line(1:80),dens5)
    if (scan_res_dens1 > 0 .or. scan_res_dens2 > 0 .or. scan_res_dens3 > 0 &
        .or. scan_res_dens4 > 0 .or. scan_res_dens5 > 0) then
        if (scan_res_dens1 > 0) then
           write (6,1002)dens1
           dtype='MP2'
        else if (scan_res_dens2 > 0) then
           write (6,1002)dens2
           dtype='MP3'
        else if (scan_res_dens3 > 0) then
           write (6,1002)dens3
           dtype='MP4'
        else if (scan_res_dens4 > 0) then
           write (6,1002)dens4
           dtype='CC '
        else if (scan_res_dens5 > 0) then
           write (6,1002)dens5
           dtype='CI '
        endif
1002 format ('Post-HF density found, will read: ',a20)
        scan_res_dens = 0
        posthf = 1
    endif
    scan_res_atm = index(line(1:80),natm)
    scan_res_stos = index(line(1:80),stos)
    if (ios == 0 .and. scan_res_atm > 0)then
       read(line(51:61),'(i11)') natoms
       if (natoms .ne. site_info%natoms) then
          write(6,*)'AH_read_fchk: natoms in fchk .ne. natoms in crd, exiting'
          write(6,*)natoms,site_info%natoms
          stop
       endif
    else if (ios == 0 .and. scan_res_stos > 0)then
       read(line(51:61),'(i11)') num_stos
       site_info%num_stos = num_stos 
     ! allocate density matrix
       elements = (num_stos*(num_stos+1))/2
       if(.not. promolfit) then
          allocate(site_info%densmat(elements), stat = allochk)
       else
          allocate(site_info%promol_dens(elements), stat = allochk)
       endif
       if (allochk .gt. 0 ) then
          write(6,*)'H_read_fchk:could not allocate densmat, exiting'
          stop
       endif
    endif !(ios==0)
  enddo !while (ios==0)
  
  rewind(inunit)

  ios = 0
  scan_res_dens = 0
  scan_res_dens1 = 0
  scan_res_dens2 = 0
  scan_res_dens3 = 0
  scan_res_dens4 = 0
  scan_res_dens5 = 0
  counter = 0
  end_of_file = .false.
  do while (ios==0 .and. .not. end_of_file)
    line = ' '  ! to avoid junk in line after some reads
    read(inunit,'(A)',iostat=ios)line

    if (posthf == 0) then
       scan_res_dens = index(line(1:80),dens)
       if (ios == 0 .and. scan_res_dens > 0)then
       do while (counter .lt. elements)
          testend = elements - counter
          if (testend .ge. 5) then
             endline = 5
          else
             endline = testend
          endif
          read(inunit,*)(tmpvec(k),k=1,endline)
          counter = counter+endline
          do j = 1, endline
             if(.not. promolfit) then
                site_info%densmat(counter-endline+j) = tmpvec(j)
             else
                site_info%promol_dens(counter-endline+j) = tmpvec(j)
             endif
          enddo
       enddo 
       endif ! ios
    else ! posthf
       scan_res_dens1 = index(line(1:80),dens1)
       scan_res_dens2 = index(line(1:80),dens2)
       scan_res_dens3 = index(line(1:80),dens3)
       scan_res_dens4 = index(line(1:80),dens4)
       scan_res_dens5 = index(line(1:80),dens5)
       if (ios == 0 .and. scan_res_dens1 > 0 .or. scan_res_dens2 > 0 .or. &
           scan_res_dens3 > 0 .or. scan_res_dens4 > 0 .or. &
           scan_res_dens5 > 0)then ! POST-HF density
       !if (scan_res_dens1 > 0) print *,'FOUND MP2',scan_res_dens1a
          do while (counter .lt. elements)
             testend = elements - counter
             if (testend .ge. 5) then
                endline = 5
             else
                endline = testend
             endif
             read(inunit,*)(tmpvec(k),k=1,endline)
             counter = counter+endline
             do j = 1, endline
                if(.not. promolfit) then
                   site_info%densmat(counter-endline+j) = tmpvec(j)
                else
                   site_info%promol_dens(counter-endline+j) = tmpvec(j)
                endif
             enddo
          enddo 
       endif ! ios
    endif ! posthf
  enddo !while (ios==0)

  if (counter .ne. elements) then
     write(6,*)'H_read_fchk:counter .ne. num_stos, exiting'
     write(6,*)counter,num_stos
     stop
  endif

  call File_close(inunit)

  end subroutine H_read_fchk
!-------------------------------------------------------------------------
