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
!--------------------------------------------------------------
subroutine AHBASE_load_auxnames(atom_auxis,parmfile)
use definition
  implicit none
  type(auxiliary_basis)::atom_auxis
  character(len=*) parmfile
  include "io.fh"
  include "basis_names.inc"

  integer parmunit,ios, ptr
  character(len=120) line,word
  integer TRDPRM_get_next_line,lenword
  integer nbasis,jbasis

  parmunit = File_Open(parmfile,"r")
  if ( parmunit < 0 )then
    write(6,*)'failed to open parmfile ',parmfile(1:Len_Trim(parmfile))
    stop
  endif
!
! read what auxiliary basis corresponds to each site type
!
  ios = 0
  nbasis = 0
  do while (ios==0)
    line = ' '  ! to avoid junk in line after some reads
    ios = TRDPRM_get_next_line(parmunit,line)
    if ( ios == 0 )then
      ptr = 1 ! first call to str_next_token needs an initial value
      call str_next_token(line,ptr,word,lenword) !get 1st token
      if ( ptr > 0 )then !check if its 'auxbasis'
        if ( word(1:lenword) == 'auxbasis' )then
           nbasis = nbasis + 1
           read(line(10:),NAMEFMT) auxname(nbasis)
        endif
      endif
    endif
  enddo !while ios==0
  call File_Close(parmunit)
  write(6,*)'loadbasis: nbasis = ',nbasis

  !ptr = DB_alloc(p_basis,DB_int,"num_auxbasis",1)
  !IH(ptr) = nbasis
  atom_auxis%num_auxbasis = nbasis

  !do jbasis = 1,nbasis
    !write(6,*)auxname(jbasis)
  !enddo

  return
end subroutine AHBASE_load_auxnames
!--------------------------------------------------------------
subroutine AHBASE_load_auxbasis(atom_auxis,auxfile)
use definition
  implicit none
  type(auxiliary_basis)::atom_auxis
  character(len=*) auxfile
  include "io.fh"
  include "basis_names.inc"

  integer auxunit,ios, allochk
  integer jbasis,nbasis,nprim,nshell,nprimshell,ordprim,order,ncoeff,k,n
  character(len=120) line,word
  integer TRDPRM_get_next_line,lenword,TRDPRM_numparms,scan_res
  double precision rnum
  logical end_of_file,end_line

  auxunit = File_Open(auxfile,"r")
  if ( auxunit < 0 )then
    write(6,*)'AHBASE_load_auxbasis: failed to open aux file ',&
               & auxfile(1:Len_Trim(auxfile))
    stop
  endif

  nprim = 0
  ncoeff = 0
  ! first pass count quantities
  nbasis = atom_auxis%num_auxbasis

  do jbasis = 1, nbasis
    rewind(auxunit)
    ios = 0
    end_of_file = .false.
    end_line = .false.
    tmp_name = auxname(jbasis)
    !print *,'tmp name ',tmp_name
    do while (ios == 0 .and. .not. end_of_file .and. .not. end_line )
      read(auxunit,'(a80)',iostat=ios) line
      if ( ios == 0 .and. line(1:3) .eq. 'END')then
        end_line = .true.
      elseif ( ios == 0 )then
        scan_res = index(line(1:80),tmp_name(1:Len_Trim(tmp_name)))
        if (scan_res > 0 )then
          read(auxunit,'(i5)')nshell
          do n = 1,nshell
            read(auxunit,*)nprimshell,ordprim
            if ( ordprim == 0 )then
              order = 1
            elseif ( ordprim == 1 )then
              order = 4
            elseif ( ordprim == -1 )then
              order = 3
            elseif ( ordprim == -2 )then ! GAC for d shells
              order = 6
            elseif ( ordprim == 2 )then
              order = 10
            elseif ( ordprim == 3 )then
              order = 20
              !print *,'NORMALIZING f shell, 20 prims'
              !print *,'order f = ',order
            else
              write(6,*)'bad ordprim value',ordprim
              stop
            endif
            nprim = nprim + nprimshell
            ncoeff = ncoeff + order*nprimshell
            !print *,'ncoeff = ',ncoeff
            do k = 1,nprimshell
              read(auxunit,*)rnum
            enddo
          enddo !n = 1,nshell
        endif !scan_res > 0
      endif !ios == 0 .and. line(1:3) .eq. 'END'
    enddo !while ios == 0 .and. .not. end_of_file .and. .not. end_line
    !print *,'end basis ',jbasis
  enddo ! jbasis = 1, nbasis

  ! next allocate
  atom_auxis%num_primitives = nprim

  atom_auxis%num_coefficients = ncoeff

  !print *,'ncoeff = ',ncoeff

  allocate(atom_auxis%num_basis_prims(nbasis), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHBASE_load_auxbasis:could not allocate num_basis_prims, exiting'
     stop
  endif

  allocate(atom_auxis%off_basis_prims(nbasis), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHBASE_load_auxbasis:could not allocate off_basis_prims, exiting'
     stop
  endif

  allocate(atom_auxis%prim_order(nprim), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHBASE_load_auxbasis:could not allocate prim_order, exiting'
     stop
  endif

  allocate(atom_auxis%prim_expon_term(nprim), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHBASE_load_auxbasis:could not allocate prim_expon_term, &
                &exiting'
     stop
  endif

  nprim = 0
  do jbasis = 1, nbasis
    rewind(auxunit)
    ios = 0
    end_of_file = .false.
    end_line = .false.
    tmp_name = auxname(jbasis)
    do while (ios == 0 .and. .not. end_of_file .and. .not. end_line )
      read(auxunit,'(a80)',iostat=ios) line
      if ( ios == 0 .and. line(1:3) .eq. 'END')then
        end_line = .true.
      elseif ( ios == 0 )then
        scan_res = index(line(1:80),tmp_name(1:Len_Trim(tmp_name)))
        if (scan_res > 0 )then
          atom_auxis%num_basis_prims(jbasis) = 0
          read(auxunit,'(i5)')nshell
          do n = 1,nshell
            read(auxunit,*)nprimshell,ordprim
            if ( ordprim == 0 )then
              order = 1
            elseif ( ordprim == 1 )then
              order = 4
            elseif ( ordprim == -1 )then
              order = 3
            elseif ( ordprim == -2 )then ! GAC for d shells
              order = 6
            elseif ( ordprim == 2 )then
              order = 10
            elseif ( ordprim == 3 )then
              !print *,'NORMALIZING f shell, 20 prims'
              order = 20
            else
              write(6,*)'bad ordprim value',ordprim
              stop
            endif
            atom_auxis%num_basis_prims(jbasis) = &
                      atom_auxis%num_basis_prims(jbasis) + nprimshell
            do k = 1,nprimshell
              nprim = nprim + 1
              atom_auxis%prim_order(nprim) = order
              !if (order == 20) print *,'order  = ',atom_auxis%prim_order(nprim)
              read(auxunit,*)atom_auxis%prim_expon_term(nprim)
            enddo
          enddo !n = 1,nshell
        endif !scan_res > 0
      endif !ios == 0 .and. line(1:3) .eq. 'END'
    enddo !(ios == 0 .and. .not. end_of_file .and. .not. end_line )
  enddo ! jbasis = 1, nbasis
  !print *,'ncoeff = ',ncoeff

  ! do offsets
  atom_auxis%off_basis_prims(1) = 0
  do jbasis = 2,nbasis
    atom_auxis%off_basis_prims(jbasis) = atom_auxis%off_basis_prims(jbasis-1)&
                                       + atom_auxis%num_basis_prims(jbasis-1)
  enddo
  call File_Close(auxunit)
  return
end subroutine AHBASE_load_auxbasis
!--------------------------------------------------------------
subroutine AHBASE_calcnorms(atom_auxis,densfile)
!subroutine AHBASE_calcnorms(atom_auxis)
!! GAC 8/05: added call to calculate number of electrons per auxiliary
use definition
  implicit none
  type(auxiliary_basis)::atom_auxis

  integer nprim,ncoeff,order,jprim,jcoeff,nx,ny,nz,n,k, allochk
  integer scan_file
  double precision pi, fac(0:14), dfac(0:18),norm, ex, neaux
  double precision zero, one, two
  character(len=*) densfile
  character(len=4) tmp_name
  logical nw_dens
  parameter(zero=0.d0, one=1.d0, two=2.d0)
  pi = 4.d0*(atan(one))

! check if it's a NWChem density file
  nw_dens = .false.
  scan_file = 0
  tmp_name = 'NW'
  scan_file = index(densfile,tmp_name)
  if (scan_file .ne. 0) nw_dens = .true.
  if(nw_dens)write(6,*)'AHBASE_calcnorms: NWChem density found, &
                        &reordering d,f'

  fac(0) = one
  fac(1) = one
  do n = 2,14
    fac(n) = dble(n)*fac(n-1)
  enddo
  dfac(0) = one
  dfac(1) = one
  dfac(2) = two
  do n = 3, 18
     dfac(n) = dble(n)*dfac(n-2)
  enddo

  nprim = atom_auxis%num_primitives
  ncoeff = atom_auxis%num_coefficients

  allocate(atom_auxis%norm_offset(nprim), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHBASE_calcnorms:could not allocate norm_offset, exiting'
     stop
  endif

  allocate(atom_auxis%prim_normalizers(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHBASE_calcnorms:could not allocate prim_normalizers, exiting'
     stop
  endif

  allocate(atom_auxis%prim_elec(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHBASE_calcnorms:could not allocate prim_elec, exiting'
     stop
  endif

  atom_auxis%norm_offset(1) = 0 
  do n = 2,nprim
    !if (atom_auxis%prim_order(n-1) .le. 10) then
       atom_auxis%norm_offset(n) = atom_auxis%norm_offset(n-1) + &
                                   atom_auxis%prim_order(n-1) 
    !else
    !   atom_auxis%norm_offset(n) = atom_auxis%norm_offset(n-1) + 10
    !endif
  enddo
  jcoeff = 0
  do jprim = 1,nprim
     order = atom_auxis%prim_order(jprim)
     if (order .eq. 1) then ! only s shell
        ex = atom_auxis%prim_expon_term(jprim)
        nx = 0
        ny = 0
        nz = 0
        norm = zero
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        jcoeff = jcoeff + 1
        atom_auxis%prim_normalizers(jcoeff) = norm
        atom_auxis%prim_elec(jcoeff) = neaux
     else if (order .eq. 3) then ! sp shell
        ! p normalization 
        nx = 1
        ny = 0
        nz = 0
        norm = zero
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        do k = 1, 3
           jcoeff = jcoeff + 1
           atom_auxis%prim_normalizers(jcoeff) = norm
           atom_auxis%prim_elec(jcoeff) = neaux
        enddo
     else if (order .eq. 4) then ! sp shell
        ! s normalization 
        ex = atom_auxis%prim_expon_term(jprim)
        nx = 0
        ny = 0
        nz = 0
        norm = zero
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        jcoeff = jcoeff + 1
        atom_auxis%prim_normalizers(jcoeff) = norm
        atom_auxis%prim_elec(jcoeff) = neaux
        ! p normalization 
        nx = 1
        ny = 0
        nz = 0
        norm = zero
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        do k = 1, 3
           jcoeff = jcoeff + 1
           atom_auxis%prim_normalizers(jcoeff) = norm
           atom_auxis%prim_elec(jcoeff) = neaux
        enddo
     else if (order .eq. 6) then ! d only shell
      if(.not.nw_dens) then
        ex = atom_auxis%prim_expon_term(jprim)
        ! d normalization 
        nx = 2
        ny = 0
        nz = 0
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        do k = 1, 3
           jcoeff = jcoeff + 1
           atom_auxis%prim_normalizers(jcoeff) = norm
           atom_auxis%prim_elec(jcoeff) = neaux
        enddo
        nx = 1
        ny = 1
        nz = 0
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        do k = 1, 3
           jcoeff = jcoeff + 1
           !RH(p_norm+jcoeff-1) = norm
           atom_auxis%prim_normalizers(jcoeff) = norm
           atom_auxis%prim_elec(jcoeff) = neaux
        enddo
      else
        ex = atom_auxis%prim_expon_term(jprim)
        ! d normalization 
        nx = 2
        ny = 0
        nz = 0
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        do k = 1, 6
           jcoeff = jcoeff + 1
           atom_auxis%prim_normalizers(jcoeff) = norm
           atom_auxis%prim_elec(jcoeff) = neaux
        enddo
      endif
     else if (order .eq. 10) then ! spd shell
        ! s normalization 
        ex = atom_auxis%prim_expon_term(jprim)
        nx = 0
        ny = 0
        nz = 0
        norm = zero
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        jcoeff = jcoeff + 1
        atom_auxis%prim_normalizers(jcoeff) = norm
        atom_auxis%prim_elec(jcoeff) = neaux
        ! p normalization 
        nx = 1
        ny = 0
        nz = 0
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        do k = 1, 3
           jcoeff = jcoeff + 1
           atom_auxis%prim_normalizers(jcoeff) = norm
           atom_auxis%prim_elec(jcoeff) = neaux
        enddo
      if(.not.nw_dens) then
        ! d normalization 
        nx = 2
        ny = 0
        nz = 0
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        do k = 1, 3
           jcoeff = jcoeff + 1
           atom_auxis%prim_normalizers(jcoeff) = norm
           atom_auxis%prim_elec(jcoeff) = neaux
        enddo
        nx = 1
        ny = 1
        nz = 0
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        do k = 1, 3
           jcoeff = jcoeff + 1
           atom_auxis%prim_normalizers(jcoeff) = norm
           atom_auxis%prim_elec(jcoeff) = neaux
        enddo
      else
        ! d normalization 
        !nx = 2
        !ny = 0
        !nz = 0
        !call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        !call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        !do k = 1, 6
        !   jcoeff = jcoeff + 1
        !   atom_auxis%prim_normalizers(jcoeff) = norm
        !   atom_auxis%prim_elec(jcoeff) = neaux
        !enddo
        nx = 2
        ny = 0
        nz = 0
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        jcoeff = jcoeff + 1
        atom_auxis%prim_normalizers(jcoeff) = norm
        atom_auxis%prim_elec(jcoeff) = neaux
        nx = 1
        ny = 1
        nz = 0
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        jcoeff = jcoeff + 1
        atom_auxis%prim_normalizers(jcoeff) = norm
        atom_auxis%prim_elec(jcoeff) = neaux
        nx = 1
        ny = 0
        nz = 1
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        jcoeff = jcoeff + 1
        atom_auxis%prim_normalizers(jcoeff) = norm
        atom_auxis%prim_elec(jcoeff) = neaux
        nx = 0
        ny = 2
        nz = 0
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        jcoeff = jcoeff + 1
        atom_auxis%prim_normalizers(jcoeff) = norm
        atom_auxis%prim_elec(jcoeff) = neaux
        nx = 0
        ny = 1
        nz = 1
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        jcoeff = jcoeff + 1
        atom_auxis%prim_normalizers(jcoeff) = norm
        atom_auxis%prim_elec(jcoeff) = neaux
        nx = 0
        ny = 0
        nz = 2
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        jcoeff = jcoeff + 1
        atom_auxis%prim_normalizers(jcoeff) = norm
        atom_auxis%prim_elec(jcoeff) = neaux
      endif
     else if (order .eq. 20) then ! f shell
       !print *,'NORMALIZING f shell, 20 prims'
        ! s normalization 
        ex = atom_auxis%prim_expon_term(jprim)
        nx = 0
        ny = 0
        nz = 0
        norm = zero
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        jcoeff = jcoeff + 1
        atom_auxis%prim_normalizers(jcoeff) = norm
        atom_auxis%prim_elec(jcoeff) = neaux
        ! p normalization 
        nx = 1
        ny = 0
        nz = 0
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        do k = 1, 3
           jcoeff = jcoeff + 1
           atom_auxis%prim_normalizers(jcoeff) = norm
           atom_auxis%prim_elec(jcoeff) = neaux
        enddo
      if(.not.nw_dens) then
        ! d normalization 
        nx = 2
        ny = 0
        nz = 0
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        do k = 1, 3
           jcoeff = jcoeff + 1
           atom_auxis%prim_normalizers(jcoeff) = norm
           atom_auxis%prim_elec(jcoeff) = neaux
        enddo
        nx = 1
        ny = 1
        nz = 0
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        do k = 1, 3
           jcoeff = jcoeff + 1
           atom_auxis%prim_normalizers(jcoeff) = norm
           atom_auxis%prim_elec(jcoeff) = neaux
        enddo
      else
        ! d normalization 
        !nx = 2
        !ny = 0
        !nz = 0
        !call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        !call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        !do k = 1, 6
        !   jcoeff = jcoeff + 1
        !   atom_auxis%prim_normalizers(jcoeff) = norm
        !   atom_auxis%prim_elec(jcoeff) = neaux
        !enddo
        nx = 2
        ny = 0
        nz = 0
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        jcoeff = jcoeff + 1
        atom_auxis%prim_normalizers(jcoeff) = norm
        atom_auxis%prim_elec(jcoeff) = neaux
        nx = 1
        ny = 1
        nz = 0
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        jcoeff = jcoeff + 1
        atom_auxis%prim_normalizers(jcoeff) = norm
        atom_auxis%prim_elec(jcoeff) = neaux
        nx = 1
        ny = 0
        nz = 1
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        jcoeff = jcoeff + 1
        atom_auxis%prim_normalizers(jcoeff) = norm
        atom_auxis%prim_elec(jcoeff) = neaux
        nx = 0
        ny = 2
        nz = 0
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        jcoeff = jcoeff + 1
        atom_auxis%prim_normalizers(jcoeff) = norm
        atom_auxis%prim_elec(jcoeff) = neaux
        nx = 0
        ny = 1
        nz = 1
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        jcoeff = jcoeff + 1
        atom_auxis%prim_normalizers(jcoeff) = norm
        atom_auxis%prim_elec(jcoeff) = neaux
        nx = 0
        ny = 0
        nz = 2
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        jcoeff = jcoeff + 1
        atom_auxis%prim_normalizers(jcoeff) = norm
        atom_auxis%prim_elec(jcoeff) = neaux
      endif
      !if(.not.nw_dens) then
        ! f^3 normalization 
        nx = 3
        ny = 0
        nz = 0
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        do k = 1, 3
           jcoeff = jcoeff + 1
           atom_auxis%prim_normalizers(jcoeff) = norm
           atom_auxis%prim_elec(jcoeff) = neaux
        enddo
        ! f^2,1 normalization 
        nx = 2
        ny = 1
        nz = 0
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        do k = 1, 6
           jcoeff = jcoeff + 1
           atom_auxis%prim_normalizers(jcoeff) = norm
           atom_auxis%prim_elec(jcoeff) = neaux
        enddo
        ! f^1,1,1 normalization 
        nx = 1
        ny = 1
        nz = 1
        call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
        call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
        jcoeff = jcoeff + 1
        atom_auxis%prim_normalizers(jcoeff) = norm
        atom_auxis%prim_elec(jcoeff) = neaux
      !else
      !  ! f normalization 
      !  nx = 3
      !  ny = 0
      !  nz = 0
      !  call normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
      !  call aux_elec(nx, ny, nz, ex, norm, neaux, pi, dfac)
      !  do k = 1, 10
      !     jcoeff = jcoeff + 1
      !     atom_auxis%prim_normalizers(jcoeff) = norm
      !     atom_auxis%prim_elec(jcoeff) = neaux
      !  enddo
      !endif
     else
      write(6,*)'AH_loadbasis: can only handle up to F auxis, exiting'
      stop
     endif ! order .eq. 1) 
  enddo !n = 1,nprim
  !do jcoeff = 1,ncoeff
  !  write(6,*)'norm = ',atom_auxis%prim_normalizers(jcoeff),&
  !                      atom_auxis%prim_elec(jcoeff)
  !!  write(6,*)'neaux = ',atom_auxis%prim_elec(jcoeff)
  !enddo

  return
end subroutine AHBASE_calcnorms
!-------------------------------------------------------------------------
subroutine normalize_aux(nx, ny, nz, ex, norm, pi, fac, dfac)
  implicit none
  integer nx,ny,nz
  double precision ex,norm,pi, fac(0:14),dfac(0:18)

  integer iax,jax,iay,jay,iaz,jaz,lc
  double precision lnorm,factor,wax,way,waz,auxfun

  lc = nx+ny+nz
  lnorm = sqrt(2.d0)*(pi/ex)**(2.5d0)/(4.d0*ex)**lc
  factor = lnorm*(fac(nx)*fac(ny)*fac(nz))**2
  norm = 0.d0

  do iax = 0,nx/2
    do jax = 0,nx/2
      wax = auxfun(nx,iax,jax,dfac,fac)
      do iay = 0,ny/2
         do jay = 0,ny/2
            way = auxfun(ny,iay,jay,dfac,fac)
            do iaz = 0,nz/2
               do jaz = 0,nz/2
                  waz = auxfun(nz,iaz,jaz,dfac,fac)
                  norm = norm + (wax*way*waz/ &
     &                  (2*(lc-iax-iay-iaz-jax-jay-jaz)+1))
               enddo
            enddo!iaz = 0,nz/2
         enddo!jay = 0,ny/2
      enddo!iay = 0,ny/2
    enddo!jax = 0,nx/2
  enddo!iax = 0,nx/2
  norm = 1.d0 / sqrt(norm*factor)
  return
end subroutine normalize_aux
!-------------------------------------------------------------------------
subroutine aux_elec(nx, ny, nz, ex, nf, neaux, pi, dfac)
  implicit none
  
  integer nx, ny, nz, dim
  double precision ex, nf, neaux, pi, tx, ty, tz
  double precision dfac(0:18), factor
  integer i, n, lc 
  double precision two, one
  parameter(one=1.d0, two=2.d0)

  ! Loop over all primitives

  neaux = 0.d0

  lc = nx+ny+nz
  factor = dsqrt(pi**3/(two**lc*ex**(lc+3)))
  
  if((mod(nx,2).eq.0).and.(mod(ny,2).eq.0).and.(mod(nz,2).eq.0)) then
     if(nx.eq.0)then
       tx = one
     else
       tx = dfac(nx-1)
     endif
     if(ny.eq.0)then
       ty = one
     else
       ty = dfac(ny-1)
     endif
            if(nz.eq.0)then
       tz = one
     else
       tz = dfac(nz-1)
     endif
     neaux = neaux+(tx*ty*tz*nf*factor)
  else
     neaux = neaux + 0.d0
  endif
  return
end subroutine aux_elec
!-------------------------------------------------------------------------
double precision function auxfun(n,i,j,dfac,fac)
  implicit none
  integer n,i,j
  double precision dfac(0:18),fac(0:14)
  
  double precision tmp
  if ((2*(n-i-j)-1) < 0) then
     tmp = 1.d0
  else
     tmp = dfac(2*(n-i-j)-1)
  endif
  auxfun = (-1.d0)**(i+j)*tmp/(fac(i)*fac(j)*fac(n-2*i)*fac(n-2*j))
  return
end function auxfun
!-------------------------------------------------------------------------
!  BASIS READ
!-------------------------------------------------------------------------
subroutine AHBASE_load_basnames(atom_basis,parmfile)
use definition
  implicit none
  type(AO_basis)::atom_basis
  character(len=*) parmfile
  include "io.fh"
  include "basis_names.inc"

  integer parmunit,ios, ptr
  character(len=120) line,word
  integer TRDPRM_get_next_line,lenword
  integer nbasis,jbasis

  parmunit = File_Open(parmfile,"r")
  if ( parmunit < 0 )then
    write(6,*)'AHBASE_load_basnames: failed to open parmfile ',&
               & parmfile(1:Len_Trim(parmfile))
    stop
  endif
!
! read what AO basis corresponds to each site type
!
  ios = 0
  nbasis = 0
  do while (ios==0)
    line = ' '  ! to avoid junk in line after some reads
    ios = TRDPRM_get_next_line(parmunit,line)
    if ( ios == 0 )then
      ptr = 1 ! first call to str_next_token needs an initial value
      call str_next_token(line,ptr,word,lenword) !get 1st token
      if ( ptr > 0 )then !check if its 'auxbasis'
        if ( word(1:lenword) == 'AObasis' )then
           nbasis = nbasis + 1
           read(line(9:),NAMEFMT) basname(nbasis)
        endif
      endif
    endif
  enddo !while ios==0
  call File_Close(parmunit)
  !write(6,*)'loadbasis: nbasis = ',nbasis

  atom_basis%num_basis = nbasis

  !do jbasis = 1,nbasis
    !write(6,*)basname(jbasis)
  !enddo

  return
end subroutine AHBASE_load_basnames
!--------------------------------------------------------------
subroutine AHBASE_load_AObasis(atom_basis,auxfile)
! GAC similar to load_basis but this one sets up structure for contracted
!     basis sets
use definition
  implicit none
  type(AO_basis)::atom_basis
  character(len=*) auxfile
  include "io.fh"
  include "basis_names.inc"

  integer auxunit,ios, allochk, times, temp, nprim2, nindex, zetind
  integer jbasis,nbasis,nprim,nshell,nprimshell,ordprim,order,ncoeff,k,n
  character(len=120) line,word
  integer TRDPRM_get_next_line,lenword,TRDPRM_numparms,scan_res,contract
  double precision rnum
  logical end_of_file,end_line

  auxunit = File_Open(auxfile,"r")
  if ( auxunit < 0 )then
    write(6,*)'AHBASE_load_AObasis: failed to open bas file ',&
               & auxfile(1:Len_Trim(auxfile))
    stop
  endif

  nprim = 0
  ncoeff = 0
  ! first pass count quantities
  nbasis = atom_basis%num_basis

  do jbasis = 1, nbasis
    rewind(auxunit)
    ios = 0
    end_of_file = .false.
    end_line = .false.
    tmp_name = basname(jbasis)
    do while (ios == 0 .and. .not. end_of_file .and. .not. end_line )
      read(auxunit,'(a80)',iostat=ios) line
      if ( ios == 0 .and. line(1:3) .eq. 'END')then
        end_line = .true.
      elseif ( ios == 0 )then
        scan_res = index(line(1:80),tmp_name(1:Len_Trim(tmp_name)))
        if (scan_res > 0 )then
          read(auxunit,'(i5)')nshell
          do n = 1,nshell
            read(auxunit,*)nprimshell,ordprim,contract
            if ( ordprim == 0 )then
              order = 1
            elseif ( ordprim == 1 )then
              order = 3
            elseif ( ordprim == 2 )then
              order = 6
            elseif ( ordprim == 3 )then
              order = 10
            else
              write(6,*)'bad ordprim value'
              stop
            endif
            nprim = nprim + contract
            ncoeff = ncoeff + order*contract
            do k = 1,contract
              read(auxunit,*)rnum
            enddo
          enddo !n = 1,nshell
        endif !scan_res > 0
      endif !ios == 0 .and. line(1:3) .eq. 'END'
    enddo !while ios == 0 .and. .not. end_of_file .and. .not. end_line
  enddo ! jbasis = 1, nbasis

  ! next allocate
  atom_basis%num_primitives = nprim

  atom_basis%num_coefficients = ncoeff

  print *,'nprim, ncoeff = ',nprim,ncoeff

  allocate(atom_basis%num_basis_prims(nbasis), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHBASE_load_AObasis:could not allocate num_basis_prims, exiting'
     stop
  endif

  allocate(atom_basis%off_basis_prims(nbasis), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHBASE_load_AObasis:could not allocate off_basis_prims, exiting'
     stop
  endif

  allocate(atom_basis%prim_order(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHBASE_load_AObasis:could not allocate prim_order, exiting'
     stop
  endif

  allocate(atom_basis%expon(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHBASE_load_AObasis:could not allocate expon, &
                &exiting'
     stop
  endif

  allocate(atom_basis%contr(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHBASE_load_AObasis:could not allocate contr, &
                &exiting'
     stop
  endif

  allocate(atom_basis%index(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHBASE_load_AObasis:could not allocate contr, &
                &exiting'
     stop
  endif

  allocate(atom_basis%orb_contr(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHBASE_load_AObasis:could not allocate orb_contr, &
                &exiting'
     stop
  endif

  nprim = 0
  nprim2 = 0
  ncoeff = 0
  do jbasis = 1, nbasis
    nindex = 0
    rewind(auxunit)
    ios = 0
    end_of_file = .false.
    end_line = .false.
    tmp_name = basname(jbasis)
    do while (ios == 0 .and. .not. end_of_file .and. .not. end_line )
      read(auxunit,'(a80)',iostat=ios) line
      if ( ios == 0 .and. line(1:3) .eq. 'END')then
        end_line = .true.
      elseif ( ios == 0 )then
        scan_res = index(line(1:80),tmp_name(1:Len_Trim(tmp_name)))
        if (scan_res > 0 )then
          read(auxunit,'(i5)')nshell
          atom_basis%num_basis_prims(jbasis) = 0
          do n = 1,nshell
            read(auxunit,*)nprimshell,ordprim,contract
            nindex = nindex + 1
            if ( ordprim == 0 )then
              order = 1
            elseif ( ordprim == 1 )then
              order = 3
            elseif ( ordprim == 2 )then 
              order = 6
            elseif ( ordprim == 3 )then
              order = 10
            else
              write(6,*)'bad ordprim value'
              stop
            endif
            atom_basis%num_basis_prims(jbasis) = &
                      atom_basis%num_basis_prims(jbasis) + contract*order
            zetind = nprim + 1
            do k = 1,contract
              nprim = nprim + 1
              nprim2 = nprim2 + 1
              atom_basis%prim_order(nprim) = order
              atom_basis%index(nprim) = nindex
              atom_basis%orb_contr(nprim) = contract
              read(auxunit,*)atom_basis%expon(nprim),&
                             atom_basis%contr(nprim)
              !print *,nprim,atom_basis%expon(nprim),&
              !     atom_basis%index(nprim),atom_basis%orb_contr(nprim),&
              !     atom_basis%prim_order(nprim)
            enddo
            call HBASE_unnorm(atom_basis%expon,atom_basis%contr,zetind,&
                              contract,ordprim)
            if (order > 1) then
               ncoeff = nprim
               temp = nprim - contract
               times = 0
               nindex = nindex + 1
               do k = 1,contract*(order-1) 
                  ncoeff = ncoeff + 1
                  times = times + 1
                  temp = temp + 1
                  atom_basis%expon(ncoeff) = &
                             atom_basis%expon(temp)
                  atom_basis%contr(ncoeff) = &
                             atom_basis%contr(temp)
                  atom_basis%orb_contr(ncoeff) = &
                             atom_basis%orb_contr(temp)
                  atom_basis%prim_order(ncoeff) = &
                             atom_basis%prim_order(temp)
                  atom_basis%index(ncoeff) = &
                             atom_basis%index(temp) + 1
                  !print *,ncoeff,atom_basis%expon(ncoeff),&
                  !     atom_basis%index(ncoeff),atom_basis%orb_contr(ncoeff),&
                  !     atom_basis%prim_order(nprim)
               enddo 
               nindex = atom_basis%index(temp) + 1
               nprim = nprim + times
            endif
          enddo !n = 1,nshell
        endif !scan_res > 0
      endif !ios == 0 .and. line(1:3) .eq. 'END'
    enddo !(ios == 0 .and. .not. end_of_file .and. .not. end_line )
  enddo ! jbasis = 1, nbasis

  if (nprim .ne. atom_basis%num_coefficients) then
     print *,'AHBASE_load_AObasis: ncoeff .ne. num_coefficients, exiting'
     stop
  endif

  atom_basis%off_basis_prims(1) = 0
  do jbasis = 2,nbasis
   atom_basis%off_basis_prims(jbasis) = atom_basis%off_basis_prims(jbasis-1)&
                                      + atom_basis%num_basis_prims(jbasis-1)
  enddo
  call File_Close(auxunit)
  return
end subroutine AHBASE_load_AObasis
!--------------------------------------------------------------
subroutine AHBASE_calcAOnorms(atom_basis,densfile)
use definition
  implicit none
  type(AO_basis)::atom_basis

  integer nprim,ncoeff,order,jcoeff,nx,ny,nz,n,k,allochk,temp,j,times,i,l,contr
  integer scan_file
  double precision pi, dfac(0:18),norm, ex, coef
  double precision zero, one, two, pifac
  double precision,allocatable::normvec(:),coefvec(:),expvec(:)
  character(len=*) densfile
  character(len=4) tmp_name
  logical gau_dens, nw_dens
  parameter(zero=0.d0, one=1.d0, two=2.d0)

! check if it's a gaussian checkpoint
  gau_dens = .false.
  scan_file = 0
  tmp_name = 'fchk'
  scan_file = index(densfile,tmp_name)
  if (scan_file .ne. 0) gau_dens = .true.
  if(gau_dens)write(6,*)'AHBASE_calcAOnorms: Gaussian fchk found, &
                         &reordering f'

! check if it's a NWChem density file
  nw_dens = .false.
  scan_file = 0
  tmp_name = 'NW'
  scan_file = index(densfile,tmp_name)
  if (scan_file .ne. 0) nw_dens = .true.
  if(nw_dens)write(6,*)'AHBASE_calcAOnorms: NWChem density found, &
                        &reordering d,f'

  pi = 4.d0*(atan(one))
  pifac = (2/pi)**0.75

  dfac(0) = one
  dfac(1) = one
  dfac(2) = two
  do n = 3, 12
     dfac(n) = dble(n)*dfac(n-2)
  enddo

  nprim = atom_basis%num_primitives
  ncoeff = atom_basis%num_coefficients

  allocate(atom_basis%norm_offset(nprim), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHBASE_calcAOnorms:could not allocate norm_offset, exiting'
     stop
  endif

  allocate(atom_basis%norm(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHBASE_calcAOnorms:could not allocate norm, exiting'
     stop
  endif
  atom_basis%norm(:)=0.d0

  atom_basis%norm_offset(1) = 0 
  do n = 2,nprim
    atom_basis%norm_offset(n) = atom_basis%norm_offset(n-1) + &
                                atom_basis%prim_order(n-1) 
  enddo
  jcoeff = 0
  times = 0
  do while (jcoeff .lt. ncoeff)
     jcoeff = jcoeff + 1
     order = atom_basis%prim_order(jcoeff)
     if (order .eq. 1) then ! only s shell
        ex = atom_basis%expon(jcoeff)
        nx = 0
        ny = 0
        nz = 0
        norm = zero
        call gto_norm(nx, ny, nz, ex, norm, pifac, dfac)
        atom_basis%norm(jcoeff) = norm
        atom_basis%contr(jcoeff) = atom_basis%contr(jcoeff)*norm
     else if (order .eq. 3) then ! p orbital
        ! p normalization 
        times = atom_basis%orb_contr(jcoeff)
        do k = 1, 3*times
           ex = atom_basis%expon(jcoeff)
           nx = 1
           ny = 0
           nz = 0
           norm = zero
           call gto_norm(nx, ny, nz, ex, norm, pifac, dfac)
           atom_basis%norm(jcoeff) = norm
           atom_basis%contr(jcoeff) = atom_basis%contr(jcoeff)*norm
           jcoeff = jcoeff + 1
        enddo
        jcoeff = jcoeff - 1
     else if (order .eq. 6) then ! d orbital
        ! d normalization 
        times = atom_basis%orb_contr(jcoeff)
        if (.not. nw_dens) then ! gaussian or cadpac or hondo density
           do k = 1, 3*times
              ex = atom_basis%expon(jcoeff)
              nx = 2
              ny = 0
              nz = 0
              call gto_norm(nx, ny, nz, ex, norm, pifac, dfac)
              atom_basis%norm(jcoeff) = norm
              atom_basis%contr(jcoeff) = atom_basis%contr(jcoeff)*norm
              jcoeff = jcoeff + 1
           enddo
           do k = 1, 3*times
              ex = atom_basis%expon(jcoeff)
              nx = 1
              ny = 1
              nz = 0
              call gto_norm(nx, ny, nz, ex, norm, pifac, dfac)
              atom_basis%norm(jcoeff) = norm
              atom_basis%contr(jcoeff) = atom_basis%contr(jcoeff)*norm
              jcoeff = jcoeff + 1
           enddo
        else ! nwchem density
           ! dxx 
!##########
!########## OJO, LO DE ABAJO ESTA MAL!!, DEBERIA SER nx=1,ny=1 LO CAMBIE
!########## PARA OBTENER MISMOS NUMEROS QUE NWCHEM!!!
!##########
           !do k = 1, times
           do k = 1, 6*times
              ex = atom_basis%expon(jcoeff)
              nx = 2
              ny = 0
              nz = 0
              call gto_norm(nx, ny, nz, ex, norm, pifac, dfac)
              atom_basis%norm(jcoeff) = norm
              atom_basis%contr(jcoeff) = atom_basis%contr(jcoeff)*norm
              jcoeff = jcoeff + 1
           enddo
        endif
        jcoeff = jcoeff - 1
     else if (order .eq. 10) then ! f orbital
        ! f normalization 
        times = atom_basis%orb_contr(jcoeff)
        if (.not. nw_dens) then ! gaussian or cadpac or hondo density
           do k = 1, 3*times
              ex = atom_basis%expon(jcoeff)
              nx = 3
              ny = 0
              nz = 0
              call gto_norm(nx, ny, nz, ex, norm, pifac, dfac)
              atom_basis%norm(jcoeff) = norm
              atom_basis%contr(jcoeff) = atom_basis%contr(jcoeff)*norm
              jcoeff = jcoeff + 1
           enddo
           do k = 1, 6*times
              ex = atom_basis%expon(jcoeff)
              nx = 2
              ny = 1
              nz = 0
              call gto_norm(nx, ny, nz, ex, norm, pifac, dfac)
              atom_basis%norm(jcoeff) = norm
              atom_basis%contr(jcoeff) = atom_basis%contr(jcoeff)*norm
              jcoeff = jcoeff + 1
           enddo
           do k = 1,times
              ex = atom_basis%expon(jcoeff)
              nx = 1
              ny = 1
              nz = 1
              call gto_norm(nx, ny, nz, ex, norm, pifac, dfac)
              atom_basis%norm(jcoeff) = norm
              atom_basis%contr(jcoeff) = atom_basis%contr(jcoeff)*norm
              jcoeff = jcoeff + 1
           enddo
        else
           ! fxxx
           do k = 1, 10*times
           !do k = 1, times
              ex = atom_basis%expon(jcoeff)
              nx = 3
              ny = 0
              nz = 0
              call gto_norm(nx, ny, nz, ex, norm, pifac, dfac)
              atom_basis%norm(jcoeff) = norm
              atom_basis%contr(jcoeff) = atom_basis%contr(jcoeff)*norm
              jcoeff = jcoeff + 1
           enddo
        endif
        jcoeff = jcoeff - 1
     endif ! order .eq. 1) 
  enddo ! do while
  atom_basis%norm(:) =1.d0
! now normalize contractions
  jcoeff = 0
  times = 0
  do while (jcoeff .lt. ncoeff)
     jcoeff = jcoeff + 1
     order = atom_basis%prim_order(jcoeff)
     contr = atom_basis%orb_contr(jcoeff)
     if (contr .gt. 1) then
        allocate(expvec(contr),normvec(contr),coefvec(contr), stat=allochk)
        if (allochk .gt.0) then
           write(6,*)'AHBASE_calcAOnorms: could not allocate vectors for&
                      & contraction normalization, exiting'
           stop
        endif
        do k = 1, contr
           expvec(k) = atom_basis%expon(jcoeff)
           !normvec(k) = atom_basis%norm(jcoeff)
           coefvec(k) = atom_basis%contr(jcoeff)
           jcoeff = jcoeff + 1
        enddo
        jcoeff = jcoeff - contr
        if (order .eq. 1) then ! only s shell
           l = 0
        else if (order .eq. 3) then ! p orbital
           l = 1
        else if (order .eq. 6) then ! d orbital
           l = 2
        else if (order .eq. 10) then ! f orbital
           l = 3
        endif ! order .eq. 1) 
        !call contr_norm(expvec, normvec, coefvec, l, contr)
        call contr_norm(expvec, coefvec, l, contr)
        do j = 1, order
           do k = 1, contr
              !atom_basis%norm(jcoeff) = normvec(k)
              atom_basis%contr(jcoeff) = coefvec(k)
           !print *,atom_basis%expon(jcoeff),expvec(k),atom_basis%norm(jcoeff)
              jcoeff = jcoeff+1
           enddo
        enddo
        jcoeff = jcoeff-1
        deallocate(expvec,normvec,coefvec)
     else  ! no contraction
        jcoeff = jcoeff+order-1
     endif
  enddo ! do while

  !do jcoeff = 1,ncoeff
  !  write(6,*)'norm = ',atom_basis%expon(jcoeff),atom_basis%norm(jcoeff),&
  !                      atom_basis%contr(jcoeff),atom_basis%index(jcoeff)
  !enddo

  return
end subroutine AHBASE_calcAOnorms
!-------------------------------------------------------------------------
      subroutine gto_norm(nx, ny, nz, ex, norm, pifac, dfac)
!
! 08/04 GAC: subroutine to normalize GTO functions (based on
!            whitepages from Fermann & Valeev).
!
      implicit none
      
      integer nprim, nx, ny, nz, dim, prim
      double precision ex, norm, pifac, dfac(0:18)
      integer i, n
      double precision pi, post, one, zero, sum_exp, double_sum_exp, tmp
      double precision denom_x, denom_y, denom_z
      parameter(one=1.d0, zero=0.d0)

      
      norm = zero  
      sum_exp = 2.d0**(nx + ny + nz)
      tmp = (((2.d0*nx)+(2.d0*ny)+(2.d0*nz))+3.d0)/4.d0
      double_sum_exp = (ex**(tmp))
      if ((2*nx)-1 .eq. -1) then
         denom_x = 1.d0
      else
         denom_x = dfac(((2*nx)-1))
      endif
      if ((2*ny)-1 .eq. -1) then
         denom_y = 1.d0
      else
         denom_y = dfac(((2*ny)-1))
      endif
      if ((2*nz)-1 .eq. -1) then
         denom_z = 1.d0
      else
         denom_z = dfac(((2*nz)-1))
      endif
      norm = (pifac*sum_exp*double_sum_exp)/&
             (dsqrt(denom_x*denom_y*denom_z))

   end subroutine gto_norm
!-------------------------------------------------------------------------
      !subroutine contr_norm(expvec, normvec, coefvec, l, contr)
      subroutine contr_norm(expvec, coefvec, l, contr)
!
! 08/04 GAC: subroutine to normalize contracted GTO functions (based on
!            whitepages from Fermann & Valeev).
!
      implicit none
      
      integer nprim, nx, ny, nz, dim, prim, i, n1, n2, llim, l, contr
      double precision pi, post, one, zero, sum_exp, double_sum_exp, tmp
      double precision denom_x, denom_y, denom_z, pi32
      double precision ex, norm, S, dfac
      double precision::expvec(contr), normvec(contr), coefvec(contr)
      parameter(one=1.d0, zero=0.d0)

      pi = 4.d0*(atan(one))
      pi32 = pi*sqrt(pi)

      norm = zero  
      do n1 = 1, contr     
         do n2 = 1, contr     
            S = pi32/(expvec(n1)+expvec(n2))**(1.5d0+l)/2d0**l
            !S = pi32/((2.d0**l)*expvec(n1)+expvec(n2))**1.5d0+l)
            norm = norm + coefvec(n1)*coefvec(n2)*S
         enddo
      enddo

      dfac = 1.d0
      llim = 2*l-1
      do n1 = llim,2,-2
         dfac=dfac*dble(n1)
      enddo
      norm = norm*dfac

      dfac = 1.d0/sqrt(norm)
      do n1 = 1,contr
         coefvec(n1) = coefvec(n1)*dfac
      enddo
         
   end subroutine contr_norm
!-------------------------------------------------------------------------
      subroutine HBASE_unnorm(zet,coef,first,n,l_max)

      use definition
      implicit none
      integer natoms, i, j, n, first, last, l_max
      double precision zet(*),coef(*),overlap(n,n),norm
      double precision tmp_s, tmp_p, tmp_d, tmp_f
      double precision fac, fac_s, fac_p, fac_d, fac_f, two_term

      double precision term,pi
      pi = 4.d0*(atan(1.d0))

      last = first + n - 1
!
! get overlap
!
      do i = first,last
        do j = first,last
          term = zet(i)+zet(j)
          two_term = term + term
          overlap(i-first+1,j-first+1) = (pi/term)*sqrt(pi)/sqrt(term) 
          if (l_max .eq. 1) then
!             print *,'p orbital found'
             overlap(i-first+1,j-first+1) = overlap(i-first+1,j-first+1)/&
                                            two_term
          elseif (l_max .eq. 2) then
!             print *,'d orbital found'
             overlap(i-first+1,j-first+1) = overlap(i-first+1,j-first+1)*&
                                            (3.d0/(two_term*two_term))
          elseif (l_max .eq. 3) then
             overlap(i-first+1,j-first+1) = overlap(i-first+1,j-first+1)*&
                                            (15.d0/&
                                            (two_term*two_term*two_term))
          endif
        enddo
      enddo
!
! calculate normalization
!
      norm = 0.d0
      do i = first,last
        do j = first,last
          norm = norm + coef(i)*coef(j)*overlap(i-first+1,j-first+1)/&
                     dsqrt(overlap(i-first+1,i-first+1)* &
                          overlap(j-first+1,j-first+1))
        enddo
      enddo
      norm = sqrt(norm)

!      print *,'norm = ',norm

      do i = first,last
         coef(i) = coef(i)/norm
      enddo
      
      end subroutine HBASE_unnorm
!-------------------------------------------------------------------------
