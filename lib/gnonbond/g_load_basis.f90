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
subroutine GBASE_load(p_basis,filename)
  implicit none
  integer p_basis
  character(len=*) filename
# include "io.fh"
# include "database.fh"

  integer inunit, ios
  integer nbasis,nprim,nshell,nprimshell,ordprim,order,j,k,n,nline,ncoeff
  integer p_numprim,p_primorder,ptr,p_loc_herm_coeff,p_primexpo,  &
          p_offprim,p_offloc_herm,&
          p_norm,p_offnorm !GAC
  character(len=120) line,word
  integer TRDPRM_get_next_line,lenword,TRDPRM_numparms
  double precision rnum
  integer naux_type,p_aux_type,st,ind

  inunit = File_Open(filename,"r")
  if ( inunit < 0 )then
    write(6,*)'failed to open parmfile ',filename(1:Len_Trim(filename))
    stop
  endif
  naux_type = TRDPRM_numparms(inunit,'auxtype')
  ptr = DB_alloc(p_basis,DB_int,"size_aux_type_table",1)
  IH(ptr) = naux_type
  p_aux_type = DB_alloc(p_basis,DB_int,"aux_type_table",2*naux_type)
  ! fill aux type table
  ios = 0
  naux_type = 0
  do while (ios==0)
    line = ' '  ! to avoid junk in line after some reads
    ios = TRDPRM_get_next_line(inunit,line)
    if ( ios == 0 )then
      ptr = 1 ! first call to str_next_token needs an initial value
      call str_next_token(line,ptr,word,lenword) !get 1st token
      if ( ptr > 0 )then !check if its 'auxbasis'
        if ( word(1:lenword) == 'auxtype' )then
           naux_type = naux_type + 1
           read(line(ptr:),*)IH(p_aux_type+2*(naux_type-1)), &
                             IH(p_aux_type+2*(naux_type-1)+1)
        endif
      endif
    endif
  enddo
  ncoeff = 0
  nbasis = 0
  nprim = 0
  ios = 0
  rewind(inunit)
  ! first pass count quantities
  do while (ios==0)
    line = ' '  ! to avoid junk in line after some reads
    ios = TRDPRM_get_next_line(inunit,line)
    if ( ios == 0 )then
      ptr = 1 ! first call to str_next_token needs an initial value
      call str_next_token(line,ptr,word,lenword) !get 1st token
      if ( ptr > 0 )then !check if its 'auxbasis'
        if ( word(1:lenword) == 'auxbasis' )then
           nbasis = nbasis + 1
           !call GBASE_getline(line,inunit)
           read(inunit,*)nshell
           do n = 1,nshell
             !call GBASE_getline(line,inunit)
             read(inunit,*)nprimshell,ordprim
             if ( ordprim == 0 )then
               order = 1
             elseif ( ordprim == 1 )then
               order = 4
             elseif ( ordprim == 2 )then
               order = 10
             else
               write(6,*)'bad ordprim value'
               !write(6,*)line
               stop
             endif
             nprim = nprim + nprimshell
             ncoeff = ncoeff + nprimshell*order
             nline = nprimshell*(1+order)
             do k = 1,nline
               !call GBASE_getline(line,inunit)
               read(inunit,*)rnum
             enddo
           enddo
        endif
      endif
    endif
  enddo
  ! next allocate
  write(6,*)'nbasis = ',nbasis
  write(6,*)'num prims = ',nprim
  write(6,*)'ncoeffs = ',ncoeff
  ptr = DB_alloc(p_basis,DB_int,"num_auxbasis",1)
  IH(ptr) = nbasis
  ptr = DB_alloc(p_basis,DB_int,"num_primitives",1)
  IH(ptr) = nprim
  ptr = DB_alloc(p_basis,DB_int,"num_coefficients",1)
  IH(ptr) = ncoeff
  p_numprim = DB_alloc(p_basis,DB_int,"num_basis_prims",nbasis)
  p_offprim = DB_alloc(p_basis,DB_int,"off_basis_prims",nbasis)
  p_primorder = DB_alloc(p_basis,DB_int,"prim_order",nprim)
  p_offnorm = DB_alloc(p_basis,DB_int,"norm_offset",nprim)  !GAC
  p_offloc_herm = DB_alloc(p_basis,DB_int,"local_hermite_coeff_offset",nprim)
  p_primexpo = DB_alloc(p_basis,DB_real,"prim_expon_term",nprim)
  p_loc_herm_coeff = DB_alloc(p_basis,DB_real,  &
                      "local_hermite_coefficients",ncoeff)
  p_norm = DB_alloc(p_basis,DB_real,"prim_normalizers",ncoeff)  !GAC
  rewind(inunit)
  ios = 0
  nbasis = 0
  nprim = 0
  ncoeff = 0
  do while (ios==0)
    line = ' '  ! to avoid junk in line after some reads
    ios = TRDPRM_get_next_line(inunit,line)
    if ( ios == 0 )then
      ptr = 1 ! first call to str_next_token needs an initial value
      call str_next_token(line,ptr,word,lenword) !get 1st token
      if ( ptr > 0 )then !check if its 'auxbasis'
        if ( word(1:lenword) == 'auxbasis' )then
           nbasis = nbasis + 1
           IH(p_numprim+nbasis-1) = 0
           read(inunit,*)nshell
           do n = 1,nshell
             read(inunit,*)nprimshell,ordprim
             IH(p_numprim+nbasis-1) = IH(p_numprim+nbasis-1) + nprimshell
             if ( ordprim == 0 )then
               order = 1
             elseif ( ordprim == 1 )then
               order = 4
             elseif ( ordprim == 2 )then
               order = 10
             else
               write(6,*)'bad ordprim value'
               !write(6,*)line
               stop
             endif
             do k = 1,nprimshell
               nprim = nprim + 1
               IH(p_primorder+nprim-1) = order
               read(inunit,*)RH(p_primexpo+nprim-1)
               do j = 1,order
                 ncoeff = ncoeff+1
                 read(inunit,*)RH(p_loc_herm_coeff+ncoeff-1)
                 RH(p_norm+ncoeff-1)=RH(p_loc_herm_coeff+ncoeff-1)
               enddo
             enddo
           enddo
        endif
      endif
    endif
  enddo
  ! do offsets
  IH(p_offprim) = 0
  do k = 2,nbasis
    IH(p_offprim+k-1) = IH(p_offprim+(k-1)-1) + IH(p_numprim+(k-1)-1)
  enddo
  IH(p_offloc_herm) = 0
  do k = 2,nprim
    IH(p_offloc_herm+k-1) = IH(p_offloc_herm+(k-1)-1) + IH(p_primorder+(k-1)-1)
  enddo
  write(6,*)'num prims: ',(IH(p_numprim+k-1),k=1,nbasis)
  !write(6,*)'exponent terms'
  !do k = 1,nprim
    !write(6,*)'k,coeff = ',k,RH(p_primexpo+k-1)
  !enddo
  !do k = 1,ncoeff
    !write(6,*)'k,herm coeff = ',k,RH(p_loc_herm_coeff+k-1)
  !enddo
call File_close(inunit)

  return
end subroutine GBASE_load
!--------------------------------------------------------------
subroutine GBASE_getline(line,inunit)
  implicit none
  character(len=*) line
  integer inunit
  integer ios
  integer TRDPRM_get_next_line
# include "io.fh"
# include "database.fh"
  line = ' '  ! to avoid junk in line after some reads
  ios = TRDPRM_get_next_line(inunit,line)
  write(6,*)'ios = ',ios
  write(6,*)'line = ',line
  if ( ios == 0 )then
     return
  else
     write(6,*)'parse error'
     stop
  endif
end subroutine GBASE_getline
