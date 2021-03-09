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
!------------------------------------------------------
subroutine GHERM_load_prim(p_site_info,p_basis)
  implicit none
  integer p_site_info,p_basis
# include "io.fh"
# include "database.fh"

  integer p_stype
  integer p_numprim,p_offprim,p_primorder,p_off_LH_coeff,p_primexpo,p_LH_coeff
  integer k,m,n,nsites,ptr,stype,nbasis,nprim,ncoeff
  integer ps_numprim,ps_offprim,ps_primexpo,ps_primorder,ps_prim_off_coeff, &
          p_G_hermite_coeff,p_L_hermite_coeff
  integer ntable,p_table,ind,p_indbasis
  integer nprimtot,ncoefftot,num,off,soff,scoff,coff,order
  character(len=120) line,word
  integer inunit,ios,lenword,TRDPRM_get_next_line
  double precision alpha,pi,factor

  pi = 3.14159265358979323846d0

  ptr = DB_Get_Pointer(p_site_info,DB_int,"numsites",1)
  nsites = IH(ptr)
  p_indbasis = DB_alloc(p_site_info,DB_int,"basis_index",nsites)
  p_stype = DB_Get_Pointer(p_site_info,DB_int,"site_type",nsites)
  ptr = DB_Get_Pointer(p_basis,DB_int,"size_aux_type_table",1)
  ntable = IH(ptr)
  p_table = DB_Get_Pointer(p_basis,DB_int,"aux_type_table",2*ntable)
  
  ! first load the basis index
  do n = 1,nsites
    stype = IH(p_stype+n-1)
    ind = -1
    do k = 1,ntable
      if ( IH(p_table+2*(k-1)) == stype )then
        ind = IH(p_table+2*(k-1)+1)
      endif
    enddo
    if ( ind < 0 )then
      write(6,*)'auxiliary index not forund for site ',n
      stop
    endif
    IH(p_indbasis+n-1) = ind
  enddo
  ptr = DB_Get_Pointer(p_basis,DB_int,"num_auxbasis",1)
  nbasis = IH(ptr)
  ptr = DB_Get_Pointer(p_basis,DB_int,"num_primitives",1)
  nprim = IH(ptr)
  ptr = DB_Get_Pointer(p_basis,DB_int,"num_coefficients",1)
  ncoeff = IH(ptr)
  p_numprim = DB_Get_Pointer(p_basis,DB_int,"num_basis_prims",nbasis)
  p_offprim = DB_Get_Pointer(p_basis,DB_int,"off_basis_prims",nbasis)
  p_primorder = DB_Get_Pointer(p_basis,DB_int,"prim_order",nprim)
  p_off_LH_coeff = DB_Get_Pointer(p_basis,DB_int,  &
                "local_hermite_coeff_offset",nprim)
  p_primexpo = DB_Get_Pointer(p_basis,DB_real,"prim_expon_term",nprim)
  p_LH_coeff = DB_Get_Pointer(p_basis,DB_real, &
               "local_hermite_coefficients",ncoeff)
  write(6,*)' next get the number of prims per site...'
  ! next get the number of prims per site and total
  ps_numprim = DB_alloc(p_site_info,DB_int,"num_site_prims",nsites)
  ps_offprim = DB_alloc(p_site_info,DB_int,"off_site_prims",nsites)
  nprimtot = 0
  do n = 1,nsites
    ind = IH(p_indbasis+n-1)
    IH(ps_numprim+n-1) = IH(p_numprim+ind-1)
    nprimtot = nprimtot + IH(p_numprim+ind-1)
  enddo
  write(6,*)'nprimtot = ',nprimtot
  ! get offsets
  IH(ps_offprim) = 0
  do n = 2,nsites
    IH(ps_offprim+n-1) = IH(ps_offprim+(n-1)-1)+IH(ps_numprim+(n-1)-1)
  enddo
  ! get site prim arrays
  ptr = DB_alloc(p_site_info,DB_int,"tot_site_prims",1)
  IH(ptr) = nprimtot
  ps_primorder = DB_alloc(p_site_info,DB_int,"prim_order",nprimtot)
  ps_prim_off_coeff = DB_alloc(p_site_info,DB_int,"coeff_offset",nprimtot)
  ps_primexpo = DB_alloc(p_site_info,DB_real,"prim_expon_term",nprimtot)
  do n = 1,nsites
    ind = IH(p_indbasis+n-1)
    num = IH(p_numprim+ind-1)
    off = IH(p_offprim+ind-1)
    soff = IH(ps_offprim+n-1)
    write(6,*)'site # ind,num,off,soff = ',n,ind,IH(ps_numprim+n-1),off,soff
    do k = 1,num
      RH(ps_primexpo+soff+k-1) = RH(p_primexpo+off+k-1)
      IH(ps_primorder+soff+k-1) = IH(p_primorder+off+k-1)
    enddo
  enddo
  write(6,*)'get coeff offsets and coeff total'
  ! get coeff offsets and coeff total
  IH(ps_prim_off_coeff) = 0
  do n = 2,nprimtot
    IH(ps_prim_off_coeff+n-1) = IH(ps_prim_off_coeff+(n-1)-1) + &
                                IH(ps_primorder+(n-1)-1)
  enddo
  ncoefftot = IH(ps_prim_off_coeff+nprimtot-1)+IH(ps_primorder+nprimtot-1)
  ptr = DB_alloc(p_site_info,DB_int,"num_hermite_coefficients",1)
  IH(ptr) = ncoefftot
  !do n = 1,nprimtot
    !write(9,*)'n,order = ',n,IH(ps_primorder+n-1)
  !enddo
  write(6,*)'ncoefftot = ',ncoefftot
  p_G_hermite_coeff = DB_alloc(p_site_info,DB_real,  &
                        "global_hermite_coeffs",ncoefftot)
  p_L_hermite_coeff = DB_alloc(p_site_info,DB_real,  &
                        "local_hermite_coeffs",ncoefftot)

  ! copy the prim local hermite coefficients
  do n = 1,nsites
    ind = IH(p_indbasis+n-1)
    num = IH(p_numprim+ind-1)
    off = IH(p_offprim+ind-1)
    soff = IH(ps_offprim+n-1)
    !write(10,*)'n,num = ',n,num,' ==============================='
    do k = 1,num
      order = IH(p_primorder+off+k-1)
      scoff = IH(ps_prim_off_coeff+soff+k-1)
      coff = IH(p_off_LH_coeff+off+k-1)
      do m = 1,order
        RH(p_L_hermite_coeff+scoff+m-1) = RH(p_LH_coeff+coff+m-1)
        !write(10,*)RH(p_L_hermite_coeff+scoff+m-1)
      enddo
      !write(10,*)'+++++++++++++++++++++++++++++++++++++'
    enddo
  enddo

  return
end subroutine GHERM_load_prim
!------------------------------------------------------
subroutine GHERM_local_to_global(p_site_info)
  implicit none
  integer p_site_info
#include "database.fh"

  integer ptr,nsites,ps_numprim,ps_offprim,ps_primorder,ps_prim_off_coeff
  integer num,off,order,scoff,dimxy,j,k,m,n,nframes,p_frame,p_frameindex
  integer p_G_hermite_coeff,p_L_hermite_coeff,nprimtot,ncoefftot
  double precision At(3,3),Mpole_xy(10,10)

  write(6,*)'get ptrs'
  ptr = DB_Get_Pointer(p_site_info,DB_int,"numsites",1)
  nsites = IH(ptr)
  ps_numprim = DB_Get_Pointer(p_site_info,DB_int,"num_site_prims",nsites)
  ps_offprim = DB_Get_Pointer(p_site_info,DB_int,"off_site_prims",nsites)
  ptr = DB_Get_Pointer(p_site_info,DB_int,"nframes",1)
  nframes = IH(ptr)
  p_frame = DB_Get_Pointer(p_site_info,DB_real,"frames",3*3*nframes)
  p_frameindex = DB_Get_Pointer(p_site_info,DB_int,"frame_index",nsites)
  ptr = DB_Get_Pointer(p_site_info,DB_int,"tot_site_prims",1)
  nprimtot = IH(ptr)
  ps_primorder = DB_Get_Pointer(p_site_info,DB_int,"prim_order",nprimtot)
  ps_prim_off_coeff = DB_Get_Pointer(p_site_info,DB_int,"coeff_offset",nprimtot)
  ncoefftot = IH(ps_prim_off_coeff+nprimtot-1)+IH(ps_primorder+nprimtot-1)
  p_G_hermite_coeff = DB_Get_Pointer(p_site_info,DB_real,  &
                        "global_hermite_coeffs",ncoefftot)
  p_L_hermite_coeff = DB_Get_Pointer(p_site_info,DB_real,  &
                        "local_hermite_coeffs",ncoefftot)
  write(6,*)'before loop'
  do n = 1,nsites
    num = IH(ps_numprim+n-1)
    off = IH(ps_offprim+n-1)
    m = IH(p_frameindex+n-1)
    write(6,*)'num,off,m = ',num,off,m
    ! transpose frame for global to local
    !call VEC3D_3x3_transpose(RH(p_frame+9*(m-1)),At)
    order = 10 ! either 1 or 10---here use 10 to get Mpole_xy
    call XFORM_MPOLE_matrix(RH(p_frame+9*(m-1)),Mpole_xy,order)
    do k = 1,num
      order = IH(ps_primorder+off+k-1)
      scoff = IH(ps_prim_off_coeff+off+k-1)
      dimxy = order
      if ( order == 1 )then ! simply copy s-orbital coeff: local=global
        RH(p_G_hermite_coeff+scoff) = RH(p_L_hermite_coeff+scoff)
      elseif ( order > 1 )then
        call XFORM_MPOLE(Mpole_xy,dimxy,RH(p_L_hermite_coeff+scoff), &
             RH(p_G_hermite_coeff+scoff),order)
      endif
    enddo
  enddo
  do n = 1,nsites
    num = IH(ps_numprim+n-1)
    off = IH(ps_offprim+n-1)
    write(12,*)'n,num = ',n,num,' ==============================='
    do k = 1,num
      order = IH(ps_primorder+off+k-1)
      scoff = IH(ps_prim_off_coeff+off+k-1)
      do m = 1,order
        write(12,*)RH(p_G_hermite_coeff+scoff+m-1)
      enddo
      write(12,*)'+++++++++++++++++++++++++++++++++++++'
    enddo
  enddo
  return
end subroutine GHERM_local_to_global
!------------------------------------------------------
