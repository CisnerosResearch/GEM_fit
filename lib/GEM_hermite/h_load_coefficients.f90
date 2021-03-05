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
subroutine AHCOEFF_load_cartesian(site_info,coeff_file,readdens,densfile)
use definition
  implicit none
  type(sites_info)::site_info
  character(len=*) coeff_file
  include "io.fh"
  integer coeff_unit,ptr,ncoeff,p_norm,k,n,nprim,coff,order,num_coefficients
  integer scan_file
  double precision alpha,pi,factor
  character(len=*) densfile
  character(len=4) tmp_name
  logical nw_dens, readdens

! check if it's a NWChem density file
  nw_dens = .false.
  scan_file = 0
  tmp_name = 'NW'
  scan_file = index(densfile,tmp_name)
  if (scan_file .ne. 0) nw_dens = .true.
  if(nw_dens)write(6,*)'AHCOEFF_load_cartesian: NWChem density found, &
                        &reordering d,f'

  pi = 4.d0*(atan(1.d0))

  ncoeff = site_info%num_coefficients
  if (.not. readdens) then
     coeff_unit = File_Open(coeff_file,"r")
     if ( coeff_unit < 0 )then
       write(6,*)'failed to open coeff_file ',coeff_file(1:Len_Trim(coeff_file))
       stop
     endif
     read(coeff_unit,*)num_coefficients
   
     if ( num_coefficients /= ncoeff )then
       write(6,*)'AHCOEFF_load_cartesian: wrong number of coeffs'
       stop
     endif

     do n = 1,ncoeff
       read(coeff_unit,*)site_info%cartesian_coeffs(n)
       !write(9,5004)site_info%cartesian_coeffs(n)
       ! normalize
       site_info%cartesian_coeffs(n) = site_info%cartesian_coeffs(n) * &
                                       site_info%cart_coeff_norms(n)
       !write(9,5005)site_info%cartesian_coeffs(n)
     enddo
  else
     do n = 1,ncoeff
       ! normalize
       site_info%cartesian_coeffs(n) = site_info%cartesian_coeffs(n) * &
                                       site_info%cart_coeff_norms(n)
       !write(9,5005)site_info%cartesian_coeffs(n)
     enddo
  endif
!5004    format (T2,'coeff1 = ',f20.10)
!5005    format (T2,'coeff2 = ',f20.10)
  nprim = site_info%tot_num_primitives

  coff = 1
  ! cart to Global Hermite
  do n = 1,nprim
    order = site_info%prim_order(n)
    alpha = site_info%prim_expo(n)

    factor = (pi/alpha)*sqrt(pi/alpha)
    if ( order == 1 )then  ! no change: cartesian = hermite
      site_info%global_hermite_coeffs(coff) = &
                site_info%cartesian_coeffs(coff)
      ! further normalize the charge dist so corresponds to H0 integrate to 1
      site_info%global_hermite_coeffs(coff) = &
                site_info%cartesian_coeffs(coff) * factor
    else if ( order == 3 )then  ! GAC : p shell
      do k = 1,3
        site_info%global_hermite_coeffs(coff+k-1) = &
                  site_info%cartesian_coeffs(coff+k-1) / (2.d0*alpha)
        site_info%global_hermite_coeffs(coff+k-1) = &
                   site_info%global_hermite_coeffs(coff+k-1) * factor
      enddo
    else if ( order == 4 )then  ! GAC : p shell
      site_info%global_hermite_coeffs(coff) = site_info%cartesian_coeffs(coff)
      do k = 2,4
        site_info%global_hermite_coeffs(coff+k-1) = &
                  site_info%cartesian_coeffs(coff+k-1) / (2.d0*alpha)
      enddo
      do k = 1,4
        site_info%global_hermite_coeffs(coff+k-1) = &
                   site_info%global_hermite_coeffs(coff+k-1) * factor
      enddo
    elseif ( order == 6 )then ! GAC : d shell
      do k = 1,6
        site_info%global_hermite_coeffs(coff+k-1) = &
                  site_info%cartesian_coeffs(coff+k-1) / (2.d0*alpha)**2
      enddo
      ! further normalize the charge dist so corresponds to H0 integrate to 1
      do k = 1,6
        site_info%global_hermite_coeffs(coff+k-1) = &
                   site_info%global_hermite_coeffs(coff+k-1) * factor
      enddo
    elseif ( order == 10 )then
      if (.not. nw_dens)then
         site_info%global_hermite_coeffs(coff)=&
             site_info%cartesian_coeffs(coff) &
            & + (site_info%cartesian_coeffs(coff+4) + &
            &    site_info%cartesian_coeffs(coff+5) + &
            &    site_info%cartesian_coeffs(coff+6)) / (2.d0*alpha)
      else
         site_info%global_hermite_coeffs(coff)=&
             site_info%cartesian_coeffs(coff) &
            & + (site_info%cartesian_coeffs(coff+4) + &
            &    site_info%cartesian_coeffs(coff+7) + &
            &    site_info%cartesian_coeffs(coff+9)) / (2.d0*alpha)
      endif
      do k = 2,4
        site_info%global_hermite_coeffs(coff+k-1) = &
                       site_info%cartesian_coeffs(coff+k-1) / (2.d0*alpha)
      enddo
      do k = 5,10
        site_info%global_hermite_coeffs(coff+k-1) = &
                       site_info%cartesian_coeffs(coff+k-1) / (2.d0*alpha)**2
      enddo
      ! further normalize the charge dist so corresponds to H0 integrate to 1
      do k = 1,10
        site_info%global_hermite_coeffs(coff+k-1) = &
                   site_info%global_hermite_coeffs(coff+k-1) * factor
      enddo
    elseif ( order == 20 )then ! GAC : d shell
      if (.not. nw_dens)then
         site_info%global_hermite_coeffs(coff)=&
             site_info%cartesian_coeffs(coff) &
            & + (site_info%cartesian_coeffs(coff+4) + &
            &    site_info%cartesian_coeffs(coff+5) + &
            &    site_info%cartesian_coeffs(coff+6)) / (2.d0*alpha)
      else
         site_info%global_hermite_coeffs(coff)=&
             site_info%cartesian_coeffs(coff) &
            & + (site_info%cartesian_coeffs(coff+4) + &
            &    site_info%cartesian_coeffs(coff+7) + &
            &    site_info%cartesian_coeffs(coff+9)) / (2.d0*alpha)
      endif
      do k = 2,4
        site_info%global_hermite_coeffs(coff+k-1) = &
                       site_info%cartesian_coeffs(coff+k-1) / (2.d0*alpha)
      enddo
      do k = 5,10
        site_info%global_hermite_coeffs(coff+k-1) = &
                       site_info%cartesian_coeffs(coff+k-1) / (2.d0*alpha)**2
      enddo
      do k = 11,20
        site_info%global_hermite_coeffs(coff+k-1) = &
                       site_info%cartesian_coeffs(coff+k-1) / (2.d0*alpha)
      enddo
      ! further normalize the charge dist so corresponds to H0 integrate to 1
      do k = 1,20
        site_info%global_hermite_coeffs(coff+k-1) = &
                   site_info%global_hermite_coeffs(coff+k-1) * factor
      enddo
    endif
    !if (order .le. 10) then
       coff = coff + order
    !else
    !   coff = coff + 10
    !endif
  enddo
  do n = 1,ncoeff
    site_info%cartesian_coeffs(n) = site_info%cartesian_coeffs(n) / &
                                    site_info%cart_coeff_norms(n)
  enddo
  !do n = 1,ncoeff
  !  write(8,5006)n,site_info%global_hermite_coeffs(n),&
  !                        site_info%cartesian_coeffs(n)
  !enddo
5006    format (T2,'n,GH = ',i3,2f23.16)

  return
end subroutine AHCOEFF_load_cartesian
!---------------------------------------------------------
subroutine AHCOEFF_golbalH_to_cart(site_info,densfile)
use definition
  implicit none
  type(sites_info)::site_info
  include "io.fh"
  integer coeff_unit,ptr,ncoeff,p_norm,k,n,nprim,coff,order,num,l,off
  integer scan_file
  double precision alpha,pi,factor
  character(len=*) densfile
  character(len=4) tmp_name
  logical nw_dens

! check if it's a NWChem density file
  nw_dens = .false.
  scan_file = 0
  tmp_name = 'NW'
  scan_file = index(densfile,tmp_name)
  if (scan_file .ne. 0) nw_dens = .true.
  if(nw_dens)write(6,*)'AHCOEFF_globalH_to_cart: NWChem density found, &
                        &reordering d,f'

  pi = 4.d0*(atan(1.d0))

  ncoeff = site_info%num_coefficients

  ! unnormalize
  !do n = 1,ncoeff
  !  site_info%global_hermite_coeffs(n) = site_info%global_hermite_coeffs(n)/&
  !                                  site_info%cart_coeff_norms(n)
  !enddo

  nprim = site_info%tot_num_primitives

  coff = 1
  ! cart to Global Hermite
  do n = 1,nprim
    order = site_info%prim_order(n)
    alpha = site_info%prim_expo(n)

    factor = (pi/alpha)*sqrt(pi/alpha)
    if ( order == 1 )then  ! no change: cartesian = hermite
      site_info%cartesian_coeffs(coff) = &
                site_info%global_hermite_coeffs(coff)
      ! further normalize the charge dist so corresponds to H0 integrate to 1
      site_info%cartesian_coeffs(coff) = &
                site_info%global_hermite_coeffs(coff) / factor
    else if ( order == 3 )then  ! GAC : p shell
      do k = 1,3
        site_info%cartesian_coeffs(coff+k-1) = &
                  site_info%global_hermite_coeffs(coff+k-1) * (2.d0*alpha)
        site_info%cartesian_coeffs(coff+k-1) = &
                   site_info%cartesian_coeffs(coff+k-1) / factor
      enddo
    else if ( order == 4 )then  ! GAC : p shell
      site_info%cartesian_coeffs(coff) = site_info%global_hermite_coeffs(coff)
      do k = 2,4
        site_info%cartesian_coeffs(coff+k-1) = &
                  site_info%global_hermite_coeffs(coff+k-1) * (2.d0*alpha)
      enddo
      do k = 1,4
        site_info%cartesian_coeffs(coff+k-1) = &
                   site_info%cartesian_coeffs(coff+k-1) / factor
      enddo
    elseif ( order == 6 )then ! GAC : d shell
      do k = 1,6
        site_info%cartesian_coeffs(coff+k-1) = &
                  site_info%global_hermite_coeffs(coff+k-1) * (2.d0*alpha)**2
      enddo
      ! further normalize the charge dist so corresponds to H0 integrate to 1
      do k = 1,6
        site_info%cartesian_coeffs(coff+k-1) = &
                   site_info%cartesian_coeffs(coff+k-1) / factor
      enddo
    elseif ( order == 10 )then
      if (.not. nw_dens)then
         site_info%cartesian_coeffs(coff)=&
             site_info%global_hermite_coeffs(coff) &
            & - (site_info%global_hermite_coeffs(coff+4) + &
            &    site_info%global_hermite_coeffs(coff+5) + &
            &    site_info%global_hermite_coeffs(coff+6)) * (2.d0*alpha)
      else
         site_info%cartesian_coeffs(coff)=&
             site_info%global_hermite_coeffs(coff) &
            & - (site_info%global_hermite_coeffs(coff+4) + &
            &    site_info%global_hermite_coeffs(coff+7) + &
            &    site_info%global_hermite_coeffs(coff+9)) * (2.d0*alpha)
      endif
      do k = 2,4
        site_info%cartesian_coeffs(coff+k-1) = &
                  site_info%global_hermite_coeffs(coff+k-1) * (2.d0*alpha)
      enddo
      do k = 5,10
        site_info%cartesian_coeffs(coff+k-1) = &
                  site_info%global_hermite_coeffs(coff+k-1) * (2.d0*alpha)**2
      enddo
      do k = 1,10
        site_info%cartesian_coeffs(coff+k-1) = &
                   site_info%cartesian_coeffs(coff+k-1) / factor
      enddo
    elseif ( order == 20 )then
      if (.not. nw_dens)then
         site_info%cartesian_coeffs(coff)=&
             site_info%global_hermite_coeffs(coff) &
            & - (site_info%global_hermite_coeffs(coff+4) + &
            &    site_info%global_hermite_coeffs(coff+5) + &
            &    site_info%global_hermite_coeffs(coff+6)) * (2.d0*alpha)
      else
         site_info%cartesian_coeffs(coff)=&
             site_info%global_hermite_coeffs(coff) &
            & - (site_info%global_hermite_coeffs(coff+4) + &
            &    site_info%global_hermite_coeffs(coff+7) + &
            &    site_info%global_hermite_coeffs(coff+9)) * (2.d0*alpha)
      endif
      do k = 2,4
        site_info%cartesian_coeffs(coff+k-1) = &
                  site_info%global_hermite_coeffs(coff+k-1) * (2.d0*alpha)
      enddo
      do k = 5,10
        site_info%cartesian_coeffs(coff+k-1) = &
                  site_info%global_hermite_coeffs(coff+k-1) * (2.d0*alpha)**2
      enddo
      do k = 11,20
        site_info%cartesian_coeffs(coff+k-1) = &
                  site_info%global_hermite_coeffs(coff+k-1) * (2.d0*alpha)
      enddo
      do k = 1,20
        site_info%cartesian_coeffs(coff+k-1) = &
                   site_info%cartesian_coeffs(coff+k-1) / factor
      enddo
    endif
    coff = coff + order
  enddo

  ! take out normalization
  do n = 1,ncoeff
    site_info%cartesian_coeffs(n) = site_info%cartesian_coeffs(n)/&
                                    site_info%cart_coeff_norms(n)
  enddo

  !do n = 1,site_info%nsites
  !  num = site_info%num_primitives(n)
  !  off = site_info%off_primitives(n)
  !  do k = 1,num
  !    order = site_info%prim_order(off+k)
  !    coff = site_info%coeff_offset(off+k)
  !    do l = 1,order
  !       print *,site_info%cartesian_coeffs(coff+l)
  !    enddo
  !  enddo
  !enddo

  return
end subroutine AHCOEFF_golbalH_to_cart
!---------------------------------------------------------
subroutine AHCOEFF_copy_norms_from_basis(site_info,atom_auxis)
use definition
  implicit none
  type(sites_info)::site_info
  type(auxiliary_basis)::atom_auxis
  integer j,k,n,nsites,off,num,ind_basis,offb,order,coff,bcoff
  integer allochk, ncoeff

  nsites = site_info%nsites

  ncoeff = site_info%num_coefficients

  ! allocate
  allocate(site_info%cart_coeff_norms(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHCOEFF_copy_norms_from_basis: could not allocate&
              & cart_coeff_norms, exiting'
     stop
  endif

  allocate(site_info%aux_elec(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHCOEFF_copy_norms_from_basis: could not allocate&
              & aux_elec, exiting'
     stop
  endif

  !print *,'size = ',size(atom_auxis%prim_normalizers)
  !print *,'size = ',size(site_info%cart_coeff_norms)
  !print *,'printing norms neaux en loadcoef'
  do n = 1,nsites
    off = site_info%off_primitives(n)
    num = site_info%num_primitives(n)
    ind_basis = site_info%basis_index(n)
    offb = atom_auxis%off_basis_prims(ind_basis)
    !print *,off,num,ind_basis,offb
    do k = 1,num
      order = site_info%prim_order(off+k)
      !print *,'inside k ',k,num,order
      !if(order.gt.10)order=10
      if(order.gt.20)then
        print *, 'Cannot handle sites with prims order>20'
        stop
      endif
      coff = site_info%coeff_offset(off+k)
      bcoff = atom_auxis%norm_offset(offb+k)
      do j = 1,order
        site_info%cart_coeff_norms(coff+j) = &
                  atom_auxis%prim_normalizers(bcoff+j)
        site_info%aux_elec(coff+j) = &
                  atom_auxis%prim_elec(bcoff+j)
        !print *,'norm = ',site_info%cart_coeff_norms(coff+j)
        !print *,'elec1 = ',site_info%aux_elec(coff+j)
         !print *,site_info%cart_coeff_norms(coff+j),&
         !        site_info%aux_elec(coff+j)
      enddo
      !print *,'finished order ',order
    enddo
    !print *,'finished num ',num
  enddo

end subroutine AHCOEFF_copy_norms_from_basis
!---------------------------------------------------------
subroutine AHCOEFF_load_local_hermite(site_info,lherm_file)
use definition
  implicit none
  type(sites_info)::site_info
  character(len=*) lherm_file
  include "io.fh"

  integer lherm_unit,ptr,ncoeff,p_lhermite,num_coefficients,jcoeff
  integer jres,numres,k

  lherm_unit = File_Open(lherm_file,"r")
  if ( lherm_unit < 0 )then
    write(6,*)'failed to open lherm_file ',lherm_file(1:Len_Trim(lherm_file))
    stop
  endif
  read(lherm_unit,*)num_coefficients
  ncoeff = site_info%num_coefficients

  numres = site_info%num_residues

  if ( numres*num_coefficients /= ncoeff )then
    write(6,*)'AHCOEFF_load_local_hermite: mismatch between num hermite&
              & coeffs',ncoeff,numres*num_coefficients
    stop
  endif
! HACK!! for now read the file in once for each residue
  jcoeff = 0
  do jres = 1,numres 
    rewind(lherm_unit)
    read(lherm_unit,*)num_coefficients
    do k = 1,num_coefficients
      jcoeff = jcoeff + 1 
      !read(lherm_unit,*)RH(p_lhermite+jcoeff-1)
      read(lherm_unit,*)site_info%local_hermite_coeffs(jcoeff)
    enddo
  enddo
  jcoeff = 0
  !do jres = 1, numres
  !  do k = 1,num_coefficients
  !    jcoeff = jcoeff + 1 
  !     write(6,5006)site_info%local_hermite_coeffs(jcoeff)
  !  enddo
  !enddo
!5006    format (T2,f20.17)
  call File_Close(lherm_unit)
  return
end subroutine AHCOEFF_load_local_hermite
!---------------------------------------------------------
subroutine AHCOEFF_mult_hermite(site_info,coeff_file)
use definition
  implicit none
  type(sites_info)::site_info
  character(len=*) coeff_file
  include "io.fh"
  integer coeff_unit,ptr,ncoeff,p_norm,k,n,nprim,coff,order,num_coefficients
  double precision alpha,pi,factor

  pi = 4.d0*(atan(1.d0))

  ncoeff = site_info%num_coefficients
  do n = 1,ncoeff
    site_info%cartesian_coeffs(n) = site_info%cartesian_coeffs(n) * &
                                    site_info%cart_coeff_norms(n)
    !write(9,5005)site_info%cartesian_coeffs(n)
  enddo
  nprim = site_info%tot_num_primitives

  coff = 1
  ! cart to Global Hermite
  do n = 1,nprim
    order = site_info%prim_order(n)
    alpha = site_info%prim_expo(n)
    factor = (pi/alpha)*sqrt(pi/alpha)
    if ( order == 1 )then  ! no change: cartesian = hermite
      site_info%global_hermite_coeffs(coff) = &
                site_info%cartesian_coeffs(coff) * factor
    else if ( order == 3 )then  ! GAC : p shell
      do k = 1,3
        site_info%global_hermite_coeffs(coff+k-1) = &
                   site_info%cartesian_coeffs(coff+k-1) * factor
      enddo
    else if ( order == 4 )then  ! GAC : p shell
      do k = 1,4
        site_info%global_hermite_coeffs(coff+k-1) = &
                   site_info%cartesian_coeffs(coff+k-1) * factor
      enddo
    elseif ( order == 6 )then ! GAC : d shell
      do k = 1,6
        site_info%global_hermite_coeffs(coff+k-1) = &
                   site_info%cartesian_coeffs(coff+k-1) * factor
      enddo
    elseif ( order == 10 )then
      do k = 1,10
        site_info%global_hermite_coeffs(coff+k-1) = &
                   site_info%cartesian_coeffs(coff+k-1) * factor
      enddo
    elseif ( order == 20 )then
      do k = 1,20
        site_info%global_hermite_coeffs(coff+k-1) = &
                   site_info%cartesian_coeffs(coff+k-1) * factor
      enddo
    endif
    coff = coff + order
  enddo

  do n = 1,ncoeff
    site_info%cartesian_coeffs(n) = site_info%cartesian_coeffs(n) / &
                                    site_info%cart_coeff_norms(n)
  enddo

  !do n = 1,ncoeff
  !  write(8,5006)n,site_info%global_hermite_coeffs(n),&
  !                        site_info%cartesian_coeffs(n)
  !enddo
5006    format (T2,'n,GH = ',i3,2f23.16)

  return
end subroutine AHCOEFF_mult_hermite
!---------------------------------------------------------
