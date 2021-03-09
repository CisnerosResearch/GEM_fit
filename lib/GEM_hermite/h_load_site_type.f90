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
!-------------------------------------------------------
subroutine AHTYPE_load_site_info(site_info,atom_auxis,atom_basis,parmfile,&
                                 readdens)
use definition
  implicit none
  type(AO_basis)::atom_basis
  type(sites_info)::site_info
  type(auxiliary_basis)::atom_auxis
  character(len=*) parmfile
  include "io.fh"
  include "site_names.inc"

  integer parmunit,AHTYPE_type_from_name,AHTYPE_basis_index
  integer j,k,n,nsites,ttype,ind_basis, allochk, natoms
  integer nbasis,nprim,nbprim,ncoeff
  integer nbasisao,nprimao,nbprimao,ncoeffao,naos
  double precision charge,mass
  integer off,offb,num
  integer coff,bcoff,order,nbcoeff
  integer offao,offbao,numao
  integer coffao,bcoffao,orderao,nbcoeffao
  logical readdens

  parmunit = File_Open(parmfile,"r")
  if ( parmunit < 0 )then
    write(6,*)'failed to open parmfile ',parmfile(1:Len_Trim(parmfile))
    stop
  endif
  nsites = site_info%nsites
  natoms = site_info%natoms

  allocate(site_info%site_type(nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*),'AHTYPE_load_site_info:could not allocate site_type, exiting'
     stop
  endif

  do n = 1,nsites
    ttype = AHTYPE_type_from_name(sitename(n),parmunit)
    if ( ttype > 0 )then
      site_info%site_type(n) = ttype
    else
      write(6,*)'AHTYPE_load_site_info: bad type for: ',sitename(n)
      stop
    endif
    !write(6,*)'n,name,type = ',n,sitename(n),site_info%site_type(n)
  enddo

  allocate(site_info%nuclear_charge(nsites), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate nuclear_charge, exiting'
     stop
  endif

  allocate(site_info%nuclear_mass(nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHTYPE_load_site_info:could not allocate nuclear_mass, exiting'
     stop
  endif

  do n = 1,nsites
    call AHTYPE_charge_and_mass(site_info%site_type(n),parmunit,charge,mass) 
    if ( charge >= 0 )then
      site_info%nuclear_charge(n) = charge
      site_info%nuclear_mass(n) = mass
    else
      write(6,*)'AHTYPE_load_site_info: bad charge,mass for: ',sitename(n)
      stop
    endif
  enddo
  ! get basis index
  nbasis = atom_auxis%num_auxbasis

  allocate(site_info%basis_index(nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHTYPE_load_site_info:could not allocate basis_index, exiting'
     stop
  endif

  allocate(site_info%ao_bas_ind(natoms), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHTYPE_load_site_info:could not allocate ao_bas_ind, exiting'
     stop
  endif

  do n = 1,nsites
    ind_basis =  AHTYPE_basis_index(site_info%site_type(n),parmunit)
    if ( ind_basis > 0 .and. ind_basis <= nbasis )then
      site_info%basis_index(n) = ind_basis
      if(n .le. natoms) site_info%ao_bas_ind(n) = ind_basis
    else
      write(6,*)'AHTYPE_load_site_info: bad basis index for: ',sitename(n)
      stop
    endif
  enddo
  call File_Close(parmunit)
  ! get the number of prims

  allocate(site_info%num_primitives(nsites), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate num_primitives, exiting'
     stop
  endif

  allocate(site_info%off_primitives(nsites), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate off_primitives, exiting'
     stop
  endif

  allocate(site_info%ao_num_prim(natoms), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate ao_num_prim, exiting'
     stop
  endif

  allocate(site_info%ao_off_prim(natoms), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate ao_off_prim, exiting'
     stop
  endif

  nprim = 0
  nprimao = 0
  do n = 1,nsites
    ind_basis = site_info%basis_index(n)
    site_info%num_primitives(n) = atom_auxis%num_basis_prims(ind_basis)
    nprim = nprim + atom_auxis%num_basis_prims(ind_basis)
    if(n .le. natoms) then 
       site_info%ao_num_prim(n) = atom_basis%num_basis_prims(ind_basis)
       nprimao = nprimao + atom_basis%num_basis_prims(ind_basis)
    endif
  enddo
  site_info%tot_num_primitives = nprim
  site_info%tot_ao_num_prim = nprimao
  write(6,*)'num primitives = ',nprim
  write(6,*)'num ao prim = ',nprimao

  site_info%off_primitives(1) = 0
  site_info%ao_off_prim(1) = 0
  do n = 2,nsites
    site_info%off_primitives(n) = site_info%off_primitives(n-1) + &
                                  site_info%num_primitives(n-1) 
    if(n .le. natoms) site_info%ao_off_prim(n) = &
                                  site_info%ao_off_prim(n-1) + &
                                  site_info%ao_num_prim(n-1) 
  enddo

  allocate(site_info%prim_expo(nprim), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate prim_expo, exiting'
     stop
  endif

  allocate(site_info%prim_order(nprim), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate prim_order, exiting'
     stop
  endif

  allocate(site_info%coeff_offset(nprim), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate coeff_offset, exiting'
     stop
  endif

  allocate(site_info%ao_expo(nprimao), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate ao_expo, exiting'
     stop
  endif

  allocate(site_info%ao_contr_coeff(nprimao), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate ao_contr_coeff, exiting'
     stop
  endif

  allocate(site_info%ao_contr_deg(nprimao), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate ao_contr_coeff, exiting'
     stop
  endif

  allocate(site_info%ao_norms(nprimao), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate ao_offset, exiting'
     stop
  endif

  allocate(site_info%ao_index(nprimao), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate ao_offset, exiting'
     stop
  endif

  allocate(site_info%ao_order(nprimao), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate ao_offset, exiting'
     stop
  endif

  nbprim = atom_auxis%num_primitives
  !nbprimao = atom_basis%num_primitives

  write(6,*)'num auxis = ',nbprim
  do n = 1,nsites
    off = site_info%off_primitives(n)
    num = site_info%num_primitives(n)
    ind_basis = site_info%basis_index(n)
    offb = atom_auxis%off_basis_prims(ind_basis) 
    do k = 1,num
      site_info%prim_order(off+k) = atom_auxis%prim_order(offb+k)
      site_info%prim_expo(off+k) = atom_auxis%prim_expon_term(offb+k)
    enddo
  enddo
  ncoeff = 0
  do n = 1,nprim
       ncoeff = ncoeff + site_info%prim_order(n)
  enddo
  write(6,*)'ncoeff = ',ncoeff
  site_info%num_coefficients = ncoeff

  site_info%coeff_offset(1) = 0
  site_info%coeff_offset(1) = 0
  do n = 2,nprim
       site_info%coeff_offset(n) = site_info%coeff_offset(n-1) + &
                                   site_info%prim_order(n-1)
  enddo
 

  ! for AOs; NOTE THAT IT INCLUDES NORMALIZATION BECAUSE AOs ARE UNROLLED
  naos = 0
  ncoeffao = 0
  do n = 1,natoms
    off = site_info%ao_off_prim(n)
    num = site_info%ao_num_prim(n)
    ind_basis = site_info%ao_bas_ind(n)
    offb = atom_basis%off_basis_prims(ind_basis) 
    do k = 1,num
      site_info%ao_index(off+k) = atom_basis%index(offb+k)
      site_info%ao_order(off+k) = atom_basis%prim_order(offb+k)
      site_info%ao_expo(off+k) = atom_basis%expon(offb+k)
      site_info%ao_contr_coeff(off+k) = atom_basis%contr(offb+k)
      site_info%ao_contr_deg(off+k) = atom_basis%orb_contr(offb+k)
      site_info%ao_norms(off+k) = atom_basis%norm(offb+k)
      ncoeffao = ncoeffao + 1
      !print *,'ao info ',n,site_info%ao_index(off+k),site_info%ao_expo(off+k),&
      !          site_info%ao_contr_coeff(off+k),site_info%ao_norms(off+k),&
      !          site_info%ao_contr_deg(off+k),site_info%ao_order(off+k)
    enddo
    naos = naos + site_info%ao_index(off+num)
  enddo
  write(6,*)'ncoeff AO = ',ncoeffao
  site_info%num_coefficients = ncoeff
  if (readdens) then
     if (naos .ne. site_info%num_stos) then
        write(6,*)'AHTYPE_load_site_info: naos .ne. num_stos from densmat,&
                   & exiting'
        write(6,*)naos,site_info%num_stos
        stop
     endif
  endif
  write(6,*)'num AOs = ',naos
  if(.not. readdens) site_info%num_stos = naos  ! GAC changed this for QMMM!!!

  allocate(site_info%cartesian_coeffs(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHTYPE_load_site_info:could not allocate cartesian_coeffs,&
                & exiting'
     stop
  endif
  allocate(site_info%global_hermite_coeffs(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
    write(6,*)'AHTYPE_load_site_info:could not allocate global_hermite_coeffs,&
               & exiting'
     stop
  endif
  allocate(site_info%local_hermite_coeffs(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate local_hermite_coeffs,&
              & exiting'
     stop
  endif

  return
end subroutine AHTYPE_load_site_info
!-------------------------------------------------------
subroutine AHTYPE_load_site_info_AO(site_info,atom_basis,parmfile)
use definition
  implicit none
  type(AO_basis)::atom_basis
  type(sites_info)::site_info
  character(len=*) parmfile
  include "io.fh"
  include "site_names.inc"

  integer parmunit,AHTYPE_type_from_name,AHTYPE_basis_index
  integer j,k,n,nsites,ttype,ind_basis, allochk, natoms
  integer nbasis,nprim,nbprim,ncoeff
  integer nbasisao,nprimao,nbprimao,ncoeffao,naos
  double precision charge,mass
  integer off,offb,num
  integer coff,bcoff,order,nbcoeff
  integer offao,offbao,numao
  integer coffao,bcoffao,orderao,nbcoeffao
  logical readdens

  parmunit = File_Open(parmfile,"r")
  if ( parmunit < 0 )then
    write(6,*)'failed to open parmfile ',parmfile(1:Len_Trim(parmfile))
    stop
  endif
  nsites = site_info%nsites
  natoms = site_info%natoms

  allocate(site_info%site_type(nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*),'AHTYPE_load_site_info:could not allocate site_type, exiting'
     stop
  endif

  do n = 1,nsites
    ttype = AHTYPE_type_from_name(sitename(n),parmunit)
    if ( ttype > 0 )then
      site_info%site_type(n) = ttype
    else
      write(6,*)'AHTYPE_load_site_info: bad type for: ',sitename(n)
      stop
    endif
    !write(6,*)'n,name,type = ',n,sitename(n),site_info%site_type(n)
  enddo

  allocate(site_info%nuclear_charge(nsites), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate nuclear_charge, exiting'
     stop
  endif

  allocate(site_info%nuclear_mass(nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHTYPE_load_site_info:could not allocate nuclear_mass, exiting'
     stop
  endif

  do n = 1,nsites
    call AHTYPE_charge_and_mass(site_info%site_type(n),parmunit,charge,mass) 
    if ( charge >= 0 )then
      site_info%nuclear_charge(n) = charge
      site_info%nuclear_mass(n) = mass
    else
      write(6,*)'AHTYPE_load_site_info: bad charge,mass for: ',sitename(n)
      stop
    endif
  enddo
  ! get basis index

  allocate(site_info%basis_index(nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHTYPE_load_site_info:could not allocate basis_index, exiting'
     stop
  endif

  allocate(site_info%ao_bas_ind(natoms), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHTYPE_load_site_info:could not allocate ao_bas_ind, exiting'
     stop
  endif

  do n = 1,nsites
    ind_basis =  AHTYPE_basis_index(site_info%site_type(n),parmunit)
    if ( ind_basis > 0 .and. ind_basis <= nbasis )then
      site_info%basis_index(n) = ind_basis
      if(n .le. natoms) site_info%ao_bas_ind(n) = ind_basis
    else
      write(6,*)'AHTYPE_load_site_info: bad basis index for: ',sitename(n)
      stop
    endif
  enddo
  call File_Close(parmunit)
  ! get the number of prims

   allocate(site_info%ao_num_prim(natoms), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate ao_num_prim, exiting'
     stop
  endif

  allocate(site_info%ao_off_prim(natoms), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate ao_off_prim, exiting'
     stop
  endif

  nprimao = 0
  do n = 1,nsites
    ind_basis = site_info%basis_index(n)
    if(n .le. natoms) then 
       site_info%ao_num_prim(n) = atom_basis%num_basis_prims(ind_basis)
       nprimao = nprimao + atom_basis%num_basis_prims(ind_basis)
    endif
  enddo
  site_info%tot_ao_num_prim = nprimao
  write(6,*)'num ao prim = ',nprimao

  site_info%off_primitives(1) = 0
  site_info%ao_off_prim(1) = 0
  do n = 2,nsites
    site_info%off_primitives(n) = site_info%off_primitives(n-1) + &
                                  site_info%num_primitives(n-1) 
    if(n .le. natoms) site_info%ao_off_prim(n) = &
                                  site_info%ao_off_prim(n-1) + &
                                  site_info%ao_num_prim(n-1) 
  enddo

  allocate(site_info%ao_expo(nprimao), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate ao_expo, exiting'
     stop
  endif

  allocate(site_info%ao_contr_coeff(nprimao), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate ao_contr_coeff, exiting'
     stop
  endif

  allocate(site_info%ao_contr_deg(nprimao), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate ao_contr_coeff, exiting'
     stop
  endif

  allocate(site_info%ao_norms(nprimao), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate ao_offset, exiting'
     stop
  endif

  allocate(site_info%ao_index(nprimao), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate ao_offset, exiting'
     stop
  endif

  allocate(site_info%ao_order(nprimao), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate ao_offset, exiting'
     stop
  endif


  ! for AOs; NOTE THAT IT INCLUDES NORMALIZATION BECAUSE AOs ARE UNROLLED
  naos = 0
  ncoeffao = 0
  do n = 1,natoms
    off = site_info%ao_off_prim(n)
    num = site_info%ao_num_prim(n)
    ind_basis = site_info%basis_index(n)
    offb = atom_basis%off_basis_prims(ind_basis) 
    do k = 1,num
      site_info%ao_index(off+k) = atom_basis%index(offb+k)
      site_info%ao_order(off+k) = atom_basis%prim_order(offb+k)
      site_info%ao_expo(off+k) = atom_basis%expon(offb+k)
      site_info%ao_contr_coeff(off+k) = atom_basis%contr(offb+k)
      site_info%ao_contr_deg(off+k) = atom_basis%orb_contr(offb+k)
      site_info%ao_norms(off+k) = atom_basis%norm(offb+k)
      ncoeffao = ncoeffao + 1
      !print *,'ao info ',site_info%ao_index(off+k),site_info%ao_expo(off+k),&
      !          site_info%ao_contr_coeff(off+k),site_info%ao_norms(off+k),&
      !          site_info%ao_contr_deg(off+k),site_info%ao_order(off+k)
    enddo
    naos = naos + site_info%ao_index(off+num)
  enddo
  write(6,*)'ncoeff AO = ',ncoeffao
  site_info%num_coefficients = ncoeff
  if (readdens) then
     if (naos .ne. site_info%num_stos) then
        write(6,*)'AHTYPE_load_site_info: naos .ne. num_stos from densmat,&
                   & exiting'
        write(6,*)naos,site_info%num_stos
        stop
     endif
  endif
  write(6,*)'num AOs = ',naos

  allocate(site_info%cartesian_coeffs(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHTYPE_load_site_info:could not allocate cartesian_coeffs,&
                & exiting'
     stop
  endif
  allocate(site_info%global_hermite_coeffs(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
    write(6,*)'AHTYPE_load_site_info:could not allocate global_hermite_coeffs,&
               & exiting'
     stop
  endif
  allocate(site_info%local_hermite_coeffs(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate local_hermite_coeffs,&
              & exiting'
     stop
  endif

  return
end subroutine AHTYPE_load_site_info_AO
!-------------------------------------------------------
subroutine AHTYPE_load_site_info_2(site_info,atom_auxis,parmfile)
use definition
  implicit none
  type(sites_info)::site_info
  type(auxiliary_basis)::atom_auxis
  character(len=*) parmfile
  include "io.fh"
  include "site_names.inc"

  integer parmunit,AHTYPE_type_from_name,AHTYPE_basis_index
  integer j,k,n,nsites,ttype,ind_basis, allochk
  integer nbasis,nprim,nbprim,ncoeff
  double precision charge,mass
  integer off,offb,num
  integer coff,bcoff,order,nbcoeff
  logical readdens

  parmunit = File_Open(parmfile,"r")
  if ( parmunit < 0 )then
    write(6,*)'failed to open parmfile ',parmfile(1:Len_Trim(parmfile))
    stop
  endif
  nsites = site_info%nsites

  allocate(site_info%site_type(nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*),'AHTYPE_load_site_info:could not allocate site_type, exiting'
     stop
  endif

  do n = 1,nsites
    ttype = AHTYPE_type_from_name(sitename(n),parmunit)
    if ( ttype > 0 )then
      site_info%site_type(n) = ttype
    else
      write(6,*)'AHTYPE_load_site_info: bad type for: ',sitename(n)
      stop
    endif
    !write(6,*)'n,name,type = ',n,sitename(n),site_info%site_type(n)
  enddo

  allocate(site_info%nuclear_charge(nsites), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate nuclear_charge, exiting'
     stop
  endif

  allocate(site_info%nuclear_mass(nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHTYPE_load_site_info:could not allocate nuclear_mass, exiting'
     stop
  endif

  do n = 1,nsites
    call AHTYPE_charge_and_mass(site_info%site_type(n),parmunit,charge,mass) 
    if ( charge >= 0 )then
      site_info%nuclear_charge(n) = charge
      site_info%nuclear_mass(n) = mass
    else
      write(6,*)'AHTYPE_load_site_info: bad charge,mass for: ',sitename(n)
      stop
    endif
  enddo
  ! get basis index
  nbasis = atom_auxis%num_auxbasis

  allocate(site_info%basis_index(nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHTYPE_load_site_info:could not allocate basis_index, exiting'
     stop
  endif

  do n = 1,nsites
    ind_basis =  AHTYPE_basis_index(site_info%site_type(n),parmunit)
    if ( ind_basis > 0 .and. ind_basis <= nbasis )then
      site_info%basis_index(n) = ind_basis
    else
      write(6,*)'AHTYPE_load_site_info: bad basis index for: ',sitename(n)
      stop
    endif
  enddo
  call File_Close(parmunit)
  ! get the number of prims

  allocate(site_info%num_primitives(nsites), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate num_primitives, exiting'
     stop
  endif

  allocate(site_info%off_primitives(nsites), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate off_primitives, exiting'
     stop
  endif

  nprim = 0
  do n = 1,nsites
    ind_basis = site_info%basis_index(n)
    site_info%num_primitives(n) = atom_auxis%num_basis_prims(ind_basis)
    nprim = nprim + atom_auxis%num_basis_prims(ind_basis)
  enddo
  site_info%tot_num_primitives = nprim
  write(6,*)'num primitives = ',nprim

  site_info%off_primitives(1) = 0
  do n = 2,nsites
    site_info%off_primitives(n) = site_info%off_primitives(n-1) + &
                                  site_info%num_primitives(n-1) 
  enddo

  allocate(site_info%prim_expo(nprim), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate prim_expo, exiting'
     stop
  endif

  allocate(site_info%prim_order(nprim), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate prim_order, exiting'
     stop
  endif

  allocate(site_info%coeff_offset(nprim), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate coeff_offset, exiting'
     stop
  endif

  nbprim = atom_auxis%num_primitives

  write(6,*)'num auxis = ',nbprim
  do n = 1,nsites
    off = site_info%off_primitives(n)
    num = site_info%num_primitives(n)
    ind_basis = site_info%basis_index(n)
    offb = atom_auxis%off_basis_prims(ind_basis) 
    do k = 1,num
      site_info%prim_order(off+k) = atom_auxis%prim_order(offb+k)
      site_info%prim_expo(off+k) = atom_auxis%prim_expon_term(offb+k)
    enddo
  enddo
  ncoeff = 0
  do n = 1,nprim
       ncoeff = ncoeff + site_info%prim_order(n)
  enddo
  write(6,*)'ncoeff = ',ncoeff
  site_info%num_coefficients = ncoeff

  site_info%coeff_offset(1) = 0
  do n = 2,nprim
    site_info%coeff_offset(n) = site_info%coeff_offset(n-1) + &
                                site_info%prim_order(n-1)
  enddo

  allocate(site_info%cartesian_coeffs(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHTYPE_load_site_info:could not allocate cartesian_coeffs,&
                & exiting'
     stop
  endif
  allocate(site_info%global_hermite_coeffs(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
    write(6,*)'AHTYPE_load_site_info:could not allocate global_hermite_coeffs,&
               & exiting'
     stop
  endif
  allocate(site_info%local_hermite_coeffs(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
   write(6,*)'AHTYPE_load_site_info:could not allocate local_hermite_coeffs,&
              & exiting'
     stop
  endif

  return
end subroutine AHTYPE_load_site_info_2
!-------------------------------------------------------
integer function AHTYPE_type_from_name(sname,parmunit)
  implicit none
  integer parmunit
  character(len=*) sname

  integer ios,ptr
  character(len=120) line,word
  integer TRDPRM_get_next_line,lenword
  rewind(parmunit)
  ios = 0
  do while (ios==0)
    line = ' '  ! to avoid junk in line after some reads
    ios = TRDPRM_get_next_line(parmunit,line)
    if ( ios == 0 )then
      ptr = 1 ! first call to str_next_token needs an initial value
      call str_next_token(line,ptr,word,lenword) !get 1st token
      if ( ptr > 0 )then !check if its 'sitetype'
        if ( word(1:lenword) == 'sitetype' )then
          call str_next_token(line,ptr,word,lenword) !get 1st token
          if ( sname(1:Len_trim(sname)) == word(1:lenword) )then
            read(line(ptr:),*)AHTYPE_type_from_name
            return
          endif
        endif
      endif !ptr > 0
    endif!ios == 0 )
  enddo !ios==0
  AHTYPE_type_from_name = -1
  return
end function AHTYPE_type_from_name
!-------------------------------------------------------
subroutine AHTYPE_charge_and_mass(stype,parmunit,charge,mass)
  implicit none
  integer stype,parmunit
  double precision charge,mass

  integer ios,ptr,ttype,tcharge
  double precision tmass
  character(len=120) line,word
  integer TRDPRM_get_next_line,lenword
  rewind(parmunit)
  ios = 0
  do while (ios==0)
    line = ' '  ! to avoid junk in line after some reads
    ios = TRDPRM_get_next_line(parmunit,line)
    if ( ios == 0 )then
      ptr = 1 ! first call to str_next_token needs an initial value
      call str_next_token(line,ptr,word,lenword) !get 1st token
      if ( ptr > 0 )then !check if its 'type'
        if ( word(1:lenword) == 'type' )then
          read(line(ptr:),*)ttype,tcharge,tmass
          if ( stype == ttype )then
            charge = tcharge
            mass = tmass
            return
          endif
        endif
      endif !ptr > 0
    endif!ios == 0 )
  enddo !ios==0

  charge = -1.d0
  mass = -1.d0
  return
end subroutine AHTYPE_charge_and_mass
!-------------------------------------------------------
integer function AHTYPE_basis_index(stype,parmunit)
  implicit none
  integer stype,parmunit

  integer ios,ptr,ttype,ind_basis
  character(len=120) line,word
  integer TRDPRM_get_next_line,lenword
  rewind(parmunit)
  ios = 0
  do while (ios==0)
    line = ' '  ! to avoid junk in line after some reads
    ios = TRDPRM_get_next_line(parmunit,line)
    if ( ios == 0 )then
      ptr = 1 ! first call to str_next_token needs an initial value
      call str_next_token(line,ptr,word,lenword) !get 1st token
      if ( ptr > 0 )then !check if its 'auxtype'
        if ( word(1:lenword) == 'auxtype' )then
          read(line(ptr:),*)ttype,ind_basis
          if ( stype == ttype )then
            AHTYPE_basis_index = ind_basis
            return
          endif
        endif
      endif !ptr > 0
    endif!ios == 0 )
  enddo !ios==0
  AHTYPE_basis_index = -1
  return
end function AHTYPE_basis_index
!-------------------------------------------------------
