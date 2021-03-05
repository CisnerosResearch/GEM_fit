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
!----------------------------------------------------------
subroutine AH_form_aux_coefs(site_info,Gmat,Jvec,aux_coefs,ncoeff,scale,&
                             constrain, debug)
  ! GAC: subroutine to invert Gmat and calculate the auxiliary coefficients
  !      in hermite local frame
use definition
  implicit none
  type(sites_info)::site_info
  integer ncoeff, nsites,i,j,k,l,lwork,info,allochk,tmpall
  integer slo1, shi1, site1, num1, off1, counter, k1, c_off1, l1, order1
  double precision Gmat(ncoeff,ncoeff)
  double precision Jvec(ncoeff)
  double precision aux_coefs(ncoeff)
  double precision,allocatable:: Jvec_bar(:), aux_coefs_bar(:),&
                                 theta(:), theta_bar(:), d(:), e(:), tau(:),&
                                 aux_charges(:), work(:), work2(:), a(:,:),&
                                 atmchg(:)
  double precision epsilon,lambda,tmp1,tmp2,num,denom,zero,scale,nelec
  double precision ecoul_aux, ecoul_J, num_elec, num_elec_atm, ecoul_intermol
  logical constrain, debug
  !parameter(epsilon = 1.d-16, zero = 0.0d0)
  parameter(epsilon = 1.d-15, zero = 0.0d0)

! first allocate all needed arrays

  allocate(a(ncoeff,ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AH_form_aux_coefs:could not allocate a, exiting'
     stop
  endif
  allocate(Jvec_bar(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AH_form_aux_coefs:could not allocate Jvec_bar, exiting'
     stop
  endif
  allocate(aux_coefs_bar(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AH_form_aux_coefs:could not allocate aux_coefs_bar, exiting'
     stop
  endif
  allocate(theta(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AH_form_aux_coefs:could not allocate theta, exiting'
     stop
  endif
  allocate(theta_bar(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AH_form_aux_coefs:could not allocate theta_bar, exiting'
     stop
  endif
  allocate(d(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AH_form_aux_coefs:could not allocate d, exiting'
     stop
  endif
  allocate(e(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AH_form_aux_coefs:could not allocate e, exiting'
     stop
  endif
  allocate(tau(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AH_form_aux_coefs:could not allocate tau, exiting'
     stop
  endif
  allocate(aux_charges(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AH_form_aux_coefs:could not allocate aux_charges, exiting'
     stop
  endif
  allocate(work2(2*ncoeff-1), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AH_form_aux_coefs:could not allocate work2, exiting'
     stop
  endif
  allocate(atmchg(site_info%nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AH_form_aux_coefs:could not allocate atmchg, exiting'
     stop
  endif

  lwork = ncoeff*32
  allocate(work(lwork), stat = allochk)
  if (allochk .gt. 0 ) then
     tmpall = allochk/ncoeff
     lwork = ncoeff*tmpall
     write(6,*)'AH_form_aux_coefs:could not allocate work in full,&
                & reducing memory allocation'
     allocate(work(lwork), stat = allochk)
     if (allochk .gt. 0 ) then
        write(6,*)'AH_form_aux_coefs:could not allocate work, exiting'
     endif
     stop
  endif

!
!!! Invert the matrix
!

! place Gmat in upper triangular matrix and neaux in theta
 
  do i=1,ncoeff
     do j = i,ncoeff
        a(i,j) = Gmat(i,j)
        if(i==j)a(i,j) = a(i,j)+scale
        !if(i==j.and.constrain)a(i,j) = a(i,j)+scale
        !print *,'gmat = ',a(i,j)
     enddo
     theta(i) = site_info%aux_elec(i)
     if(debug)write(12,*),'theta = ',theta(i)
  enddo

! decompose matrix into tridiagonal form
 
  call dsytrd('U', ncoeff, a, ncoeff, d, e, tau, work, lwork, info)
  if (info .ne. 0)then
     write (*,5001) -info
5001 format (/,T2,'Argument ',i3,' had illegal value in dsytrd, exiting'/)
     stop
  endif
  !print *,'work(1) = ',work(1), ncoeff

! form orthogonal matrix Q
 
  call dorgtr('U', ncoeff, a, ncoeff, tau, work, lwork, info)
  if (info .ne. 0)then
     write (*,5002) -info
5002 format (/,T2,'Argument ',i3,' had illegal value in dorgtr, exiting'/)
     stop
  endif

! diagonalize matrix (calculate eigenvectors [stored in a] _and_
!                               eigenvalues [stored in d])

  call dsteqr('V', ncoeff, d, e, a, ncoeff, work2, info)
  if (info .lt. 0)then
     write (*,5003) -info
5003 format (/,T2,'Argument ',i3,' had illegal value in dsteqr, exiting'/)
     stop
  else if (info .gt. 0)then
     write (*,5004) info
5004 format (/,T2,'dsteqr failed to find eigenvalue for ',i3,', exiting'/)
     stop
  endif

! quench and invert eigenvalues, and form y^-, theta^- and theta^-t
 
  if (debug) then
     do i=1,ncoeff
           write(12,*)'eigval = ',d(i)
     enddo
  endif
  do i = 1, ncoeff
     if (d(i) .lt. epsilon) then
        d(i) = zero
     else
        d(i) = 1/d(i)
     endif
  enddo 

  do i = 1, ncoeff
     Jvec_bar(i) = zero
     theta_bar(i) = zero
     do j = 1,ncoeff
        Jvec_bar(i) = Jvec_bar(i) + a(j,i)*Jvec(j)
        theta_bar(i) = theta_bar(i) + a(j,i)*theta(j)
     enddo
        !print *,Jvec_bar(i),theta_bar(i),d(i)
  enddo

! calculate lambda (Lagrange multiplier)
 
  nelec = zero
  do i = 1, site_info%natoms
     nelec = nelec + site_info%nuclear_charge(i)
  enddo
  nelec = nelec - site_info%mol_chg
  lambda = zero
  num = zero
  denom = zero
  do i = 1, ncoeff
     num = num + d(i)*theta_bar(i)*Jvec_bar(i)
     denom = denom + d(i)*theta_bar(i)*theta_bar(i)
  enddo
  num = nelec - num
  print *,'nelec = ',nelec,num
  lambda = num/denom

! calculate aux_coefs_bar (x^-) enforcing charge constraint (or not)
 
  if (constrain) then
     do i = 1, ncoeff
        aux_coefs_bar(i) = d(i)*Jvec_bar(i) + lambda*d(i)*theta_bar(i)
     enddo
  else
     do i = 1, ncoeff
        aux_coefs_bar(i) = d(i)*Jvec_bar(i)
     enddo
  endif

! calculate aux_coefs (x)
 
  do i = 1, ncoeff
      aux_coefs(i) = zero
      do j = 1,ncoeff
         aux_coefs(i) = aux_coefs(i) + a(i,j)*aux_coefs_bar(j)
      enddo
      if(debug)write(12,*),'aux coef = ',aux_coefs(i)
  enddo

! calculate electrostatic energy 

  ecoul_aux = zero
  ecoul_J = zero
  do i = 1, ncoeff
     aux_charges(i) = zero
     do j = 1, ncoeff
        aux_charges(i) = aux_charges(i) + Gmat(i,j)*aux_coefs(j)
     enddo
     ecoul_aux = ecoul_aux + aux_coefs(i)*aux_charges(i)
     ecoul_J = ecoul_J + Jvec(i)*aux_coefs(i)
  enddo
  print *,'electrostatic energy from aux = ',ecoul_aux/2.d0

! calculate number of electrons in auxiliary basis

  num_elec = 0.d0
  do i = 1, ncoeff
     num_elec = num_elec + theta(i)*aux_coefs(i)
  enddo
  print *,'Number of electrons in Aux = ',num_elec

! calculate number of electrons per site
5005    format (T2,'Site ',i3,' has ',f5.2,' electrons')
  atmchg(:) = 0.0d0
  counter = 0
  slo1 = site_info%residue_start(1)
  shi1 = site_info%residue_end(1)
  do site1 = slo1,shi1
    num1 = site_info%num_primitives(site1)
    off1 = site_info%off_primitives(site1)
    do k1 = 1,num1
       order1 = site_info%prim_order(off1+k1)
       do l1 = 1, order1
          counter = counter + 1
          atmchg(site1) = atmchg(site1) + theta(counter)*aux_coefs(counter)
       enddo
    enddo
    write (*,5005) site1, atmchg(site1)
  enddo
  !write (*,5005) site1, atmchg(site1)
       

!  num_elec_atm = 0.d0
!  curr_atm = sto_auxis(1)%atom_ctr_of_shell
!  do i = 1, ncoeff
!     if (sto_auxis(i)%atom_ctr_of_shell .eq. curr_atm) then
!        num_elec_atm = num_elec_atm + theta(i)*aux_coefs(i)
!     else
!!        write (*,5005) curr_atm, num_elec_atm
!        atmchg(curr_atm) = num_elec_atm
!        num_elec_atm = 0.d0
!        curr_atm = sto_auxis(i)%atom_ctr_of_shell
!        num_elec_atm = num_elec_atm + theta(i)*aux_coefs(i)
!     endif
!  enddo
!        curr_atm = sto_auxis(num_aux_stos)%atom_ctr_of_shell
!        atmchg(curr_atm) = num_elec_atm
!!        write (*,5006) curr_atm, num_elec_atm
!!5006    format (T2,'Atom ',i3,' has ',f5.2,' electrons')

  site_info%cartesian_coeffs(:) = aux_coefs(:)
  
  deallocate(a,Jvec_bar,aux_coefs_bar,theta,theta_bar,d,e,tau,aux_charges,&
             work2,work,atmchg)

  return
end subroutine AH_form_aux_coefs
!---------------------------------------------------------
subroutine AH_form_aux_coefs2(site_info,Gmat,Jvec,aux_coefs,ncoeff,scale,&
                             constrain, debug)
  ! GAC: subroutine to invert Gmat and calculate the auxiliary coefficients
  !      in hermite local frame
use definition
  implicit none
  type(sites_info)::site_info
  integer ncoeff, nsites,i,j,k,l,lwork,info,allochk,tmpall
  double precision Gmat(ncoeff,ncoeff)
  double precision Jvec(ncoeff)
  double precision aux_coefs(ncoeff)
  double precision,allocatable:: Jvec_bar(:), aux_coefs_bar(:),&
                                 theta(:), theta_bar(:), d(:), e(:), tau(:),&
                                 aux_charges(:), work(:), work2(:), a(:,:)
  double precision epsilon,lambda,tmp1,tmp2,num,denom,zero,scale,nelec
  double precision ecoul_aux, ecoul_J, num_elec, num_elec_atm, ecoul_intermol
  logical constrain, debug
  parameter(epsilon = 1.d-8, zero = 0.0d0)

! first allocate all needed arrays

  allocate(a(ncoeff,ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AH_form_aux_coefs:could not allocate a, exiting'
     stop
  endif
  allocate(Jvec_bar(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AH_form_aux_coefs:could not allocate Jvec_bar, exiting'
     stop
  endif
  allocate(aux_coefs_bar(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AH_form_aux_coefs:could not allocate aux_coefs_bar, exiting'
     stop
  endif
  allocate(theta(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AH_form_aux_coefs:could not allocate theta, exiting'
     stop
  endif
  allocate(theta_bar(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AH_form_aux_coefs:could not allocate theta_bar, exiting'
     stop
  endif
  allocate(d(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AH_form_aux_coefs:could not allocate d, exiting'
     stop
  endif
  allocate(e(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AH_form_aux_coefs:could not allocate e, exiting'
     stop
  endif
  allocate(tau(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AH_form_aux_coefs:could not allocate tau, exiting'
     stop
  endif
  allocate(aux_charges(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AH_form_aux_coefs:could not allocate aux_charges, exiting'
     stop
  endif
  allocate(work2(2*ncoeff-1), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AH_form_aux_coefs:could not allocate work2, exiting'
     stop
  endif

  lwork = ncoeff*32
  allocate(work(lwork), stat = allochk)
  if (allochk .gt. 0 ) then
     tmpall = allochk/ncoeff
     lwork = ncoeff*tmpall
     write(6,*)'AH_form_aux_coefs:could not allocate work in full,&
                & reducing memory allocation'
     allocate(work(lwork), stat = allochk)
     if (allochk .gt. 0 ) then
        write(6,*)'AH_form_aux_coefs:could not allocate work, exiting'
     endif
     stop
  endif

!
!!! Invert the matrix
!

! place Gmat in upper triangular matrix and neaux in theta
 
  do i=1,ncoeff
     do j = i,ncoeff
        a(i,j) = Gmat(i,j)
        if(i==j)a(i,j) = a(i,j)+scale
        !if(i==j.and.constrain)a(i,j) = a(i,j)+scale
        !print *,'gmat = ',a(i,j)
     enddo
     theta(i) = site_info%aux_elec(i)
     print *,'theta = ',theta(i)
  enddo

! decompose matrix into tridiagonal form
 
  call dsytrd('U', ncoeff, a, ncoeff, d, e, tau, work, lwork, info)
  if (info .ne. 0)then
     write (*,5001) -info
5001 format (/,T2,'Argument ',i3,' had illegal value in dsytrd, exiting'/)
     stop
  endif
  !print *,'work(1) = ',work(1), ncoeff

! form orthogonal matrix Q
 
  call dorgtr('U', ncoeff, a, ncoeff, tau, work, lwork, info)
  if (info .ne. 0)then
     write (*,5002) -info
5002 format (/,T2,'Argument ',i3,' had illegal value in dorgtr, exiting'/)
     stop
  endif

! diagonalize matrix (calculate eigenvectors [stored in a] _and_
!                               eigenvalues [stored in d])

  call dsteqr('V', ncoeff, d, e, a, ncoeff, work2, info)
  if (info .lt. 0)then
     write (*,5003) -info
5003 format (/,T2,'Argument ',i3,' had illegal value in dsteqr, exiting'/)
     stop
  else if (info .gt. 0)then
     write (*,5004) info
5004 format (/,T2,'dsteqr failed to find eigenvalue for ',i3,', exiting'/)
     stop
  endif

! quench and invert eigenvalues, and form y^-, theta^- and theta^-t
 
  if (debug) then
     do i=1,ncoeff
           write(12,*)'eigval PROMOL = ',d(i)
     enddo
  endif
  do i = 1, ncoeff
     if (d(i) .lt. epsilon) then
        d(i) = zero
     else
        d(i) = 1/d(i)
     endif
  enddo 

  do i = 1, ncoeff
     Jvec_bar(i) = zero
     theta_bar(i) = zero
     do j = 1,ncoeff
        Jvec_bar(i) = Jvec_bar(i) + a(j,i)*Jvec(j)
        theta_bar(i) = theta_bar(i) + a(j,i)*theta(j)
     enddo
  enddo

! calculate lambda (Lagrange multiplier)
 
  nelec = zero
  do i = 1, site_info%natoms
     nelec = nelec + site_info%nuclear_charge(i)
  enddo
  nelec = zero
  lambda = zero
  num = zero
  denom = zero
  do i = 1, ncoeff
     num = num + d(i)*theta_bar(i)*Jvec_bar(i)
     denom = denom + d(i)*theta_bar(i)*theta_bar(i)
  enddo
  num = nelec - num
  print *,'nelec = ',nelec,num
  lambda = num/denom

! calculate aux_coefs_bar (x^-) enforcing charge constraint (or not)
 
  if (constrain) then
     do i = 1, ncoeff
        aux_coefs_bar(i) = d(i)*Jvec_bar(i) + lambda*d(i)*theta_bar(i)
     enddo
  else
     do i = 1, ncoeff
        aux_coefs_bar(i) = d(i)*Jvec_bar(i)
     enddo
  endif

! calculate aux_coefs (x)
 
  do i = 1, ncoeff
      aux_coefs(i) = zero
      do j = 1,ncoeff
         aux_coefs(i) = aux_coefs(i) + a(i,j)*aux_coefs_bar(j)
      enddo
      if(debug)write(12,*),'aux coef = ',aux_coefs(i)
  enddo

! calculate electrostatic energy 

  ecoul_aux = zero
  ecoul_J = zero
  do i = 1, ncoeff
     aux_charges(i) = zero
     do j = 1, ncoeff
        aux_charges(i) = aux_charges(i) + Gmat(i,j)*aux_coefs(j)
     enddo
     ecoul_aux = ecoul_aux + aux_coefs(i)*aux_charges(i)
     ecoul_J = ecoul_J + Jvec(i)*aux_coefs(i)
  enddo
  print *,'electrostatic energy from aux = ',ecoul_aux/2.d0

! calculate number of electrons in auxiliary basis

  num_elec = 0.d0
  do i = 1, ncoeff
     num_elec = num_elec + theta(i)*aux_coefs(i)
  enddo
  print *,'Number of electrons in Aux = ',num_elec

  site_info%cartesian_coeffs(:) = aux_coefs(:)
  
  deallocate(a,Jvec_bar,aux_coefs_bar,theta,theta_bar,d,e,tau,aux_charges,&
             work2,work)

  return
end subroutine AH_form_aux_coefs2
!----------------------------------------------------------
