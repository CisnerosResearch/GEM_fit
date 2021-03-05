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
subroutine H_form_J_vec_C(site_info,auxis,basis,Jvec,ncoeff,num_stos,&
                          inttype,beta,debug,cube,promolfit,twocentfit,&
                          densfile)
  ! GAC: subroutine to calculate <k||l> and <k|l> in local hermite.
use definition
  implicit none
  type(sites_info)::site_info
  type(aux_orbitals)::auxis(ncoeff)
  type(ao_orbitals)::basis(num_stos)

  integer site1,allochk,inttype
  integer k1,num1,off1,c_off1,npa,npb,nx1,ny1,nz1
  integer ncoeff,num_stos,nx2,ny2,nz2,natoms,lim1lo,lim1high
  integer ires,jres,n,i,j,k,deg1,deg2,i2,j2,counter,nx3,ny3,nz3,l,m,init,fin
  integer,allocatable::limitvec(:)
  double precision alphap,R_ab,brak,pre,temp,px,py,pz,norma,temp2
  double precision dx,dy,dz,Jvec(ncoeff),NUCvec(ncoeff)
  double precision x1,x2,y1,y2,z1,z2,x3,y3,z3
  double precision boyspar,result,beta,beta2,boyspar2,&
                   expon3,norm3,result2
  double precision,allocatable::expon1(:),expon2(:),norm1(:),norm2(:),&
                   contr1(:),contr2(:)
  double precision D(0:3,0:3,0:6),E(0:3,0:3,0:6),F(0:3,0:3,0:6),&
                   D1(0:3,0:3),E1(0:3,0:3),F1(0:3,0:3),&
                   DDS(0:3,0:3),EDS(0:3,0:3),FDS(0:3,0:3),&
                   densmat(num_stos,num_stos) 
  character(len=*) densfile
  logical normalize, debug, cube, promolfit, twocentfit

  Jvec(:) = 0.0d0
  NUCvec(:) = 0.0d0
  if(inttype == 1)beta2 = beta*beta
  
  counter = 0
  do i = 1, num_stos
     do j = 1, i
        counter = counter+1
        if (.not. promolfit) then
           densmat(j,i) = site_info%densmat(counter)
        else
           densmat(j,i) = site_info%promol_dens(counter)
        endif
        densmat(i,j) = densmat(j,i)
     enddo
  enddo

  if (.not. cube .and. .not. promolfit) then
     call H_form_basis(site_info,basis,num_stos,densfile)
     if (twocentfit) then
        natoms = site_info%natoms
        allocate(limitvec(natoms), stat=allochk)
        if (allochk .gt. 0) then
           write(6,*)'H_form_J_vec_C: could not allocate limitvec, exiting'
           stop
        endif
        call H_basis_limits(site_info,basis,num_stos,natoms,limitvec)
     endif
  endif

! para imprimir los coefficientes de 2 centros
  if(twocentfit)then
  open(61,file="blocked_mat",status="unknown")
  write(61,*)num_stos
  counter=1
  do i = 1, num_stos
     if (i .eq. 1) then
        lim1lo = 1
        lim1high = limitvec(counter)
     else if ((i .eq. lim1high+1) .and. (i .lt. num_stos)) then
        counter = counter+1
        lim1lo = lim1high+1
        lim1high = limitvec(counter)
     else if (i .eq. num_stos) then
        if (counter .ne. natoms) then
           print *,'H_exact_esp: something wrong with fit limits, &
                    &exiting'
           stop
        endif
     endif
     do j = 1, i
        !if ((j .ge. lim1lo) .or. (j .le. lim1high))then
        if ((j .lt. lim1lo) .or. (j .gt. lim1high))then
           !write(61,*) 0.0d0
           write(61,*) densmat(i,j)
        else   
           write(61,*) 0.0d0
           !write(61,*) densmat(i,j)
        endif
     enddo
  enddo
  close(61)
  endif

  counter = 1
  !do i = 1, num_stos
  do i = 1, num_stos
     if (twocentfit) then
        if (i .eq. 1) then
           lim1lo = 1
           lim1high = limitvec(counter)
        else if ((i .eq. lim1high+1) .and. i .lt. num_stos) then
           counter = counter+1
           lim1lo = lim1high+1
           lim1high = limitvec(counter)
        else if (i .eq. num_stos) then
           if (counter .ne. natoms) then
              print *,'H_form_J_vec_C: something wrong with fit limits, &
                       &exiting'
              stop
           endif
        endif
     endif
     npa = basis(i)%deg_contr  
     allocate(expon1(npa),contr1(npa),norm1(npa), stat = allochk)
     if (allochk .gt. 0 ) then
        write(6,*)'H_form_Jvec:could not allocate expon1, contr1, norm1;&
                   & exiting'
        stop
     endif     
     expon1(:) = basis(i)%expo(:)
     contr1(:) = basis(i)%contr(:)
     norm1(:) = basis(i)%norm(:)
     nx1 = basis(i)%x  
     ny1 = basis(i)%y 
     nz1 = basis(i)%z 
     x1 = basis(i)%coords(1)  
     y1 = basis(i)%coords(2) 
     z1 = basis(i)%coords(3) 
     if (twocentfit) then

        do j = 1, lim1lo-1
         if ((j .lt. lim1lo) .or. (j .gt. lim1high))then
           !write(60,*),'doing pair ',i,j
           npb = basis(j)%deg_contr  
           allocate(expon2(npb),contr2(npb),norm2(npb), stat = allochk)
           if (allochk .gt. 0 ) then
              write(6,*)'H_form_Jvec:could not allocate expon2, contr2, norm2;&
                         & exiting'
              stop
           endif     
           expon2(:) = basis(j)%expo(:)
           contr2(:) = basis(j)%contr(:)
           norm2(:) = basis(j)%norm(:)
           nx2 = basis(j)%x  
           ny2 = basis(j)%y 
           nz2 = basis(j)%z 
           x2 = basis(j)%coords(1)  
           y2 = basis(j)%coords(2) 
           z2 = basis(j)%coords(3) 
           deg1 = nx1+ny1+nz1+nx2+ny2+nz2 
           R_ab = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
           do k = 1, ncoeff
              expon3 = auxis(k)%expo
              nx3 = auxis(k)%x  
              ny3 = auxis(k)%y 
              nz3 = auxis(k)%z 
              norma = auxis(k)%norm
              deg2 = nx3+ny3+nz3
              call H_form_exp_coef_noS(nx3,x3,expon3,DDS)
              call H_form_exp_coef_noS(ny3,y3,expon3,EDS)
              call H_form_exp_coef_noS(nz3,z3,expon3,FDS)
              temp = 0.0d0
              do i2 = 1, npa
                 do j2 = 1, npb
                    alphap = expon1(i2) + expon2(j2)
                    px = (expon1(i2)*x1 + expon2(j2)*x2)/alphap
                    py = (expon1(i2)*y1 + expon2(j2)*y2)/alphap
                    pz = (expon1(i2)*z1 + expon2(j2)*z2)/alphap
                    call H_form_exp_coef(nx1,nx2,x1,x2,px,alphap,D)
                    call H_form_exp_coef(ny1,ny2,y1,y2,py,alphap,E)
                    call H_form_exp_coef(nz1,nz2,z1,z2,pz,alphap,F)
                    brak = exp(-1.d0*expon1(i2)*expon2(j2)*R_ab/alphap)
                    dx = px - auxis(k)%coords_G(1) ! because J needs to be in
                    dy = py - auxis(k)%coords_G(2) ! global coordinates
                    dz = pz - auxis(k)%coords_G(3)
                    boyspar = alphap*expon3 / (alphap + expon3)
                    if (inttype == 0) then
                       pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*brak*norma
                       call GN_MD_REC_J_C(deg1,deg2,boyspar,dx,dy,dz,result,&
                                     nx1,ny1,nz1,nx2,ny2,nz2,nx3,ny3,nz3,&
                                     alphap,expon3,D,E,F,DDS,EDS,FDS,0)
                    else if (inttype == 1) then
                       pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*brak*norma
                       boyspar2 = (beta2*alphap*expon3)/(beta2*(alphap+expon3)+&
                                alphap*expon3)
                       call GN_MD_REC_erfc_J_C(deg1,deg2,boyspar,dx,dy,dz,&
                                        result,nx1,ny1,nz1,nx2,ny2,nz2,nx3,ny3,&
                                        nz3,alphap,expon3,boyspar2,D,E,F,DDS,&
                                        EDS,FDS)
                    else if (inttype == 2) then
                       pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*brak*norma
                       call GN_MD_REC_OVERLAP_J_C(deg1,deg2,boyspar,dx,dy,dz,&
                                        result,nx1,ny1,nz1,nx2,ny2,nz2,nx3,ny3,&
                                        nz3,alphap,expon3,D,E,F,DDS,EDS,FDS)
                    endif
                    temp = temp+result*pre
                    !if(debug)write(12,*)result,result*pre
                 enddo ! j2 = 1, npb
              enddo ! i2 = 1, npa
              if (i .eq. j) then
                 Jvec(k) = Jvec(k) + temp*densmat(i,j)
              else
                 Jvec(k) = Jvec(k) + temp*densmat(i,j)*2.0d0
              endif
           enddo ! k = 1, ncoeff
           deallocate(expon2,contr2,norm2)
         endif ! i for lim1lo lim1high
        enddo ! j = 1, lim1lo

     else !not two cent fit, fit the whole densmat
        do j = 1, i
           npb = basis(j)%deg_contr  
           allocate(expon2(npb),contr2(npb),norm2(npb), stat = allochk)
           if (allochk .gt. 0 ) then
              write(6,*)'H_form_Jvec:could not allocate expon2, contr2, norm2;&
                         & exiting'
              stop
           endif     
           expon2(:) = basis(j)%expo(:)
           contr2(:) = basis(j)%contr(:)
           norm2(:) = basis(j)%norm(:)
           nx2 = basis(j)%x  
           ny2 = basis(j)%y 
           nz2 = basis(j)%z 
           x2 = basis(j)%coords(1)  
           y2 = basis(j)%coords(2) 
           z2 = basis(j)%coords(3) 
           deg1 = nx1+ny1+nz1+nx2+ny2+nz2 
           R_ab = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
           do k = 1, ncoeff
              expon3 = auxis(k)%expo
              nx3 = auxis(k)%x  
              ny3 = auxis(k)%y 
              nz3 = auxis(k)%z 
              norma = auxis(k)%norm
              deg2 = nx3+ny3+nz3
              call H_form_exp_coef_noS(nx3,x3,expon3,DDS)
              call H_form_exp_coef_noS(ny3,y3,expon3,EDS)
              call H_form_exp_coef_noS(nz3,z3,expon3,FDS)
              temp = 0.0d0
              do i2 = 1, npa
                 do j2 = 1, npb
                    alphap = expon1(i2) + expon2(j2)
                    px = (expon1(i2)*x1 + expon2(j2)*x2)/alphap
                    py = (expon1(i2)*y1 + expon2(j2)*y2)/alphap
                    pz = (expon1(i2)*z1 + expon2(j2)*z2)/alphap
                    call H_form_exp_coef(nx1,nx2,x1,x2,px,alphap,D)
                    call H_form_exp_coef(ny1,ny2,y1,y2,py,alphap,E)
                    call H_form_exp_coef(nz1,nz2,z1,z2,pz,alphap,F)
                    brak = exp(-1.d0*expon1(i2)*expon2(j2)*R_ab/alphap)
                    dx = px - auxis(k)%coords_G(1) ! because J needs to be in
                    dy = py - auxis(k)%coords_G(2) ! global coordinates
                    dz = pz - auxis(k)%coords_G(3)
                    boyspar = alphap*expon3 / (alphap + expon3)
                    if (inttype == 0) then
                      pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*brak*norma
                       call GN_MD_REC_J_C(deg1,deg2,boyspar,dx,dy,dz,result,&
                                     nx1,ny1,nz1,nx2,ny2,nz2,nx3,ny3,nz3,&
                                     alphap,expon3,D,E,F,DDS,EDS,FDS,0)
                    else if (inttype == 1) then
                      pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*brak*norma
                       boyspar2 = (beta2*alphap*expon3)/(beta2*(alphap+expon3)+&
                                alphap*expon3)
                       call GN_MD_REC_erfc_J_C(deg1,deg2,boyspar,dx,dy,dz,&
                                        result,nx1,ny1,nz1,nx2,ny2,nz2,nx3,ny3,&
                                        nz3,alphap,expon3,boyspar2,D,E,F,DDS,&
                                        EDS,FDS)
                    else if (inttype == 2) then
                      pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*brak*norma
                       call GN_MD_REC_OVERLAP_J_C(deg1,deg2,boyspar,dx,dy,dz,&
                                        result,nx1,ny1,nz1,nx2,ny2,nz2,nx3,ny3,&
                                        nz3,alphap,expon3,D,E,F,DDS,EDS,FDS)
                    endif
                    temp = temp+result*pre
                    !if(debug)write(12,*)result,result*pre
                 enddo ! j2 = 1, npb
              enddo ! i2 = 1, npa
              if (i .eq. j) then
                 Jvec(k) = Jvec(k) + temp*densmat(i,j)
              else
                 Jvec(k) = Jvec(k) + temp*densmat(i,j)*2.0d0
              endif
           enddo ! k = 1, ncoeff
           deallocate(expon2,contr2,norm2)
        enddo ! j = 1, num_stos
     endif ! twocentfit
     deallocate(expon1,contr1,norm1)
  enddo ! i = 1, num_stos

!!!! NUCLEAR ESP FOR TOM
!  DDS(:,:) = 1.0d0
!  EDS(:,:) = 1.0d0
!  FDS(:,:) = 1.0d0
!  do i = 1, ncoeff
!     expon3 = auxis(i)%expo
!     x3 = auxis(i)%x  
!     y3 = auxis(i)%y 
!     z3 = auxis(i)%z 
!     deg1 = x3+y3+z3 
!     norm3 = auxis(i)%norm
!     call H_form_exp_coef_noS(x3,auxis(i)%coords_G(1),expon3,D1)
!     call H_form_exp_coef_noS(y3,auxis(i)%coords_G(2),expon3,E1)
!     call H_form_exp_coef_noS(z3,auxis(i)%coords_G(3),expon3,F1)
!     boyspar = expon3
!     temp2 = 0.0d0
!     do j = 1, site_info%natoms
!        deg2 = 0
!        dx = auxis(i)%coords_G(1) - site_info%site_crds(3*(j-1)+1)
!        dy = auxis(i)%coords_G(2) - site_info%site_crds(3*(j-1)+2)
!        dz = auxis(i)%coords_G(3) - site_info%site_crds(3*(j-1)+3)
!        call GN_MD_REC_C(deg1,deg2,boyspar,dx,dy,dz,result2,x3,y3,z3,&
!                         0,0,0,expon3,3.14159265358979323846d0,D1,E1,F1,&
!                         DDS,EDS,FDS,0)
!        temp2 = temp2 + result2*site_info%nuclear_charge(j)
!     enddo
!     NUCvec(i) = NUCvec(i) + temp2*norm3
!  enddo
!!!! NUCLEAR ESP FOR TOM
!
!  open(75,file='Jvec',status='unknown')
!  open(76,file='ESPvec',status='unknown')
!  write(75,*)ncoeff
!  write(76,*)ncoeff
!  do i = 1, ncoeff
!     write(75,*)Jvec(i)
!     write(76,*)Jvec(i)+NUCvec(i)
!  enddo
!  close(75)
!  close(76)
     
  if (debug) then
     do i = 1, ncoeff
        write(12,1001)i,Jvec(i)
     enddo
  endif
1001    format ('Jvec (',i4,') = ',f27.20)

  return
end subroutine H_form_J_vec_C
!---------------------------------------------------------
