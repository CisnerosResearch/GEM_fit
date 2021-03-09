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
subroutine H_exact_esp(site_info,basis,cx,cy,cz,tot_esp,limitvec,&
                       natoms,num_stos,cube,twocentfit)
use definition
  implicit none
  type(sites_info)::site_info
  type(ao_orbitals)::basis(num_stos)

  integer site1,allochk,inttype,natoms
  integer k1,num1,off1,c_off1,npa,npb,nx1,ny1,nz1
  integer ncoeff,num_stos,nx2,ny2,nz2,limitvec(natoms)
  integer ires,jres,n,i,j,k,deg1,deg2,i2,j2,counter,nx3,ny3,nz3,l,m
  integer n_esp_points,slo1,shi1,lim1lo,lim1high
  double precision alphap,R_ab,brak,pre,temp,px,py,pz,norma,tot_esp,tot_N
  double precision dx,dy,dz,cx,cy,cz,dist,nuc_chg1
  double precision x1,x2,y1,y2,z1,z2,x3,y3,z3
  double precision boyspar,result,beta,beta2,boyspar2,&
                   expon3,norm3
  double precision,allocatable::expon1(:),expon2(:),norm1(:),norm2(:),&
                   contr1(:),contr2(:)
  double precision D(0:3,0:3,0:6),E(0:3,0:3,0:6),F(0:3,0:3,0:6),&
                   DDS(0:3,0:3),EDS(0:3,0:3),FDS(0:3,0:3),&
                   densmat(num_stos,num_stos) 
  double precision,parameter::AtoBohr = 0.529177249d0,tiny=1.d-12
  logical normalize, debug, cube, off_site, twocentfit

  counter = 0
  do i = 1, num_stos
     do j = 1, i
        counter = counter+1
        densmat(j,i) = site_info%densmat(counter)
        densmat(i,j) = densmat(j,i)
     enddo
  enddo

  nx3 = 0
  ny3 = 0
  nz3 = 0
  deg2 = 0
  DDS(:,:) = 1.0d0
  EDS(:,:) = 1.0d0
  FDS(:,:) = 1.0d0
  ! Just in case we want to calculate ESP at a single point
  inquire(file = 'esp_point', exist = off_site)
  if (off_site .and. .not. cube) then
     open(13,file='esp_point',status='unknown')
     write (*,*) 'H_exact_esp: found esp_point, will calculate ESP at points'
     read (13,*),n_esp_points
  else
     n_esp_points = 1
  endif
  do l = 1, n_esp_points
     tot_esp = 0.0d0
     tot_N = 0.0d0
     if (off_site .and. .not. cube) then
        read (13,*),cx, cy, cz
        cx=cx/AtoBohr
        cy=cy/AtoBohr 
        cz=cz/AtoBohr
     endif
     counter = 1
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
                 print *,'H_exact_esp: something wrong with fit limits, &
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
              npb = basis(j)%deg_contr  
              allocate(expon2(npb),contr2(npb),norm2(npb), stat = allochk)
              if (allochk .gt. 0 ) then
                 write(6,*)'H_nuc_elec:could not allocate expon2, contr2, &
                            &norm2; exiting'
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
                    boyspar = alphap
                    pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*brak
                    dx = px - cx
                    dy = py - cy
                    dz = pz - cz
                    expon3 = 3.14159265358979323846d0
                    call GN_MD_REC_J_C(deg1,deg2,boyspar,dx,dy,dz,result,&
                                       nx1,ny1,nz1,nx2,ny2,nz2,nx3,ny3,nz3,&
                                       alphap,expon3,D,E,F,DDS,EDS,FDS,0)
                    temp = temp+result*pre
                 enddo ! j2 = 1, npb
              enddo ! i2 = 1, npa
              if (i .eq. j) then
                 tot_esp = tot_esp - temp*densmat(i,j)
              else
                 tot_esp = tot_esp - temp*densmat(i,j)*2.0d0
              endif
              deallocate(expon2,contr2,norm2)
             endif ! only do one center part
           enddo ! j = 1, i
        else
           do j = 1, i
              npb = basis(j)%deg_contr  
              allocate(expon2(npb),contr2(npb),norm2(npb), stat = allochk)
              if (allochk .gt. 0 ) then
                 write(6,*)'H_nuc_elec:could not allocate expon2, contr2, &
                            &norm2; exiting'
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
                    boyspar = alphap
                    pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*brak
                    dx = px - cx
                    dy = py - cy
                    dz = pz - cz
                    expon3 = 3.14159265358979323846d0
                    call GN_MD_REC_J_C(deg1,deg2,boyspar,dx,dy,dz,result,&
                                       nx1,ny1,nz1,nx2,ny2,nz2,nx3,ny3,nz3,&
                                       alphap,expon3,D,E,F,DDS,EDS,FDS,0)
                    temp = temp+result*pre
                 enddo ! j2 = 1, npb
              enddo ! i2 = 1, npa
              if (i .eq. j) then
                 tot_esp = tot_esp - temp*densmat(i,j)
              else
                 tot_esp = tot_esp - temp*densmat(i,j)*2.0d0
              endif
              deallocate(expon2,contr2,norm2)
           enddo ! j = 1, i
        endif ! for twocentfit integrals only
        deallocate(expon1,contr1,norm1)
     enddo ! i = 1, num_stos
   ! if (.not. twocentfit) then !HABRA QUE QUITAR ESP DE LOS NUCLEOS?
   ! now do nuclear ESP
     slo1 = site_info%residue_start(1)
     shi1 = site_info%residue_end(1)
     do site1 = slo1,shi1
       dx = site_info%site_crds(3*(site1-1)+1) - cx
       if (abs(dx) .lt. tiny) dx = 0.0d0
       dy = site_info%site_crds(3*(site1-1)+2) - cy
       if (abs(dy) .lt. tiny) dy = 0.0d0
       dz = site_info%site_crds(3*(site1-1)+3) - cz
       if (abs(dz) .lt. tiny) dz = 0.0d0
       dist = dsqrt(dx*dx+dy*dy+dz*dz)
       nuc_chg1 = site_info%nuclear_charge(site1)
       if(site1 .le. natoms) then
         if (dist .gt. tiny) tot_N = tot_N + nuc_chg1/dist
       endif
     enddo
     if (off_site)write(6,1001) cx*AtoBohr,cy*AtoBohr,cz*AtoBohr,tot_esp
     if (off_site)write(6,1002) cx*AtoBohr,cy*AtoBohr,cz*AtoBohr,tot_N
     tot_esp = tot_esp+tot_N
   ! endif ! for twocentfit
     if (.not. cube)write(6,1003) cx*AtoBohr,cy*AtoBohr,cz*AtoBohr,tot_esp
     if (.not. cube)write(6,1004)
  enddo ! n_esp_points
1001 format (1X,'EXACT ELECTRONIC ESP at',3(1X,F10.6),' = ',F20.15)
1002 format (1X,'EXACT NUCLEAR ESP at',3(1X,F10.6),' = ',F20.15)
1003 format (1X,'EXACT ESP at',3(1X,F10.6),' = ',F20.15)
1004 format (' ')
  
  close(13)
     
  return
end subroutine H_exact_esp
!---------------------------------------------------------
subroutine H_exact_fld(site_info,basis,cx,cy,cz,tot_fld,limitvec,&
                       natoms,num_stos,cube,twocentfit)
use definition
  implicit none
  type(sites_info)::site_info
  type(ao_orbitals)::basis(num_stos)

  integer site1,allochk,inttype,natoms
  integer k1,num1,off1,c_off1,npa,npb,nx1,ny1,nz1
  integer ncoeff,num_stos,nx2,ny2,nz2,limitvec(natoms)
  integer ires,jres,n,i,j,k,deg1,deg2,i2,j2,counter,nx3,ny3,nz3,l,m
  integer n_esp_points,slo1,shi1,fldint,lim1lo,lim1high
  double precision alphap,R_ab,brak,pre,temp,px,py,pz,norma,tot_fld(3),tot_N
  double precision dx,dy,dz,cx,cy,cz,dist,nuc_chg1,tot_E
  double precision x1,x2,y1,y2,z1,z2,x3,y3,z3,xga,yga,zga,ar3
  double precision boyspar,result,beta,beta2,boyspar2,&
                   expon3,norm3
  double precision,allocatable::expon1(:),expon2(:),norm1(:),norm2(:),&
                   contr1(:),contr2(:)
  double precision D(0:3,0:3,0:6),E(0:3,0:3,0:6),F(0:3,0:3,0:6),&
                   DDS(0:3,0:3),EDS(0:3,0:3),FDS(0:3,0:3),&
                   densmat(num_stos,num_stos) 
  double precision,parameter::AtoBohr = 0.529177249d0,tiny=1.d-15
  logical normalize, debug, cube, off_site, twocentfit
  character*1 int_type

  counter = 0
  do i = 1, num_stos
     do j = 1, i
        counter = counter+1
        densmat(j,i) = site_info%densmat(counter)
        densmat(i,j) = densmat(j,i)
     enddo
  enddo

  slo1 = site_info%residue_start(1)
  shi1 = site_info%residue_end(1)
  nx3 = 0
  ny3 = 0
  nz3 = 0
  deg2 = 0
  tot_fld(:) = 0.0d0
  DDS(:,:) = 1.0d0
  EDS(:,:) = 1.0d0
  FDS(:,:) = 1.0d0
  ! Just in case we want to calculate ESP at a single point
  inquire(file = 'esp_point', exist = off_site)
  if (off_site .and. .not. cube) then
     open(13,file='esp_point',status='unknown')
     read (13,*),n_esp_points
  else
     n_esp_points = 1
  endif
  do l = 1, n_esp_points
     if (off_site .and. .not. cube) then
        read (13,*),cx, cy, cz
        cx=cx/AtoBohr
        cy=cy/AtoBohr 
        cz=cz/AtoBohr
     endif
     do fldint = 1, 3
        tot_E = 0.0d0
        tot_N = 0.0d0
        counter = 1
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
                    print *,'H_exact_esp: something wrong with fit limits, &
                             &exiting'
                    stop
                 endif
              endif
           endif
           npa = basis(i)%deg_contr  
           allocate(expon1(npa),contr1(npa),norm1(npa), stat = allochk)
           if (allochk .gt. 0 ) then
              write(6,*)'H_form_Jvec:could not allocate expon1, contr1,&
                         & norm1; exiting'
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
           if (twocentfit) then ! do only one center part
             do j = 1, lim1lo-1
               if ((j .lt. lim1lo) .or. (j .gt. lim1high))then
                 npb = basis(j)%deg_contr  
                 allocate(expon2(npb),contr2(npb),norm2(npb), stat = allochk)
                 if (allochk .gt. 0 ) then
                    write(6,*)'H_nuc_elec:could not allocate expon2, contr2,&
                               & norm2; exiting'
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
                       boyspar = alphap
                       pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*brak
                       dx = px - cx
                       dy = py - cy
                       dz = pz - cz
                       expon3 = 3.14159265358979323846d0
                       call GN_MD_REC_J_C(deg1,deg2,boyspar,dx,dy,dz,result,&
                                          nx1,ny1,nz1,nx2,ny2,nz2,nx3,ny3,nz3,&
                                          alphap,expon3,D,E,F,DDS,EDS,FDS,&
                                          fldint)
                       temp = temp+result*pre
                    enddo ! j2 = 1, npb
                 enddo ! i2 = 1, npa
                 if (i .eq. j) then
                    tot_E = tot_E + temp*densmat(i,j)
                 else
                    tot_E = tot_E + temp*densmat(i,j)*2.0d0
                 endif
                 deallocate(expon2,contr2,norm2)
              endif ! lim1lo and lim1high
             enddo ! j = 1, i
           else !do whole fit
              do j = 1, i
                 npb = basis(j)%deg_contr  
                 allocate(expon2(npb),contr2(npb),norm2(npb), stat = allochk)
                 if (allochk .gt. 0 ) then
                    write(6,*)'H_nuc_elec:could not allocate expon2, contr2,&
                               & norm2; exiting'
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
                       boyspar = alphap
                       pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*brak
                       dx = px - cx
                       dy = py - cy
                       dz = pz - cz
                       expon3 = 3.14159265358979323846d0
                       call GN_MD_REC_J_C(deg1,deg2,boyspar,dx,dy,dz,result,&
                                          nx1,ny1,nz1,nx2,ny2,nz2,nx3,ny3,nz3,&
                                          alphap,expon3,D,E,F,DDS,EDS,FDS,&
                                          fldint)
                       temp = temp+result*pre
                    enddo ! j2 = 1, npb
                 enddo ! i2 = 1, npa
                 if (i .eq. j) then
                    tot_E = tot_E + temp*densmat(i,j)
                 else
                    tot_E = tot_E + temp*densmat(i,j)*2.0d0
                 endif
                 deallocate(expon2,contr2,norm2)
              enddo ! j = 1, num_stos
           endif ! for twocent fit only
           deallocate(expon1,contr1,norm1)
        enddo ! i = 1, num_stos

      ! now do nuclear FLD
      ! if (.not. twocentfit) then ! QUITAR CAMPO NUCLEAR TAMBIEN?
        do site1 = slo1,shi1
          dx = site_info%site_crds(3*(site1-1)+1) - cx
          if (abs(dx) .lt. tiny) dx = 0.0d0
          dy = site_info%site_crds(3*(site1-1)+2) - cy
          if (abs(dy) .lt. tiny) dy = 0.0d0
          dz = site_info%site_crds(3*(site1-1)+3) - cz
          if (abs(dz) .lt. tiny) dz = 0.0d0
          dist = dx*dx + dy*dy + dz*dz
          nuc_chg1 = site_info%nuclear_charge(site1)
          if (dist .gt. tiny) then
             dist = 1.0d0/dist
             xga = cx - site_info%site_crds(3*(site1-1)+1)
             yga = cy - site_info%site_crds(3*(site1-1)+2)
             zga = cz - site_info%site_crds(3*(site1-1)+3)
             ar3 = nuc_chg1 * dsqrt(dist) * dist
             if(fldint .eq. 1) then
               tot_N = tot_N + ar3*xga
             else if(fldint .eq. 2) then
               tot_N = tot_N + ar3*yga
             else if(fldint .eq. 3) then
               tot_N = tot_N + ar3*zga
             endif
          endif
        enddo
        ! endif ! for twocentfit
        tot_fld(fldint) = (tot_E+tot_N)
        if (.not. cube) then
           if (fldint .eq. 1) then
              int_type = 'X'
           elseif (fldint .eq. 2) then
              int_type = 'Y'
           elseif (fldint .eq. 3) then
              int_type = 'Z'
           endif
           write(6,1001) int_type,cx*AtoBohr,cy*AtoBohr,&
                         cz*AtoBohr,tot_fld(fldint)
           write(6,1002) int_type,cx*AtoBohr,cy*AtoBohr,&
                         cz*AtoBohr,tot_E
           write(6,1003) int_type,cx*AtoBohr,cy*AtoBohr,&
                         cz*AtoBohr,tot_N
           write(6,1004) 
        endif
     enddo ! fldint
  enddo ! n_esp_points
1001 format (1X,'EXACT Electric Field (',a1,') at',3(1X,F10.6),' = ',F11.6)
1002 format (1X,'EXACT Electronic EF (',a1,') at',3(1X,F10.6),' = ',F11.6)
1003 format (1X,'EXACT Nuclear EF (',a1,') at',3(1X,F10.6),' = ',F11.6)
1004 format (' ')
  
  close(13)
     
  return
end subroutine H_exact_fld
!---------------------------------------------------------
subroutine H_exact_dens(site_info,basis,cx,cy,cz,tot_dens,limitvec,&
                        natoms,num_stos,cube,twocentfit)
use definition
  implicit none
  type(sites_info)::site_info
  type(ao_orbitals)::basis(num_stos)

  integer site1,allochk,inttype,natoms
  integer k1,num1,off1,c_off1,npa,npb,nx1,ny1,nz1
  integer ncoeff,num_stos,nx2,ny2,nz2,limitvec(natoms)
  integer ires,jres,n,i,j,k,deg1,deg2,i2,j2,counter,nx3,ny3,nz3,l,m
  integer n_esp_points,slo1,shi1,lim1lo,lim1high
  double precision alphap,R_ab,brak,pre,temp,px,py,pz,norma,tot_dens
  double precision dx,dy,dz,cx,cy,cz,dist,nuc_chg1
  double precision x1,x2,y1,y2,z1,z2,x3,y3,z3
  double precision boyspar,result,beta,beta2,boyspar2,&
                   expon3,norm3,gaua,gaub,R1,R2
  double precision,allocatable::expon1(:),expon2(:),norm1(:),norm2(:),&
                   contr1(:),contr2(:)
  double precision D(0:3,0:3),E(0:3,0:3),F(0:3,0:3),&
                   DDS(0:3,0:3),EDS(0:3,0:3),FDS(0:3,0:3),&
                   densmat(num_stos,num_stos) 
  double precision,parameter::AtoBohr = 0.529177249d0,tiny=1.d-12
  logical normalize, debug, cube, off_site, twocentfit

  counter = 0
  do i = 1, num_stos
     do j = 1, i
        counter = counter+1
        densmat(j,i) = site_info%densmat(counter)
        densmat(i,j) = densmat(j,i)
     enddo
  enddo

  ! Just in case we want to calculate DENS at a single point
  inquire(file = 'esp_point', exist = off_site)
  if (off_site .and. .not. cube) then
     open(13,file='esp_point',status='unknown')
     read (13,*),n_esp_points
  else
     n_esp_points = 1
  endif
  do l = 1, n_esp_points
     tot_dens = 0.0d0
     if (off_site .and. .not. cube) then
        read (13,*),cx, cy, cz
        cx=cx/AtoBohr
        cy=cy/AtoBohr 
        cz=cz/AtoBohr
     endif
     counter = 1
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
                 print *,'H_exact_esp: something wrong with fit limits, &
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
        x1 = cx - basis(i)%coords(1)
        y1 = cy - basis(i)%coords(2)
        z1 = cz - basis(i)%coords(3)
        R1 = x1**2 + y1**2 + z1**2
        if (twocentfit) then
          do j = 1, lim1lo-1
            if ((j .lt. lim1lo) .or. (j .gt. lim1high))then
              npb = basis(j)%deg_contr  
              allocate(expon2(npb),contr2(npb),norm2(npb), stat = allochk)
              if (allochk .gt. 0 ) then
                 write(6,*)'H_nuc_elec:could not allocate expon2, contr2,&
                            &norm2; exiting'
                 stop
              endif     
              expon2(:) = basis(j)%expo(:)
              contr2(:) = basis(j)%contr(:)
              norm2(:) = basis(j)%norm(:)
              nx2 = basis(j)%x  
              ny2 = basis(j)%y 
              nz2 = basis(j)%z 
              x2 = cx - basis(j)%coords(1)
              y2 = cy - basis(j)%coords(2)
              z2 = cz - basis(j)%coords(3)
              R2 = x2**2 + y2**2 + z2**2
              temp = 0.0d0
              do i2 = 1, npa
                 do j2 = 1, npb
                    pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)
                    gaua = (x1**nx1)*(y1**ny1)*(z1**nz1)*&
                           exp(-expon1(i2)*R1)
                    gaub = (x2**nx2)*(y2**ny2)*(z2**nz2)*&
                           exp(-expon2(j2)*R2)
                    temp = temp+gaua*gaub*pre
                 enddo ! j2 = 1, npb
              enddo ! i2 = 1, npa
              if (i .eq. j) then
                 tot_dens = tot_dens + temp*densmat(i,j)
              else
                 tot_dens = tot_dens + temp*densmat(i,j)*2.0d0
              endif
              deallocate(expon2,contr2,norm2)
            endif ! do only 1 center part
          enddo ! j = 1, i
        else
           do j = 1, i
              npb = basis(j)%deg_contr  
              allocate(expon2(npb),contr2(npb),norm2(npb), stat = allochk)
              if (allochk .gt. 0 ) then
                 write(6,*)'H_nuc_elec:could not allocate expon2, contr2,&
                            &norm2; exiting'
                 stop
              endif     
              expon2(:) = basis(j)%expo(:)
              contr2(:) = basis(j)%contr(:)
              norm2(:) = basis(j)%norm(:)
              nx2 = basis(j)%x  
              ny2 = basis(j)%y 
              nz2 = basis(j)%z 
              x2 = cx - basis(j)%coords(1)
              y2 = cy - basis(j)%coords(2)
              z2 = cz - basis(j)%coords(3)
              R2 = x2**2 + y2**2 + z2**2
              temp = 0.0d0
              do i2 = 1, npa
                 do j2 = 1, npb
                    pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)
                    gaua = (x1**nx1)*(y1**ny1)*(z1**nz1)*&
                           exp(-expon1(i2)*R1)
                    gaub = (x2**nx2)*(y2**ny2)*(z2**nz2)*&
                           exp(-expon2(j2)*R2)
                    temp = temp+gaua*gaub*pre
                 enddo ! j2 = 1, npb
              enddo ! i2 = 1, npa
              if (i .eq. j) then
                 tot_dens = tot_dens + temp*densmat(i,j)
              else
                 tot_dens = tot_dens + temp*densmat(i,j)*2.0d0
              endif
              deallocate(expon2,contr2,norm2)
           enddo ! j = 1, num_stos
        endif ! for twocentfit
        deallocate(expon1,contr1,norm1)
     enddo ! i = 1, num_stos
     if (.not. cube)write(6,1003) cx*AtoBohr,cy*AtoBohr,cz*AtoBohr,tot_dens
  enddo ! n_esp_points
1003 format (1X,'EXACT DENS at',3(1X,F10.6),' = ',F20.15)
  
  close(13)
     
  return
end subroutine H_exact_dens
!----------------------------------------------------------
subroutine H_exact_esp2(site_info,basis,cx,cy,cz,tot_esp,num_stos,cube)
use definition
  implicit none
  type(sites_info)::site_info
  type(ao_orbitals)::basis(num_stos)

  integer site1,allochk,inttype,natoms
  integer k1,num1,off1,c_off1,npa,npb,nx1,ny1,nz1
  integer ncoeff,num_stos,nx2,ny2,nz2
  integer ires,jres,n,i,j,k,deg1,deg2,i2,j2,counter,nx3,ny3,nz3,l,m
  integer n_esp_points,slo1,shi1
  double precision alphap,R_ab,brak,pre,temp,px,py,pz,norma,tot_esp,tot_N
  double precision dx,dy,dz,cx,cy,cz,dist,nuc_chg1
  double precision x1,x2,y1,y2,z1,z2,x3,y3,z3
  double precision boyspar,result,beta,beta2,boyspar2,&
                   expon3,norm3
  double precision,allocatable::expon1(:),expon2(:),norm1(:),norm2(:),&
                   contr1(:),contr2(:)
  double precision D(0:3,0:3,0:6),E(0:3,0:3,0:6),F(0:3,0:3,0:6),&
                   DDS(0:3,0:3),EDS(0:3,0:3),FDS(0:3,0:3),&
                   densmat(num_stos,num_stos) 
  double precision,parameter::AtoBohr = 0.529177249d0,tiny=1.d-12
  logical normalize, debug, cube

  natoms = site_info%natoms
  nx3 = 0
  ny3 = 0
  nz3 = 0
  deg2 = 0
  DDS(:,:) = 1.0d0
  EDS(:,:) = 1.0d0
  FDS(:,:) = 1.0d0
  tot_esp = 0.0d0
  tot_N = 0.0d0

  do i = 1, num_stos
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
     do j = i, num_stos
        npb = basis(j)%deg_contr  
        allocate(expon2(npb),contr2(npb),norm2(npb), stat = allochk)
        if (allochk .gt. 0 ) then
           write(6,*)'H_nuc_elec:could not allocate expon2, contr2, norm2;&
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
              boyspar = alphap
              pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*brak
              dx = px - cx
              dy = py - cy
              dz = pz - cz
              expon3 = 3.14159265358979323846d0
              call GN_MD_REC_J_C(deg1,deg2,boyspar,dx,dy,dz,result,&
                                 nx1,ny1,nz1,nx2,ny2,nz2,nx3,ny3,nz3,&
                                 alphap,expon3,D,E,F,DDS,EDS,FDS,0)
              temp = temp+result*pre
           enddo ! j2 = 1, npb
        enddo ! i2 = 1, npa
        if (i .eq. j) then
           tot_esp = tot_esp - temp
        else
           tot_esp = tot_esp - temp*2.0d0
        endif
        deallocate(expon2,contr2,norm2)
     enddo ! j = 1, num_stos
     deallocate(expon1,contr1,norm1)
  enddo ! i = 1, num_stos

! now do nuclear ESP
  slo1 = site_info%residue_start(1)
  shi1 = site_info%residue_end(1)
  do site1 = slo1,shi1
    dx = site_info%site_crds(3*(site1-1)+1) - cx
    if (abs(dx) .lt. tiny) dx = 0.0d0
    dy = site_info%site_crds(3*(site1-1)+2) - cy
    if (abs(dy) .lt. tiny) dy = 0.0d0
    dz = site_info%site_crds(3*(site1-1)+3) - cz
    if (abs(dz) .lt. tiny) dz = 0.0d0
    dist = dsqrt(dx*dx+dy*dy+dz*dz)
    nuc_chg1 = site_info%nuclear_charge(site1)
    if(site1 .le. natoms) then
      if (dist .gt. tiny) tot_N = tot_N + nuc_chg1/dist
    endif
  enddo
  tot_esp = tot_esp+tot_N
  if (.not. cube)write(6,1003) cx*AtoBohr,cy*AtoBohr,cz*AtoBohr,tot_esp
1003 format (1X,'EXACT ESP at',3(1X,F10.6),' = ',F20.15)

  return
end subroutine H_exact_esp2
!---------------------------------------------------------
subroutine H_exact_fld2(site_info,basis,cx,cy,cz,tot_fld,&
                        num_stos,cube)
use definition
  implicit none
  type(sites_info)::site_info
  type(ao_orbitals)::basis(num_stos)

  integer site1,allochk,inttype,natoms
  integer k1,num1,off1,c_off1,npa,npb,nx1,ny1,nz1
  integer ncoeff,num_stos,nx2,ny2,nz2
  integer ires,jres,n,i,j,k,deg1,deg2,i2,j2,counter,nx3,ny3,nz3,l,m
  integer n_esp_points,slo1,shi1,fldint,lim1lo,lim1high
  double precision alphap,R_ab,brak,pre,temp,px,py,pz,norma,tot_fld(3),tot_N
  double precision dx,dy,dz,cx,cy,cz,dist,nuc_chg1,tot_E
  double precision x1,x2,y1,y2,z1,z2,x3,y3,z3,xga,yga,zga,ar3
  double precision boyspar,result,beta,beta2,boyspar2,&
                   expon3,norm3
  double precision,allocatable::expon1(:),expon2(:),norm1(:),norm2(:),&
                   contr1(:),contr2(:)
  double precision D(0:3,0:3,0:6),E(0:3,0:3,0:6),F(0:3,0:3,0:6),&
                   DDS(0:3,0:3),EDS(0:3,0:3),FDS(0:3,0:3),&
                   densmat(num_stos,num_stos) 
  double precision,parameter::AtoBohr = 0.529177249d0,tiny=1.d-15
  logical normalize, debug, cube, off_site
  character*1 int_type

  counter = 0
  do i = 1, num_stos
     do j = 1, i
        counter = counter+1
        densmat(j,i) = site_info%densmat(counter)
        densmat(i,j) = densmat(j,i)
     enddo
  enddo

  slo1 = site_info%residue_start(1)
  shi1 = site_info%residue_end(1)
  nx3 = 0
  ny3 = 0
  nz3 = 0
  deg2 = 0
  tot_fld(:) = 0.0d0
  DDS(:,:) = 1.0d0
  EDS(:,:) = 1.0d0
  FDS(:,:) = 1.0d0
  ! Just in case we want to calculate ESP at a single point
  inquire(file = 'esp_point', exist = off_site)
  if (off_site .and. .not. cube) then
     open(13,file='esp_point',status='unknown')
     read (13,*),n_esp_points
  else
     n_esp_points = 1
  endif
  do l = 1, n_esp_points
     if (off_site .and. .not. cube) then
        read (13,*),cx, cy, cz
        cx=cx/AtoBohr
        cy=cy/AtoBohr 
        cz=cz/AtoBohr
     endif
     do fldint = 1, 3
        tot_E = 0.0d0
        tot_N = 0.0d0
        counter = 1
        do i = 1, num_stos
           npa = basis(i)%deg_contr  
           allocate(expon1(npa),contr1(npa),norm1(npa), stat = allochk)
           if (allochk .gt. 0 ) then
              write(6,*)'H_form_Jvec:could not allocate expon1, contr1,&
                         & norm1; exiting'
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
           do j = 1, i
              npb = basis(j)%deg_contr  
              allocate(expon2(npb),contr2(npb),norm2(npb), stat = allochk)
              if (allochk .gt. 0 ) then
                 write(6,*)'H_nuc_elec:could not allocate expon2, contr2,&
                            & norm2; exiting'
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
                    boyspar = alphap
                    pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*brak
                    dx = px - cx
                    dy = py - cy
                    dz = pz - cz
                    expon3 = 3.14159265358979323846d0
                    call GN_MD_REC_J_C(deg1,deg2,boyspar,dx,dy,dz,result,&
                                       nx1,ny1,nz1,nx2,ny2,nz2,nx3,ny3,nz3,&
                                       alphap,expon3,D,E,F,DDS,EDS,FDS,&
                                       fldint)
                    temp = temp+result*pre
                 enddo ! j2 = 1, npb
              enddo ! i2 = 1, npa
              if (i .eq. j) then
                 tot_E = tot_E + temp*densmat(i,j)
              else
                 tot_E = tot_E + temp*densmat(i,j)*2.0d0
              endif
              deallocate(expon2,contr2,norm2)
           enddo ! j = 1, num_stos
           deallocate(expon1,contr1,norm1)
        enddo ! i = 1, num_stos

      ! now do nuclear FLD
        do site1 = slo1,shi1
          dx = site_info%site_crds(3*(site1-1)+1) - cx
          if (abs(dx) .lt. tiny) dx = 0.0d0
          dy = site_info%site_crds(3*(site1-1)+2) - cy
          if (abs(dy) .lt. tiny) dy = 0.0d0
          dz = site_info%site_crds(3*(site1-1)+3) - cz
          if (abs(dz) .lt. tiny) dz = 0.0d0
          dist = dx*dx + dy*dy + dz*dz
          nuc_chg1 = site_info%nuclear_charge(site1)
          if (dist .gt. tiny) then
             dist = 1.0d0/dist
             xga = cx - site_info%site_crds(3*(site1-1)+1)
             yga = cy - site_info%site_crds(3*(site1-1)+2)
             zga = cz - site_info%site_crds(3*(site1-1)+3)
             !xga = site_info%site_crds(3*(site1-1)+1)-cx
             !yga = site_info%site_crds(3*(site1-1)+2)-cy
             !zga = site_info%site_crds(3*(site1-1)+3)-cz
             ar3 = nuc_chg1 * dsqrt(dist) * dist
             if(fldint .eq. 1) then
               tot_N = tot_N + ar3*xga
             else if(fldint .eq. 2) then
               tot_N = tot_N + ar3*yga
             else if(fldint .eq. 3) then
               tot_N = tot_N + ar3*zga
             endif
          endif
        enddo
        tot_fld(fldint) = (tot_E+tot_N)
        if (.not. cube) then
           if (fldint .eq. 1) then
              int_type = 'X'
           elseif (fldint .eq. 2) then
              int_type = 'Y'
           elseif (fldint .eq. 3) then
              int_type = 'Z'
           endif
           write(6,1001) int_type,cx*AtoBohr,cy*AtoBohr,&
                         cz*AtoBohr,tot_fld(fldint)
        endif
     enddo ! fldint
  enddo ! n_esp_points
1001 format (1X,'EXACT Electric Field (',a1,') at',3(1X,F10.6),' = ',F11.6)
  
  close(13)
     
  return
end subroutine H_exact_fld2
!---------------------------------------------------------
!----------------------------------------------------------
subroutine H_exact_esp3(site_info,basis,cx,cy,cz,tot_esp,&
                        num_stos,cube)
use definition
  implicit none
  type(sites_info)::site_info
  type(ao_orbitals)::basis(num_stos)

  integer site1,allochk,inttype,natoms
  integer k1,num1,off1,c_off1,npa,npb,nx1,ny1,nz1
  integer ncoeff,num_stos,nx2,ny2,nz2
  integer ires,jres,n,i,j,k,deg1,deg2,i2,j2,counter,nx3,ny3,nz3,l,m
  integer n_esp_points,slo1,shi1,lim1lo,lim1high
  double precision alphap,R_ab,brak,pre,temp,px,py,pz,norma,tot_esp,tot_N
  double precision dx,dy,dz,cx,cy,cz,dist,nuc_chg1
  double precision x1,x2,y1,y2,z1,z2,x3,y3,z3
  double precision boyspar,result,beta,beta2,boyspar2,&
                   expon3,norm3
  double precision,allocatable::expon1(:),expon2(:),norm1(:),norm2(:),&
                   contr1(:),contr2(:)
  double precision D(0:3,0:3,0:6),E(0:3,0:3,0:6),F(0:3,0:3,0:6),&
                   DDS(0:3,0:3),EDS(0:3,0:3),FDS(0:3,0:3),&
                   densmat(num_stos,num_stos) 
  double precision,parameter::AtoBohr = 0.529177249d0,tiny=1.d-12
  logical normalize, debug, cube, off_site

  counter = 0
  do i = 1, num_stos
     do j = 1, i
        counter = counter+1
        densmat(j,i) = site_info%densmat(counter)
        densmat(i,j) = densmat(j,i)
     enddo
  enddo
 
  natoms = site_info%natoms
  nx3 = 0
  ny3 = 0
  nz3 = 0
  deg2 = 0
  DDS(:,:) = 1.0d0
  EDS(:,:) = 1.0d0
  FDS(:,:) = 1.0d0
  ! Just in case we want to calculate ESP at a single point
  inquire(file = 'esp_point', exist = off_site)
  if (off_site .and. .not. cube) then
     open(13,file='esp_point',status='unknown')
     read (13,*),n_esp_points
  else
     n_esp_points = 1
  endif
  do l = 1, n_esp_points
     tot_esp = 0.0d0
     tot_N = 0.0d0
     if (off_site .and. .not. cube) then
        read (13,*),cx, cy, cz
        cx=cx/AtoBohr
        cy=cy/AtoBohr 
        cz=cz/AtoBohr
     endif
     counter = 1
     do i = 1, num_stos
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
           do j = 1, i
              npb = basis(j)%deg_contr  
              allocate(expon2(npb),contr2(npb),norm2(npb), stat = allochk)
              if (allochk .gt. 0 ) then
                 write(6,*)'H_nuc_elec:could not allocate expon2, contr2, &
                            &norm2; exiting'
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
                    boyspar = alphap
                    pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*brak
                    dx = px - cx
                    dy = py - cy
                    dz = pz - cz
                    expon3 = 3.14159265358979323846d0
                    call GN_MD_REC_J_C(deg1,deg2,boyspar,dx,dy,dz,result,&
                                       nx1,ny1,nz1,nx2,ny2,nz2,nx3,ny3,nz3,&
                                       alphap,expon3,D,E,F,DDS,EDS,FDS,0)
                    temp = temp+result*pre
                 enddo ! j2 = 1, npb
              enddo ! i2 = 1, npa
              if (i .eq. j) then
                 tot_esp = tot_esp - temp*densmat(i,j)
              else
                 tot_esp = tot_esp - temp*densmat(i,j)*2.0d0
              endif
              deallocate(expon2,contr2,norm2)
           enddo ! j = 1, i
        deallocate(expon1,contr1,norm1)
     enddo ! i = 1, num_stos
   ! now do nuclear ESP
     slo1 = site_info%residue_start(1)
     shi1 = site_info%residue_end(1)
     do site1 = slo1,shi1
       dx = site_info%site_crds(3*(site1-1)+1) - cx
       if (abs(dx) .lt. tiny) dx = 0.0d0
       dy = site_info%site_crds(3*(site1-1)+2) - cy
       if (abs(dy) .lt. tiny) dy = 0.0d0
       dz = site_info%site_crds(3*(site1-1)+3) - cz
       if (abs(dz) .lt. tiny) dz = 0.0d0
       dist = dsqrt(dx*dx+dy*dy+dz*dz)
       nuc_chg1 = site_info%nuclear_charge(site1)
       if(site1 .le. natoms) then
         if (dist .gt. tiny) tot_N = tot_N + nuc_chg1/dist
       endif
     enddo
     tot_esp = tot_esp+tot_N
   ! endif ! for twocentfit
     if (.not. cube)write(6,1003) cx*AtoBohr,cy*AtoBohr,cz*AtoBohr,tot_esp
  enddo ! n_esp_points
1003 format (1X,'EXACT ESP at',3(1X,F10.6),' = ',F20.10)
  
  close(13)
     
  return
end subroutine H_exact_esp3
!---------------------------------------------------------
subroutine H_exact_dens2(site_info,basis,cx,cy,cz,tot_dens,num_stos,cube)
use definition
  implicit none
  type(sites_info)::site_info
  type(ao_orbitals)::basis(num_stos)

  integer site1,allochk,inttype
  integer k1,num1,off1,c_off1,npa,npb,nx1,ny1,nz1
  integer ncoeff,num_stos,nx2,ny2,nz2
  integer ires,jres,n,i,j,k,deg1,deg2,i2,j2,counter,nx3,ny3,nz3,l,m
  integer n_esp_points,slo1,shi1,lim1lo,lim1high
  double precision alphap,R_ab,brak,pre,temp,px,py,pz,norma,tot_dens
  double precision dx,dy,dz,cx,cy,cz,dist,nuc_chg1
  double precision x1,x2,y1,y2,z1,z2,x3,y3,z3
  double precision boyspar,result,beta,beta2,boyspar2,&
                   expon3,norm3,gaua,gaub,R1,R2
  double precision,allocatable::expon1(:),expon2(:),norm1(:),norm2(:),&
                   contr1(:),contr2(:)
  double precision D(0:3,0:3),E(0:3,0:3),F(0:3,0:3),&
                   DDS(0:3,0:3),EDS(0:3,0:3),FDS(0:3,0:3),&
                   densmat(num_stos,num_stos) 
  double precision,parameter::AtoBohr = 0.529177249d0,tiny=1.d-12
  logical normalize, debug, cube, off_site, twocentfit

  counter = 0
  do i = 1, num_stos
     do j = 1, i
        counter = counter+1
        densmat(j,i) = site_info%densmat(counter)
        densmat(i,j) = densmat(j,i)
     enddo
  enddo

  ! Just in case we want to calculate DENS at a single point
  inquire(file = 'esp_point', exist = off_site)
  if (off_site .and. .not. cube) then
     open(13,file='esp_point',status='unknown')
     read (13,*),n_esp_points
  else
     n_esp_points = 1
  endif
  do l = 1, n_esp_points
     tot_dens = 0.0d0
     if (off_site .and. .not. cube) then
        read (13,*),cx, cy, cz
        cx=cx/AtoBohr
        cy=cy/AtoBohr 
        cz=cz/AtoBohr
     endif
     counter = 1
     do i = 1, num_stos
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
        x1 = cx - basis(i)%coords(1)
        y1 = cy - basis(i)%coords(2)
        z1 = cz - basis(i)%coords(3)
        R1 = x1**2 + y1**2 + z1**2
        do j = 1, i
           npb = basis(j)%deg_contr  
           allocate(expon2(npb),contr2(npb),norm2(npb), stat = allochk)
           if (allochk .gt. 0 ) then
              write(6,*)'H_nuc_elec:could not allocate expon2, contr2,&
                         &norm2; exiting'
              stop
           endif     
           expon2(:) = basis(j)%expo(:)
           contr2(:) = basis(j)%contr(:)
           norm2(:) = basis(j)%norm(:)
           nx2 = basis(j)%x  
           ny2 = basis(j)%y 
           nz2 = basis(j)%z 
           x2 = cx - basis(j)%coords(1)
           y2 = cy - basis(j)%coords(2)
           z2 = cz - basis(j)%coords(3)
           R2 = x2**2 + y2**2 + z2**2
           temp = 0.0d0
           do i2 = 1, npa
              do j2 = 1, npb
                 pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)
                 gaua = (x1**nx1)*(y1**ny1)*(z1**nz1)*&
                        exp(-expon1(i2)*R1)
                 gaub = (x2**nx2)*(y2**ny2)*(z2**nz2)*&
                        exp(-expon2(j2)*R2)
                 temp = temp+gaua*gaub*pre
              enddo ! j2 = 1, npb
           enddo ! i2 = 1, npa
           if (i .eq. j) then
              tot_dens = tot_dens + temp*densmat(i,j)
           else
              tot_dens = tot_dens + temp*densmat(i,j)*2.0d0
           endif
           deallocate(expon2,contr2,norm2)
        enddo ! j = 1, num_stos
        deallocate(expon1,contr1,norm1)
     enddo ! i = 1, num_stos
     if (.not. cube)write(6,1003) cx*AtoBohr,cy*AtoBohr,cz*AtoBohr,tot_dens
  enddo ! n_esp_points
1003 format (1X,'EXACT DENS at',3(1X,F10.6),' = ',F20.15)
  
  close(13)
     
  return
end subroutine H_exact_dens2
!----------------------------------------------------------
