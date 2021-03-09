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
!---------------------------------------------------------
subroutine H_calc_esp_CART(site_info,auxis,cx,cy,cz,total,ncoeff,cartfit,cube)
  ! GAC: subroutine to calculate ESP with cartesian gaussians.
use definition
  implicit none
  type(sites_info)::site_info
  type(aux_orbitals)::auxis(ncoeff)

  integer site1,i,x1,y1,z1,deg1,deg2,ncoeff
  integer k1,num1,off1,c_off1,n_esp_points,l
  integer order1,order2,allochk
  integer slo1,shi1,ires,natoms,nsites,fldint
  double precision prim_prim_EE,tot_esp,tot_E,tot_N,total,result,norm1
  double precision dx,dy,dz,cx,cy,cz,dist,tmp_vec(10)
  double precision expon1,expon2,boyspar,nuc_chg1,nuc_chg2,ene,sys_ene
  !double precision,allocatable::elec_esp(:),nuc_esp(:)
  double precision,allocatable::nuc_esp(:)
  double precision D(0:3,0:3),E(0:3,0:3),F(0:3,0:3),&
                   DDS(0:3,0:3),EDS(0:3,0:3),FDS(0:3,0:3)
  double precision,parameter::AtoBohr = 0.529177249d0,hartree=627.51d0,&
                   tiny=1.d-12
  logical off_site,cartfit,cube

  natoms = site_info%natoms
  slo1 = site_info%residue_start(1)
  shi1 = site_info%residue_end(1)
  fldint = 0

  !allocate(elec_esp(nsites),nuc_esp(natoms), stat = allochk)
  allocate(nuc_esp(natoms), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'H_calc_esp:could not allocate elec_esp, nuc_esp; exiting'
     stop
  endif

  n_esp_points = 1
  ! Just in case we want to calculate ESP at a single point
  inquire(file = 'esp_point', exist = off_site)
  if (off_site .and. .not. cube) then
     open(13,file='esp_point',status='unknown')
     read (13,*),n_esp_points
  else
     n_esp_points = 1
  endif

  DDS(:,:) = 1.0d0
  EDS(:,:) = 1.0d0
  FDS(:,:) = 1.0d0
  do l = 1, n_esp_points

  if (off_site .and. .not. cube) then
     read (13,*),cx, cy, cz
     cx=cx/AtoBohr
     cy=cy/AtoBohr 
     cz=cz/AtoBohr
  endif

  !elec_esp(:) = 0.0d0
  nuc_esp(:) = 0.0d0
  tot_N = 0.0d0
  tot_E = 0.0d0
  do i = 1, ncoeff
     expon1 = auxis(i)%expo
     x1 = auxis(i)%x  
     y1 = auxis(i)%y 
     z1 = auxis(i)%z 
     deg1 = x1+y1+z1 
     norm1 = auxis(i)%norm
     call H_form_exp_coef_noS(x1,auxis(i)%coords_G(1),expon1,D)
     call H_form_exp_coef_noS(y1,auxis(i)%coords_G(2),expon1,E)
     call H_form_exp_coef_noS(z1,auxis(i)%coords_G(3),expon1,F)
     dx = auxis(i)%coords_G(1) - cx
     dy = auxis(i)%coords_G(2) - cy
     dz = auxis(i)%coords_G(3) - cz
     boyspar = expon1
     deg2 = 0
     expon2 = 3.14159265358979323846d0
     call GN_MD_REC_C(deg1,deg2,boyspar,dx,dy,dz,result,x1,y1,z1,&
                      0,0,0,expon1,expon2,D,E,F,DDS,EDS,FDS,fldint)
     !elec_esp(site1) = elec_esp(site1) - result*norm1*&
     !                                    site_info%cartesian_coeffs(i)
     if (cartfit) then
        tot_E = tot_E - result*norm1*site_info%cartesian_coeffs(i)
     else
        tot_E = tot_E - result*site_info%global_hermite_coeffs(i)*&
                        (expon1/expon2)*sqrt(expon1/expon2)
     endif
  enddo
  !if(n_esp_points .eq. 1.and..not.cube)write(6,*),'total electronic = ',tot_E
 
! now do nuclear ESP
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
      if (dist .gt. tiny) nuc_esp(site1) = nuc_chg1 / dist
      tot_N = tot_N + nuc_esp(site1)
    endif
  enddo
  
  total = (tot_E + tot_N) 
  !do site1 = 1, nsites
  !   write(6,1001) site1,elec_esp(site1)
  !enddo
  if (.not. cube) then
     !do site1 = 1, natoms
     !   if(n_esp_points .eq. 1)write(6,1002) site1,nuc_esp(site1)
     !enddo
     write(6,1003) cx*AtoBohr,cy*AtoBohr,cz*AtoBohr,total
  endif

  enddo ! n_esp_points

  close(13)

1001 format (' Electronic ESP of Atom ',i3,' is ',f14.6)
1002 format (' Nuclear ESP of Atom ',i3,' is ',f14.6)
1003 format (1X,'Total ESP at',3(1X,F10.6),' = ',F20.10)

  return
end subroutine H_calc_esp_CART
!----------------------------------------------------------
