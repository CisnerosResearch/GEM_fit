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
subroutine H_esp_aux(auxis,cx,cy,cz,tot_E,coeff)
  ! GAC: subroutine to calculate ESP with cartesian gaussians.
use definition
  implicit none
  type(aux_orbitals)::auxis(*)

  integer site1,i,x1,y1,z1,deg1,deg2,coeff
  integer k1,num1,off1,c_off1,l
  integer order1,order2,allochk
  integer fldint
  double precision tot_E,total,result,norm1
  double precision dx,dy,dz,cx,cy,cz,dist
  double precision expon1,expon2,boyspar,nuc_chg1,nuc_chg2,ene,sys_ene
  double precision D(0:3,0:3),E(0:3,0:3),F(0:3,0:3),&
                   DDS(0:3,0:3),EDS(0:3,0:3),FDS(0:3,0:3)
  logical off_site,cartfit

  fldint = 0

  DDS(:,:) = 1.0d0
  EDS(:,:) = 1.0d0
  FDS(:,:) = 1.0d0

  tot_E = 0.0d0
  i = coeff
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
  tot_E = tot_E - result*norm1
 
  return
end subroutine H_esp_aux
!----------------------------------------------------------
subroutine H_esp_aux2(auxis,cx,cy,cz,tot_E,coeff)
  ! GAC: subroutine to calculate ESP with cartesian gaussians.
use definition
  implicit none
  type(aux_orbitals)::auxis(*)

  integer site1,i,x1,y1,z1,deg1,deg2,coeff
  integer k1,num1,off1,c_off1,l
  integer order1,order2,allochk
  integer fldint
  double precision tot_E,total,result,norm1
  double precision dx,dy,dz,cx,cy,cz,dist
  double precision expon1,expon2,boyspar,nuc_chg1,nuc_chg2,ene,sys_ene
  !double precision D(0:3,0:3),E(0:3,0:3),F(0:3,0:3),&
  !                 DDS(0:3,0:3),EDS(0:3,0:3),FDS(0:3,0:3)
  logical off_site,cartfit

  !fldint = 0

  !D(:,:) = 1.0d0
  !E(:,:) = 1.0d0
  !F(:,:) = 1.0d0
  !DDS(:,:) = 1.0d0
  !EDS(:,:) = 1.0d0
  !FDS(:,:) = 1.0d0

  tot_E = 0.0d0
  i = coeff
  expon1 = auxis(i)%expo
  x1 = auxis(i)%x  
  y1 = auxis(i)%y 
  z1 = auxis(i)%z 
  deg1 = x1+y1+z1 
  !norm1 = auxis(i)%norm
  !call H_form_exp_coef_noS(x1,auxis(i)%coords_G(1),expon1,D)
  !call H_form_exp_coef_noS(y1,auxis(i)%coords_G(2),expon1,E)
  !call H_form_exp_coef_noS(z1,auxis(i)%coords_G(3),expon1,F)
  dx = auxis(i)%coords_G(1) - cx
  dy = auxis(i)%coords_G(2) - cy
  dz = auxis(i)%coords_G(3) - cz
  boyspar = expon1
  deg2 = 0
  !expon2 = 3.14159265358979323846d0
  !call GN_MD_REC_C(deg1,deg2,boyspar,dx,dy,dz,result,x1,y1,z1,&
  !               0,0,0,expon1,expon2,D,E,F,DDS,EDS,FDS,fldint)
  call GN_MD_REC(deg2,deg1,boyspar,dx,dy,dz,result,0,0,0,&
                 x1,y1,z1)
  tot_E = tot_E + result!*norm1
 
  return
end subroutine H_esp_aux2
!----------------------------------------------------------
subroutine H_esp_aux3(site_info,site_info2,tot_NE)
  ! GAC: subroutine to calculate ESP with hermite gaussians.
use definition
  implicit none
  type(sites_info)::site_info,site_info2

  integer site1,site2
  integer k1,num1,off1,c_off1,k2,num2,off2,c_off2
  integer order1,order2,i
  integer num_res,slo1,shi1,slo2,shi2,ires,jres,num_res2
  double precision prim_prim_EE,tot_NN,tot_EE,tot_NE
  double precision site_site_EE,site_site_NN,site_site_NE,ss,sd,ds,dd
  double precision sp, ps, pp, pd, dp, sd2, d2s, d2d2, d2d, dd2 ! GAC 
  double precision dx,dy,dz,tmp_vec(100)
  double precision expon1,expon2,boyspar,nuc_chg1,nuc_chg2,hartree,ene,sys_ene

  slo1 = site_info%residue_start(1)
  shi1 = site_info%natoms

  tot_NE = 0.d0
  do jres = 1,site_info2%num_residues
     slo2 = site_info2%residue_start(jres)
     shi2 = site_info2%residue_end(jres)
     do site1 = slo1,shi1
       !num1 = site_info%num_primitives(site1)
       !off1 = site_info%off_primitives(site1)
       do site2 = slo2,shi2
         num2 = site_info2%num_primitives(site2)
         off2 = site_info2%off_primitives(site2)
         dx = site_info2%site_crds(3*(site2-1)+1) - &
                        site_info%site_crds(3*(site1-1)+1)
         dy = site_info2%site_crds(3*(site2-1)+2) - &
                        site_info%site_crds(3*(site1-1)+2)
         dz = site_info2%site_crds(3*(site2-1)+3) - &
                        site_info%site_crds(3*(site1-1)+3)
         nuc_chg1 = site_info%nuclear_charge(site1)
         site_site_NE = 0.d0
         do k2 = 1,num2
           order2 = site_info2%prim_order(off2+k2)
           expon2 = site_info2%prim_expo(off2+k2)
           c_off2 = site_info2%coeff_offset(off2+k2)
           boyspar = expon2
           if ( order2 == 1 )then
             call GN_MCMUR_DAV_ss(boyspar,dx,dy,dz,  &
                 nuc_chg1,site_info2%global_hermite_coeffs(c_off2+1),ene,&
                   tmp_vec)
           elseif ( order2 == 4 )then ! GAC
             call GN_MCMUR_DAV_sp(boyspar,dx,dy,dz,  &
                 nuc_chg1,site_info2%global_hermite_coeffs(c_off2+1),ene,&
                   tmp_vec)
           elseif ( order2 == 6 )then ! GAC
             call GN_MCMUR_DAV_sda(boyspar,dx,dy,dz,  &
                 nuc_chg1,site_info2%global_hermite_coeffs(c_off2+1),ene,&
                   tmp_vec)
           elseif ( order2 == 10 )then
             call GN_MCMUR_DAV_sd(boyspar,dx,dy,dz,  &
                 nuc_chg1,site_info2%global_hermite_coeffs(c_off2+1),ene,&
                   tmp_vec)
           endif
           site_site_NE = site_site_NE + ene
         enddo
         tot_NE = tot_NE + site_site_NE
       enddo!site2 = slo2,shi2
     enddo!site1 = slo1,slo2
  enddo ! jres = 1, nres
    
  return
end subroutine H_esp_aux3
!----------------------------------------------------------
