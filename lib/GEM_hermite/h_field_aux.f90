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
subroutine H_fld_aux(auxis,cx,cy,cz,total,coeff)
  ! GAC: subroutine to calculate ESP with cartesian gaussians.
use definition
  implicit none
  type(aux_orbitals)::auxis(*)

  integer site1,i,x1,y1,z1,deg1,deg2,coeff,dim
  integer k1,num1,off1,c_off1,n_esp_points,l
  integer order1,order2,allochk
  integer shi1,fldint
  double precision tot_E,result,norm1
  double precision dx,dy,dz,cx,cy,cz,dist,total(3)
  double precision expon1,expon2,boyspar
  double precision D(0:3,0:3),E(0:3,0:3),F(0:3,0:3),&
                   DDS(0:3,0:3),EDS(0:3,0:3),FDS(0:3,0:3)
  double precision,allocatable::nuc_fld(:,:)
  logical off_site,cartfit
  character*1 int_type

  total(:) = 0.0d0
  DDS(:,:) = 1.0d0
  EDS(:,:) = 1.0d0
  FDS(:,:) = 1.0d0

  i = coeff
  do fldint = 1, 3
     tot_E = 0.0d0
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
     tot_E = tot_E + result*norm1
     total(fldint) = tot_E
  enddo !fldint
  return
end subroutine H_fld_aux
!----------------------------------------------------------
