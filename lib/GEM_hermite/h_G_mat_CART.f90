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
subroutine AH_form_G_mat_CART(site_info,auxis,Gmat,ncoeff,inttype,beta,&
                              debug,cube,densfile)
  ! GAC: subroutine to calculate <k||l> and <k|l> in local hermite.
use definition
  implicit none
  type(sites_info)::site_info
  type(aux_orbitals)::auxis(ncoeff)

  integer site1,allochk,inttype
  integer k1,num1,off1,c_off1 
  integer ncoeff,x1,x2,y1,y2,z1,z2
  integer ires,jres,n,i,j,deg1,deg2,fldint
  double precision norm1,norm2 
  double precision dx,dy,dz,Gmat(ncoeff,ncoeff)
  double precision expon1,expon2,boyspar,result,beta,beta2,boyspar2
  double precision D(0:3,0:3),E(0:3,0:3),F(0:3,0:3),&
                   DDS(0:3,0:3),EDS(0:3,0:3),FDS(0:3,0:3)
  character(len=*) densfile
  logical debug, cube

  Gmat(:,:) = 0.d0
  if(inttype == 1)beta2 = beta*beta
  if(debug)write(12,*)'---------------- Gmat ----------------------'

  fldint = 0
  if(.not. cube)call AH_form_auxis(site_info,auxis,ncoeff,densfile)

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
     do j = 1, i
        expon2 = auxis(j)%expo
        x2 = auxis(j)%x  
        y2 = auxis(j)%y 
        z2 = auxis(j)%z 
        deg2 = x2+y2+z2 
        norm2 = auxis(j)%norm
        call H_form_exp_coef_noS(x2,auxis(j)%coords_G(1),expon2,DDS)
        call H_form_exp_coef_noS(y2,auxis(j)%coords_G(2),expon2,EDS)
        call H_form_exp_coef_noS(z2,auxis(j)%coords_G(3),expon2,FDS)
        dx = auxis(i)%coords_G(1) - auxis(j)%coords_G(1)
        dy = auxis(i)%coords_G(2) - auxis(j)%coords_G(2)
        dz = auxis(i)%coords_G(3) - auxis(j)%coords_G(3)
        boyspar = expon1*expon2 / (expon1 + expon2)
        if (inttype == 0) then
           call GN_MD_REC_C(deg1,deg2,boyspar,dx,dy,dz,result,x1,y1,z1,&
                            x2,y2,z2,expon1,expon2,D,E,F,DDS,EDS,FDS,fldint)
        elseif (inttype == 1) then
           boyspar2 = (beta2*expon1*expon2)/(beta2*(expon1+expon2)+&
                                             expon1*expon2)
           call GN_MD_REC_erfc_C(deg1,deg2,boyspar,dx,dy,dz,result,&
                              x1,y1,z1,x2,y2,z2,expon1,expon2,boyspar2,&
                              D,E,F,DDS,EDS,FDS)
        elseif (inttype == 2) then
           call GN_MD_REC_OVERLAP_C(deg1,deg2,boyspar,dx,dy,dz,result,&
                                 x1,y1,z1,x2,y2,z2,expon1,expon2,D,E,F,&
                                 DDS,EDS,FDS)
        endif
        Gmat(i,j) = result*norm1*norm2
        if(debug)write(12,*)result,result*norm1*norm2
     enddo
  enddo
  
! symmetrize
  do i = 1, ncoeff
     do j = 1, ncoeff
        Gmat(i,j) = Gmat(j,i)
        if (debug) write(12,1001)i,j,Gmat(i,j)
     enddo
  enddo
1001    format ('Gmat (',i4,',',i4,') = ',f27.20)

  return
end subroutine AH_form_G_mat_CART
!----------------------------------------------------------
subroutine AH_form_auxis(site_info,auxis,ncoeff,densfile)
  ! GAC: subroutine to form ordered auxiliary orbitals for G matrix
use definition
  implicit none
  integer ncoeff
  type(sites_info)::site_info
  type(aux_orbitals)::auxis(ncoeff)

  integer site1,allochk, order1, scan_file
  integer k1,num1,off1,c_off1 !,t1,t2
  integer counter, res, c1, c2
  integer slo1,shi1,ires,jres,n,i,j
  double precision expon1,expon2
  character(len=*) densfile
  character(len=4) tmp_name
  logical nw_dens

! check if it's a NWChem density file
  nw_dens = .false.
  scan_file = 0
  tmp_name = 'NW'
  scan_file = index(densfile,tmp_name)
  if (scan_file .ne. 0) nw_dens = .true.
  if(nw_dens)write(6,*)'AH_form_auxis: NWChem density found, reordering d,f'

  res = site_info%residue_num(1)
  slo1 = site_info%residue_start(res)
  shi1 = site_info%residue_end(res)
  counter = 0
  do site1 = slo1,shi1
    num1 = site_info%num_primitives(site1)
    off1 = site_info%off_primitives(site1)
    do k1 = 1,num1
      order1 = site_info%prim_order(off1+k1)
      expon1 = site_info%prim_expo(off1+k1)
      c_off1 = site_info%coeff_offset(off1+k1)
      if (order1 == 1) then ! s shell
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 0
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+1)
      elseif (order1 == 4) then ! sp shell
         ! *** s orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 0
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+1)
         ! *** px orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 1
         auxis(counter)%y = 0
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+2)
         ! *** py orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 1
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+3)
         ! *** pz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 0
         auxis(counter)%z = 1
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+4)
      elseif (order1 == 6) then ! d only shell
       if (.not.nw_dens) then
         ! *** dxx orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 2
         auxis(counter)%y = 0
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+1)
         ! *** dyy orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 2
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+2)
         ! *** dzz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 0
         auxis(counter)%z = 2
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+3)
         ! *** dxy orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 1
         auxis(counter)%y = 1
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+4)
         ! *** dxz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 1
         auxis(counter)%y = 0
         auxis(counter)%z = 1
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+5)
         ! *** dyz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 1
         auxis(counter)%z = 1
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+6)
       else
         ! *** dxx orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 2
         auxis(counter)%y = 0
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+1)
         ! *** dxy orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 1
         auxis(counter)%y = 1
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+4)
         ! *** dxz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 1
         auxis(counter)%y = 0
         auxis(counter)%z = 1
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+5)
         ! *** dyy orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 2
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+2)
         ! *** dyz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 1
         auxis(counter)%z = 1
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+6)
         ! *** dzz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 0
         auxis(counter)%z = 2
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+3)
       endif
      elseif (order1 == 10) then ! spd shell
         ! *** s orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 0
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+1)
         ! *** px orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 1
         auxis(counter)%y = 0
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+2)
         ! *** py orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 1
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+3)
         ! *** pz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 0
         auxis(counter)%z = 1
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+4)
       if(.not.nw_dens)then
         ! *** dxx orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 2
         auxis(counter)%y = 0
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+5)
         ! *** dyy orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 2
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+6)
         ! *** dzz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 0
         auxis(counter)%z = 2
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+7)
         ! *** dxy orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 1
         auxis(counter)%y = 1
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+8)
         ! *** dxz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 1
         auxis(counter)%y = 0
         auxis(counter)%z = 1
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+9)
         ! *** dyz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 1
         auxis(counter)%z = 1
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+10)
       else
         ! *** dxx orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 2
         auxis(counter)%y = 0
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+5)
         ! *** dxy orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 1
         auxis(counter)%y = 1
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+6)
         ! *** dxz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 1
         auxis(counter)%y = 0
         auxis(counter)%z = 1
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+7)
         ! *** dyy orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 2
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+8)
         ! *** dyz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 1
         auxis(counter)%z = 1
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+9)
         ! *** dzz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 0
         auxis(counter)%z = 2
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+10)
       endif
      elseif (order1 == 20) then ! f shell
         ! *** s orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 0
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+1)
         ! *** px orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 1
         auxis(counter)%y = 0
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+2)
         ! *** py orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 1
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+3)
         ! *** pz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 0
         auxis(counter)%z = 1
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+4)
       if(.not.nw_dens)then
         ! *** dxx orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 2
         auxis(counter)%y = 0
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+5)
         ! *** dyy orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 2
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+6)
         ! *** dzz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 0
         auxis(counter)%z = 2
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+7)
         ! *** dxy orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 1
         auxis(counter)%y = 1
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+8)
         ! *** dxz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 1
         auxis(counter)%y = 0
         auxis(counter)%z = 1
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+9)
         ! *** dyz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 1
         auxis(counter)%z = 1
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+10)
       else
         ! *** dxx orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 2
         auxis(counter)%y = 0
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+5)
         ! *** dxy orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 1
         auxis(counter)%y = 1
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+6)
         ! *** dxz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 1
         auxis(counter)%y = 0
         auxis(counter)%z = 1
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+7)
         ! *** dyy orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 2
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+8)
         ! *** dyz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 1
         auxis(counter)%z = 1
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+9)
         ! *** dzz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 0
         auxis(counter)%z = 2
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+10)
       endif
       if (.not.nw_dens) then
         ! *** fxxx orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 3
         auxis(counter)%y = 0
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+1)
         ! *** fyyy orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 3
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+2)
         ! *** fzzz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 0
         auxis(counter)%z = 3
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+3)
         ! *** fxxy orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 2
         auxis(counter)%y = 1
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+4)
         ! *** fxxz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 2
         auxis(counter)%y = 0
         auxis(counter)%z = 1
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+5)
         ! *** fxyy orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 1
         auxis(counter)%y = 2
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+6)
         ! *** fyyz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 2
         auxis(counter)%z = 1
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+7)
         ! *** fxzz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 1
         auxis(counter)%y = 0
         auxis(counter)%z = 2
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+8)
         ! *** fyzz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 1
         auxis(counter)%z = 2
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+9)
         ! *** fxyz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 1
         auxis(counter)%y = 1
         auxis(counter)%z = 1
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+10)
       else
         ! *** fxxx orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 3
         auxis(counter)%y = 0
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+1)
         ! *** fxxy orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 2
         auxis(counter)%y = 1
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+4)
         ! *** fxxz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 2
         auxis(counter)%y = 0
         auxis(counter)%z = 1
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+5)
         ! *** fxyy orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 1
         auxis(counter)%y = 2
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+6)
         ! *** fxyz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 1
         auxis(counter)%y = 1
         auxis(counter)%z = 1
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+10)
         ! *** fxzz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 1
         auxis(counter)%y = 0
         auxis(counter)%z = 2
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+8)
         ! *** fyyy orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 3
         auxis(counter)%z = 0
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+2)
         ! *** fyyz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 2
         auxis(counter)%z = 1
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+7)
         ! *** fyzz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 1
         auxis(counter)%z = 2
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+9)
         ! *** fzzz orbital*** 
         counter = counter + 1
         auxis(counter)%expo = expon1
         auxis(counter)%x = 0
         auxis(counter)%y = 0
         auxis(counter)%z = 3
         auxis(counter)%coords(1) = site_info%local_crds(3*(site1-1)+1)
         auxis(counter)%coords(2) = site_info%local_crds(3*(site1-1)+2)
         auxis(counter)%coords(3) = site_info%local_crds(3*(site1-1)+3)
         auxis(counter)%coords_G(1) = site_info%site_crds(3*(site1-1)+1)
         auxis(counter)%coords_G(2) = site_info%site_crds(3*(site1-1)+2)
         auxis(counter)%coords_G(3) = site_info%site_crds(3*(site1-1)+3)
         auxis(counter)%norm = site_info%cart_coeff_norms(c_off1+3)
       endif
      endif ! select order
    enddo ! k1
  enddo ! site1
  if(counter .ne. ncoeff) then
    write(6,*)'AH_form_auxis: something wrong while forming auxis, exiting'
    write(6,*)counter,ncoeff
    stop
  endif
  return
end subroutine AH_form_auxis
