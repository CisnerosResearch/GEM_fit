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
subroutine AHSITESITE_interact(site_info,sys_ene)
use definition
  implicit none
  type(sites_info)::site_info

  integer site1,site2
  integer k1,num1,off1,c_off1,k2,num2,off2,c_off2
  integer order1,order2
  integer num_res,slo1,shi1,slo2,shi2,ires,jres
  double precision prim_prim_EE,tot_NN,tot_EE,tot_NE
  double precision site_site_EE,site_site_NN,site_site_NE,ss,sd,ds,dd
  double precision sp, ps, pp, pd, dp, sd2, d2s, d2d2, d2d, dd2, sp1, ps1, pp1
  double precision dx,dy,dz,tmp_vec(100)
  double precision expon1,expon2,boyspar,nuc_chg1,nuc_chg2,hartree,ene,sys_ene

  hartree = 627.51d0

  num_res = site_info%num_residues

  sys_ene = 0.d0
  do ires = 1,num_res-1
    slo1 = site_info%residue_start(ires)
    shi1 = site_info%residue_end(ires)
    do jres = ires+1,num_res
      slo2 = site_info%residue_start(jres)
      shi2 = site_info%residue_end(jres)
      tot_NN = 0.d0
      tot_EE = 0.d0
      tot_NE = 0.d0
      do site1 = slo1,shi1
        num1 = site_info%num_primitives(site1)
        off1 = site_info%off_primitives(site1)
        do site2 = slo2,shi2
          num2 = site_info%num_primitives(site2)
          off2 = site_info%off_primitives(site2)
          dx = site_info%site_crds(3*(site2-1)+1) - &
                         site_info%site_crds(3*(site1-1)+1)
          dy = site_info%site_crds(3*(site2-1)+2) - &
                         site_info%site_crds(3*(site1-1)+2)
          dz = site_info%site_crds(3*(site2-1)+3) - &
                         site_info%site_crds(3*(site1-1)+3)
          nuc_chg1 = site_info%nuclear_charge(site1)
          nuc_chg2 = site_info%nuclear_charge(site2)
          site_site_NN = nuc_chg1*nuc_chg2 / dsqrt(dx*dx+dy*dy+dz*dz)
          site_site_EE = 0.d0
          ss = 0.d0
          sd = 0.d0
          ds = 0.d0
          dd = 0.d0
          ! GAC
          sp = 0.d0
          sp1 = 0.d0
          ps = 0.d0
          ps1 = 0.d0
          pp = 0.d0
          pp1 = 0.d0
          pd = 0.d0
          dp = 0.d0
          sd2 = 0.d0
          d2s = 0.d0
          d2d2 = 0.d0
          d2d = 0.d0
          dd2 = 0.d0
          do k1 = 1,num1
            order1 = site_info%prim_order(off1+k1)
            expon1 = site_info%prim_expo(off1+k1)
            c_off1 = site_info%coeff_offset(off1+k1)
            do k2 = 1,num2
              order2 = site_info%prim_order(off2+k2)
              expon2 = site_info%prim_expo(off2+k2)
              c_off2 = site_info%coeff_offset(off2+k2)
              if ( order1 == 1 .and. order2 == 1 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_ss(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    ss = ss + prim_prim_EE
              elseif ( order1 == 1 .and. order2 == 3 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_sp1(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    sp1 = sp1 + prim_prim_EE
              elseif ( order1 == 1 .and. order2 == 4 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_sp(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    sp = sp + prim_prim_EE
              elseif ( order1 == 1 .and. order2 == 6 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_sda(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    sd2 = sd2 + prim_prim_EE
              elseif ( order1 == 1 .and. order2 == 10 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_sd(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    sd = sd + prim_prim_EE
              elseif ( order1 == 3 .and. order2 == 1 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_ps1(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    ps1 = ps1 + prim_prim_EE
              elseif ( order1 == 4 .and. order2 == 1 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_ps(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    ps = ps + prim_prim_EE
              elseif ( order1 == 3 .and. order2 == 3 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_pp1(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    pp1 = pp1 + prim_prim_EE
              elseif ( order1 == 4 .and. order2 == 4 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_pp(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    pp = pp + prim_prim_EE
              elseif ( order1 == 6 .and. order2 == 1 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_das(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    d2s = d2s + prim_prim_EE
              elseif ( order1 == 10 .and. order2 == 1 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_ds(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    ds = ds + prim_prim_EE
              elseif ( order1 == 6 .and. order2 == 6 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_dada(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    d2d2 = d2d2 + prim_prim_EE
              elseif ( order1 == 6 .and. order2 == 10 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_dad(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    d2d = d2d + prim_prim_EE
              elseif ( order1 == 10 .and. order2 == 6 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_dda(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    d2d = dd2 + prim_prim_EE
              elseif ( order1 == 10 .and. order2 == 10 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_dd(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    dd = dd + prim_prim_EE
              ! GAC : consider p shells only, with d and spd shells
              elseif ( (order1 == 3 .and. order2 .ge. 6) .or. &
                       (order1 .ge. 6 .and. order2 == 3) ) then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_RECp(order1,order2,boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
              ! GAC : consider sp shells only, with spd shells
              elseif ( (order1 == 4 .and. order2 == 10) .or. &
                       (order1 == 10 .and. order2 == 4) ) then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_RECp(order1,order2,boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    if(order1 == 4)pd = pd + prim_prim_EE
                    if(order2 == 4)dp = dp + prim_prim_EE
              else ! catch all for F functions
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_REC(order1,order2,boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
              endif
              site_site_EE = site_site_EE + prim_prim_EE
              !print *,'ss_EE = ',site_site_EE
            enddo !k2 = 1,num2
          enddo!k1 = 1,num1
          site_site_NE = 0.d0
          do k1 = 1,num1
            order1 = site_info%prim_order(off1+k1)
            expon1 = site_info%prim_expo(off1+k1)
            c_off1 = site_info%coeff_offset(off1+k1)
            boyspar = expon1
            if ( order1 == 1 )then
              call GN_MCMUR_DAV_ss(boyspar,dx,dy,dz,  &
                  site_info%global_hermite_coeffs(c_off1+1),nuc_chg2,ene,&
                    tmp_vec)
            elseif ( order1 == 3 )then ! GAC
              call GN_MCMUR_DAV_ps1(boyspar,dx,dy,dz,  &
                  site_info%global_hermite_coeffs(c_off1+1),nuc_chg2,ene,&
                    tmp_vec)
            elseif ( order1 == 4 )then ! GAC
              call GN_MCMUR_DAV_ps(boyspar,dx,dy,dz,  &
                  site_info%global_hermite_coeffs(c_off1+1),nuc_chg2,ene,&
                    tmp_vec)
            elseif ( order1 == 6 )then ! GAC
              call GN_MCMUR_DAV_das(boyspar,dx,dy,dz,  &
                  site_info%global_hermite_coeffs(c_off1+1),nuc_chg2,ene,&
                    tmp_vec)
            elseif ( order1 == 10 )then
              call GN_MCMUR_DAV_ds(boyspar,dx,dy,dz,  &
                  site_info%global_hermite_coeffs(c_off1+1),nuc_chg2,ene,&
                    tmp_vec)
            elseif ( order1 == 20 ) then
              call GN_MCMUR_DAV_REC(order1,1,boyspar,dx,dy,dz,  &
                  site_info%global_hermite_coeffs(c_off1+1),nuc_chg2,ene,&
                    tmp_vec)
            endif
            site_site_NE = site_site_NE + ene
          enddo
          do k2 = 1,num2
            order2 = site_info%prim_order(off2+k2)
            expon2 = site_info%prim_expo(off2+k2)
            c_off2 = site_info%coeff_offset(off2+k2)
            boyspar = expon2
            if ( order2 == 1 )then
              call GN_MCMUR_DAV_ss(boyspar,dx,dy,dz,  &
                  nuc_chg1,site_info%global_hermite_coeffs(c_off2+1),ene,&
                    tmp_vec)
            elseif ( order2 == 3 )then ! GAC
              call GN_MCMUR_DAV_sp(boyspar,dx,dy,dz,  &
                  nuc_chg1,site_info%global_hermite_coeffs(c_off2+1),ene,&
                    tmp_vec)
            elseif ( order2 == 4 )then ! GAC
              call GN_MCMUR_DAV_sp(boyspar,dx,dy,dz,  &
                  nuc_chg1,site_info%global_hermite_coeffs(c_off2+1),ene,&
                    tmp_vec)
            elseif ( order2 == 6 )then ! GAC
              call GN_MCMUR_DAV_sda(boyspar,dx,dy,dz,  &
                  nuc_chg1,site_info%global_hermite_coeffs(c_off2+1),ene,&
                    tmp_vec)
            elseif ( order2 == 10 )then
              call GN_MCMUR_DAV_sd(boyspar,dx,dy,dz,  &
                  nuc_chg1,site_info%global_hermite_coeffs(c_off2+1),ene,&
                    tmp_vec)
            elseif ( order2 == 20 )then
              call GN_MCMUR_DAV_REC(1,order2,boyspar,dx,dy,dz,  &
                  nuc_chg1,site_info%global_hermite_coeffs(c_off2+1),ene,&
                    tmp_vec)
            endif
            site_site_NE = site_site_NE + ene
          enddo
          !write(37,*)'site1,site2,site_site_NN,site_site_EE,site_site_NE = ',&
                 !site1,site2,site_site_NN,site_site_EE,site_site_NE
          !write(37,*)'=====>ss,sd,ds,dd,sp1,ps1,sp,ps,dp,pd = ',&
          !                   ss,sd,ds,dd,sp,ps,dp,pd
          tot_EE = tot_EE + site_site_EE
          tot_NN = tot_NN + site_site_NN
          tot_NE = tot_NE + site_site_NE
          !print *,'e_pp = ',tot_NN + tot_EE - tot_NE
        enddo!site2 = slo2,shi2
      enddo!site1 = slo1,slo2
      if(site_info%num_residues .lt. 3) then
         write(6,*)'ires,jres = ',ires,jres
         write(6,*)'tot_NN,tot_EE,tot_NE = ',tot_NN,tot_EE,tot_NE
        !write(37,*)'tot_NN,tot_EE,tot_NE = ',tot_NN,tot_EE,tot_NE
         write(6,*)'tot respair ene = ',tot_NN + tot_EE - tot_NE, &
          '   in kcals = ',hartree*(tot_NN + tot_EE - tot_NE)
      endif
      sys_ene = sys_ene + hartree*(tot_NN + tot_EE - tot_NE)
    enddo !jres = ires+1,num_res
  enddo! ires = 1,num_res-1
  write(6,*)'total system energy = ',sys_ene
  sys_ene = sys_ene/hartree
  return
end subroutine AHSITESITE_interact
!----------------------------------------------------------
subroutine AHSITESITE_interact2(site_info,site_info2,sys_ene)
use definition
  implicit none
  type(sites_info)::site_info,site_info2

  integer site1,site2
  integer k1,num1,off1,c_off1,k2,num2,off2,c_off2
  integer order1,order2,i,count1,count2
  integer slo1,shi1,slo2,shi2,ires,jres
  double precision prim_prim_EE,tot_NN,tot_EE,tot_NE
  double precision site_site_EE,site_site_NN,site_site_NE,ss,sd,ds,dd
  double precision sp, ps, pp, pd, dp, sd2, d2s, d2d2, d2d, dd2 ! GAC 
  double precision dx,dy,dz,tmp_vec(100)
  double precision expon1,expon2,boyspar,nuc_chg1,nuc_chg2,hartree,ene,sys_ene

  hartree = 627.51d0

  ! FOR JP FIT
  !jres=0
  !do ires = 1,site_info%num_residues
  !  slo1 = site_info%residue_start(ires)
  !  shi1 = site_info%residue_end(ires)
  !    do site1 = slo1,shi1
  !      num1 = site_info%num_primitives(site1)
  !      off1 = site_info%off_primitives(site1)
  !        do k1 = 1,num1
  !          order1 = site_info%prim_order(off1+k1)
  !          expon1 = site_info%prim_expo(off1+k1)
  !          c_off1 = site_info%coeff_offset(off1+k1)
  !          if ( order1 == 1 )then
  !            jres= jres+1
  !            write(66,*),jres,site1,expon1,order1
  !          else
  !            do k2 = 1, 10 !1=only s; 4 = sp
  !               jres= jres+1
  !              write(66,*),jres,site1,expon1,order1
  !            enddo
  !          endif
  !          !if ( order1 == 10 )then
  !          !  do k2 = 1, 9 !1=only s; 4 = sp
  !          !   site_info%global_hermite_coeffs(c_off1+k2+1) = 0.d0
  !          !   site_info%local_hermite_coeffs(c_off1+k2+1) = 0.d0
  !          !  enddo
  !          !endif
  !        enddo
  !    enddo
  !enddo

  count1=0
  count2=0
  sys_ene = 0.d0
  do ires = 1,site_info%num_residues
    slo1 = site_info%residue_start(ires)
    shi1 = site_info%residue_end(ires)
    do jres = 1,site_info2%num_residues
      slo2 = site_info2%residue_start(jres)
      shi2 = site_info2%residue_end(jres)
      tot_NN = 0.d0
      tot_EE = 0.d0
      tot_NE = 0.d0
      do site1 = slo1,shi1
        num1 = site_info%num_primitives(site1)
        off1 = site_info%off_primitives(site1)
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
          nuc_chg2 = site_info2%nuclear_charge(site2)
          if (dsqrt(dx*dx+dy*dy+dz*dz)>1.d-10) then
             site_site_NN = nuc_chg1*nuc_chg2 / dsqrt(dx*dx+dy*dy+dz*dz)
          else
             site_site_NN = 0.d0
          endif
          site_site_EE = 0.d0
          ss = 0.d0
          sd = 0.d0
          ds = 0.d0
          dd = 0.d0
          ! GAC
          sp = 0.d0
          ps = 0.d0
          pp = 0.d0
          pd = 0.d0
          dp = 0.d0
          sd2 = 0.d0
          d2s = 0.d0
          d2d2 = 0.d0
          d2d = 0.d0
          dd2 = 0.d0
          do k1 = 1,num1
            order1 = site_info%prim_order(off1+k1)
            expon1 = site_info%prim_expo(off1+k1)
            c_off1 = site_info%coeff_offset(off1+k1)
            do k2 = 1,num2
              order2 = site_info2%prim_order(off2+k2)
              expon2 = site_info2%prim_expo(off2+k2)
              c_off2 = site_info2%coeff_offset(off2+k2)
              if ( order1 == 1 .and. order2 == 1 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_ss(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    ss = ss + prim_prim_EE
              !elseif ( order1 == 1 .and. order2 == 3 )then
              !  boyspar = expon1*expon2 / (expon1 + expon2)
              !  call GN_MCMUR_DAV_sp1(boyspar,dx,dy,dz,  &
              !      site_info%global_hermite_coeffs(c_off1+1),  &
              !      site_info2%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
              !      tmp_vec)
              elseif ( order1 == 1 .and. order2 == 4 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_sp(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    sp = sp + prim_prim_EE
              elseif ( order1 == 1 .and. order2 == 6 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_sda(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    sd2 = sd2 + prim_prim_EE
              elseif ( order1 == 1 .and. order2 == 10 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_sd(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    sd = sd + prim_prim_EE
              !elseif ( order1 == 3 .and. order2 == 1 )then
              !  boyspar = expon1*expon2 / (expon1 + expon2)
              !  call GN_MCMUR_DAV_ps1(boyspar,dx,dy,dz,  &
              !      site_info%global_hermite_coeffs(c_off1+1),  &
              !      site_info2%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
              !      tmp_vec)
              elseif ( order1 == 4 .and. order2 == 1 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_ps(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    ps = ps + prim_prim_EE
              elseif ( order1 == 3 .and. order2 == 3 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_pp1(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
              elseif ( order1 == 4 .and. order2 == 4 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_pp(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    pp = pp + prim_prim_EE
              elseif ( order1 == 6 .and. order2 == 1 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_das(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    d2s = d2s + prim_prim_EE
              elseif ( order1 == 10 .and. order2 == 1 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_ds(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    ds = ds + prim_prim_EE
              elseif ( order1 == 6 .and. order2 == 6 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_dada(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    d2d2 = d2d2 + prim_prim_EE
              elseif ( order1 == 6 .and. order2 == 10 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_dad(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    d2d = d2d + prim_prim_EE
              elseif ( order1 == 10 .and. order2 == 6 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_dda(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    d2d = dd2 + prim_prim_EE
              elseif ( order1 == 10 .and. order2 == 10 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_dd(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    dd = dd + prim_prim_EE
              ! GAC : consider p shells only on A with anything 
              elseif ( order1 == 3 .and. order2 .ne. 3) then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_RECp2(order1,order2,boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
              ! GAC : consider p shells only on B with anything 
              elseif (order2 == 3 .and. order1 .ne. 3) then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_RECp2a(order1,order2,boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
              ! GAC : consider sp shells only, with spd shells
              elseif ( (order1 == 4 .and. order2 == 10) .or. &
                       (order1 == 10 .and. order2 == 4) ) then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_RECp(order1,order2,boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    if(order1 == 4)pd = pd + prim_prim_EE
                    if(order2 == 4)dp = dp + prim_prim_EE
              else ! catch all for F functions
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_REC(order1,order2,boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
              endif
              site_site_EE = site_site_EE + prim_prim_EE
            enddo !k2 = 1,num2
          enddo!k1 = 1,num1
          site_site_NE = 0.d0
          do k1 = 1,num1
            order1 = site_info%prim_order(off1+k1)
            expon1 = site_info%prim_expo(off1+k1)
            c_off1 = site_info%coeff_offset(off1+k1)
            boyspar = expon1
            if ( order1 == 1 )then
              call GN_MCMUR_DAV_ss(boyspar,dx,dy,dz,  &
                  site_info%global_hermite_coeffs(c_off1+1),nuc_chg2,ene,&
                    tmp_vec)
            elseif ( order1 == 3 )then ! GAC
              call GN_MCMUR_DAV_ps1(boyspar,dx,dy,dz,  &
                  site_info%global_hermite_coeffs(c_off1+1),nuc_chg2,ene,&
                    tmp_vec)
            elseif ( order1 == 4 )then ! GAC
              call GN_MCMUR_DAV_ps(boyspar,dx,dy,dz,  &
                  site_info%global_hermite_coeffs(c_off1+1),nuc_chg2,ene,&
                    tmp_vec)
            elseif ( order1 == 6 )then ! GAC
              call GN_MCMUR_DAV_das(boyspar,dx,dy,dz,  &
                  site_info%global_hermite_coeffs(c_off1+1),nuc_chg2,ene,&
                    tmp_vec)
            elseif ( order1 == 10 )then
              call GN_MCMUR_DAV_ds(boyspar,dx,dy,dz,  &
                  site_info%global_hermite_coeffs(c_off1+1),nuc_chg2,ene,&
                    tmp_vec)
            elseif ( order1 == 20 ) then
              call GN_MCMUR_DAV_REC(order1,1,boyspar,dx,dy,dz,  &
                  site_info%global_hermite_coeffs(c_off1+1),nuc_chg2,ene,&
                    tmp_vec)
            endif
            site_site_NE = site_site_NE + ene
          enddo
          do k2 = 1,num2
            order2 = site_info2%prim_order(off2+k2)
            expon2 = site_info2%prim_expo(off2+k2)
            c_off2 = site_info2%coeff_offset(off2+k2)
            boyspar = expon2
            if ( order2 == 1 )then
              call GN_MCMUR_DAV_ss(boyspar,dx,dy,dz,  &
                  nuc_chg1,site_info2%global_hermite_coeffs(c_off2+1),ene,&
                    tmp_vec)
            elseif ( order2 == 3 )then ! GAC
              call GN_MCMUR_DAV_sp1(boyspar,dx,dy,dz,  &
                  nuc_chg1,site_info2%global_hermite_coeffs(c_off2+1),ene,&
                    tmp_vec)
                  !do k1 = 1,3
                  !   print *,site_info2%global_hermite_coeffs(c_off2+k1)
                  !enddo
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
            elseif ( order2 == 10 )then
              call GN_MCMUR_DAV_REC(1,order2,boyspar,dx,dy,dz,  &
                  nuc_chg1,site_info2%global_hermite_coeffs(c_off2+1),ene,&
                    tmp_vec)
            endif
            site_site_NE = site_site_NE + ene
          enddo
          !write(37,*)'site1,site2,site_site_NN,site_site_EE,site_site_NE = ',&
                 !site1,site2,site_site_NN,site_site_EE,site_site_NE
          !write(37,*)'=====>ss,sd,ds,dd,sp,ps,dp,pd = ',&
          !                   ss,sd,ds,dd,sp,ps,dp,pd
          tot_EE = tot_EE + site_site_EE
          tot_NN = tot_NN + site_site_NN
          tot_NE = tot_NE + site_site_NE
        enddo!site2 = slo2,shi2
      enddo!site1 = slo1,slo2
      if((site_info%num_residues .eq. 1) .and. &
         (site_info2%num_residues .eq. 1)) then
         !write(6,*)'ires,jres = ',ires,jres
         write(6,*)'tot_NN,tot_EE,tot_NE = ',tot_NN,tot_EE,tot_NE
        !write(37,*)'tot_NN,tot_EE,tot_NE = ',tot_NN,tot_EE,tot_NE
         write(6,*)'tot respair ene = ',tot_NN + tot_EE - tot_NE, &
          '   in kcals = ',hartree*(tot_NN + tot_EE - tot_NE)
      endif
      sys_ene = sys_ene + hartree*(tot_NN + tot_EE - tot_NE)
    enddo !jres = 1,site_info2%num_res
  enddo! ires = 1,site_info%num_res
  write(6,*)'total system energy = ',sys_ene
  sys_ene = sys_ene/hartree
  return
end subroutine AHSITESITE_interact2
!----------------------------------------------------------
subroutine AHSITESITE_interact_CHG(site_info,sys_ene)
use definition
  implicit none
  type(sites_info)::site_info

  integer site1,site2
  integer k1,num1,off1,c_off1,k2,num2,off2,c_off2, l
  integer order1,order2
  integer num_res,slo1,shi1,slo2,shi2,ires,jres
  double precision prim_prim_EE,tot_NN,tot_EE,tot_NE,coulomb
  double precision site_site_EE,site_site_NN,site_site_NE,ss,sd,ds,dd
  double precision sp, ps, pp, pd, dp, sd2, d2s, d2d2, d2d, dd2, sp1, ps1, pp1
  double precision dx,dy,dz,bohr,tmp_vec(100)
  double precision expon1,expon2,boyspar,nuc_chg1,nuc_chg2,hartree,ene,sys_ene

  hartree = 627.51d0
  coulomb = 332.05382d0 ! THIS IS THE COULOMB CONSTANT FROM TINKER
  bohr = 0.529177249d0

  print *,' SITE SITE WITH NUCLEAR CHARGE ADDED TO MONOPOLE'
  print *,' MAKE SURE THAT THE MULTIPOLES ARE IN ANGTROM, NOT a.u.!'
  print *,' AND THAT THEY DO NOT INCLUDE 2/3 FACTOR FROM STONE!'

  num_res = site_info%num_residues

  sys_ene = 0.d0
  do ires = 1,num_res-1
    slo1 = site_info%residue_start(ires)
    shi1 = site_info%residue_end(ires)
    do jres = ires+1,num_res
      slo2 = site_info%residue_start(jres)
      shi2 = site_info%residue_end(jres)
      tot_NN = 0.d0
      tot_EE = 0.d0
      tot_NE = 0.d0
      do site1 = slo1,shi1
        num1 = site_info%num_primitives(site1)
        off1 = site_info%off_primitives(site1)
        do site2 = slo2,shi2
          num2 = site_info%num_primitives(site2)
          off2 = site_info%off_primitives(site2)
          dx = (site_info%site_crds(3*(site2-1)+1) - &
                         site_info%site_crds(3*(site1-1)+1))*bohr
          dy = (site_info%site_crds(3*(site2-1)+2) - &
                         site_info%site_crds(3*(site1-1)+2))*bohr
          dz = (site_info%site_crds(3*(site2-1)+3) - &
                         site_info%site_crds(3*(site1-1)+3))*bohr
          !print *,'dx,dy,dz',dx,dy,dz
          nuc_chg1 = site_info%nuclear_charge(site1)
          nuc_chg2 = site_info%nuclear_charge(site2)
          site_site_EE = 0.d0
          ss = 0.d0
          sd = 0.d0
          ds = 0.d0
          dd = 0.d0
          ! GAC
          sp = 0.d0
          sp1 = 0.d0
          ps = 0.d0
          ps1 = 0.d0
          pp = 0.d0
          pp1 = 0.d0
          pd = 0.d0
          dp = 0.d0
          sd2 = 0.d0
          d2s = 0.d0
          d2d2 = 0.d0
          d2d = 0.d0
          dd2 = 0.d0
          do k1 = 1,num1
            order1 = site_info%prim_order(off1+k1)
            expon1 = site_info%prim_expo(off1+k1)
            c_off1 = site_info%coeff_offset(off1+k1)
            do k2 = 1,num2
              order2 = site_info%prim_order(off2+k2)
              expon2 = site_info%prim_expo(off2+k2)
              c_off2 = site_info%coeff_offset(off2+k2)
              if ( order1 == 1 .and. order2 == 1 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_ss(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    ss = ss + prim_prim_EE
              elseif ( order1 == 1 .and. order2 == 3 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_sp1(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    sp1 = sp1 + prim_prim_EE
              elseif ( order1 == 1 .and. order2 == 4 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_sp(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    sp = sp + prim_prim_EE
              elseif ( order1 == 1 .and. order2 == 6 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_sda(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    sd2 = sd2 + prim_prim_EE
              elseif ( order1 == 1 .and. order2 == 10 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_sd(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    sd = sd + prim_prim_EE
              elseif ( order1 == 3 .and. order2 == 1 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_ps1(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    ps1 = ps1 + prim_prim_EE
              elseif ( order1 == 4 .and. order2 == 1 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_ps(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    ps = ps + prim_prim_EE
              elseif ( order1 == 3 .and. order2 == 3 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_pp1(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    pp1 = pp1 + prim_prim_EE
              elseif ( order1 == 4 .and. order2 == 4 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_pp(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    pp = pp + prim_prim_EE
              elseif ( order1 == 6 .and. order2 == 1 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_das(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    d2s = d2s + prim_prim_EE
              elseif ( order1 == 10 .and. order2 == 1 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_ds(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    ds = ds + prim_prim_EE
              elseif ( order1 == 6 .and. order2 == 6 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_dada(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    d2d2 = d2d2 + prim_prim_EE
              elseif ( order1 == 6 .and. order2 == 10 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_dad(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    d2d = d2d + prim_prim_EE
              elseif ( order1 == 10 .and. order2 == 6 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_dda(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    d2d = dd2 + prim_prim_EE
              elseif ( order1 == 10 .and. order2 == 10 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_dd(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                !print *,'e_pp = ',prim_prim_EE 
                !do l = 1,10
                !   print *, site_info%global_hermite_coeffs(c_off1+1+l-1)
                !enddo
                !do l = 1,10
                !   print *, site_info%global_hermite_coeffs(c_off2+1+l-1)
                !enddo
                    dd = dd + prim_prim_EE
              ! GAC : consider p shells only, with d and spd shells
              elseif ( (order1 == 3 .and. order2 .ge. 6) .or. &
                       (order1 .ge. 6 .and. order2 == 3) ) then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_RECp(order1,order2,boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
              ! GAC : consider sp shells only, with spd shells
              elseif ( (order1 == 4 .and. order2 == 10) .or. &
                       (order1 == 10 .and. order2 == 4) ) then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_RECp(order1,order2,boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
                    if(order1 == 4)pd = pd + prim_prim_EE
                    if(order2 == 4)dp = dp + prim_prim_EE
              else ! catch all for F functions
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_REC(order1,order2,boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),prim_prim_EE,&
                    tmp_vec)
              endif
              site_site_EE = site_site_EE + prim_prim_EE
            enddo !k2 = 1,num2
          enddo!k1 = 1,num1
          tot_EE = tot_EE + site_site_EE
        enddo!site2 = slo2,shi2
      enddo!site1 = slo1,slo2
      if(site_info%num_residues .lt. 3) then
         write(6,*)'ires,jres = ',ires,jres
         write(6,*)'tot_NN,tot_EE,tot_NE = ',tot_NN,tot_EE,tot_NE
        !write(37,*)'tot_NN,tot_EE,tot_NE = ',tot_NN,tot_EE,tot_NE
         write(6,*)'tot respair ene = ',tot_NN + tot_EE - tot_NE, &
          '   in kcals = ',coulomb*(tot_NN + tot_EE - tot_NE)
      endif
      sys_ene = sys_ene + coulomb*(tot_NN + tot_EE - tot_NE)
    enddo !jres = ires+1,num_res
  enddo! ires = 1,num_res-1
  write(6,*)'total system energy = ',sys_ene
  sys_ene = sys_ene/coulomb
  return
end subroutine AHSITESITE_interact_CHG
!----------------------------------------------------------
