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
subroutine GSITESITE_setup(p_site_info,p_site_site)
  implicit none
  integer p_site_info,p_site_site
# include "database.fh"
  integer ptr,p_mcmd_ss,p_mcmd_sd,p_mcmd_ds,p_mcmd_dd,deg1,deg2,itoptab
  integer maxnum_s,maxnum_d,topdeg
  p_mcmd_ss = DB_alloc(p_site_site,DB_node,"ss_mcmurchie_davidson",30)
  p_mcmd_ds = DB_alloc(p_site_site,DB_node,"ds_mcmurchie_davidson",30)
  p_mcmd_sd = DB_alloc(p_site_site,DB_node,"sd_mcmurchie_davidson",30)
  p_mcmd_dd = DB_alloc(p_site_site,DB_node,"dd_mcmurchie_davidson",30)
  write(6,*)'going into GSITESITE_get_max_prims'
  call GSITESITE_get_max_prims(p_site_info,maxnum_s,maxnum_d)
  write(6,*)'back from GSITESITE_get_max_prims'
  ptr = DB_alloc(p_site_site,DB_int,"maximum_num_s_prims",1)
  IH(ptr) = maxnum_s
  ptr = DB_alloc(p_site_site,DB_int,"maximum_num_d_prims",1)
  IH(ptr) = maxnum_d
  ptr = DB_alloc(p_site_site,DB_real,"exponents_s",maxnum_s)
  ptr = DB_alloc(p_site_site,DB_real,"exponents_d",maxnum_d)
  ptr = DB_alloc(p_site_site,DB_real,"coefficients_s",maxnum_s)
  ptr = DB_alloc(p_site_site,DB_real,"coefficients_d",maxnum_d*10)
  ptr = DB_alloc(p_site_site,DB_real,"scratch",maxnum_s+maxnum_d)
  ! ss
  deg1 = 0
  deg2 = 0
  topdeg = deg1 + deg2 + 1 !one more for grad
  call GMcMur_Dav_init(p_mcmd_ss,topdeg,itoptab,deg1,deg2)
  ptr = DB_alloc(p_site_site,DB_int,"size_ss_mcmd_table",1)
  IH(ptr) = itoptab
  ptr = DB_alloc(p_site_site,DB_real,"boys_function_ss",maxnum_s*(topdeg+1))
  ptr = DB_alloc(p_site_site,DB_real,"mcmd_recur_ss",maxnum_s*itoptab)
  ptr = DB_alloc(p_site_site,DB_real,"boys_par_ss",maxnum_s)
  ! sd
  deg1 = 0
  deg2 = 2
  topdeg = deg1 + deg2 + 1 !one more for grad
  call GMcMur_Dav_init(p_mcmd_sd,topdeg,itoptab,deg1,deg2)
  ptr = DB_alloc(p_site_site,DB_int,"size_sd_mcmd_table",1)
  IH(ptr) = itoptab
  ptr = DB_alloc(p_site_site,DB_real,"boys_function_sd",maxnum_d*(topdeg+1))
  ptr = DB_alloc(p_site_site,DB_real,"mcmd_recur_sd",maxnum_d*itoptab)
  ptr = DB_alloc(p_site_site,DB_real,"boys_par_sd",maxnum_d)
  !  ds
  deg1 = 2
  deg2 = 0
  topdeg = deg1 + deg2 + 1 !one more for grad
  call GMcMur_Dav_init(p_mcmd_ds,topdeg,itoptab,deg1,deg2)
  ptr = DB_alloc(p_site_site,DB_int,"size_ds_mcmd_table",1)
  IH(ptr) = itoptab
  ptr = DB_alloc(p_site_site,DB_real,"boys_function_ds",maxnum_s*(topdeg+1))
  ptr = DB_alloc(p_site_site,DB_real,"mcmd_recur_ds",maxnum_s*itoptab)
  ptr = DB_alloc(p_site_site,DB_real,"boys_par_ds",maxnum_s)
  ! dd
  deg1 = 2
  deg2 = 2
  topdeg = deg1 + deg2 + 1 !one more for grad
  call GMcMur_Dav_init(p_mcmd_dd,topdeg,itoptab,deg1,deg2)
  ptr = DB_alloc(p_site_site,DB_int,"size_dd_mcmd_table",1)
  IH(ptr) = itoptab
  ptr = DB_alloc(p_site_site,DB_real,"boys_function_dd",maxnum_d*(topdeg+1))
  ptr = DB_alloc(p_site_site,DB_real,"mcmd_recur_dd",maxnum_d*itoptab)
  ptr = DB_alloc(p_site_site,DB_real,"boys_par_dd",maxnum_d)


end subroutine GSITESITE_setup
!---------------------------------------------------------
subroutine GSITESITE_get_max_prims(p_site_info,maxnum_s,maxnum_d)
  implicit none
  integer p_site_info,maxnum_s,maxnum_d
# include "database.fh"
  integer ptr,nsites,k,n,off,num,order,nprim,ps_numprim,ps_offprim,ps_primorder
  integer num_s,num_d,p_nums,p_numd

  ptr = DB_Get_Pointer(p_site_info,DB_int,"numsites",1)
  nsites = IH(ptr)
  ps_numprim = DB_Get_Pointer(p_site_info,DB_int,"num_site_prims",nsites)
  ps_offprim = DB_Get_Pointer(p_site_info,DB_int,"off_site_prims",nsites)
  ptr = DB_Get_Pointer(p_site_info,DB_int,"tot_site_prims",1)
  nprim = IH(ptr)
  ps_primorder = DB_Get_Pointer(p_site_info,DB_int,"prim_order",nprim)
  p_nums = DB_alloc(p_site_info,DB_int,"num_s_prims",nsites)
  p_numd = DB_alloc(p_site_info,DB_int,"num_d_prims",nsites)

  maxnum_s = -1
  maxnum_d = -1
  do n = 1,nsites
    num_s = 0
    num_d = 0
    num = IH(ps_numprim+n-1)
    off = IH(ps_offprim+n-1)
    do k = 1,num
      order = IH(ps_primorder+off+k-1)
      if ( order == 1 )then
        num_s = num_s + 1
      elseif ( order == 10 )then
        num_d = num_d + 1
      endif
    enddo
    IH(p_nums+n-1) = num_s
    IH(p_numd+n-1) = num_d
    if ( num_s > maxnum_s )maxnum_s = num_s
    if ( num_d > maxnum_d )maxnum_d = num_d
  enddo
  write(6,*)'maxnum_s,maxnum_d = ',maxnum_s,maxnum_d

  return
end subroutine GSITESITE_get_max_prims
!---------------------------------------------------------
subroutine GSITESITE_interact(p_site_info,p_site_site)
  implicit none
  integer p_site_info,p_site_site
# include "database.fh"
  integer ptr,p_sitecrds,ps_numprim,ps_offprim,ps_primorder,  &
          ps_primexpo,ps_prim_off_coeff,p_G_hermite_coeff
  integer nsites,site1,site2,nprim,nums,numd,p_nums,p_numd
  integer k1,num1,off1,order1,c_off1,num2,off2
  double precision expon1,ene_ss,ene_sd,ene_ds,ene_dd,ene
  integer p_exps,p_expd,p_coeff_s,p_coeff_d,p_scratch
  integer itoptab,topdeg,ncoeff,maxnum_s,maxnum_d
  integer p_mcmd_ss,p_boyspar_ss,p_boystab_ss,p_recurtab_ss
  integer p_mcmd_sd,p_boyspar_sd,p_boystab_sd,p_recurtab_sd
  integer p_mcmd_ds,p_boyspar_ds,p_boystab_ds,p_recurtab_ds
  integer p_mcmd_dd,p_boyspar_dd,p_boystab_dd,p_recurtab_dd
! GAC
  integer p_numsites 
! GAC
  double precision dx,dy,dz
  double precision tot_ene_ss,tot_ene_sd,tot_ene_ds,tot_ene_dd,tot_ene

  write(6,*)'inside GSITESITE_interact'
  p_mcmd_ss = DB_Get_Node(p_site_site,"ss_mcmurchie_davidson")
  p_mcmd_ds = DB_Get_Node(p_site_site,"ds_mcmurchie_davidson")
  p_mcmd_sd = DB_Get_Node(p_site_site,"sd_mcmurchie_davidson")
  p_mcmd_dd = DB_Get_Node(p_site_site,"dd_mcmurchie_davidson")
  ptr = DB_Get_Pointer(p_site_info,DB_int,"numsites",1)
  nsites = IH(ptr)
  p_sitecrds = DB_Get_Pointer(p_site_info,DB_real,"site_crds",3*nsites)
  ps_numprim = DB_Get_Pointer(p_site_info,DB_int,"num_site_prims",nsites)
  ps_offprim = DB_Get_Pointer(p_site_info,DB_int,"off_site_prims",nsites)
  ptr = DB_Get_Pointer(p_site_info,DB_int,"tot_site_prims",1)
  nprim = IH(ptr)
  ps_primorder = DB_Get_Pointer(p_site_info,DB_int,"prim_order",nprim)
  ps_primexpo = DB_Get_Pointer(p_site_info,DB_real,"prim_expon_term",nprim)
  ps_prim_off_coeff = DB_Get_Pointer(p_site_info,DB_int,"coeff_offset",nprim)
  p_nums = DB_Get_Pointer(p_site_info,DB_int,"num_s_prims",nsites)
  p_numd = DB_Get_Pointer(p_site_info,DB_int,"num_d_prims",nsites)
  ptr = DB_Get_Pointer(p_site_info,DB_int,"num_hermite_coefficients",1)
  ncoeff = IH(ptr)
  p_G_hermite_coeff = DB_Get_Pointer(p_site_info,DB_real,  &
                        "global_hermite_coeffs",ncoeff)
  ptr = DB_Get_Pointer(p_site_site,DB_int,"maximum_num_s_prims",1)
  maxnum_s = IH(ptr)
  ptr = DB_Get_Pointer(p_site_site,DB_int,"maximum_num_d_prims",1)
  maxnum_d = IH(ptr)
  p_exps = DB_Get_Pointer(p_site_site,DB_real,"exponents_s",maxnum_s)
  p_expd = DB_Get_Pointer(p_site_site,DB_real,"exponents_d",maxnum_d)
  p_coeff_s = DB_Get_Pointer(p_site_site,DB_real,"coefficients_s",maxnum_s)
  p_coeff_d = DB_Get_Pointer(p_site_site,DB_real,"coefficients_d",maxnum_d*10)
  p_scratch = DB_Get_Pointer(p_site_site,DB_real,"scratch",maxnum_s+maxnum_d)

  ptr = DB_Get_Pointer(p_site_site,DB_int,"size_ss_mcmd_table",1)
  itoptab = IH(ptr)
  topdeg = 1
  p_boystab_ss = DB_Get_Pointer(p_site_site,DB_real, &
           "boys_function_ss",maxnum_s*(topdeg+1))
  p_recurtab_ss = DB_Get_Pointer(p_site_site,DB_real,  &
           "mcmd_recur_ss",maxnum_s*itoptab)
  p_boyspar_ss = DB_Get_Pointer(p_site_site,DB_real,"boys_par_ss",maxnum_s)
  ptr = DB_Get_Pointer(p_site_site,DB_int,"size_sd_mcmd_table",1)
  itoptab = IH(ptr)
  topdeg = 3
  p_boystab_sd = DB_Get_Pointer(p_site_site,DB_real, &
           "boys_function_sd",maxnum_d*(topdeg+1))
  p_recurtab_sd = DB_Get_Pointer(p_site_site,DB_real,  &
           "mcmd_recur_sd",maxnum_d*itoptab)
  p_boyspar_sd = DB_Get_Pointer(p_site_site,DB_real,"boys_par_sd",maxnum_d)
  ptr = DB_Get_Pointer(p_site_site,DB_int,"size_ds_mcmd_table",1)
  itoptab = IH(ptr)
  topdeg = 3
  p_boystab_ds = DB_Get_Pointer(p_site_site,DB_real, &
           "boys_function_ds",maxnum_s*(topdeg+1))
  p_recurtab_ds = DB_Get_Pointer(p_site_site,DB_real,  &
           "mcmd_recur_ds",maxnum_s*itoptab)
  p_boyspar_ds = DB_Get_Pointer(p_site_site,DB_real,"boys_par_ds",maxnum_s)
  ptr = DB_Get_Pointer(p_site_site,DB_int,"size_dd_mcmd_table",1)
  itoptab = IH(ptr)
  topdeg = 5
  p_boystab_dd = DB_Get_Pointer(p_site_site,DB_real, &
           "boys_function_dd",maxnum_d*(topdeg+1))
  p_recurtab_dd = DB_Get_Pointer(p_site_site,DB_real,  &
           "mcmd_recur_dd",maxnum_d*itoptab)
  p_boyspar_dd = DB_Get_Pointer(p_site_site,DB_real,"boys_par_dd",maxnum_d)
  write(6,*)'done allocing'
! GAC: retrieve number of sites in molecule to calculate G
  ptr = DB_Get_Pointer(p_site_info,DB_int,"numsites",1)
  p_numsites = IH(ptr)
  print *,'GAC: p_numsites = ',p_numsites
! GAC: retrieve number of sites in molecule to calculate G

  tot_ene_ss = 0.d0
  tot_ene_sd = 0.d0
  tot_ene_ds = 0.d0
  tot_ene_dd = 0.d0
  do site1 = 1,5
    num1 = IH(ps_numprim+site1-1)
    off1 = IH(ps_offprim+site1-1)
    do site2 = 6,10
      num2 = IH(ps_numprim+site2-1)
      off2 = IH(ps_offprim+site2-1)
      nums = IH(p_nums+site2-1)
      numd = IH(p_numd+site2-1)
      call GSITESITE_fillsite_coeffs_expons(nums,numd,num2,off2, &
                    RH(p_G_hermite_coeff),RH(ps_primexpo),IH(ps_primorder), &
                    IH(ps_prim_off_coeff),RH(p_coeff_s),RH(p_exps), &
                    RH(p_coeff_d),RH(p_expd))
      write(6,*)'nums,numd = ',nums,numd
      ene_ss = 0.d0
      ene_sd = 0.d0
      ene_ds = 0.d0
      ene_dd = 0.d0
      dx = RH(p_sitecrds+3*(site2-1)) - RH(p_sitecrds+3*(site1-1))
      dy = RH(p_sitecrds+3*(site2-1)+1) - RH(p_sitecrds+3*(site1-1)+1)
      dz = RH(p_sitecrds+3*(site2-1)+2) - RH(p_sitecrds+3*(site1-1)+2)
      do k1 = 1,num1
        order1 = IH(ps_primorder+off1+k1-1)
        expon1 = RH(ps_primexpo+off1+k1-1)
        c_off1 = IH(ps_prim_off_coeff+off1+k1-1)
        write(6,*)'order1 = ',order1
        if ( order1 == 1 )then
          call GSITESITE_prim_group_prims_ene(p_mcmd_ss,  &
             nums,RH(p_coeff_s),RH(p_exps), &
             expon1,RH(p_G_hermite_coeff+c_off1),dx,dy,dz, &
             RH(p_boyspar_ss),RH(p_boystab_ss),RH(p_recurtab_ss), &
             RH(p_scratch),ene)
             ene_ss = ene_ss + ene
          call GSITESITE_prim_group_prims_ene(p_mcmd_sd,  &
             numd,RH(p_coeff_d),RH(p_expd), &
             expon1,RH(p_G_hermite_coeff+c_off1),dx,dy,dz, &
             RH(p_boyspar_sd),RH(p_boystab_sd),RH(p_recurtab_sd), &
             RH(p_scratch),ene)
             ene_sd = ene_sd + ene
        elseif ( order1 == 10 )then
          call GSITESITE_prim_group_prims_ene(p_mcmd_ds,  &
             nums,RH(p_coeff_s),RH(p_exps), &
             expon1,RH(p_G_hermite_coeff+c_off1),dx,dy,dz, &
             RH(p_boyspar_ds),RH(p_boystab_ds),RH(p_recurtab_ds), &
             RH(p_scratch),ene)
             ene_ds = ene_ds + ene
          call GSITESITE_prim_group_prims_ene(p_mcmd_dd,  &
             numd,RH(p_coeff_d),RH(p_expd), &
             expon1,RH(p_G_hermite_coeff+c_off1),dx,dy,dz, &
             RH(p_boyspar_dd),RH(p_boystab_dd),RH(p_recurtab_dd), &
             RH(p_scratch),ene)
             ene_dd = ene_dd + ene
        endif
      enddo
      write(16,*)'site1,site2,ene_ss,ene_sd,ene_ds,ene_dd = ',  &
                  site1,site2,ene_ss,ene_sd,ene_ds,ene_dd
      tot_ene_ss = tot_ene_ss + ene_ss
      tot_ene_sd = tot_ene_sd + ene_sd
      tot_ene_ds = tot_ene_ds + ene_ds
      tot_ene_dd = tot_ene_dd + ene_dd
    enddo
  enddo
  tot_ene = tot_ene_ss + tot_ene_sd + tot_ene_ds + tot_ene_dd
  write(6,*)'tot_ene_ss,tot_ene_sd,tot_ene_ds,tot_ene_dd,tot_ene = ', &
             tot_ene_ss,tot_ene_sd,tot_ene_ds,tot_ene_dd,tot_ene
  return
end subroutine GSITESITE_interact
!---------------------------------------------------------
subroutine GSITESITE_prim_group_prims_ene(p_mcmd,  &
             numgroup,coeff_group,expon_group, &
             expon_prim,coeff_prim,dx,dy,dz, &
             boyspar,boys_table,recur_table,scratch,ene)
  ! calculate energy between a single prim of site1 and the collected prims
  ! of a given type (s or d) of site2
  implicit none
  integer p_mcmd,numgroup
  double precision coeff_group(numgroup,*),expon_group(numgroup)
  double precision expon_prim,coeff_prim(*),dx,dy,dz
  double precision boyspar(numgroup),boys_table(numgroup,*),&
                   recur_table(numgroup,*),scratch(numgroup),ene,& !GAC
                   lambda(numgroup)

  if ( numgroup == 0 )then
    ene = 0.d0
    return
  endif
  write(6,*)'inside GSITESITE_prim_group_prims_ene, numgroup = ',numgroup
  write(6,*)'dx,dy,dz = ',dx,dy,dz
  !call GSITESITE_calc_boyspar(numgroup,expon_prim,expon_group,boyspar)
  !call GMCMUR_DAV_recur(p_mcmd,numgroup,boyspar,dx,dy,dz, &
  !                      recur_table,boys_table,scratch)
  call GSITESITE_calc_boyspar(numgroup,expon_prim,expon_group,boyspar,lambda)
  call GMCMUR_DAV_recur(p_mcmd,numgroup,boyspar,dx,dy,dz, &
                        recur_table,boys_table,scratch)
  call GMCMUR_DAV_tensor(p_mcmd,numgroup,recur_table,  &
                         coeff_prim,coeff_group,ene,lambda)
  return
end subroutine GSITESITE_prim_group_prims_ene
!---------------------------------------------------------
subroutine GSITESITE_fillsite_coeffs_expons(nums,numd,numprim,offprim, &
                    coeff,expon,order, &
                    off_coeff,coeff_s,expon_s,coeff_d,expon_d)
  implicit none
  integer nums,numd,numprim,offprim,order(*),off_coeff(*)
  double precision coeff(*),expon(*),coeff_s(nums),expon_s(nums), &
                    coeff_d(numd,10),expon_d(numd)

  integer k,offc,ord,ms,md,j
  double precision expo,boys
  ms = 0
  md = 0
  do k = 1,numprim
    ord = order(offprim+k)
    if ( ord == 1 )then
      ms = ms + 1
      expon_s(ms) = expon(offprim+k)
      offc = off_coeff(offprim+k)
      coeff_s(ms) = coeff(offc+1)
    elseif ( ord == 10 )then
      md = md + 1
      expon_d(md) = expon(offprim+k)
      offc = off_coeff(offprim+k)
      do j = 1,10
        coeff_d(md,j) = coeff(offc+j)
      enddo
    endif 
  enddo
  return
end subroutine GSITESITE_fillsite_coeffs_expons
!---------------------------------------------------------
subroutine GSITESITE_calc_boyspar(mlist,expo1,expo_array,boyspar,lambda)
  implicit none
  integer mlist, i
  double precision expo1,expo_array(mlist),boyspar(mlist),lambda(mlist) ! GAC
  double precision pi, pi52

  integer m

  pi = 4.d0*(atan(1.d0))
  pi52 = pi
  do i = 1, 4
     pi52 = pi52*pi
  enddo
  pi52 = sqrt(pi52)

  do m = 1,mlist
    boyspar(m) = expo1*expo_array(m) / (expo1 + expo_array(m))
    lambda(m) = (pi52*2.d0)/(expo1*expo_array(m)*&
   &                          sqrt(expo1+expo_array(m)))
    !lambda(m) = (pi*2.d0)/(expo1+expo_array(m))
  enddo
  return
end subroutine GSITESITE_calc_boyspar
!---------------------------------------------------------
