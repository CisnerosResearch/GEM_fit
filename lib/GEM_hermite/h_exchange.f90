!----------------------------------------------------------
subroutine H_exchange_CART(site_info,auxis,ncoeff,k_parm,densfile)
  ! GAC: subroutine to calculate <k||l> and <k|l> in local hermite.
use definition
  implicit none
  type(sites_info)::site_info
  type(aux_orbitals)::auxis(ncoeff)

  integer site1,allochk,inttype
  integer k1,num1,off1,c_off1 
  integer ncoeff,x1,x2,y1,y2,z1,z2
  integer ires,jres,n,i,j,deg1,deg2,fldint
  double precision norm1,norm2,exch_res,k_parm,Sij
  double precision dx,dy,dz,aux_chg(ncoeff),aux_coefs(ncoeff)
  double precision expon1,expon2,boyspar,result,beta,beta2,boyspar2
  double precision D(0:3,0:3),E(0:3,0:3),F(0:3,0:3),&
                   DDS(0:3,0:3),EDS(0:3,0:3),FDS(0:3,0:3)
  character(len=*) densfile

  exch_res = 0.0d0
  !aux_coefs(:) = site_info%global_hermite_coeffs(:)
  aux_coefs(:) = site_info%cartesian_coeffs(:)

  print *,'before form auxis'
  call AH_form_auxis(site_info,auxis,ncoeff,densfile)
  print *,'after form auxis'

  do i = 1, ncoeff
     aux_chg(i) = 0.0d0
     expon1 = auxis(i)%expo
     x1 = auxis(i)%x  
     y1 = auxis(i)%y 
     z1 = auxis(i)%z 
     deg1 = x1+y1+z1 
     norm1 = auxis(i)%norm
     call H_form_exp_coef_noS(x1,auxis(i)%coords(1),expon1,D)
     call H_form_exp_coef_noS(y1,auxis(i)%coords(2),expon1,E)
     call H_form_exp_coef_noS(z1,auxis(i)%coords(3),expon1,F)
     do j = 1, i
        expon2 = auxis(j)%expo
        x2 = auxis(j)%x  
        y2 = auxis(j)%y 
        z2 = auxis(j)%z 
        deg2 = x2+y2+z2 
        norm2 = auxis(j)%norm
        call H_form_exp_coef_noS(x2,auxis(j)%coords(1),expon2,DDS)
        call H_form_exp_coef_noS(y2,auxis(j)%coords(2),expon2,EDS)
        call H_form_exp_coef_noS(z2,auxis(j)%coords(3),expon2,FDS)
        dx = auxis(i)%coords(1) - auxis(j)%coords(1)
        dy = auxis(i)%coords(2) - auxis(j)%coords(2)
        dz = auxis(i)%coords(3) - auxis(j)%coords(3)
        boyspar = expon1*expon2 / (expon1 + expon2)
        call GN_MD_REC_OVERLAP_C(deg1,deg2,boyspar,dx,dy,dz,result,&
                                 x1,y1,z1,x2,y2,z2,expon1,expon2,D,E,F,&
                                 DDS,EDS,FDS)
        aux_chg(i) = aux_chg(i) + result*norm1*norm2*aux_coefs(j)*k_parm
     enddo
     exch_res = exch_res + aux_chg(i)*aux_coefs(i)
  enddo
  
  write(6,*)'total exchange = ',exch_res
  return
end subroutine H_exchange_CART
!----------------------------------------------------------
subroutine H_exchange(site_info,k_parm,tot_exch)
use definition
  implicit none
  type(sites_info)::site_info

  integer site1,site2
  integer k1,num1,off1,c_off1,k2,num2,off2,c_off2
  integer order1,order2
  integer num_res,slo1,shi1,slo2,shi2,ires,jres
  double precision tot_exch,Sij,k_parm,tmp_exch
  double precision dx,dy,dz,tmp_vec(100)
  double precision expon1,expon2,boyspar,hartree,ene,sys_ene

  hartree = 627.51d0
  tot_exch = 0.0d0

  num_res = site_info%num_residues

  do ires = 1,num_res-1
    slo1 = site_info%residue_start(ires)
    shi1 = site_info%residue_end(ires)
    do jres = ires+1,num_res
      slo2 = site_info%residue_start(jres)
      shi2 = site_info%residue_end(jres)
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
          tmp_exch = 0.0d0
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
                call GN_MCMUR_DAV_S_ss(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              elseif ( order1 == 1 .and. order2 == 4 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_S_sp(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              elseif ( order1 == 1 .and. order2 == 6 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_S_sda(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              elseif ( order1 == 1 .and. order2 == 10 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_S_sd(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              elseif ( order1 == 4 .and. order2 == 1 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_S_ps(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              elseif ( order1 == 4 .and. order2 == 4 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_S_pp(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              elseif ( order1 == 6 .and. order2 == 1 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_S_das(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              elseif ( order1 == 10 .and. order2 == 1 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_S_ds(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              elseif ( order1 == 6 .and. order2 == 6 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_S_dada(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              elseif ( order1 == 6 .and. order2 == 10 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_S_dad(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              elseif ( order1 == 10 .and. order2 == 6 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_S_dda(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              elseif ( order1 == 10 .and. order2 == 10 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_S_dd(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              ! GAC : consider p shells only, with spd shells
              elseif ( (order1 == 4 .and. order2 == 10) .or. &
                       (order1 == 10 .and. order2 == 4) ) then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_S_RECp(order1,order2,boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              endif
              !tmp_exch = tmp_exch + k_parm*Sij
              tmp_exch = tmp_exch + Sij
            enddo !k2 = 1,num2
          enddo!k1 = 1,num1
          tot_exch = tot_exch + tmp_exch
        enddo!site2 = slo2,shi2
      enddo!site1 = slo1,slo2
    enddo !jres = ires+1,num_res
  enddo! ires = 1,num_res-1
  !tot_exch = k_parm*(tot_exch**.98)
  tot_exch = k_parm*tot_exch
  write(6,*)'total exchange (kcal/mol) = ',tot_exch*hartree
  return
end subroutine H_exchange
!----------------------------------------------------------
subroutine H_exchange2(site_info,site_info2,k_parm,tot_exch)
use definition
  implicit none
  type(sites_info)::site_info,site_info2

  integer site1,site2
  integer k1,num1,off1,c_off1,k2,num2,off2,c_off2
  integer order1,order2
  integer slo1,shi1,slo2,shi2,ires,jres
  double precision Sij,tot_exch,k_parm,tmp_exch
  double precision dx,dy,dz,tmp_vec(100)
  double precision expon1,expon2,boyspar,nuc_chg1,nuc_chg2,hartree,ene,sys_ene

  hartree = 627.51d0

  tot_exch = 0.d0
  do ires = 1,site_info%num_residues
    slo1 = site_info%residue_start(ires)
    shi1 = site_info%residue_end(ires)
    do jres = 1,site_info2%num_residues
      slo2 = site_info2%residue_start(jres)
      shi2 = site_info2%residue_end(jres)
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
          tmp_exch = 0.d0
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
                call GN_MCMUR_DAV_S_ss(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              elseif ( order1 == 1 .and. order2 == 4 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_S_sp(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              elseif ( order1 == 1 .and. order2 == 6 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_S_sda(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              elseif ( order1 == 1 .and. order2 == 10 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_S_sd(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              elseif ( order1 == 4 .and. order2 == 1 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_S_ps(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              elseif ( order1 == 4 .and. order2 == 4 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_S_pp(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              elseif ( order1 == 6 .and. order2 == 1 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_S_das(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              elseif ( order1 == 10 .and. order2 == 1 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_S_ds(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              elseif ( order1 == 6 .and. order2 == 6 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_S_dada(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              elseif ( order1 == 6 .and. order2 == 10 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_S_dad(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              elseif ( order1 == 10 .and. order2 == 6 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_S_dda(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              elseif ( order1 == 10 .and. order2 == 10 )then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_S_dd(boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              ! GAC : consider p shells only, with spd shells
              elseif ( (order1 == 4 .and. order2 == 10) .or. &
                       (order1 == 10 .and. order2 == 4) ) then
                boyspar = expon1*expon2 / (expon1 + expon2)
                call GN_MCMUR_DAV_S_RECp(order1,order2,boyspar,dx,dy,dz,  &
                    site_info%global_hermite_coeffs(c_off1+1),  &
                    site_info2%global_hermite_coeffs(c_off2+1),Sij,&
                    tmp_vec)
              endif
              tmp_exch = tmp_exch + Sij*k_parm
            enddo !k2 = 1,num2
          enddo!k1 = 1,num1
        tot_exch = tot_exch + tmp_exch
        enddo!site2 = slo2,shi2
      enddo!site1 = slo1,slo2
    enddo !jres = 1,site_info2%num_res
  enddo! ires = 1,site_info%num_res
  write(6,*)'total exchange (kcal/mol) = ',tot_exch*hartree
  return
end subroutine H_exchange2
!----------------------------------------------------------
subroutine H_qmmm_exch(site_info,site_info2,basis,qmmmexch,k_parm,num_stos)
use definition
  implicit none
  type(sites_info)::site_info
  type(sites_info)::site_info2
  type(ao_orbitals)::basis(num_stos)

  integer site1,allochk,inttype
  integer k1,num1,off1,c_off,npa,npb,nx1,ny1,nz1
  integer ncoeff,num_stos,nx2,ny2,nz2,slo1,shi1,order1
  integer ires,jres,n,i,j,k,deg1,deg2,i2,j2,counter,nx3,ny3,nz3,l,m
  double precision alphap,R_ab,brak,pre,temp,px,py,pz,norma,esp_E
  double precision dx,dy,dz,dist,cx,cy,cz,k_parm,qmmmexch
  double precision x1,x2,y1,y2,z1,z2,x3,y3,z3,nuc_chg1,nuc_chg2
  double precision boyspar,result,beta,beta2,boyspar2,&
                   expon3,norm3,nuc_chg_E, nuc_el_E, nuc_nuc_E, el_nuc_E
  double precision Vmat(num_stos,num_stos),densmat(num_stos,num_stos)
  double precision,allocatable::expon1(:),expon2(:),norm1(:),norm2(:),&
                   contr1(:),contr2(:)
  double precision D(0:3,0:3,0:6),E(0:3,0:3,0:6),F(0:3,0:3,0:6)
  double precision,parameter::AtoBohr = 0.529177249d0,tiny=1.d-12

  call HLOAD_dens(site_info,'gaudens.dat',.false.)
  counter = 0
  do i = 1, num_stos
     do j = 1, i
        counter = counter+1
        densmat(j,i) = site_info%densmat(counter)
        densmat(i,j) = densmat(j,i)
     enddo
  enddo

  qmmmexch = 0.0d0
  expon3 = 3.14159265358979323846d0
! calculate H_core for electron-aux
  do i = 1, num_stos
     npa = basis(i)%deg_contr  
     allocate(expon1(npa),contr1(npa),norm1(npa), stat = allochk)
     if (allochk .gt. 0 ) then
        write(6,*)'H_V_mat2:could not allocate expon1, contr1, norm1;&
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
           write(6,*)'H_V_mat2:could not allocate expon2, contr2, norm2;&
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
              do ires = 1,site_info2%num_residues
                 slo1 = site_info2%residue_start(ires)
                 shi1 = site_info2%residue_end(ires)
                 do k = slo1, shi1
                    num1 = site_info2%num_primitives(k)
                    off1 = site_info2%off_primitives(k)
                    do k1 = 1, num1
                       order1 = site_info2%prim_order(off1+k1)
                       expon3 = site_info2%prim_expo(off1+k1)
                       c_off = site_info2%coeff_offset(off1+k1)
                       if (order1 ==1) then
                          nx3 = 0
                          ny3 = 0
                          nz3 = 0 
                          deg2 = nx3+ny3+nz3
                          dx = px - site_info2%site_crds(3*(k-1)+1) 
                          dy = py - site_info2%site_crds(3*(k-1)+2)
                          dz = pz - site_info2%site_crds(3*(k-1)+3)
                          boyspar = alphap*expon3 / (alphap + expon3)
                          pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*&
                                brak*site_info2%global_hermite_coeffs(c_off+1)
                          call GN_MD_REC_OVERLAP_J(deg1,deg2,boyspar,dx,dy,&
                                                   dz,result,nx1,ny1,nz1,nx2,&
                                                   ny2,nz2,nx3,ny3,nz3,&
                                                   alphap,expon3,D,E,F)
                          temp = temp + result*pre
                       else if (order1 == 4) then
                          nx3 = 0
                          ny3 = 0
                          nz3 = 0 
                          deg2 = nx3+ny3+nz3
                          dx = px - site_info2%site_crds(3*(k-1)+1) 
                          dy = py - site_info2%site_crds(3*(k-1)+2)
                          dz = pz - site_info2%site_crds(3*(k-1)+3)
                          boyspar = alphap*expon3 / (alphap + expon3)
                          pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*&
                                brak*site_info2%global_hermite_coeffs(c_off+1)
                          call GN_MD_REC_OVERLAP_J(deg1,deg2,boyspar,dx,dy,&
                                                   dz,result,nx1,ny1,nz1,nx2,&
                                                   ny2,nz2,nx3,ny3,nz3,&
                                                   alphap,expon3,D,E,F)
                          temp = temp + result*pre
                          nx3 = 1
                          ny3 = 0
                          nz3 = 0 
                          deg2 = nx3+ny3+nz3
                          dx = px - site_info2%site_crds(3*(k-1)+1) 
                          dy = py - site_info2%site_crds(3*(k-1)+2)
                          dz = pz - site_info2%site_crds(3*(k-1)+3)
                          boyspar = alphap*expon3 / (alphap + expon3)
                          pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*&
                                brak*site_info2%global_hermite_coeffs(c_off+2)
                          call GN_MD_REC_OVERLAP_J(deg1,deg2,boyspar,dx,dy,&
                                                   dz,result,nx1,ny1,nz1,nx2,&
                                                   ny2,nz2,nx3,ny3,nz3,&
                                                   alphap,expon3,D,E,F)
                          temp = temp + result*pre
                          nx3 = 0
                          ny3 = 1
                          nz3 = 0 
                          deg2 = nx3+ny3+nz3
                          dx = px - site_info2%site_crds(3*(k-1)+1) 
                          dy = py - site_info2%site_crds(3*(k-1)+2)
                          dz = pz - site_info2%site_crds(3*(k-1)+3)
                          boyspar = alphap*expon3 / (alphap + expon3)
                          pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*&
                                brak*site_info2%global_hermite_coeffs(c_off+3)
                          call GN_MD_REC_OVERLAP_J(deg1,deg2,boyspar,dx,dy,&
                                                   dz,result,nx1,ny1,nz1,nx2,&
                                                   ny2,nz2,nx3,ny3,nz3,&
                                                   alphap,expon3,D,E,F)
                          temp = temp + result*pre
                          nx3 = 0
                          ny3 = 0
                          nz3 = 1 
                          deg2 = nx3+ny3+nz3
                          dx = px - site_info2%site_crds(3*(k-1)+1) 
                          dy = py - site_info2%site_crds(3*(k-1)+2)
                          dz = pz - site_info2%site_crds(3*(k-1)+3)
                          boyspar = alphap*expon3 / (alphap + expon3)
                          pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*&
                                brak*site_info2%global_hermite_coeffs(c_off+4)
                          call GN_MD_REC_OVERLAP_J(deg1,deg2,boyspar,dx,dy,&
                                                   dz,result,nx1,ny1,nz1,nx2,&
                                                   ny2,nz2,nx3,ny3,nz3,&
                                                   alphap,expon3,D,E,F)
                          temp = temp + result*pre
                       else if (order1 == 6) then
                          nx3 = 2
                          ny3 = 0
                          nz3 = 0 
                          deg2 = nx3+ny3+nz3
                          dx = px - site_info2%site_crds(3*(k-1)+1) 
                          dy = py - site_info2%site_crds(3*(k-1)+2)
                          dz = pz - site_info2%site_crds(3*(k-1)+3)
                          boyspar = alphap*expon3 / (alphap + expon3)
                          pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*&
                                brak*site_info2%global_hermite_coeffs(c_off+1)
                          call GN_MD_REC_OVERLAP_J(deg1,deg2,boyspar,dx,dy,&
                                                   dz,result,nx1,ny1,nz1,nx2,&
                                                   ny2,nz2,nx3,ny3,nz3,&
                                                   alphap,expon3,D,E,F)
                          temp = temp + result*pre
                          nx3 = 0
                          ny3 = 2
                          nz3 = 0 
                          deg2 = nx3+ny3+nz3
                          dx = px - site_info2%site_crds(3*(k-1)+1) 
                          dy = py - site_info2%site_crds(3*(k-1)+2)
                          dz = pz - site_info2%site_crds(3*(k-1)+3)
                          boyspar = alphap*expon3 / (alphap + expon3)
                          pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*&
                                brak*site_info2%global_hermite_coeffs(c_off+2)
                          call GN_MD_REC_OVERLAP_J(deg1,deg2,boyspar,dx,dy,&
                                                   dz,result,nx1,ny1,nz1,nx2,&
                                                   ny2,nz2,nx3,ny3,nz3,&
                                                   alphap,expon3,D,E,F)
                          temp = temp + result*pre
                          nx3 = 0
                          ny3 = 0
                          nz3 = 2 
                          deg2 = nx3+ny3+nz3
                          dx = px - site_info2%site_crds(3*(k-1)+1) 
                          dy = py - site_info2%site_crds(3*(k-1)+2)
                          dz = pz - site_info2%site_crds(3*(k-1)+3)
                          boyspar = alphap*expon3 / (alphap + expon3)
                          pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*&
                                brak*site_info2%global_hermite_coeffs(c_off+3)
                          call GN_MD_REC_OVERLAP_J(deg1,deg2,boyspar,dx,dy,&
                                                   dz,result,nx1,ny1,nz1,nx2,&
                                                   ny2,nz2,nx3,ny3,nz3,&
                                                   alphap,expon3,D,E,F)
                          temp = temp + result*pre
                          nx3 = 1
                          ny3 = 1
                          nz3 = 0 
                          deg2 = nx3+ny3+nz3
                          dx = px - site_info2%site_crds(3*(k-1)+1) 
                          dy = py - site_info2%site_crds(3*(k-1)+2)
                          dz = pz - site_info2%site_crds(3*(k-1)+3)
                          boyspar = alphap*expon3 / (alphap + expon3)
                          pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*&
                                brak*site_info2%global_hermite_coeffs(c_off+4)
                          call GN_MD_REC_OVERLAP_J(deg1,deg2,boyspar,dx,dy,&
                                                   dz,result,nx1,ny1,nz1,nx2,&
                                                   ny2,nz2,nx3,ny3,nz3,&
                                                   alphap,expon3,D,E,F)
                          temp = temp + result*pre
                          nx3 = 1
                          ny3 = 0
                          nz3 = 1 
                          deg2 = nx3+ny3+nz3
                          dx = px - site_info2%site_crds(3*(k-1)+1) 
                          dy = py - site_info2%site_crds(3*(k-1)+2)
                          dz = pz - site_info2%site_crds(3*(k-1)+3)
                          boyspar = alphap*expon3 / (alphap + expon3)
                          pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*&
                                brak*site_info2%global_hermite_coeffs(c_off+5)
                          call GN_MD_REC_OVERLAP_J(deg1,deg2,boyspar,dx,dy,&
                                                   dz,result,nx1,ny1,nz1,nx2,&
                                                   ny2,nz2,nx3,ny3,nz3,&
                                                   alphap,expon3,D,E,F)
                          temp = temp + result*pre
                          nx3 = 0
                          ny3 = 1
                          nz3 = 1 
                          deg2 = nx3+ny3+nz3
                          dx = px - site_info2%site_crds(3*(k-1)+1) 
                          dy = py - site_info2%site_crds(3*(k-1)+2)
                          dz = pz - site_info2%site_crds(3*(k-1)+3)
                          boyspar = alphap*expon3 / (alphap + expon3)
                          pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*&
                                brak*site_info2%global_hermite_coeffs(c_off+6)
                          call GN_MD_REC_OVERLAP_J(deg1,deg2,boyspar,dx,dy,&
                                                   dz,result,nx1,ny1,nz1,nx2,&
                                                   ny2,nz2,nx3,ny3,nz3,&
                                                   alphap,expon3,D,E,F)
                          temp = temp + result*pre
                       else if (order1 == 10) then
                          nx3 = 0
                          ny3 = 0
                          nz3 = 0 
                          deg2 = nx3+ny3+nz3
                          dx = px - site_info2%site_crds(3*(k-1)+1) 
                          dy = py - site_info2%site_crds(3*(k-1)+2)
                          dz = pz - site_info2%site_crds(3*(k-1)+3)
                          boyspar = alphap*expon3 / (alphap + expon3)
                          pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*&
                                brak*site_info2%global_hermite_coeffs(c_off+1)
                          call GN_MD_REC_OVERLAP_J(deg1,deg2,boyspar,dx,dy,&
                                                   dz,result,nx1,ny1,nz1,nx2,&
                                                   ny2,nz2,nx3,ny3,nz3,&
                                                   alphap,expon3,D,E,F)
                          temp = temp + result*pre
                          nx3 = 1
                          ny3 = 0
                          nz3 = 0 
                          deg2 = nx3+ny3+nz3
                          dx = px - site_info2%site_crds(3*(k-1)+1) 
                          dy = py - site_info2%site_crds(3*(k-1)+2)
                          dz = pz - site_info2%site_crds(3*(k-1)+3)
                          boyspar = alphap*expon3 / (alphap + expon3)
                          pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*&
                                brak*site_info2%global_hermite_coeffs(c_off+2)
                          call GN_MD_REC_OVERLAP_J(deg1,deg2,boyspar,dx,dy,&
                                                   dz,result,nx1,ny1,nz1,nx2,&
                                                   ny2,nz2,nx3,ny3,nz3,&
                                                   alphap,expon3,D,E,F)
                          temp = temp + result*pre
                          nx3 = 0
                          ny3 = 1
                          nz3 = 0 
                          deg2 = nx3+ny3+nz3
                          dx = px - site_info2%site_crds(3*(k-1)+1) 
                          dy = py - site_info2%site_crds(3*(k-1)+2)
                          dz = pz - site_info2%site_crds(3*(k-1)+3)
                          boyspar = alphap*expon3 / (alphap + expon3)
                          pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*&
                                brak*site_info2%global_hermite_coeffs(c_off+3)
                          call GN_MD_REC_OVERLAP_J(deg1,deg2,boyspar,dx,dy,&
                                                   dz,result,nx1,ny1,nz1,nx2,&
                                                   ny2,nz2,nx3,ny3,nz3,&
                                                   alphap,expon3,D,E,F)
                          temp = temp + result*pre
                          nx3 = 0
                          ny3 = 0
                          nz3 = 1 
                          deg2 = nx3+ny3+nz3
                          dx = px - site_info2%site_crds(3*(k-1)+1) 
                          dy = py - site_info2%site_crds(3*(k-1)+2)
                          dz = pz - site_info2%site_crds(3*(k-1)+3)
                          boyspar = alphap*expon3 / (alphap + expon3)
                          pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*&
                                brak*site_info2%global_hermite_coeffs(c_off+4)
                          call GN_MD_REC_OVERLAP_J(deg1,deg2,boyspar,dx,dy,&
                                                   dz,result,nx1,ny1,nz1,nx2,&
                                                   ny2,nz2,nx3,ny3,nz3,&
                                                   alphap,expon3,D,E,F)
                          temp = temp + result*pre
                          nx3 = 2
                          ny3 = 0
                          nz3 = 0 
                          deg2 = nx3+ny3+nz3
                          dx = px - site_info2%site_crds(3*(k-1)+1) 
                          dy = py - site_info2%site_crds(3*(k-1)+2)
                          dz = pz - site_info2%site_crds(3*(k-1)+3)
                          boyspar = alphap*expon3 / (alphap + expon3)
                          pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*&
                                brak*site_info2%global_hermite_coeffs(c_off+5)
                          call GN_MD_REC_OVERLAP_J(deg1,deg2,boyspar,dx,dy,&
                                                   dz,result,nx1,ny1,nz1,nx2,&
                                                   ny2,nz2,nx3,ny3,nz3,&
                                                   alphap,expon3,D,E,F)
                          temp = temp + result*pre
                          nx3 = 0
                          ny3 = 2
                          nz3 = 0 
                          deg2 = nx3+ny3+nz3
                          dx = px - site_info2%site_crds(3*(k-1)+1) 
                          dy = py - site_info2%site_crds(3*(k-1)+2)
                          dz = pz - site_info2%site_crds(3*(k-1)+3)
                          boyspar = alphap*expon3 / (alphap + expon3)
                          pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*&
                                brak*site_info2%global_hermite_coeffs(c_off+6)
                          call GN_MD_REC_OVERLAP_J(deg1,deg2,boyspar,dx,dy,&
                                                   dz,result,nx1,ny1,nz1,nx2,&
                                                   ny2,nz2,nx3,ny3,nz3,&
                                                   alphap,expon3,D,E,F)
                          temp = temp + result*pre
                          nx3 = 0
                          ny3 = 0
                          nz3 = 2 
                          deg2 = nx3+ny3+nz3
                          dx = px - site_info2%site_crds(3*(k-1)+1) 
                          dy = py - site_info2%site_crds(3*(k-1)+2)
                          dz = pz - site_info2%site_crds(3*(k-1)+3)
                          boyspar = alphap*expon3 / (alphap + expon3)
                          pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*&
                                brak*site_info2%global_hermite_coeffs(c_off+7)
                          call GN_MD_REC_OVERLAP_J(deg1,deg2,boyspar,dx,dy,&
                                                   dz,result,nx1,ny1,nz1,nx2,&
                                                   ny2,nz2,nx3,ny3,nz3,&
                                                   alphap,expon3,D,E,F)
                          temp = temp + result*pre
                          nx3 = 1
                          ny3 = 1
                          nz3 = 0 
                          deg2 = nx3+ny3+nz3
                          dx = px - site_info2%site_crds(3*(k-1)+1) 
                          dy = py - site_info2%site_crds(3*(k-1)+2)
                          dz = pz - site_info2%site_crds(3*(k-1)+3)
                          boyspar = alphap*expon3 / (alphap + expon3)
                          pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*&
                                brak*site_info2%global_hermite_coeffs(c_off+8)
                          call GN_MD_REC_OVERLAP_J(deg1,deg2,boyspar,dx,dy,&
                                                   dz,result,nx1,ny1,nz1,nx2,&
                                                   ny2,nz2,nx3,ny3,nz3,&
                                                   alphap,expon3,D,E,F)
                          temp = temp + result*pre
                          nx3 = 1
                          ny3 = 0
                          nz3 = 1 
                          deg2 = nx3+ny3+nz3
                          dx = px - site_info2%site_crds(3*(k-1)+1) 
                          dy = py - site_info2%site_crds(3*(k-1)+2)
                          dz = pz - site_info2%site_crds(3*(k-1)+3)
                          boyspar = alphap*expon3 / (alphap + expon3)
                          pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*&
                                brak*site_info2%global_hermite_coeffs(c_off+9)
                          call GN_MD_REC_OVERLAP_J(deg1,deg2,boyspar,dx,dy,&
                                                   dz,result,nx1,ny1,nz1,nx2,&
                                                   ny2,nz2,nx3,ny3,nz3,&
                                                   alphap,expon3,D,E,F)
                          temp = temp + result*pre
                          nx3 = 0
                          ny3 = 1
                          nz3 = 1 
                          deg2 = nx3+ny3+nz3
                          dx = px - site_info2%site_crds(3*(k-1)+1) 
                          dy = py - site_info2%site_crds(3*(k-1)+2)
                          dz = pz - site_info2%site_crds(3*(k-1)+3)
                          boyspar = alphap*expon3 / (alphap + expon3)
                          pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*&
                                brak*site_info2%global_hermite_coeffs(c_off+10)
                          call GN_MD_REC_OVERLAP_J(deg1,deg2,boyspar,dx,dy,&
                                                   dz,result,nx1,ny1,nz1,nx2,&
                                                   ny2,nz2,nx3,ny3,nz3,&
                                                   alphap,expon3,D,E,F)
                          temp = temp + result*pre
                      endif
                    enddo! k1 = 1, num1
                 enddo ! k = 1, ncoeff
              enddo ! ires 1, num_res
           enddo ! j2 = 1, npb
        enddo ! i2 = 1, npa
        if (i .eq. j) then
           qmmmexch = qmmmexch + temp*densmat(i,j)
        else
           qmmmexch = qmmmexch + temp*densmat(i,j)*2.0d0
        endif
        deallocate(expon2,contr2,norm2)
     enddo ! j = i, num_stos
     deallocate(expon1,contr1,norm1)
  enddo ! i = 1, num_stos
  qmmmexch = qmmmexch*k_parm

end subroutine H_qmmm_exch
!----------------------------------------------------------
