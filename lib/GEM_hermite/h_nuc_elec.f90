!---------------------------------------------------------
subroutine H_nuc_elec_AUX(site_info,auxis,ncoeff,cartfit)
  ! GAC: subroutine to calculate ESP with cartesian gaussians.
use definition
  implicit none
  type(sites_info)::site_info
  type(aux_orbitals)::auxis(ncoeff)

  integer site1,i,x1,y1,z1,deg1,deg2,ncoeff
  integer k1,num1,off1,c_off1,n_esp_points,l
  integer order1,order2,allochk
  integer slo1,shi1,ires,nsites,fldint
  double precision prim_prim_EE,tot_esp,tot_E,tot_N,total,result,norm1
  double precision dx,dy,dz,cx,cy,cz,dist,tmp_vec(10),tot_NE
  double precision expon1,expon2,boyspar,nuc_chg1,nuc_chg2,ene,sys_ene
  double precision D(0:3,0:3),E(0:3,0:3),F(0:3,0:3),&
                   DDS(0:3,0:3),EDS(0:3,0:3),FDS(0:3,0:3)
  double precision,parameter::AtoBohr = 0.529177249d0,hartree=627.51d0,&
                   tiny=1.d-12
  logical off_site,cartfit

  fldint = 0
  tot_NE = 0.0d0
  DDS(:,:) = 1.0d0
  EDS(:,:) = 1.0d0
  FDS(:,:) = 1.0d0
  slo1 = site_info%residue_start(1)
  shi1 = site_info%natoms
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
     do site1 = slo1, shi1
        dx = auxis(i)%coords_G(1) - site_info%site_crds(3*(site1-1)+1)
        if (abs(dx) .lt. tiny) dx = 0.0d0
        dy = auxis(i)%coords_G(2) - site_info%site_crds(3*(site1-1)+2)
        if (abs(dy) .lt. tiny) dy = 0.0d0
        dz = auxis(i)%coords_G(3) - site_info%site_crds(3*(site1-1)+3)
        if (abs(dz) .lt. tiny) dz = 0.0d0
        boyspar = expon1
        deg2 = 0.0d0
        expon2 = 3.14159265358979323846d0
        nuc_chg1 = site_info%nuclear_charge(site1)
        call GN_MD_REC_C(deg1,deg2,boyspar,dx,dy,dz,result,x1,y1,z1,&
                         0,0,0,expon1,expon2,D,E,F,DDS,EDS,FDS,fldint)
        if (cartfit) then
           tot_NE = tot_NE - result*norm1*site_info%cartesian_coeffs(i)*&
                             nuc_chg1
        else
           tot_NE = tot_NE - result*site_info%global_hermite_coeffs(i)*&
                             (expon1/expon2)*sqrt(expon1/expon2)*nuc_chg1
        endif
     enddo ! site1
  enddo ! ncoeff
  write(6,*),'Total nuc-elec from AUXILIARY = ',tot_NE
 
  return
end subroutine H_nuc_elec_AUX
!----------------------------------------------------------
subroutine H_nuc_elec(site_info,basis,num_stos)
use definition
  implicit none
  type(sites_info)::site_info
  type(ao_orbitals)::basis(num_stos)

  integer site1,allochk,inttype
  integer k1,num1,off1,c_off1,npa,npb,nx1,ny1,nz1
  integer ncoeff,num_stos,nx2,ny2,nz2
  integer ires,jres,n,i,j,k,deg1,deg2,i2,j2,counter,nx3,ny3,nz3,l,m
  integer slo1,shi1
  double precision alphap,R_ab,brak,pre,temp,px,py,pz,norma,tot_NE
  double precision dx,dy,dz
  double precision x1,x2,y1,y2,z1,z2,x3,y3,z3
  double precision boyspar,result,beta,beta2,boyspar2,&
                   expon3,norm3
  double precision,allocatable::expon1(:),expon2(:),norm1(:),norm2(:),&
                   contr1(:),contr2(:)
  double precision D(0:3,0:3,0:6),E(0:3,0:3,0:6),F(0:3,0:3,0:6),&
                   DDS(0:3,0:3),EDS(0:3,0:3),FDS(0:3,0:3),&
                   densmat(num_stos,num_stos) 
  logical normalize, debug

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
  tot_NE = 0.0d0
  slo1 = site_info%residue_start(1)
  shi1 = site_info%natoms
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
              do site1 = slo1, shi1
                 pre = contr1(i2)*contr2(j2)*norm1(i2)*norm2(j2)*brak*&
                       site_info%nuclear_charge(site1)
                 dx = px - site_info%site_crds(3*(site1-1)+1)
                 dy = py - site_info%site_crds(3*(site1-1)+2)
                 dz = pz - site_info%site_crds(3*(site1-1)+3)
                 expon3 = 3.14159265358979323846d0
                 call GN_MD_REC_J_C(deg1,deg2,boyspar,dx,dy,dz,result,&
                                    nx1,ny1,nz1,nx2,ny2,nz2,nx3,ny3,nz3,&
                                    alphap,expon3,D,E,F,DDS,EDS,FDS,0)
                 temp = temp+result*pre
               enddo ! site1
           enddo ! j2 = 1, npb
        enddo ! i2 = 1, npa
        if (i .eq. j) then
           tot_NE = tot_NE - temp*densmat(i,j)
        else
           tot_NE = tot_NE - temp*densmat(i,j)*2.0d0
        endif
        deallocate(expon2,contr2,norm2)
     enddo ! j = 1, num_stos
     deallocate(expon1,contr1,norm1)
  enddo ! i = 1, num_stos
  write(6,*),'Total nuc-elec EXACT = ',tot_NE
     
  return
end subroutine H_nuc_elec
!---------------------------------------------------------
