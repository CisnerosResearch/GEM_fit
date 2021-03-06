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
!------------------------------------------------------
subroutine CH_COEFF_mpoles_NEW(atom_auxis,site_info,filename,filename2)
use definition
  implicit none
  !include "definition.fh"
  integer p_site_info
  type(auxiliary_basis)::atom_auxis
  type(sites_info)::site_info
  character(len=*) filename,filename2
  character(len=45) filename3,filename4
  integer ptr,nsites,ps_numprim,ps_offprim,ps_primorder,  &
          ps_prim_off_coeff,ncoefftot,p_L_hermite_coeff,ps_atomic_num, &
          p_M_hermite_coeff,ps_primexpo
  double precision lpole(10),gpole(3),finaldip,fac,pi,sqpi,origin(3),bohr
  integer num,off,order,scoff,i,j,k,m,n,allochk
  integer p_chg,p_dip,p_sec_mom
  integer p_Lchg,p_Ldip,p_Lsec_mom
  integer imol,num_mols,nlo,nhi,p_molptr,p_molsitecrd
  double precision,allocatable::Lchg(:),Ldip(:),Lsec_mom(:),trace(:),Gdip(:)
  double precision,allocatable::Lsec_mom_TR(:)
  double precision,allocatable::nucdip(:),elecdip(:),chgcontr(:),&
                                atmdip(:),atmchg(:),totdip(:)
  double precision alpha

  filename3=trim(filename)//'.TINK'
  filename4=trim(filename)//'.AMB'
  !open(unit=14,file=filename,status='unknown')
  open(unit=15,file=filename2,status='unknown')
  open(unit=16,file=filename3,status='unknown')
  !open(unit=17,file=filename4,status='unknown')

  nsites = site_info%nsites

  allocate(Lchg(nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'CH_COEFF_mpoles_NEW:could not allocate Lchg, exiting'
     stop
  endif
  allocate(Ldip(3*nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'CH_COEFF_mpoles_NEW:could not allocate Ldip, exiting'
     stop
  endif
  allocate(Gdip(3*nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'CH_COEFF_mpoles_NEW:could not allocate Gdip, exiting'
     stop
  endif
  allocate(Lsec_mom(6*nsites),Lsec_mom_TR(6*nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'CH_COEFF_mpoles_NEW: could not allocate Lsec_mom, exiting'
     stop
  endif
  allocate(trace(nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'CH_COEFF_mpoles_NEW: could not allocate trace, exiting'
     stop
  endif
  allocate(nucdip(3), elecdip(3), chgcontr(3), totdip(3),&
 &         atmchg(nsites), atmdip(3*nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'CH_COEFF_mpoles_NEW:could not allocate nucdip or elecdip&
    & or atmchg or totdip or chgcontr or atmdip, exiting'
     stop
  endif

  trace(:) = 0.0d0
  nucdip(:) = 0.0d0
  elecdip(:) = 0.0d0
  chgcontr(:) = 0.0d0
  atmdip(:) = 0.0d0
  totdip(:) = 0.0d0
  origin(:) = 0.0d0
  atmchg(:) = 0.0d0
  pi = 4.d0*(atan(1.0d0))
  sqpi = dsqrt(pi)
  fac = 0.529177249d0*4.8032068d-10*1.d10
  bohr = 0.529177249d0

  num_mols = 1

  do n = 1,nsites
    num = site_info%num_primitives(n)
    off = site_info%off_primitives(n)

    gpole(:) = 0.d0
    do k = 1,10
      lpole(k) = 0.d0
    enddo
    do k = 1,num
      order = site_info%prim_order(off+k)

      scoff = site_info%coeff_offset(off+k)

      alpha = site_info%prim_expo(off+k)

      lpole(1) = lpole(1) + site_info%local_hermite_coeffs(scoff+1)

      atmchg(n) = atmchg(n) + site_info%aux_elec(scoff+1)*&
                              site_info%cartesian_coeffs(scoff+1)
      
      if ( order == 4 )then
        do j = 2,4
          lpole(j) = lpole(j) + site_info%local_hermite_coeffs(scoff+j)
          gpole(j-1) = gpole(j-1) + site_info%cartesian_coeffs(scoff+j)*&
                                   (site_info%cart_coeff_norms(scoff+j)*&
                                    pi*sqpi)/(2.d0*alpha**2.5)
          atmchg(n) = atmchg(n) + site_info%aux_elec(scoff+j)*&
                                  site_info%cartesian_coeffs(scoff+j)
        enddo
      else if ( order == 6 )then
        do j = 2,4
          lpole(j) = lpole(j) + 2.d0*site_info%local_hermite_coeffs(scoff+j)
          atmchg(n) = atmchg(n) + site_info%aux_elec(scoff+j)*&
                                  site_info%cartesian_coeffs(scoff+j)
        enddo
        do j = 5,7
          lpole(j) = lpole(j) + site_info%local_hermite_coeffs(scoff+j)
          atmchg(n) = atmchg(n) + site_info%aux_elec(scoff+j)*&
                                  site_info%cartesian_coeffs(scoff+j)
        enddo
      else if ( order == 10 )then
        do j = 2,4
          lpole(j) = lpole(j) + site_info%local_hermite_coeffs(scoff+j)
          gpole(j-1) = gpole(j-1) + site_info%cartesian_coeffs(scoff+j)*&
                                   (site_info%cart_coeff_norms(scoff+j)*&
                                    pi*sqpi)/(2.d0*alpha**2.5)
          atmchg(n) = atmchg(n) + site_info%aux_elec(scoff+j)*&
                                  site_info%cartesian_coeffs(scoff+j)
        enddo
        do j = 5,7
          lpole(j) = lpole(j) + 2.d0*site_info%local_hermite_coeffs(scoff+j)
          atmchg(n) = atmchg(n) + site_info%aux_elec(scoff+j)*&
                                  site_info%cartesian_coeffs(scoff+j)
        enddo
        do j = 8,10
          lpole(j) = lpole(j) + site_info%local_hermite_coeffs(scoff+j)
          atmchg(n) = atmchg(n) + site_info%aux_elec(scoff+j)*&
                                  site_info%cartesian_coeffs(scoff+j)
        enddo
      endif
      ! next add charge contribution to xx,yy,zz moments
      do j = 5,7
        lpole(j) = lpole(j) + site_info%local_hermite_coeffs(scoff+1)&
                              / (2.d0*alpha)
      enddo
    enddo
    !RH(p_Lchg+n-1) = lpole(1)
    Lchg(n) = lpole(1)
    do j = 1,3
      !RH(p_Ldip+3*(n-1)+j-1) = lpole(j+1)
      Ldip(3*(n-1)+j) = lpole(j+1)
      Gdip(3*(n-1)+j) = gpole(j)
    enddo
    do j = 1,6
      !RH(p_Lsec_mom+6*(n-1)+j-1) = lpole(j+4)
      Lsec_mom(6*(n-1)+j) = lpole(j+4)
      Lsec_mom_TR(6*(n-1)+j) = lpole(j+4)
      if (j .le. 3) then
         trace(n) = trace(n) + Lsec_mom(6*(n-1)+j)
      endif
    enddo
    trace(n) = (1.d0/3.d0)*trace(n)
  enddo
  !Lsec_mom_TR(:)=Lsec_mom(:)

  print *,'Adding trace to quadrupoles'
  do n = 1, nsites
    do j = 1,3
      Lsec_mom(6*(n-1)+j) = 0.5*(Lsec_mom(6*(n-1)+j)-trace(n))
    enddo
    !print *,'atmchg = ',atmchg(n)
  enddo

  ! calculating atomic dipoles for polarizabilities and total dipole
100 format(20x,'Dipole moment in Debye (w.r.t. origin)')
110 format(21X,'x',15X,'y',15X,'z'&
     &/1X,'Electronic',3(2X,F14.8)&
     &/1X,'Nuclear',3X,3(2X,F14.8)&
     &/1X,'Elec Chg',2X,3(2X,F14.8)&
     &/1X,'Total',5X,3(2X,F14.8)&
     &/1X,'Dipole Moment',F15.8)
102 format(20x,'Atomic Dipoles in Debye'&
          /21X,'x',15X,'y',15X,'z')
112 format(1X,'site ',i3,2x,3(2X,F14.8))

  do j = 1,3
     do n = 1,nsites
        atmdip(3*(n-1)+j) = site_info%nuclear_charge(n)*&   !nuc chg contr.
                            (site_info%site_crds(3*(n-1)+j)-origin(j))&
                            - Gdip(3*(n-1)+j) -&            !elec dip
                            atmchg(n)*&                     !chg contr
                            (site_info%site_crds(3*(n-1)+j)-origin(j))
        chgcontr(j) = chgcontr(j) + atmchg(n)*&
                                    (site_info%site_crds(3*(n-1)+j)&
                                     -origin(j))
        nucdip(j) = nucdip(j) + site_info%nuclear_charge(n)*&
                                (site_info%site_crds(3*(n-1)+j)-origin(j))
        elecdip(j) = elecdip(j) + Gdip(3*(n-1)+j)
     enddo
     totdip(j) = totdip(j) + nucdip(j) - elecdip(j) - chgcontr(j)
  enddo
  elecdip(:) = elecdip(:)*fac
  nucdip(:) = nucdip(:)*fac
  chgcontr(:) = chgcontr(:)*fac
  totdip(:) = totdip(:)*fac
  atmdip(:) = atmdip(:)*fac
     
  finaldip = dsqrt(totdip(1)*totdip(1)+totdip(2)*totdip(2)+&
                   totdip(3)*totdip(3))

  write(*,102)
  do n = 1, nsites
     write(*,112) n,(atmdip(3*(n-1)+j),j=1,3)
  enddo

  print *,''
  write(*,100)
  write(*,110) elecdip,nucdip,chgcontr,totdip,finaldip


  !write(14,*)'=============================================================='
  !write(14,*)filename(1:70)
  !!write(14,*)'=============================================================='
  write(16,114)site_info%natoms
  do imol = 1,num_mols
    nlo = 1
    nhi = nsites
    !!write(14,*)'molecule ',imol,' =========Site props========================='
    write(15,*)nsites*10
    do n = nlo,nhi
      !!write(14,*)'---------- site: ',n,'  ---------------------------'
      !write(14,*)'charge,dip = ',RH(ps_atomic_num+n-1)-RH(p_Lchg+n-1), &
      !                       (-RH(p_Ldip+3*(n-1)+j-1),j=1,3)
      !write(14,*)'sec mom = ',(-RH(p_Lsec_mom+6*(n-1)+j-1),j=1,3)
      !write(14,*)'sec mom = ',(-RH(p_Lsec_mom+6*(n-1)+j-1),j=4,6)
      !!write(14,*)'charge,dip = ',site_info%nuclear_charge(n)-Lchg(n), &
      !!                       (-Ldip(3*(n-1)+j),j=1,3)
      !!write(14,*)'sec mom (TR) = ',(-Lsec_mom_TR(6*(n-1)+j),j=1,3)
      !!write(14,*)'sec mom (TR) = ',(-Lsec_mom_TR(6*(n-1)+j),j=4,6)
      ! GAC for mpole analysis
      !if(n.eq.6)write(15,*)'-------------------------------------'
      !write(15,*)Lchg(n)!-site_info%nuclear_charge(n)
      write(15,*)Lchg(n)
      write(15,*)Ldip(3*(n-1)+1)
      write(15,*)Ldip(3*(n-1)+2)
      write(15,*)Ldip(3*(n-1)+3)
      write(15,*)Lsec_mom(6*(n-1)+1)
      write(15,*)Lsec_mom(6*(n-1)+2)
      write(15,*)Lsec_mom(6*(n-1)+3)
      write(15,*)Lsec_mom(6*(n-1)+4)
      write(15,*)Lsec_mom(6*(n-1)+5)
      write(15,*)Lsec_mom(6*(n-1)+6)
! FOR POLEDIT AMOEBA
      write(16,114)int(site_info%nuclear_charge(n))
      write(16,113)(site_info%site_crds(3*(n-1)+j)*bohr,j=1,3)
      write(16,111)(site_info%nuclear_charge(n)-Lchg(n))
      write(16,111)-Ldip(3*(n-1)+1)!*bohr
      write(16,111)-Ldip(3*(n-1)+2)!*bohr
      write(16,111)-Ldip(3*(n-1)+3)!*bohr
      !! PARA MULTIPOLOS EN TINKER, NECESITA 3/2 para QUADS y 2 para diag
      write(16,111)-(Lsec_mom(6*(n-1)+1)*3.0d0)
      write(16,111)-(Lsec_mom(6*(n-1)+2)*3.0d0)
      write(16,111)-(Lsec_mom(6*(n-1)+3)*3.0d0)
      write(16,111)-(Lsec_mom(6*(n-1)+4)*1.5d0)
      write(16,111)-(Lsec_mom(6*(n-1)+5)*1.5d0)
      write(16,111)-(Lsec_mom(6*(n-1)+6)*1.5d0)
      !! PARA MULTIPOLOS EN GEM, NO USAR 2/3 de STONE!
      !write(16,111)-(Lsec_mom(6*(n-1)+1)*bohr*bohr)
      !write(16,111)-(Lsec_mom(6*(n-1)+2)*bohr*bohr)
      !write(16,111)-(Lsec_mom(6*(n-1)+3)*bohr*bohr)
      !write(16,111)-(Lsec_mom(6*(n-1)+4)*bohr*bohr)
      !write(16,111)-(Lsec_mom(6*(n-1)+5)*bohr*bohr)
      !write(16,111)-(Lsec_mom(6*(n-1)+6)*bohr*bohr)
!! FOR AMBER AMOEBA implementation direct
!      write(17,*)(site_info%nuclear_charge(n)-Lchg(n))
!      write(17,*)-Ldip(3*(n-1)+1)
!      write(17,*)-Ldip(3*(n-1)+2)
!      write(17,*)-Ldip(3*(n-1)+3)
!      !! PARA MULTIPOLOS EN AMBER, NECESITA 3/2 para TODOS LOS QUADS
!      write(17,*)-(Lsec_mom(6*(n-1)+1)*1.5d0)
!      write(17,*)-(Lsec_mom(6*(n-1)+2)*1.5d0)
!      write(17,*)-(Lsec_mom(6*(n-1)+3)*1.5d0)
!      write(17,*)-(Lsec_mom(6*(n-1)+4)*1.5d0)
!      write(17,*)-(Lsec_mom(6*(n-1)+5)*1.5d0)
!      write(17,*)-(Lsec_mom(6*(n-1)+6)*1.5d0)
    enddo
  enddo
111 format(F18.10)
113 format (3(1X,F14.8))
114 format(i3)
  !close(unit=14)
  close(unit=15)
  close(unit=16)
  !close(unit=17)
  deallocate(Lchg,Ldip,Lsec_mom)
  return
end subroutine CH_COEFF_mpoles_NEW
!------------------------------------------------------
subroutine CH_COEFF_mpoles_GLOBAL(atom_auxis,site_info,filename2)
use definition
  implicit none
  !include "definition.fh"
  integer p_site_info
  type(auxiliary_basis)::atom_auxis
  type(sites_info)::site_info
  character(len=*) filename2
  integer ptr,nsites,ps_numprim,ps_offprim,ps_primorder,  &
          ps_prim_off_coeff,ncoefftot,p_L_hermite_coeff,ps_atomic_num, &
          p_M_hermite_coeff,ps_primexpo
  double precision lpole(10),gpole(3),finaldip,fac,pi,sqpi,origin(3)
  integer num,off,order,scoff,i,j,k,m,n,allochk
  integer p_chg,p_dip,p_sec_mom
  integer p_Lchg,p_Ldip,p_Lsec_mom
  integer imol,num_mols,nlo,nhi,p_molptr,p_molsitecrd
  double precision,allocatable::Lchg(:),Ldip(:),Lsec_mom(:),trace(:),Gdip(:)
  double precision,allocatable::nucdip(:),elecdip(:),chgcontr(:),&
                                atmdip(:),atmchg(:),totdip(:)
  double precision alpha

  open(unit=16,file=filename2,status='unknown')

  nsites = site_info%nsites

  allocate(Lchg(nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'CH_COEFF_mpoles_NEW:could not allocate Lchg, exiting'
     stop
  endif
  allocate(Ldip(3*nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'CH_COEFF_mpoles_NEW:could not allocate Ldip, exiting'
     stop
  endif
  allocate(Gdip(3*nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'CH_COEFF_mpoles_NEW:could not allocate Gdip, exiting'
     stop
  endif
  allocate(Lsec_mom(6*nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'CH_COEFF_mpoles_NEW: could not allocate Lsec_mom, exiting'
     stop
  endif
  allocate(trace(nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'CH_COEFF_mpoles_NEW: could not allocate trace, exiting'
     stop
  endif
  allocate(nucdip(3), elecdip(3), chgcontr(3), totdip(3),&
 &         atmchg(nsites), atmdip(3*nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'CH_COEFF_mpoles_NEW:could not allocate nucdip or elecdip&
    & or atmchg or totdip or chgcontr or atmdip, exiting'
     stop
  endif

  trace(:) = 0.0d0
  nucdip(:) = 0.0d0
  elecdip(:) = 0.0d0
  chgcontr(:) = 0.0d0
  atmdip(:) = 0.0d0
  totdip(:) = 0.0d0
  origin(:) = 0.0d0
  atmchg(:) = 0.0d0
  pi = 4.d0*(atan(1.0d0))
  sqpi = dsqrt(pi)
  fac = 0.529177249d0*4.8032068d-10*1.d10

  num_mols = 1

  do n = 1,nsites
    num = site_info%num_primitives(n)
    off = site_info%off_primitives(n)

    gpole(:) = 0.d0
    do k = 1,10
      lpole(k) = 0.d0
    enddo
    do k = 1,num
      order = site_info%prim_order(off+k)

      scoff = site_info%coeff_offset(off+k)

      alpha = site_info%prim_expo(off+k)

      lpole(1) = lpole(1) + site_info%global_hermite_coeffs(scoff+1)

      atmchg(n) = atmchg(n) + site_info%aux_elec(scoff+1)*&
                              site_info%cartesian_coeffs(scoff+1)
      
      if ( order == 4 )then
        do j = 2,4
          lpole(j) = lpole(j) + site_info%global_hermite_coeffs(scoff+j)
          gpole(j-1) = gpole(j-1) + site_info%cartesian_coeffs(scoff+j)*&
                                   (site_info%cart_coeff_norms(scoff+j)*&
                                    pi*sqpi)/(2.d0*alpha**2.5)
          atmchg(n) = atmchg(n) + site_info%aux_elec(scoff+j)*&
                                  site_info%cartesian_coeffs(scoff+j)
        enddo
      else if ( order == 6 )then
        do j = 2,4
          lpole(j) = lpole(j) + 2.d0*site_info%global_hermite_coeffs(scoff+j)
          atmchg(n) = atmchg(n) + site_info%aux_elec(scoff+j)*&
                                  site_info%cartesian_coeffs(scoff+j)
        enddo
        do j = 5,7
          lpole(j) = lpole(j) + site_info%global_hermite_coeffs(scoff+j)
          atmchg(n) = atmchg(n) + site_info%aux_elec(scoff+j)*&
                                  site_info%cartesian_coeffs(scoff+j)
        enddo
      else if ( order == 10 )then
        do j = 2,4
          lpole(j) = lpole(j) + site_info%global_hermite_coeffs(scoff+j)
          gpole(j-1) = gpole(j-1) + site_info%cartesian_coeffs(scoff+j)*&
                                   (site_info%cart_coeff_norms(scoff+j)*&
                                    pi*sqpi)/(2.d0*alpha**2.5)
          atmchg(n) = atmchg(n) + site_info%aux_elec(scoff+j)*&
                                  site_info%cartesian_coeffs(scoff+j)
        enddo
        do j = 5,7
          lpole(j) = lpole(j) + 2.d0*site_info%global_hermite_coeffs(scoff+j)
          atmchg(n) = atmchg(n) + site_info%aux_elec(scoff+j)*&
                                  site_info%cartesian_coeffs(scoff+j)
        enddo
        do j = 8,10
          lpole(j) = lpole(j) + site_info%global_hermite_coeffs(scoff+j)
          atmchg(n) = atmchg(n) + site_info%aux_elec(scoff+j)*&
                                  site_info%cartesian_coeffs(scoff+j)
        enddo
      endif
      ! next add charge contribution to xx,yy,zz moments
      do j = 5,7
        lpole(j) = lpole(j) + site_info%global_hermite_coeffs(scoff+1)&
                              / (2.d0*alpha)
      enddo
    enddo
    !RH(p_Lchg+n-1) = lpole(1)
    Lchg(n) = lpole(1)
    do j = 1,3
      !RH(p_Ldip+3*(n-1)+j-1) = lpole(j+1)
      Ldip(3*(n-1)+j) = lpole(j+1)
      Gdip(3*(n-1)+j) = gpole(j)
    enddo
    do j = 1,6
      !RH(p_Lsec_mom+6*(n-1)+j-1) = lpole(j+4)
      Lsec_mom(6*(n-1)+j) = lpole(j+4)
      if (j .le. 3) then
         trace(n) = trace(n) + Lsec_mom(6*(n-1)+j)
      endif
    enddo
    trace(n) = (1.d0/3.d0)*trace(n)
  enddo

  !print *,'Adding trace to quadrupoles'
  do n = 1, nsites
    do j = 1,3
      Lsec_mom(6*(n-1)+j) = 0.5*(Lsec_mom(6*(n-1)+j)-trace(n))
    enddo
    !print *,'atmchg = ',atmchg(n)
  enddo

  ! calculating atomic dipoles for polarizabilities and total dipole
100 format(20x,'Dipole moment in Debye (w.r.t. origin)')
110 format(21X,'x',15X,'y',15X,'z'&
     &/1X,'Electronic',3(2X,F14.8)&
     &/1X,'Nuclear',3X,3(2X,F14.8)&
     &/1X,'Elec Chg',2X,3(2X,F14.8)&
     &/1X,'Total',5X,3(2X,F14.8)&
     &/1X,'Dipole Moment',F15.8)
102 format(20x,'Atomic Dipoles in Debye'&
          /21X,'x',15X,'y',15X,'z')
112 format(1X,'site ',i3,2x,3(2X,F14.8))

  do j = 1,3
     do n = 1,nsites
        atmdip(3*(n-1)+j) = site_info%nuclear_charge(n)*&   !nuc chg contr.
                            (site_info%site_crds(3*(n-1)+j)-origin(j))&
                            - Gdip(3*(n-1)+j) -&            !elec dip
                            atmchg(n)*&                     !chg contr
                            (site_info%site_crds(3*(n-1)+j)-origin(j))
        chgcontr(j) = chgcontr(j) + atmchg(n)*&
                                    (site_info%site_crds(3*(n-1)+j)&
                                     -origin(j))
        nucdip(j) = nucdip(j) + site_info%nuclear_charge(n)*&
                                (site_info%site_crds(3*(n-1)+j)-origin(j))
        elecdip(j) = elecdip(j) + Gdip(3*(n-1)+j)
     enddo
     totdip(j) = totdip(j) + nucdip(j) - elecdip(j) - chgcontr(j)
  enddo
  elecdip(:) = elecdip(:)*fac
  nucdip(:) = nucdip(:)*fac
  chgcontr(:) = chgcontr(:)*fac
  totdip(:) = totdip(:)*fac
  atmdip(:) = atmdip(:)*fac
     
  finaldip = dsqrt(totdip(1)*totdip(1)+totdip(2)*totdip(2)+&
                   totdip(3)*totdip(3))

  !write(*,102)
  !do n = 1, nsites
  !   write(*,112) n,(atmdip(3*(n-1)+j),j=1,3)
  !enddo

  !print *,''
  !write(*,100)
  !write(*,110) elecdip,nucdip,chgcontr,totdip,finaldip


  do imol = 1,num_mols
    nlo = 1
    nhi = nsites
    write(16,*)nsites*10
    do n = nlo,nhi
      ! GAC for mpole analysis
      !if(n.eq.6)write(16,*)'-------------------------------------'
      !write(16,*)Lchg(n)!-site_info%nuclear_charge(n)
      write(16,*)Lchg(n)
      write(16,*)Ldip(3*(n-1)+1)
      write(16,*)Ldip(3*(n-1)+2)
      write(16,*)Ldip(3*(n-1)+3)
      write(16,*)Lsec_mom(6*(n-1)+1)
      write(16,*)Lsec_mom(6*(n-1)+2)
      write(16,*)Lsec_mom(6*(n-1)+3)
      write(16,*)Lsec_mom(6*(n-1)+4)
      write(16,*)Lsec_mom(6*(n-1)+5)
      write(16,*)Lsec_mom(6*(n-1)+6)
    enddo
  enddo
  close(unit=16)
  deallocate(Lchg,Ldip,Lsec_mom)
  return
end subroutine CH_COEFF_mpoles_GLOBAL
