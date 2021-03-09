!-----------------------------------------------------------------------------
subroutine GridGEN(site_info,radnum,angnum)

use definition

  IMPLICIT NONE

  type(sites_info)::site_info
  INTEGER   :: i,j,maxpts,radnum,angnum,nattypes,allochk
  integer,allocatable,dimension(:)::intchg
  double precision :: bohr
  double precision BSrad(105)
  character*2,allocatable,dimension(:):: atmnames
  include "site_names.inc"
!
!
!     Table of Bragg-Slater Atomic Radii (Angstroms)
!
!     Bragg-Slater radii: J.C. Slater, Symmetry and Energy Bands in Crystals,
!                         Dover, N.Y. 1972, page 55.
!                         The radii of noble gas atoms are set to be equal 
!                         to the radii of the corresponding halogen atoms.
!                         The radius of At is set to be equal to the radius of
!                         Po. 
!                         The radius of Fr is set to be equal to the radius of
!                         Cs. 
!
!                  H    He   Li   Be    B    C    N    O    F   Ne   
      Data BSrad/0.35,0.35,1.45,1.05,0.85,0.70,0.65,0.60,0.50,0.50,&
!                  Na   Mg   Al   Si    P    S   Cl   Ar    K   Ca
                 1.80,1.50,1.25,1.10,1.00,1.00,1.00,1.00,2.20,1.80,&
!                  Sc   Ti    V   Cr   Mn   Fe   Co   Ni   Cu   Zn
                 1.60,1.40,1.35,1.40,1.40,1.40,1.35,1.35,1.35,1.35,&
!                  Ga   Ge   As   Se   Br   Kr   Rb   Sr    Y   Zr
                 1.30,1.25,1.15,1.15,1.15,1.15,2.35,2.00,1.80,1.55,&
!                  Nb   Mo   Tc   Ru   Rh   Pd   Ag   Cd   In   Sn
                 1.45,1.45,1.35,1.30,1.35,1.40,1.60,1.55,1.55,1.45,&
!                  Sb   Te    I   Xe   Cs   Ba   La   Ce   Pr   Nd
                 1.45,1.40,1.40,1.40,2.60,2.15,1.95,1.85,1.85,1.85,&
!                  Pm   Sm   Eu   Gd   Tb   Dy   Ho   Er   Tm   Yb
                 1.85,1.85,1.85,1.80,1.75,1.75,1.75,1.75,1.75,1.75,&
!                  Lu   Hf   Ta    W   Re   Os   Ir   Pt   Au   Hg
                 1.75,1.55,1.45,1.35,1.35,1.30,1.35,1.35,1.35,1.50,&
!                  Tl   Pb   Bi   Po   At   Rn   Fr   Ra   Ac   Th
                 1.90,1.80,1.60,1.90,1.90,1.90,2.60,2.15,1.95,1.80,&
!                  Pa    U   Np   Pu   Am   Cm   Bk   Cf   Es   Fm
                 1.80,1.75,1.75,1.75,1.75,1.75,1.75,1.75,1.75,1.75,&
!                  Md   No   Lr  Unq  Unp
                 1.75,1.75,1.75,1.55,1.55/

  bohr=0.529177249d0

  write(6,*)'Generating grid with',angnum,' angular and ',radnum,' radial pts'

!
! write input file for grid generation
!
  ! have to write names of atoms, need to get from sitenames
  allocate(atmnames(site_info%natoms), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'h_Grid_Gen:could not allocate atmnames, exiting'
     stop
  endif

  call fix_atm_names(site_info,atmnames,site_info%natoms,&
                     site_info%nsites,nattypes,intchg,1)

  open(15,file='grid.input',status='unknown')

  !Write molecular coords
  write(15,*)'@Geometry'
  write(15,105)site_info%natoms
105 format (5X,'n_atoms ',i4)
  do i = 1, site_info%natoms
    write(15,'(10X,a2,2x,3f15.9)')atmnames(i)(1:2),&
         site_info%site_crds(3*(i-1)+1)*bohr,&
         site_info%site_crds(3*(i-1)+2)*bohr,&
         site_info%site_crds(3*(i-1)+3)*bohr
  enddo

  allocate(intchg(nattypes), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'h_Grid_Gen:could not allocate atmnames, exiting'
     stop
  endif

  call fix_atm_names(site_info,atmnames,site_info%natoms,&
                     site_info%nsites,nattypes,intchg,0)

  ! Write Atom Types
  write(15,*)'@Element'
  write(15,103)nattypes
  do i = 1, nattypes
    write(15,'(10X,a2)')atmnames(i)
  enddo

  ! Write Nuclear Charges
  write(15,*)'@Z'
  write(15,103)nattypes
  do i = 1, nattypes
    !intchg(i) = int(site_info%nuclear_charge(i))
    write(15,104)atmnames(i),intchg(i)
  enddo
104 format (10X,a2,10X,i3,'.0')

  close (15)
!
! write data file for grids
!
  open(15,file='Grid_Parm.dat',status='unknown')

  ! Write Grid parameters
  write(15,*)'@Instruct'
  write(15,101)angnum
101 format (5X,'Lebedev_Ang_Type ',i4,' 1')
  write(15,102)radnum
102 format (5X,'N_Radial_Pt ',i4)

  ! Write Atom Types
  write(15,*)'@Element'
  write(15,103)nattypes
103 format (5X,'n_element ',i4)
  do i = 1, nattypes
    write(15,'(10X,a2)')atmnames(i)
  enddo

  ! Write Slater-Bragg Radii
  write(15,*)'@Bragg_Radii'
  write(15,103)nattypes
  do i = 1, nattypes
    write(15,'(10X,a2,2xf6.4)')atmnames(i),BSRad(intchg(i))
  enddo

  close (15)

  
END SUBROUTINE GridGEN

!-----------------------------------------------------------------------------
subroutine fix_atm_names(site_info,atmnames,natoms,nsites,nattypes,&
                         intchg,full)
use definition

  implicit none
  include "site_names.inc"

  type(sites_info)::site_info
  integer  i,j,natoms,nsites,nattypes,there(nsites),intchg(nsites),&
           full,nums,ntyptmp,esta,jj
  double precision chgtmp,chg(nsites)
  character*2 atmnames(natoms)

  ! figure out how many different atoms we have
  nattypes=maxval(site_info%site_type,nsites)

  there(1) = 0
  there(2:nsites) = 1
  chg(1) = site_info%nuclear_charge(1)   
  if (site_info%nuclear_charge(1).eq.0.0d0) then
      write (6,*)'Charge on first site is 0, there is something wrong'
      stop
  endif
  
  chg(:) = site_info%nuclear_charge(:)   

  i = 1
  j = 2
  ntyptmp = 1
  do while (ntyptmp .lt. nattypes)
     if (i .eq. 1) then
        do while (j .lt. site_info%nsites)
           if (chg(i) .eq. site_info%nuclear_charge(j) )then
              there(j) = 1 
              j = j + 1
           else
              there(j) = 0
              ntyptmp = ntyptmp + 1
              i = j
              j = site_info%nsites
           endif
        enddo
     else
        j = i + 1
        do while (j .lt. site_info%nsites)
           esta = 0
           do jj = 1, i
              if (chg(jj) .eq. site_info%nuclear_charge(j) )then
                 esta = 1
                 goto 101
              endif
           enddo
101        continue
           if (esta == 1) then
              there(j) = 1 
              j = j + 1
           else
              there(j) = 0
              ntyptmp = ntyptmp + 1
              i = j
              j = site_info%nsites
           endif
        enddo
     endif
  enddo

  if ((ntyptmp .eq. nattypes) .and. i .lt. site_info%nsites) then
     do j = i+1,nsites
        there(j) = 1
     enddo
  endif

  ntyptmp = 1
  do i = 2, site_info%nsites
     if (there(i) .eq. 0 )then
        ntyptmp = ntyptmp + 1
        chg(ntyptmp) = site_info%nuclear_charge(i)
     endif
  enddo

  ! get atom names for each different atom
  if (full == 1) then 
     nums = site_info%natoms
     chg(:) = site_info%nuclear_charge(:)
     do i = 1, nums
        if (chg(i)==1.d0 .or. chg(i)==5.d0 .or. chg(i)==6.d0 .or.&
            chg(i)==7.d0 .or. chg(i)==8.d0 .or. chg(i)==9.d0 .or.&
            chg(i)==15.d0 .or. chg(i)==16.d0 .or. chg(i)==19.d0.or.&
            chg(i)==23.d0 .or. chg(i)==39.d0 .or. chg(i)==53.d0.or.&
            chg(i)==74.d0) then
            atmnames(i)(1:1)=sitename(i)(1:1)
            atmnames(i)(2:2)=' '
        else
            atmnames(i)(1:2)=sitename(i)(1:2)
        endif
     enddo
  else
     nums = 0
     do i = 1, site_info%natoms
        if (there(i) .eq. 0) then
           if (chg(i)==1.d0 .or. chg(i)==5.d0 .or. chg(i)==6.d0 .or.&
               chg(i)==7.d0 .or. chg(i)==8.d0 .or. chg(i)==9.d0 .or.&
               chg(i)==15.d0 .or. chg(i)==16.d0 .or. chg(i)==19.d0.or.&
               chg(i)==23.d0 .or. chg(i)==39.d0 .or. chg(i)==53.d0.or.&
               chg(i)==74.d0) then
               nums = nums + 1
               intchg(nums) = int(site_info%nuclear_charge(i))
               atmnames(nums)(1:1)=sitename(i)(1:1)
               atmnames(nums)(2:2)=' '
           else
               nums = nums + 1
               intchg(nums) = int(site_info%nuclear_charge(i))
               atmnames(nums)(1:2)=sitename(i)(1:2)
           endif
        endif
     enddo
  endif

         

end subroutine fix_atm_names
!-----------------------------------------------------------------------------
