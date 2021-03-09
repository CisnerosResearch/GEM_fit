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
program GEM_numfit

use definition
use cubes

  implicit none

  include "io.fh"
  type(AO_basis)::atom_basis
  type(auxiliary_basis)::atom_auxis
  type(sites_info)::site_info
  type(aux_orbitals),allocatable::auxis(:)
  type(ao_orbitals),allocatable::basis(:)
  type(exact_cubes)::ex_cubes
  type(auxiliary_cubes)::aux_cubes
  integer argc,arg,l,allochk,i,iargc,ncoeff,num_stos,inttype,cubetype
  integer npts,ntheta,nang,nphi,nrad,count,ntheta2,npts2,natoms,rdstat
  integer maxpts,radnum,angnum,grdsize
  double precision rmin, rmax, delta_r, delta_ang, circum, delta_phi, theta
  double precision rad, angle, rmin2, rmax2, delta_r2, pi, cutoff
  double precision beta, scale, tot_esp, tot_fld(3), exact_esp, exact_fld(3)
  double precision tx, ty, tz, wt
  double precision,allocatable,dimension(:)::Jvec,aux_coefs,JvecS,wght
  double precision,allocatable,dimension(:,:)::Gmat,GmatS,rad_grid
  character(len=PATH_MAX) argv
  character(len=3)dtype
  character(len=80)parmfile,auxfile,crdfile,densfile,basfile,coeff_file,&
                   lhermitefile,mpole_file,mpole_file2,controlfile,filehead,&
                   promolfile,mpole_file3,wfnfile,cmd
  logical constrain, readdens, readbas, debug, cartfit, cube, printcube,&
          twocentfit, promolfit, sphere, file_there, read_grid, rd_grd2,&
          HLYfit,gemesp

  controlfile = ''
  arg=1
  argc=IArgC()
  do while (arg <= argc)
    call GetArg(arg,argv)
    if(argv == '-control') then
      arg=arg+1
      call GetArg(arg,controlfile)
    endif
    arg = arg + 1
  enddo
  if (controlfile == '')then
       write(6,*)'usage: GEM_fit -control controlfile'
       stop
  endif
  cube = .true.

  maxpts = 1000000 ! hardcoded, need a better way to handle this
 
  call HCRD_read_control_NEW(controlfile,parmfile,auxfile,basfile,crdfile,&
                         densfile,promolfile,coeff_file,lhermitefile,&
                         mpole_file,mpole_file2,mpole_file3,&
                         wfnfile,radnum,angnum,grdsize,&
                         constrain,readdens,readbas,debug,cartfit,inttype,&
                         beta,scale,cutoff,printcube,gemesp,&
                         promolfit,twocentfit,cubetype,sphere,HLYfit)
  cartfit = .true.
  print *,densfile
  call AHCRD_read_file(site_info,crdfile)
  call AHBASE_load_auxnames(atom_auxis,parmfile)
  call AHBASE_load_auxbasis(atom_auxis,auxfile)
  call AHBASE_calcnorms(atom_auxis,densfile)
  if (readdens) call HLOAD_dens(site_info,densfile,dtype,.false.)
  if (readbas) then
     call AHBASE_load_basnames(atom_basis,parmfile)
     call AHBASE_load_AObasis(atom_basis,basfile)
     call AHBASE_calcAOnorms(atom_basis,densfile)
  endif
  call AHTYPE_load_site_info(site_info,atom_auxis,atom_basis,parmfile,readdens)
  call AHCOEFF_copy_norms_from_basis(site_info,atom_auxis)
  call AHFRAME_get_molframe(site_info,parmfile)
  call AHFRAME_load_deflist(site_info,parmfile,readdens,debug)
  call AHFRAME_build_frames(site_info)

!! ALLOCATE G matrix, and Jvec
  ncoeff = site_info%num_coefficients
  num_stos = site_info%num_stos
  allocate(Gmat(ncoeff,ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'GEM_fit:could not allocate Gmat, exiting'
     if (debug)close(12)
     stop
  endif
  allocate(Jvec(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'GEM_fit:could not allocate Jvec, exiting'
     if (debug)close(12)
     stop
  endif
  if (sphere) then
     allocate(GmatS(ncoeff,ncoeff), stat = allochk)
     if (allochk .gt. 0 ) then
        write(6,*)'GEM_fit:could not allocate GmatS, exiting'
        if (debug)close(12)
        stop
     endif
     GmatS(:,:)=0.d0
     allocate(JvecS(ncoeff), stat = allochk)
     if (allochk .gt. 0 ) then
        write(6,*)'GEM_fit:could not allocate JvecS, exiting'
        if (debug)close(12)
        stop
     endif
     JvecS(:)=0.d0
  endif ! sphere
  allocate(aux_coefs(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'GEM_fit:could not allocate aux_coefs, exiting'
     if (debug)close(12)
     stop
  endif
  allocate(auxis(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'GEM_fit:could not allocate auxis, exiting'
     if (debug)close(12)
     stop
  endif
  allocate(basis(num_stos), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'GEM_fit:could not allocate basis, exiting'
     if (debug)close(12)
     stop
  endif

! Calculate exact cubes (or spheres) of ESP 

  if (.not. sphere) then
     write(6,*)' GEM_fit: Performing fit from rectangular cube'
     call H_cube_gen_exact(site_info,basis,ex_cubes,num_stos,cube,printcube,&
                           twocentfit,densfile)
  else
      write(6,*)' GEM_fit: Performing fit from spherical grid'

      ! Determine size of grid

      if ((radnum .eq. 0) .and. (angnum .eq. 0)) then
         if (grdsize .eq. 1) then
            !radnum = 21
            radnum = 21
            angnum = 23
            write (6,*),' Using extra coarse grids'
         else if (grdsize .eq. 2) then
            radnum = 35
            angnum = 29
            write (6,*),' Using coarse grids'
         else if (grdsize .eq. 3) then
            radnum = 49
            angnum = 41
            write (6,*),' Using medium grids'
         else if (grdsize .eq. 4) then
            radnum = 70
            angnum = 47
            write (6,*),' Using fine grids'
         else if (grdsize .eq. 5) then
            radnum = 100
            angnum = 53
            write (6,*),' Using extra fine grids'
         endif
      endif
      if ((radnum .lt. 0).or.(radnum .gt. 200)) then
         write(6,*)'Number of radial grid points out of range, setting to 21'
         radnum = 21
      endif
      if ((angnum .ne. 110).and.(angnum .ne. 194).and.(angnum .ne. 302)&
           .and.(angnum .ne. 590).and.(angnum .ne. 770).and.(angnum &
           .ne. 974).and.(angnum.ne.0).and.(grdsize.lt.1)) then
         angnum = 194
         write(6,*)'Number of angular grid points out of range, setting&
                    & to 194 (extra coarse)',radnum,angnum
      endif

      ! Generate grid points if grids not there 
      rd_grd2= .false.
      inquire(file = 'grd_pts.xyz', exist = rd_grd2)
      if (rd_grd2) then

         write (6,*)'Found Grid File, reading in points'
         open(13,file='grd_pts.xyz',status='unknown')
         rdstat = 0
         npts = -1
         do while(rdstat==0)
            read(13,*,iostat=rdstat) tx, ty, tz, wt
            npts = npts + 1
         enddo
         rewind(13)
         !read(13,*)npts
         write(6,*)'GEM_fit: Using Becke grids to fit ',npts,' points'

         allocate(rad_grid(npts,3),wght(npts), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'GEM_G_mat:could not allocate rad_grid, exiting'
            stop
         endif
         rad_grid(:,:)=0.d0
         wght(:) = 1.d0
         do i = 1, npts
            read(13,*)rad_grid(i,1),rad_grid(i,2),rad_grid(i,3),wght(i)
            rad_grid(i,1) = rad_grid(i,1)/0.529177249d0
            rad_grid(i,2) = rad_grid(i,2)/0.529177249d0
            rad_grid(i,3) = rad_grid(i,3)/0.529177249d0
         enddo
         close(13)

      else

         ! Generate grids (code from Paul Ayers)
         !call wfnBECKE(rad_grid,wght,npts,maxpts,radnum,angnum,wfnfile)

         !generate inputs and call Grid program from D. Elking
         call GridGEN(site_info,radnum,angnum)
         write(cmd,*)'$GEMNF_HOME/grid_gen.x grid.input'
         call system(cmd)
         !clean up files for grid generation
         write(cmd,*)'/bin/rm grid.input Grid_Parm.dat grid.debug'
         call system(cmd)
 
         ! read grids from file
         open(13,file='grd_pts.xyz',status='unknown')
         rdstat = 0
         npts = -1
         do while(rdstat==0)
            read(13,*,iostat=rdstat) tx, ty, tz, wt
            npts = npts + 1
         enddo
         rewind(13)
         !read (13,*)npts
         if (debug)write (12,*)npts
         write (6,*)'GEM_fit: Grid has ',npts,' points'
  
         if (npts .gt. maxpts) then
            write (6,*)'Maximum Grid points exceeded, exiting'
            stop
         endif

         ! ALLOCATE rad_grid with correct npts
         allocate(rad_grid(npts,3),wght(npts), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'GEM_fit:could not allocate rad_grid, wght, exiting'
            stop
         endif
         rad_grid(:,:)=0.d0
         wght(:)=0.d0

         do i = 1,npts
            read(13,*)rad_grid(i,1),rad_grid(i,2),rad_grid(i,3),wght(i)
            rad_grid(i,1) = rad_grid(i,1)/0.529177249d0
            rad_grid(i,2) = rad_grid(i,2)/0.529177249d0
            rad_grid(i,3) = rad_grid(i,3)/0.529177249d0
            if (debug)write (12,*)rad_grid(i,1),rad_grid(i,2),&
                                  rad_grid(i,3), wght(i)
         enddo
         close(13)
         !print *,'stopping after reading grid', npts
         !stop

      endif

      ! Calculate (or read) ESP at each grid point
      call H_spherical_grid2(site_info,basis,ex_cubes,rad_grid,wght,npts,&
                             num_stos,densfile,dtype,printcube,debug,&
                             cube,gemesp)

  endif

! Calculate Gmat and Jvec with cubes

  if (.not. sphere) then
     call H_cube_aux(site_info,auxis,aux_cubes,ex_cubes,Gmat,Jvec,cutoff,&
                     rad_grid,wght,1,ncoeff,debug,cubetype,sphere,.false.,&
                     .false.,densfile,HLYfit)
  else
     call H_cube_aux(site_info,auxis,aux_cubes,ex_cubes,GmatS,JvecS,cutoff,&
                     rad_grid,wght,npts,ncoeff,debug,cubetype,sphere,&
                     promolfit,.false.,densfile,HLYfit)
     if (promolfit) then
        write (6,*) 'GEM_fit: promol keyword found, calculating analytic &
                     &Jvec and Gmat'
        call AH_form_G_mat_CART(site_info,auxis,Gmat,ncoeff,inttype,beta,&
                                debug,cube)
        call H_form_J_vec_C(site_info,auxis,basis,Jvec,ncoeff,num_stos,&
                            inttype,beta,debug,cube,.false.,twocentfit,densfile)
        Gmat(:,:) = Gmat(:,:) - GmatS(:,:)
        Jvec(:) = Jvec(:) - JvecS(:)
     else
        Gmat(:,:) = GmatS(:,:)
        Jvec(:) = JvecS(:)
     endif
  endif

! Invert G matrix and form coefficients
  call AH_form_aux_coefs(site_info,Gmat,Jvec,aux_coefs,ncoeff,scale,&
                         constrain,debug)
  deallocate(aux_coefs,Jvec,Gmat)

!! NOW TRANSFORM CARTESIAN COEFS TO LOCAL HERMITE
  call AHCOEFF_load_cartesian(site_info,coeff_file,readdens,densfile)
  call AHFRAME_rotate_hermites(site_info,lhermitefile)
  call CH_COEFF_mpoles_NEW(atom_auxis,site_info,mpole_file,mpole_file2)
  !call CH_COEFF_mpoles_GLOBAL(atom_auxis,site_info,mpole_file3)
  call H_calc_esp_CART(site_info,auxis,0.0d0,0.0d0,0.0d0,tot_esp,&
                       ncoeff,.false.,.false.)
  !call H_calc_fld_CART(site_info,auxis,0.0d0,0.0d0,0.0d0,tot_fld,&
  !                     ncoeff,cartfit,cube)

  !! Calculate auxiliary cubes and take difference
  !call H_cube_aux_coef(site_info,auxis,aux_cubes,ex_cubes,ncoeff,.true.)

  deallocate(auxis)

  !if (debug)close(12)

end program GEM_numfit
