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
subroutine H_cube_gen_exact(site_info,basis,ex_cubes,num_stos,cube,printcube,&
                            twocentfit,densfile)
use definition
use cubes
  implicit none
  type(sites_info)::site_info
  type(ao_orbitals)::basis(num_stos)
  type(exact_cubes)::ex_cubes

  integer i,j,k,l
  integer n_cub_atoms, num_a, num_b, num_c, num_stos, allochk, istat
  integer tmp_cub_atm, tmp_num_a, natoms
  integer,allocatable::cub_atom_num(:),limitvec(:)
  integer,parameter::size = 250
  double precision esp,dens,field(3),x_orig,y_orig,z_orig,x1,y1,z1,&
                   x2,y2,z2,x3,y3,z3,xin,yin,zin
  double precision tmp_x, tmp_y, tmp_z, tmp_x1, tmp_y1, tmp_z1
  double precision densmat(num_stos,num_stos) 
  double precision,allocatable::cub_chg(:),cub_coord(:,:)
  double precision,parameter::AtoBohr = 0.529177249d0,small=1.d-5
  character(len=*) densfile
  logical cube, file_there, printcube, cub1, cub2, cub3, cub4, cub5, &
          basalloc, twocentfit
  character*1 int_type

  cub1 = .false.
  cub2 = .false.
  cub3 = .false.
  cub4 = .false.
  cub5 = .false.

! Create basis info
  call H_form_basis(site_info,basis,num_stos,densfile)
  natoms = site_info%natoms
  allocate(limitvec(natoms), stat=allochk)
  if (allochk .gt. 0) then
     write(6,*)'H_cube_gen_exact: could not allocate limitvec, exiting'
     stop
  endif
  call H_basis_limits(site_info,basis,num_stos,natoms,limitvec)

! Read cube info
  inquire(file = 'cube_file', exist = file_there)
  if (.not. file_there) then
     write(6,*)'H_cube_gen_exact: Could not find cube info file &
                &(cube_file), exiting'
     stop
  endif
  open(13,file='cube_file',status='unknown')
  read (13,'(i5,3f12.6)',iostat=istat) n_cub_atoms,x_orig, y_orig, z_orig
  read (13,'(i5,3f12.6)',iostat=istat) num_a, x1, y1, z1
  read (13,'(i5,3f12.6)',iostat=istat) num_b, x2, y2, z2
  read (13,'(i5,3f12.6)',iostat=istat) num_c, x3, y3, z3
  ex_cubes%n_cub_atoms = n_cub_atoms
  ex_cubes%num_a = num_a
  ex_cubes%num_b = num_b
  ex_cubes%num_c = num_c
  ex_cubes%x_orig = x_orig
  ex_cubes%y_orig = y_orig
  ex_cubes%z_orig = z_orig
  ex_cubes%x1 = x1
  ex_cubes%y1 = y1
  ex_cubes%z1 = z1
  ex_cubes%x2 = x2
  ex_cubes%y2 = y2
  ex_cubes%z2 = z2
  ex_cubes%x3 = x3
  ex_cubes%y3 = y3
  ex_cubes%z3 = z3
  !!! Check if number of atoms in info file matches
  if (n_cub_atoms .ne. site_info%natoms) then
     write(6,*)'H_cube_gen_exact: Number of atoms in cube_file &
                &does not match, exiting'
     stop
  endif
  !!! allocate local vars
  allocate(cub_atom_num(n_cub_atoms),cub_chg(n_cub_atoms), stat=allochk)
  if (allochk .gt. 0) then
     write(6,*)'H_cube_gen_exact: could not allocate cub_atom_num, &
                &cub_chg; exiting'
     stop
  endif
  allocate(cub_coord(3,n_cub_atoms), stat=allochk)
  if (allochk .gt. 0) then
     write(6,*)'H_cube_gen_exact: could not allocate cub_coord; exiting'
     stop
  endif
  do i = 1, n_cub_atoms
     read (13,'(i5,4f12.6)',iostat=istat) cub_atom_num(i),&
 &                                       cub_chg(i),&
 &                                       cub_coord(1,i),&
 &                                       cub_coord(2,i),&
 &                                       cub_coord(3,i)
  enddo
  if (istat < 0) then
     close(13)
     write(6,*)'H_cube_gen_exact: Error reading cube_file; exiting'
     stop
  endif
  close(13)

! make sure coords are in au
  if (num_a .lt. 0) then
     do i = 1, 3
        do j = 1, n_cub_atoms
           cub_coord(i,j) = cub_coord(i,j)/AtoBohr
        enddo
     enddo
  endif

! check if we have same coords for molecule and cube, and # of grid pts
  do i = 1, n_cub_atoms
     if((site_info%site_crds(3*(i-1)+1)-cub_coord(1,i) .gt. small) .or. &
        (site_info%site_crds(3*(i-1)+2)-cub_coord(2,i) .gt. small) .or. &
        (site_info%site_crds(3*(i-1)+3)-cub_coord(3,i) .gt. small))then
        write(6,*)'H_cube_gen_exact: Coordinates from molecule and cube &
                   &do not match, exiting'
        stop
     endif
  enddo
  if (num_a .gt. size .or. num_b .gt. size .or. num_c .gt. size) then 
     write(6,*),'H_cube_gen_exact: Size of grid in cube_file too big;&
                 & exiting'
     stop
  endif

! allocate exact cubes
  allocate(ex_cubes%esp_cube(num_c,num_b,num_a),&
           ex_cubes%x_cube(num_c,num_b,num_a),&
           ex_cubes%y_cube(num_c,num_b,num_a),&
           ex_cubes%z_cube(num_c,num_b,num_a),&
           ex_cubes%dens_cube(num_c,num_b,num_a),&
           stat=allochk)
  if (allochk .gt. 0) then
     write(6,*)'H_cube_gen_exact: could not allocate cubes; exiting'
     stop
  endif

! if cube files are already there return

  inquire(file = 'exact_DENS.cube', exist = cub5)
  if (cub5) then
     write(6,*)'H_cube_gen_exact: Cube files found, skipping cube calculation'
     do i = 2, 2
           open(13,file='exact_ESP.cube',status='unknown')
        read (13,*) 
        read (13,*) 
        read (13,'(i5,3f12.6)',iostat=istat) tmp_cub_atm,tmp_x,tmp_y,tmp_z
        read (13,'(i5,3f12.6)',iostat=istat) tmp_num_a, tmp_x1, tmp_y1, tmp_z1
        read (13,*) 
        read (13,*) 
        if (tmp_cub_atm .ne. n_cub_atoms .or. tmp_x .ne. x_orig .or. &
            tmp_y .ne. y_orig .or. tmp_z .ne. z_orig .or. tmp_num_a .ne. &
            num_a .or. tmp_x1 .ne. x1 .or. tmp_y1 .ne. y1 .or. tmp_z1 .ne. &
            z1) then
            write (6,*)'H_cube_gen_exact: "cube_file" does not match exact &
                         &cubes, check input files; exiting'
            close(13)
            stop
        endif
        do j = 1, n_cub_atoms
           read (13,*) 
        enddo
        do l = 1, num_a
           do j = 1, num_b
                 read(13, '(6E13.5)') (ex_cubes%esp_cube(k,j,l),k=1,num_c)
           enddo ! num_a
        enddo ! num_b
     close(13)
     enddo ! i
     return
  endif

! calculate cubes
  do i = 1, num_a
     do j = 1, num_b
        do k = 1, num_c
          ! Calculate each point starting from origin
           xin = x_orig + (i-1)*x1 + (j-1)*x2 + (k-1)*x3
           yin = y_orig + (i-1)*y1 + (j-1)*y2 + (k-1)*y3
           zin = z_orig + (i-1)*z1 + (j-1)*z2 + (k-1)*z3
           call H_exact_esp(site_info,basis,xin,yin,zin,esp,limitvec,&
                            natoms,num_stos,cube,twocentfit)
           ex_cubes%esp_cube(k,j,i) = esp
       enddo
     enddo
  enddo
  
! print cubes to file
  if (printcube) then
     do i = 2, 2
        if (i .eq. 1) then
           open(13,file='exact_ESP.cube',status='unknown')
        else if (i .eq. 2) then
           open(13,file='exact_DENS.cube',status='unknown')
        else if (i .eq. 3) then
           open(13,file='exact_FLD_X.cube',status='unknown')
        else if (i .eq. 4) then
           open(13,file='exact_FLD_Y.cube',status='unknown')
        else if (i .eq. 5) then
           open(13,file='exact_FLD_Z.cube',status='unknown')
        endif
        write (13,*)'title'
        write (13,*)'SCF Total Density'
        write (13,'(i5,3f12.6)') n_cub_atoms, x_orig, y_orig, z_orig
        write (13,'(i5,3f12.6)') num_a, x1, y1, z1
        write (13,'(i5,3f12.6)') num_b, x2, y2, z2
        write (13,'(i5,3f12.6)') num_c, x3, y3, z3
        do j = 1, n_cub_atoms
           write (13,'(i5,4f12.6)') cub_atom_num(j),&
                                    cub_chg(j),&
                                    cub_coord(1,j),&
                                    cub_coord(2,j),&
                                    cub_coord(3,j)
        enddo
        do l = 1, num_a
           do j = 1, num_b
              if (i .eq. 1) then
                 write(13, '(6E13.5)') (ex_cubes%esp_cube(k,j,l),k=1,num_c)
              elseif (i .eq. 2) then
                 write(13, '(6E13.5)') (ex_cubes%dens_cube(k,j,l),k=1,num_c)
              elseif (i .eq. 3) then
                 write(13, '(6E13.5)') (ex_cubes%x_cube(k,j,l),k=1,num_c)
              elseif (i .eq. 4) then
                 write(13, '(6E13.5)') (ex_cubes%y_cube(k,j,l),k=1,num_c)
              elseif (i .eq. 5) then
                 write(13, '(6E13.5)') (ex_cubes%z_cube(k,j,l),k=1,num_c)
              endif
           enddo ! num_a
        enddo ! num_b
        close(13)
     enddo ! i
  endif ! printcube

end subroutine H_cube_gen_exact
!----------------------------------------------------------
