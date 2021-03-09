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
subroutine H_cube_aux_coef(site_info,auxis,aux_cubes,ex_cubes,ncoeff,cube)
use definition
use cubes
  implicit none
  type(sites_info)::site_info
  type(aux_orbitals)::auxis(ncoeff)
  type(auxiliary_cubes)::aux_cubes
  type(exact_cubes)::ex_cubes

  integer i,j,k,l,m,n
  integer n_cub_atoms, num_a, num_b, num_c, num_stos, allochk, istat
  integer tmp_cub_atm, tmp_num_a, ncoeff
  integer,allocatable::cub_atom_num(:)
  integer,parameter::size = 250
  double precision esp,field(3),x_orig,y_orig,z_orig,x1,y1,z1,&
                   x2,y2,z2,x3,y3,z3,xin,yin,zin,dens
  double precision tmp_x,tmp_y,tmp_z,tmp_x1,tmp_y1,tmp_z1,dist
  double precision,allocatable::cub_chg(:),cub_coord(:,:)
  double precision,allocatable::diff_esp(:,:,:),&
                                diff_dens(:,:,:),&
                                diff_fldx(:,:,:),diff_fldy(:,:,:),&
                                diff_fldz(:,:,:)
  double precision,parameter::AtoBohr = 0.529177249d0,small=1.d-6,&
                              cutoff = 0.0d0
                              !cutoff = 2.0d0 ! for exch-coul diff dens
  logical cube, file_there, skip
  character*1 int_type

! get cube info
  n_cub_atoms = ex_cubes%n_cub_atoms
  num_a = ex_cubes%num_a
  num_b = ex_cubes%num_b
  num_c = ex_cubes%num_c
  x_orig = ex_cubes%x_orig
  y_orig = ex_cubes%y_orig
  z_orig = ex_cubes%z_orig
  x1 = ex_cubes%x1
  y1 = ex_cubes%y1
  z1 = ex_cubes%z1
  x2 = ex_cubes%x2
  y2 = ex_cubes%y2
  z2 = ex_cubes%z2
  x3 = ex_cubes%x3
  y3 = ex_cubes%y3
  z3 = ex_cubes%z3
  !!! Check if number of atoms in info file matches
  if (n_cub_atoms .ne. site_info%natoms) then
     write(6,*)'H_cube_aux: Number of atoms in cube_file &
                &does not match, exiting'
     stop
  endif

 ! get rest of info for cube
  inquire(file = 'cube_file', exist = file_there)
  if (.not. file_there) then
     write(6,*)'H_cube_aux_coef: Could not find cube info file &
                &(cube_file), exiting'
     stop
  endif
  open(13,file='cube_file',status='unknown')
  read (13,'(i5,3f12.6)',iostat=istat) 
  read (13,'(i5,3f12.6)',iostat=istat)
  read (13,'(i5,3f12.6)',iostat=istat)
  read (13,'(i5,3f12.6)',iostat=istat)
  !!! allocate local vars
  allocate(cub_atom_num(n_cub_atoms),cub_chg(n_cub_atoms), stat=allochk)
  if (allochk .gt. 0) then
     write(6,*)'H_cube_aux_coef: could not allocate cub_atom_num, &
                &cub_chg; exiting'
     stop
  endif
  allocate(cub_coord(3,n_cub_atoms), stat=allochk)
  if (allochk .gt. 0) then
     write(6,*)'H_cube_aux_coef: could not allocate cub_coord; exiting'
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
     write(6,*)'H_cube_aux_coef: Error reading cube_file; exiting'
     stop
  endif
  close(13)

! allocate auxiliary cubes
  allocate(aux_cubes%esp_cube(1,num_c,num_b,num_a),&
           stat=allochk)
  if (allochk .gt. 0) then
     write(6,*)'H_cube_aux: could not allocate cubes; exiting'
     stop
  endif

! allocate temporary cubes for differences
  allocate(diff_esp(num_c,num_b,num_a),&
           stat=allochk)
  if (allochk .gt. 0) then
     write(6,*)'H_cube_aux: could not allocate difference cubes; exiting'
     stop
  endif

! calculate cubes
  do i = 1, num_a
     do j = 1, num_b
        do k = 1, num_c
          ! Calculate each point starting from origin
           xin = x_orig + (i-1)*x1 + (j-1)*x2 + (k-1)*x3
           yin = y_orig + (i-1)*y1 + (j-1)*y2 + (k-1)*y3
           zin = z_orig + (i-1)*z1 + (j-1)*z2 + (k-1)*z3
           do n = 1, n_cub_atoms
              dist = sqrt( (site_info%site_crds(3*(n-1)+1)-xin)**2 + &
                           (site_info%site_crds(3*(n-1)+2)-yin)**2 + &
                           (site_info%site_crds(3*(n-1)+3)-zin)**2 ) 
              if (dist .lt. cutoff) skip = .true.
           enddo

           call H_calc_esp_CART(site_info,auxis,xin,yin,zin,esp,&
                                ncoeff,.true.,cube)
           aux_cubes%esp_cube(1,k,j,i) = esp
           diff_esp(k,j,i) = abs(esp) - abs(ex_cubes%esp_cube(k,j,i))

       enddo
     enddo
  enddo

! Print difference cubes
     open(13,file='diff_ESP.cube',status='unknown')

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
              write(13, '(6E13.5)') (diff_esp(k,j,l),k=1,num_c)
        enddo ! num_a
     enddo ! num_b
     close(13)

! Print auxiliary cubes

        open(13,file='aux_ESP.cube',status='unknown')

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
              write(13, '(6E13.5)') (aux_cubes%esp_cube(1,k,j,l),k=1,num_c)
        enddo ! num_a
     enddo ! num_b
     close(13)

! deallocate arrays
  !deallocate(aux_cubes%esp_cube,aux_cubes%x_cube,aux_cubes%y_cube,&
  !           aux_cubes%z_cube,aux_cubes%dens_cube,&
  !           diff_fld,diff_esp,diff_dens)
  deallocate(diff_esp,aux_cubes%esp_cube)

end subroutine H_cube_aux_coef
!----------------------------------------------------------
