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
subroutine H_spherical_grid2(site_info,basis,ex_cubes,rad_grid,wght,npts,&
                             num_stos,densfile,dtype,printcube,debug,&
                             cube,gemesp)
use definition
use cubes
  implicit none
  type(sites_info)::site_info
  type(ao_orbitals)::basis(num_stos)
  type(exact_cubes)::ex_cubes

  integer i,j,k,l,m,p,q, natoms,num_stos, istat, allochk, rdstat
  integer npts, ntheta, nang, nphi, nrad, count, ntheta2, npts2, count2
  double precision rmin, rmax, delta_r, delta_ang, circum, delta_phi, theta
  double precision rad, angle, rmin2, rmax2, delta_r2, pi, xin, yin, zin, phi
  double precision esp, dens, radx, rady, radz, field(3), d1, d2, d3
  double precision rad_grid(npts,3), wght(npts), tx,ty,tz,wt
  double precision,allocatable::ang_mat(:,:),phi_vec(:),theta_vec(:)
  double precision,parameter::AtoBohr = 0.529177249d0,small=1.d-6
  character(len=*) densfile
  character(len=3) dtype
  character(len=150) cmd
  logical cube, file_there, printcube, cub1, cub2, cub3, cub4, cub5, &
          basalloc,twocentfit,debug,read_sphere,pt_there,rdsph2,rdsph1,gemesp
  !character*1 int_type

  cub1 = .false.

! Create basis info
  call H_form_basis(site_info,basis,num_stos,densfile)

! allocate exact cubes
  allocate(ex_cubes%esp_cube(npts,1,1),&
           stat=allochk)
  if (allochk .gt. 0) then
     write(6,*)'H_cube_gen_exact: could not allocate cubes; exiting'
     stop
  endif

  ! Check to see if grided property exists, if it does read it in and return
  inquire(file = 'sphere_ESP.cube', exist = cub1)
  if (cub1) then
     write(6,*)'H_spherical_grid: ESP Grid file found, skipping ESP calculation'
        open(13,file='sphere_ESP.cube',status='unknown')
         rdstat = 0
         npts2 = -1
         do while(rdstat==0)
            read(13,*,iostat=rdstat) tx, ty, tz, wt
            npts2 = npts2 + 1
         enddo
         rewind(13)
        !read (13,*) npts2
        if (npts .ne. npts2) then
            write (6,*)'H_spherical_grid: "sphere_file" does not match exact &
                         &grid, check input files; exiting'
            close(13)
            stop
        endif
        do l = 1, npts
           !read(13,*) ex_cubes%esp_cube(l,1,1)
           read(13,*) d1,d2,d3,ex_cubes%esp_cube(l,1,1)
        enddo ! npts
     close(13)
     return
  endif

! calculate ESP on each grid point
  if (gemesp) then ! calculate ESP with GEM integrals

     !print *,'USING GEM FOR ESP',cube
     do i = 1, npts
       ! Calculate each point from radial grid
        xin = rad_grid(i,1)
        yin = rad_grid(i,2)
        zin = rad_grid(i,3)
        call H_exact_esp3(site_info,basis,xin,yin,zin,esp,num_stos,cube)
        ex_cubes%esp_cube(i,1,1) = esp
     enddo
  
   ! print cubes to file
1001 format (4(1X,F16.10))
     open(13,file='sphere_ESP.cube',status='unknown')
     do l = 1, npts
        write(13,1001) rad_grid(l,1),rad_grid(l,2),rad_grid(l,3),&
                    ex_cubes%esp_cube(l,1,1)
     enddo ! npts
     close(13)

  else ! USE CUBEGEN (default)

     write(cmd,*)'cubegen 0 potential=',trim(dtype),' ',trim(densfile),' sphere_ESP.cube -5 n < grd_pts2.xyz'
     write(6,*) cmd
     call system(cmd)
     open(13,file='sphere_ESP.cube',status='unknown')
     do l = 1, npts
        read(13,*) tx,ty,tz,ex_cubes%esp_cube(l,1,1)
     enddo ! npts
     close(13)
 
  endif

end subroutine H_spherical_grid2
!----------------------------------------------------------
