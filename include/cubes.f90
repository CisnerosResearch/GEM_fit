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
module cubes

implicit none

  type exact_cubes

       integer n_cub_atoms, num_a, num_b, num_c
  
       double precision x_orig,y_orig,z_orig,x1,y1,z1,x2,y2,z2,x3,y3,z3

       double precision,allocatable:: esp_cube(:,:,:), x_cube(:,:,:),&
                                      y_cube(:,:,:), z_cube(:,:,:), &
                                      dens_cube(:,:,:)

  end type exact_cubes

  type auxiliary_cubes

       double precision,allocatable:: esp_cube(:,:,:,:), x_cube(:,:,:,:),&
                                      y_cube(:,:,:,:), z_cube(:,:,:,:),&
                                      dens_cube(:,:,:,:)

  end type auxiliary_cubes

end module cubes 
