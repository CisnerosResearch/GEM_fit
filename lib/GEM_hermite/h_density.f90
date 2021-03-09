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
subroutine H_dens_aux(site_info,auxis,xin,yin,zin,dens,ncoeff)
  ! GAC: subroutine to calculate density at given point from auxiliarys.
use definition
  implicit none
  type(sites_info)::site_info
  type(aux_orbitals)::auxis(ncoeff)

  integer k1,num1,off1,c_off1 
  integer ncoeff,x1,x2,y1,y2,z1,z2
  integer ires,jres,n,i,j,deg1,deg2,fldint
  double precision norm1,norm2,dens_int,dens,xin,yin,zin,R2
  double precision dx,dy,dz,aux_coefs(ncoeff)
  double precision expon1

  dens = 0.0d0
  aux_coefs(:) = site_info%cartesian_coeffs(:)

  dens = 0.0d0
  do i = 1, ncoeff
     expon1 = auxis(i)%expo
     x1 = auxis(i)%x
     y1 = auxis(i)%y
     z1 = auxis(i)%z
     !dx = auxis(i)%coords(1) - xin
     !dy = auxis(i)%coords(2) - yin
     !dz = auxis(i)%coords(3) - zin
     dx = xin - auxis(i)%coords(1)
     dy = yin - auxis(i)%coords(2)
     dz = zin - auxis(i)%coords(3)
     R2 = dx*dx + dy*dy + dz*dz
     norm1 = auxis(i)%norm
     dens_int = (dx**x1)*(dy**y1)*(dz**z1)*exp(-expon1*R2)
     dens = dens + dens_int*norm1*aux_coefs(i)
  enddo
  return
end subroutine H_dens_aux
!----------------------------------------------------------
