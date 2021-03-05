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
!---------------------------------------------------------
subroutine H_esp_nuc(site_info,cx,cy,cz,total)
  ! GAC: subroutine to calculate ESP with cartesian gaussians.
use definition
  implicit none
  type(sites_info)::site_info

  integer site1,i,x1,y1,z1
  integer slo1,shi1,ires,natoms
  double precision tot_N,total
  double precision dx,dy,dz,cx,cy,cz,dist
  double precision nuc_chg1
  double precision,parameter:: tiny=1.d-12

  natoms = site_info%natoms
  slo1 = site_info%residue_start(1)
  shi1 = site_info%residue_end(1)

! do nuclear ESP
  tot_N = 0.0d0
  do site1 = slo1,shi1
    dx = site_info%site_crds(3*(site1-1)+1) - cx
    if (abs(dx) .lt. tiny) dx = 0.0d0
    dy = site_info%site_crds(3*(site1-1)+2) - cy
    if (abs(dy) .lt. tiny) dy = 0.0d0
    dz = site_info%site_crds(3*(site1-1)+3) - cz
    if (abs(dz) .lt. tiny) dz = 0.0d0
    dist = dsqrt(dx*dx+dy*dy+dz*dz)
    nuc_chg1 = site_info%nuclear_charge(site1)
    if(site1 .le. natoms .and. dist .gt. tiny) then
      tot_N = tot_N + nuc_chg1 / dist
    endif
  enddo
  
  total = tot_N 

  return
end subroutine H_esp_nuc
!----------------------------------------------------------
