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
!-----------------------------------------------------------------
subroutine GEOM_torsion(crd_abcd,gradphi_abcd,cosphi,sinphi)
  implicit none
  double precision crd_abcd(12),gradphi_abcd(12),cosphi,sinphi
! given coords of points a,b,c,d this routine calculates
! cosine and sine of torsion phi
! as well as gradient of phi with respect to coords of a,b,c,d
! units are in radians

  double precision rab(3),rcb(3),rdc(3),ucb(3),upab(3),updc(3),rcross(3), &
         upabc(3),upbcd(3),siz,sizcb,dotp_ab_cb,sizpab,dotp_dc_cb, &
         sizpdc,S(3)
  double precision VEC3D_unitvec,VEC3D_dotprod,VEC3D_unitperpto_unitvec,VEC3D_crossprod
  integer m
  do m = 1,3
    rab(m) = crd_abcd(m) -   crd_abcd(m+3)
    rcb(m) = crd_abcd(m+6) - crd_abcd(m+3)
    rdc(m) = crd_abcd(m+9) - crd_abcd(m+6)
  enddo
  sizcb = VEC3D_unitvec(rcb,ucb)
  dotp_ab_cb = VEC3D_dotprod(rab,ucb)
! upab is unit vector along component rab perp to ucb
  sizpab = VEC3D_unitperpto_unitvec(rab,ucb,upab)
  dotp_dc_cb = VEC3D_dotprod(rdc,ucb)
! updc is unit vector along component rdc perp to ucb
  sizpdc = VEC3D_unitperpto_unitvec(rdc,ucb,updc)
! cosine of phi is given by dot product of upab and updc
  cosphi = VEC3D_dotprod(upab,updc)
! sine of phi is given by dot product of ucb and upab x updc
  siz = VEC3D_crossprod(upab,updc,rcross)
  sinphi = VEC3D_dotprod(rcross,ucb)
! gradient of phi wrt ra is perp to abc plane---movement of ra by dr perp
! to abc plane results in dphi of dr/sizpab
! perp to abc given by upab x ucb  (these are orthogonal unit vectors)
  siz = VEC3D_crossprod(upab,ucb,upabc)
! grad of phi wrt rd is perp to bcd plane--calc sim to grad phi wrt ra
! perp given by updc x ucb or ucb x updc
  siz = VEC3D_crossprod(ucb,updc,upbcd)
! now have enough for gradphi for a and d
  do m = 1,3
    gradphi_abcd(m) = upabc(m) / sizpab
    gradphi_abcd(9+m) = upbcd(m) / sizpdc
  enddo
! following chap 5 of thesis of Bekker we have grad phi wrt b = -grad phi wrt a
! plus some vec S and rad phi wrt c = -grad phi wrt d - S
! S is perp to rcb; using simple torque rule and identity for 
! triple cross product he derives S (eqn 5.20)
  do m = 1,3
    S(m) = (dotp_ab_cb/sizcb)*gradphi_abcd(m) + &
           (dotp_dc_cb/sizcb)*gradphi_abcd(m+9)
    gradphi_abcd(m+3) = S(m) - gradphi_abcd(m)
    gradphi_abcd(m+6) = -S(m) - gradphi_abcd(m+9)
  enddo

  return
end subroutine GEOM_torsion
!-------------------------------------------------------------
