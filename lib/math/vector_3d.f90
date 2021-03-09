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
!------------------------------------------------------------
subroutine VEC3D_get_perp_to_vecs(v1,v2,p,dp_dv1_p,dp_dv2_p)
  implicit none
  double precision v1(3),v2(3),p(3),dp_dv1_p(3),dp_dv2_p(3)
! given two vectors get the unit cross product p = v1 x v2 / |v1 x v2|
! a change dv1 in v1 within v1v2 plane has no effect on p
! only a change dv1_p along p can change p---equiv to rotation of siz
!  dv1_p / |v1| about axis p x v1/|v1| : change in p is therefore
! dp = (-dv1_p/|v1|)*v1/|v1| thus
! dp_dv1_p = -v1 / (v1*v1)
  double precision siz,u1(3),u2(3),u1p(3),u2p(3),dotp
  
  siz = v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3)
  u1(1) = v1(1) / siz
  u1(2) = v1(2) / siz
  u1(3) = v1(3) / siz
  siz = v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3)
  u2(1) = v2(1) / siz
  u2(2) = v2(2) / siz
  u2(3) = v2(3) / siz
! p = u1 x u2 normed since u1 not perp to u2 necessarily
  p(1) = u1(2)*u2(3) - u1(3)*u2(2)
  p(2) = u1(3)*u2(1) - u1(1)*u2(3)
  p(3) = u1(1)*u2(2) - u1(2)*u2(1)
  siz = sqrt(p(1)*p(1)+p(2)*p(2)+p(3)*p(3))
  p(1) = p(1)/siz
  p(2) = p(2)/siz
  p(3) = p(3)/siz
! u1p is unit vec perp to u1 and p : p x u1
  u1p(1) = p(2)*u1(3) - p(3)*u1(2)
  u1p(2) = p(3)*u1(1) - p(1)*u1(3)
  u1p(3) = p(1)*u1(2) - p(2)*u1(1)
! u2p is unit vec perp to u2 and p : u2 x p
  u2p(1) = u2(2)*p(3) - u2(3)*p(2)
  u2p(2) = u2(3)*p(1) - u2(1)*p(3)
  u2p(3) = u2(1)*p(2) - u2(2)*p(1)
! change in v1 along p gives dp along -u2p
  dotp = v1(1)*u2p(1)+v1(2)*u2p(2)+v1(3)*u2p(3)
  dp_dv1_p(1) = -u2p(1)/dotp
  dp_dv1_p(2) = -u2p(2)/dotp
  dp_dv1_p(3) = -u2p(3)/dotp
! change in v2 along p gives dp along -u1p
  dotp = v2(1)*u1p(1)+v2(2)*u1p(2)+v2(3)*u1p(3)
  dp_dv2_p(1) = -u1p(1)/dotp
  dp_dv2_p(2) = -u1p(2)/dotp
  dp_dv2_p(3) = -u1p(3)/dotp
  return
end subroutine VEC3D_get_perp_to_vecs
!----------------------------------------------------------
double precision function VEC3D_dotprod(v1,v2)
implicit none
  double precision v1(3),v2(3)
  VEC3D_dotprod = v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
  return
end function VEC3D_dotprod
!----------------------------------------------------------
double precision function VEC3D_unitvec(v,u)
implicit none
! changes v into unit vector along v--returns size of v
  double precision v(3),u(3),siz
  siz = sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
  u(1) = v(1)/siz
  u(2) = v(2)/siz
  u(3) = v(3)/siz
  VEC3D_unitvec = siz
  return
end function VEC3D_unitvec
!----------------------------------------------------------
double precision function VEC3D_crossprod(v1,v2,v3)
implicit none
! v3 = v1 x v2;  returns size of v3
  double precision v1(3),v2(3),v3(3)
  v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
  v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
  v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
  VEC3D_crossprod = sqrt(v3(1)*v3(1)+v3(2)*v3(2)+v3(3)*v3(3))
  return
end function VEC3D_crossprod
!----------------------------------------------------------
double precision function VEC3D_project_unitvec(v,u,w)
implicit none
! projects v onto unit vector u to get w--returns length of w
  double precision v(3),u(3),w(3),dot
  dot = v(1)*u(1)+v(2)*u(2)+v(3)*u(3)
  w(1) = dot*u(1)
  w(2) = dot*u(2)
  w(3) = dot*u(3)
  VEC3D_project_unitvec = dot
  return
end function VEC3D_project_unitvec
!----------------------------------------------------------
double precision function VEC3D_perpto_unitvec(v,u,w)
implicit none
! removes component of v along unit vector u --returns length of new v
  double precision v(3),u(3),w(3),siz,dot
  dot = v(1)*u(1)+v(2)*u(2)+v(3)*u(3)
  w(1) = v(1) - dot*u(1)
  w(2) = v(2) - dot*u(2)
  w(3) = v(3) - dot*u(3)
  siz = sqrt(w(1)*w(1)+w(2)*w(2)+w(3)*w(3))
  VEC3D_perpto_unitvec = siz
  return
end function VEC3D_perpto_unitvec
!----------------------------------------------------------
double precision function VEC3D_unitperpto_unitvec(v,u,w)
implicit none
! removes component of v along unit vector u --returns length of new w
! normalizes resulting w
  double precision v(3),u(3),w(3),siz,dot
  dot = v(1)*u(1)+v(2)*u(2)+v(3)*u(3)
  w(1) = v(1) - dot*u(1)
  w(2) = v(2) - dot*u(2)
  w(3) = v(3) - dot*u(3)
  siz = sqrt(w(1)*w(1)+w(2)*w(2)+w(3)*w(3))
  w(1) = w(1) / siz
  w(2) = w(2) / siz
  w(3) = w(3) / siz
  VEC3D_unitperpto_unitvec = siz
  return
end function VEC3D_unitperpto_unitvec
!---------------------------------------------------------
subroutine dot(v1,v2,result)
  implicit none
  double precision v1(3),v2(3),result
  result = v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
  return
end subroutine dot
!--------------------------------------------------------
subroutine cross(v1,v2,v12)
  implicit none
!
!    v12 is cross product of v1 and v2
!
  double precision v1(3),v2(3),v12(3)
  v12(1) = v1(2)*v2(3)-v1(3)*v2(2)
  v12(2) = v1(3)*v2(1)-v1(1)*v2(3)
  v12(3) = v1(1)*v2(2)-v1(2)*v2(1)
  return
end subroutine cross
!--------------------------------------------------------
subroutine VEC3D_3x3_transpose(A,At)
  implicit none
  double precision A(3,3),At(3,3)

  integer i,j
  do i = 1,3
    do j = 1,3
      At(i,j) = A(j,i)
    enddo
  enddo

  return
end subroutine VEC3D_3x3_transpose
!--------------------------------------------------------
