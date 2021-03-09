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
! computes the cross product of two vectors in 3D.
! M = X x Y
subroutine vec3_cross(M,X,Y)
  implicit none
  double precision M(3), X(3), Y(3)
  M(1) = X(2) * Y(3) - X(3) * Y(2)
  M(2) = X(3) * Y(1) - X(1) * Y(3)
  M(3) = X(1) * Y(2) - X(2) * Y(1)
end

! Inverts a 3x3 matrix:  A = inverse of B
subroutine mat3_inv(a,b)
  implicit none
  double precision a(9),b(9),c

  c = (b(5)*b(9)-b(6)*b(8))*b(1) &
     +(b(6)*b(7)-b(4)*b(9))*b(2) &
     +(b(4)*b(8)-b(5)*b(7))*b(3)
  a(1)=(b(5)*b(9)-b(6)*b(8))/c
  a(4)=(b(6)*b(7)-b(4)*b(9))/c
  a(7)=(b(4)*b(8)-b(5)*b(7))/c
  a(2)=(b(8)*b(3)-b(9)*b(2))/c
  a(5)=(b(9)*b(1)-b(7)*b(3))/c
  a(8)=(b(7)*b(2)-b(8)*b(1))/c
  a(3)=(b(2)*b(6)-b(3)*b(5))/c
  a(6)=(b(3)*b(4)-b(1)*b(6))/c
  a(9)=(b(1)*b(5)-b(2)*b(4))/c
end subroutine mat3_inv

subroutine mat3_transpose(M,Mtrans)
  implicit none
  double precision M(3,3),Mtrans(3,3)

  Mtrans(1,1)=M(1,1)
  Mtrans(1,2)=M(2,1)
  Mtrans(1,3)=M(3,1)

  Mtrans(2,1)=M(1,2)
  Mtrans(2,2)=M(2,2)
  Mtrans(2,3)=M(3,2)

  Mtrans(3,1)=M(1,3)
  Mtrans(3,2)=M(2,3)
  Mtrans(3,3)=M(3,3)
end subroutine mat3_transpose

subroutine mat3_trans_inplace(M)
  implicit none
  double precision M(3,3), n
  n=M(2,1); M(2,1)=M(1,2); M(1,2)=n
  n=M(3,1); M(3,1)=M(1,3); M(1,3)=n
  n=M(3,2); M(3,2)=M(2,3); M(2,3)=n
end subroutine mat3_trans_inplace

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Return phase in degrees from a complex number
function phase(x)
  implicit none
  complex*16 x
  include "math.fh"
  real*8 phase
  real*8 a,b

  !a=real(x)
  !b=aimag(x)
  a=real(real(x))
  b=real(aimag(x))

  phase=atan2(b,a)*M_R2D

end function phase

