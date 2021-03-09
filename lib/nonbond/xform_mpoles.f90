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
!----------------------------------------------------------------
subroutine XFORM_MPOLE_matrix(A_xy,Mpole_xy,order)
! A_xy is matrix of coefficients expanding y in terms of x
! i.e. y_i = A_i1 * x_1 + A_i2 * x_2 + A_i3 * x_3
! Mpole_xy is resulting matrix getting expansion of multipoles
! with basis y in terms of those with basis x
! order is order of multipoles 1 for charge, 4 for charge-dipoles
! 10 for charge,dipole,quadrupole, 20 for up to octupole and
! finally 35 for up to hexadecapole
 
! Mpole 1 is charge... Mpole 10 is x_3,x_3 quadrupolar coefficients etc.
  implicit none
  integer order
  double precision A_xy(3,3), Mpole_xy(order,order)
  integer i,j,k,l,m,n,p,ind1,ind2,jj,kk
  integer qind1(6),qind2(6),oind1(10),oind2(10),oind3(10)
  integer hind1(15),hind2(15),hind3(15),hind4(15)
  common/MxfInd/qind1,qind2,oind1,oind2,oind3,hind1,hind2,hind3,hind4
  data qind1   /1,2,3,1,1,2/
  data qind2   /1,2,3,2,3,3/
  data oind1   /1,2,3,1,1,2,2,3,3,1/
  data oind2   /1,2,3,1,1,2,2,3,3,2/
  data oind3   /1,2,3,2,3,1,3,1,2,3/
  data hind1   /1,2,3,1,1,2,2,3,3,1,1,2,1,2,3/
  data hind2   /1,2,3,1,1,2,2,3,3,1,1,2,1,2,3/
  data hind3   /1,2,3,1,1,2,2,3,3,2,3,3,2,1,1/
  data hind4   /1,2,3,2,3,1,3,1,2,2,3,3,3,3,2/
      

! CHARGE case
  Mpole_xy(1,1) = 1.d0
  if ( order .eq. 1 )return
  do j = 2,order
   Mpole_xy(1,j) = 0.d0
   Mpole_xy(j,1) = 0.d0
  enddo
! DIPOLES
! d_yk = A_xy(k,1)*d_x1 + A_xy(k,2)*d_x2 + A_xy(k,3)*d_x3
! D'_j
  do j = 2,4
    do k = 2,4
      Mpole_xy(j,k) = A_xy(j-1,k-1)
    enddo
  enddo
  if ( order .eq. 4 )return
  do j = 2,4
    do k = 5,order
      Mpole_xy(j,k) = 0.d0
      Mpole_xy(k,j) = 0.d0
    enddo
  enddo
! QUADRUPOLES
! q_ykyl = sum_i,j (A_xy(k,i)*A_xy(l,j)+A_xy(k,j)*A_xy(l,i))*(q_xixj+q_xjxi)
! Mp(5) = q_y1y1,..,Mp(7) = q_y3y3,Mp(8) = q_y1y2+q_y2y1,.,Mp(10)=q_y2y3+q_y3y2
! Q'_kk
  do ind1 = 1,3
    k = qind1(ind1)
    do ind2 = 1,6
      i = qind1(ind2)
      j = qind2(ind2)
      Mpole_xy(ind1+4,ind2+4) = A_xy(k,i)*A_xy(k,j)
    enddo
  enddo
! Q'_kl
  do ind1 = 4,6
    k = qind1(ind1)
    l = qind2(ind1)
    do ind2 = 1,6
      i = qind1(ind2)
      j = qind2(ind2)
      Mpole_xy(ind1+4,ind2+4) = A_xy(k,i)*A_xy(l,j) + A_xy(k,j)*A_xy(l,i)
    enddo
  enddo
  if ( order .eq. 10 )return
  do j = 5,10
    do k = 11,order
      Mpole_xy(k,j) = 0.d0
      Mpole_xy(j,k) = 0.d0
    enddo
  enddo
! OCTUPOLES
! O'_lll
  do ind1 = 1,3
    l = oind1(ind1)
    do ind2 = 1,10
      i = oind1(ind2)
      j = oind2(ind2)
      k = oind3(ind2)
      Mpole_xy(ind1+10,ind2+10) = A_xy(l,i)*A_xy(l,j)*A_xy(l,k)
    enddo
  enddo
! O'_llm
  do ind1 = 4,9
    l = oind1(ind1)
    m = oind3(ind1)
    do ind2 = 1,9
      i = oind1(ind2)
      k = oind3(ind2)
      Mpole_xy(ind1+10,ind2+10) =     &
           A_xy(l,i)*A_xy(l,i)*A_xy(m,k) +   &
           2.d0 * A_xy(l,i)*A_xy(m,i)*A_xy(l,k)
    enddo
    Mpole_xy(ind1+10,20) =   &
           A_xy(l,1)*A_xy(l,2)*A_xy(m,3) +   &
           A_xy(l,1)*A_xy(m,2)*A_xy(l,3) +   &
           A_xy(m,1)*A_xy(l,2)*A_xy(l,3)
  enddo
! O'_123
  Mpole_xy(20,11) = 6.d0*A_xy(1,1)*A_xy(2,1)*A_xy(3,1)
  Mpole_xy(20,12) = 6.d0*A_xy(1,2)*A_xy(2,2)*A_xy(3,2)
  Mpole_xy(20,13) = 6.d0*A_xy(1,3)*A_xy(2,3)*A_xy(3,3)
  do ind2 = 4,9
    i = oind1(ind2)
    k = oind3(ind2)
    Mpole_xy(20,10+ind2) = 2.d0*   &
           (  A_xy(1,i)*A_xy(2,i)*A_xy(3,k) +   &
              A_xy(1,i)*A_xy(2,k)*A_xy(3,i) +   &
              A_xy(1,k)*A_xy(2,i)*A_xy(3,i) )
  enddo
  Mpole_xy(20,20) =    &
              A_xy(1,1)*A_xy(2,2)*A_xy(3,3) +   &
              A_xy(1,1)*A_xy(3,2)*A_xy(2,3) +   &
              A_xy(2,1)*A_xy(1,2)*A_xy(3,3) +   &
              A_xy(2,1)*A_xy(3,2)*A_xy(1,3) +   &
              A_xy(3,1)*A_xy(1,2)*A_xy(2,3) +   &
              A_xy(3,1)*A_xy(2,2)*A_xy(1,3) 
  if ( order .eq. 20 )return
  do j = 11,20
    do k = 21,order
      Mpole_xy(k,j) = 0.d0
      Mpole_xy(j,k) = 0.d0
    enddo
  enddo
! HEXADECAPOLES
! H'_mmmm
  do ind1 = 1,3
    m = hind1(ind1)
    do ind2 = 1,15
      i = hind1(ind2)
      j = hind2(ind2)
      k = hind3(ind2)
      l = hind4(ind2)
      Mpole_xy(ind1+20,ind2+20) = A_xy(m,i)*A_xy(m,j)*A_xy(m,k)*A_xy(m,l)
    enddo
  enddo
! H'_mmmn
  do ind1 = 4,9
    m = hind1(ind1)
    n = hind4(ind1)
! can put iiii & iiij together
    do ind2 = 1,9
      i = hind1(ind2)
      j = hind4(ind2)
      Mpole_xy(ind1+20,ind2+20) =   &
               A_xy(m,i)*A_xy(m,i)*A_xy(m,i)*A_xy(n,j) + 3.d0*  &
               A_xy(m,i)*A_xy(m,i)*A_xy(n,i)*A_xy(m,j) 
    enddo
! can put iijj & iijk together
    do ind2 = 10,15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
      Mpole_xy(ind1+20,ind2+20) =    &
               A_xy(m,i)*A_xy(m,i)*A_xy(m,j)*A_xy(n,k) +        &
               A_xy(m,i)*A_xy(m,i)*A_xy(m,k)*A_xy(n,j) + 2.d0*  &
               A_xy(m,i)*A_xy(m,j)*A_xy(m,k)*A_xy(n,i) 
    enddo
  enddo
! H'_mmnn
  do ind1 = 10,12
    m = hind1(ind1)
    n = hind3(ind1)
! can put iiii & iiij together
    do ind2 = 1,9
      i = hind1(ind2)
      j = hind4(ind2)
      Mpole_xy(ind1+20,ind2+20) = 3.d0 * (   &
               A_xy(m,i)*A_xy(m,i)*A_xy(n,i)*A_xy(n,j) +    &
               A_xy(m,i)*A_xy(m,j)*A_xy(n,i)*A_xy(n,i)  )
    enddo
! can put iijj & iijk together
    do ind2 = 10,15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
      Mpole_xy(ind1+20,ind2+20) =   &
               A_xy(m,i)*A_xy(m,i)*A_xy(n,j)*A_xy(n,k) +         &
               A_xy(m,j)*A_xy(m,k)*A_xy(n,i)*A_xy(n,i) +  2.d0*  &
               A_xy(m,i)*A_xy(m,j)*A_xy(n,i)*A_xy(n,k) +  2.d0*  &
               A_xy(m,i)*A_xy(m,k)*A_xy(n,i)*A_xy(n,j) 
    enddo
  enddo
! H'_mmnp
  do ind1 = 13,15
    m = hind1(ind1)
    n = hind3(ind1)
    p = hind4(ind1)
! can put iiii & iiij together
    do ind2 = 1,9
      i = hind1(ind2)
      j = hind4(ind2)
      Mpole_xy(ind1+20,ind2+20) =     &
              3.d0*A_xy(m,i)*A_xy(m,i)*A_xy(n,i)*A_xy(p,j) +     &
              3.d0*A_xy(m,i)*A_xy(m,i)*A_xy(n,j)*A_xy(p,i) +     &
              6.d0*A_xy(m,i)*A_xy(m,j)*A_xy(n,i)*A_xy(p,i) 
    enddo
! can put iijj & iijk together
    do ind2 = 10,15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
      Mpole_xy(ind1+20,ind2+20) =     &
               A_xy(m,i)*A_xy(m,i)*A_xy(n,j)*A_xy(p,k) +     &
               A_xy(m,i)*A_xy(m,i)*A_xy(n,k)*A_xy(p,j) + 2.d0 * (    &
               A_xy(m,i)*A_xy(m,j)*A_xy(n,i)*A_xy(p,k) +     &
               A_xy(m,i)*A_xy(m,j)*A_xy(n,k)*A_xy(p,i) +     &
               A_xy(m,i)*A_xy(m,k)*A_xy(n,i)*A_xy(p,j) +     &
               A_xy(m,i)*A_xy(m,k)*A_xy(n,j)*A_xy(p,i) +     &
               A_xy(m,j)*A_xy(m,k)*A_xy(n,i)*A_xy(p,i) )
    enddo
  enddo
  return
end subroutine XFORM_MPOLE_matrix
!------------------------------------------------------------
subroutine XFORM_MPOLE_deriv_matrix(A_xy,DA_xy,DMpole_xy,order)
! calculate derivative of Mpole_xy as DMpole_xy
! A_xy is matrix of coefficients expanding y in terms of x
! i.e. y_i = A_i1 * x_1 + A_i2 * x_2 + A_i3 * x_3
! DA_xy is deriv of A_xy wrt to some parameter
! DMpole_xy is deriv of Mpole_xy wrt to the same parameter
! order is order of multipoles 1 for charge, 4 for charge-dipoles
! 10 for charge,dipole,quadrupole, 20 for up to octupole and
! finally 35 for up to hexadecapole
 
! Mpole 1 is charge... Mpole 10 is x_3,x_3 quadrupolar coefficients etc.
  implicit none
  integer order
  double precision A_xy(3,3)
  double precision DA_xy(3,3), DMpole_xy(order,order)
  integer i,j,k,l,m,n,p,ind1,ind2,jj,kk
  integer qind1(6),qind2(6),oind1(10),oind2(10),oind3(10)
  integer hind1(15),hind2(15),hind3(15),hind4(15)
  common/MxfInd/qind1,qind2,oind1,oind2,oind3,hind1,hind2,hind3,hind4

! CHARGE case
!     Mpole_xy(1,1) = 1.d0
  DMpole_xy(1,1) = 0.d0
  if ( order .eq. 1 )return
  do j = 2,order
    DMpole_xy(1,j) = 0.d0
    DMpole_xy(j,1) = 0.d0
  enddo
! DIPOLES
! d_yk = A_xy(k,1)*d_x1 + A_xy(k,2)*d_x2 + A_xy(k,3)*d_x3
! D'_j
  do j = 2,4
    do k = 2,4
!     Mpole_xy(j,k) = A_xy(j-1,k-1)
      DMpole_xy(j,k) = DA_xy(j-1,k-1)
    enddo
  enddo
  if ( order .eq. 4 )return
  do j = 2,4
    do k = 5,order
      DMpole_xy(j,k) = 0.d0
      DMpole_xy(k,j) = 0.d0
    enddo
  enddo
! QUADRUPOLES
! q_ykyl = sum_i,j (A_xy(k,i)*A_xy(l,j)+A_xy(k,j)*A_xy(l,i))*(q_xixj+q_xjxi)
! Mp(5) = q_y1y1,..,Mp(7) = q_y3y3,Mp(8) = q_y1y2+q_y2y1,.,Mp(10)=q_y2y3+q_y3y2
! Q'_kk
  do ind1 = 1,3
    k = qind1(ind1)
    do ind2 = 1,6
      i = qind1(ind2)
      j = qind2(ind2)
!     Mpole_xy(ind1+4,ind2+4) = A_xy(k,i)*A_xy(k,j)
      DMpole_xy(ind1+4,ind2+4) = DA_xy(k,i)*A_xy(k,j) + A_xy(k,i)*DA_xy(k,j)
    enddo
  enddo
! Q'_kl
  do ind1 = 4,6
    k = qind1(ind1)
    l = qind2(ind1)
    do ind2 = 1,6
      i = qind1(ind2)
      j = qind2(ind2)
!     Mpole_xy(ind1+4,ind2+4) =   &
!            A_xy(k,i)*A_xy(l,j) + A_xy(k,j)*A_xy(l,i)
      DMpole_xy(ind1+4,ind2+4) =     &
            DA_xy(k,i)*A_xy(l,j) + DA_xy(k,j)*A_xy(l,i) +    &
             A_xy(k,i)*DA_xy(l,j) + A_xy(k,j)*DA_xy(l,i)
    enddo
  enddo
  if ( order .eq. 10 )return
  do j = 5,10
    do k = 11,order
      DMpole_xy(k,j) = 0.d0
      DMpole_xy(j,k) = 0.d0
    enddo
  enddo
! OCTUPOLES
! O'_lll
  do ind1 = 1,3
    l = oind1(ind1)
    do ind2 = 1,10
      i = oind1(ind2)
      j = oind2(ind2)
      k = oind3(ind2)
!     Mpole_xy(ind1+10,ind2+10) = A_xy(l,i)*A_xy(l,j)*A_xy(l,k)
      DMpole_xy(ind1+10,ind2+10) =    &
                DA_xy(l,i)*A_xy(l,j)*A_xy(l,k) +    &
                A_xy(l,i)*DA_xy(l,j)*A_xy(l,k) +    &
                A_xy(l,i)*A_xy(l,j)*DA_xy(l,k) 
    enddo
  enddo
! O'_llm
  do ind1 = 4,9
    l = oind1(ind1)
    m = oind3(ind1)
    do ind2 = 1,9
      i = oind1(ind2)
      k = oind3(ind2)
!     Mpole_xy(ind1+10,ind2+10) =   &
!          A_xy(l,i)*A_xy(l,i)*A_xy(m,k) +    &
!          2.d0 * A_xy(l,i)*A_xy(m,i)*A_xy(l,k)
      DMpole_xy(ind1+10,ind2+10) =   &
           DA_xy(l,i)*A_xy(l,i)*A_xy(m,k) +     &
           2.d0 * DA_xy(l,i)*A_xy(m,i)*A_xy(l,k)
      DMpole_xy(ind1+10,ind2+10) = DMpole_xy(ind1+10,ind2+10) +   &
           A_xy(l,i)*DA_xy(l,i)*A_xy(m,k) +   &
           2.d0 * A_xy(l,i)*DA_xy(m,i)*A_xy(l,k)
      DMpole_xy(ind1+10,ind2+10) = DMpole_xy(ind1+10,ind2+10) +    &
           A_xy(l,i)*A_xy(l,i)*DA_xy(m,k) +    &
           2.d0 * A_xy(l,i)*A_xy(m,i)*DA_xy(l,k)
    enddo
!   Mpole_xy(ind1+10,20) =
!          A_xy(l,1)*A_xy(l,2)*A_xy(m,3) +
!          A_xy(l,1)*A_xy(m,2)*A_xy(l,3) +
!          A_xy(m,1)*A_xy(l,2)*A_xy(l,3)
    DMpole_xy(ind1+10,20) =    &
           DA_xy(l,1)*A_xy(l,2)*A_xy(m,3) +    &
           DA_xy(l,1)*A_xy(m,2)*A_xy(l,3) +    &
           DA_xy(m,1)*A_xy(l,2)*A_xy(l,3)
    DMpole_xy(ind1+10,20) = DMpole_xy(ind1+10,20) +   &
           A_xy(l,1)*DA_xy(l,2)*A_xy(m,3) +   &
           A_xy(l,1)*DA_xy(m,2)*A_xy(l,3) +   &
           A_xy(m,1)*DA_xy(l,2)*A_xy(l,3)
    DMpole_xy(ind1+10,20) = DMpole_xy(ind1+10,20) +   &
           A_xy(l,1)*A_xy(l,2)*DA_xy(m,3) +   &
           A_xy(l,1)*A_xy(m,2)*DA_xy(l,3) +   &
           A_xy(m,1)*A_xy(l,2)*DA_xy(l,3)
  enddo
! O'_123
! Mpole_xy(20,11) = 6.d0*A_xy(1,1)*A_xy(2,1)*A_xy(3,1)
  DMpole_xy(20,11) = 6.d0*DA_xy(1,1)*A_xy(2,1)*A_xy(3,1)
  DMpole_xy(20,11) = DMpole_xy(20,11) + 6.d0*A_xy(1,1)*DA_xy(2,1)*A_xy(3,1)
  DMpole_xy(20,11) = DMpole_xy(20,11) + 6.d0*A_xy(1,1)*A_xy(2,1)*DA_xy(3,1)
! Mpole_xy(20,12) = 6.d0*A_xy(1,2)*A_xy(2,2)*A_xy(3,2)
  DMpole_xy(20,12) = 6.d0*DA_xy(1,2)*A_xy(2,2)*A_xy(3,2)
  DMpole_xy(20,12) = DMpole_xy(20,12) + 6.d0*A_xy(1,2)*DA_xy(2,2)*A_xy(3,2)
  DMpole_xy(20,12) = DMpole_xy(20,12) + 6.d0*A_xy(1,2)*A_xy(2,2)*DA_xy(3,2)
! Mpole_xy(20,13) = 6.d0*A_xy(1,3)*A_xy(2,3)*A_xy(3,3)
  DMpole_xy(20,13) = 6.d0*DA_xy(1,3)*A_xy(2,3)*A_xy(3,3)
  DMpole_xy(20,13) = DMpole_xy(20,13) + 6.d0*A_xy(1,3)*DA_xy(2,3)*A_xy(3,3)
  DMpole_xy(20,13) = DMpole_xy(20,13) + 6.d0*A_xy(1,3)*A_xy(2,3)*DA_xy(3,3)
  do ind2 = 4,9
    i = oind1(ind2)
    k = oind3(ind2)
!   Mpole_xy(20,10+ind2) = 2.d0*
!          (  A_xy(1,i)*A_xy(2,i)*A_xy(3,k) +
!             A_xy(1,i)*A_xy(2,k)*A_xy(3,i) +
!             A_xy(1,k)*A_xy(2,i)*A_xy(3,i) )
    DMpole_xy(20,10+ind2) = 2.d0*   &
           (  DA_xy(1,i)*A_xy(2,i)*A_xy(3,k) +    &
              DA_xy(1,i)*A_xy(2,k)*A_xy(3,i) +    &
              DA_xy(1,k)*A_xy(2,i)*A_xy(3,i) )
    DMpole_xy(20,10+ind2) = DMpole_xy(20,10+ind2) + 2.d0*   &
           (  A_xy(1,i)*DA_xy(2,i)*A_xy(3,k) +    &
              A_xy(1,i)*DA_xy(2,k)*A_xy(3,i) +    &
              A_xy(1,k)*DA_xy(2,i)*A_xy(3,i) )
    DMpole_xy(20,10+ind2) = DMpole_xy(20,10+ind2) + 2.d0*    &
           (  A_xy(1,i)*A_xy(2,i)*DA_xy(3,k) +    &
              A_xy(1,i)*A_xy(2,k)*DA_xy(3,i) +   &
              A_xy(1,k)*A_xy(2,i)*DA_xy(3,i) )
  enddo
!     Mpole_xy(20,20) =
!             A_xy(1,1)*A_xy(2,2)*A_xy(3,3) +
!             A_xy(1,1)*A_xy(3,2)*A_xy(2,3) +
!             A_xy(2,1)*A_xy(1,2)*A_xy(3,3) +
!             A_xy(2,1)*A_xy(3,2)*A_xy(1,3) +
!             A_xy(3,1)*A_xy(1,2)*A_xy(2,3) +
!             A_xy(3,1)*A_xy(2,2)*A_xy(1,3)
  DMpole_xy(20,20) =    &
              DA_xy(1,1)*A_xy(2,2)*A_xy(3,3) +   &
              DA_xy(1,1)*A_xy(3,2)*A_xy(2,3) +   &
              DA_xy(2,1)*A_xy(1,2)*A_xy(3,3) +   &
              DA_xy(2,1)*A_xy(3,2)*A_xy(1,3) +   &
              DA_xy(3,1)*A_xy(1,2)*A_xy(2,3) +   &
              DA_xy(3,1)*A_xy(2,2)*A_xy(1,3)
  DMpole_xy(20,20) = DMpole_xy(20,20) +    &
              A_xy(1,1)*DA_xy(2,2)*A_xy(3,3) +   &
              A_xy(1,1)*DA_xy(3,2)*A_xy(2,3) +   &
              A_xy(2,1)*DA_xy(1,2)*A_xy(3,3) +   &
              A_xy(2,1)*DA_xy(3,2)*A_xy(1,3) +   &
              A_xy(3,1)*DA_xy(1,2)*A_xy(2,3) +   &
              A_xy(3,1)*DA_xy(2,2)*A_xy(1,3)
  DMpole_xy(20,20) = DMpole_xy(20,20) +   &
              A_xy(1,1)*A_xy(2,2)*DA_xy(3,3) +   &
              A_xy(1,1)*A_xy(3,2)*DA_xy(2,3) +   &
              A_xy(2,1)*A_xy(1,2)*DA_xy(3,3) +   &
              A_xy(2,1)*A_xy(3,2)*DA_xy(1,3) +   &
              A_xy(3,1)*A_xy(1,2)*DA_xy(2,3) +   &
              A_xy(3,1)*A_xy(2,2)*DA_xy(1,3)
  if ( order .eq. 20 )return
  do j = 11,20
    do k = 21,order
      DMpole_xy(k,j) = 0.d0
      DMpole_xy(j,k) = 0.d0
    enddo
  enddo

! HEXADECAPOLES
! H'_mmmm
  do ind1 = 1,3
    m = hind1(ind1)
    do ind2 = 1,15
      i = hind1(ind2)
      j = hind2(ind2)
      k = hind3(ind2)
      l = hind4(ind2)
!         Mpole_xy(ind1+20,ind2+20) =
!              A_xy(m,i)*A_xy(m,j)*A_xy(m,k)*A_xy(m,l)
      DMpole_xy(ind1+20,ind2+20) =   &
               DA_xy(m,i)*A_xy(m,j)*A_xy(m,k)*A_xy(m,l)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +  &
               A_xy(m,i)*DA_xy(m,j)*A_xy(m,k)*A_xy(m,l)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +  &
               A_xy(m,i)*A_xy(m,j)*DA_xy(m,k)*A_xy(m,l)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +  &
               A_xy(m,i)*A_xy(m,j)*A_xy(m,k)*DA_xy(m,l)
    enddo
  enddo
! H'_mmmn
  do ind1 = 4,9
    m = hind1(ind1)
    n = hind4(ind1)
! can put iiii & iiij together
    do ind2 = 1,9
      i = hind1(ind2)
      j = hind4(ind2)
!         Mpole_xy(ind1+20,ind2+20) =
!              A_xy(m,i)*A_xy(m,i)*A_xy(m,i)*A_xy(n,j) + 3.d0*
!              A_xy(m,i)*A_xy(m,i)*A_xy(n,i)*A_xy(m,j)
      DMpole_xy(ind1+20,ind2+20) =   &
               DA_xy(m,i)*A_xy(m,i)*A_xy(m,i)*A_xy(n,j) + 3.d0*   &
               DA_xy(m,i)*A_xy(m,i)*A_xy(n,i)*A_xy(m,j)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +  &
               A_xy(m,i)*DA_xy(m,i)*A_xy(m,i)*A_xy(n,j) + 3.d0*  &
               A_xy(m,i)*DA_xy(m,i)*A_xy(n,i)*A_xy(m,j)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +  &
               A_xy(m,i)*A_xy(m,i)*DA_xy(m,i)*A_xy(n,j) + 3.d0*  &
               A_xy(m,i)*A_xy(m,i)*DA_xy(n,i)*A_xy(m,j)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +  &
               A_xy(m,i)*A_xy(m,i)*A_xy(m,i)*DA_xy(n,j) + 3.d0*  &
               A_xy(m,i)*A_xy(m,i)*A_xy(n,i)*DA_xy(m,j)
    enddo
! can put iijj & iijk together
    do ind2 = 10,15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
!         Mpole_xy(ind1+20,ind2+20) =
!              A_xy(m,i)*A_xy(m,i)*A_xy(m,j)*A_xy(n,k) +
!              A_xy(m,i)*A_xy(m,i)*A_xy(m,k)*A_xy(n,j) + 2.d0*
!              A_xy(m,i)*A_xy(m,j)*A_xy(m,k)*A_xy(n,i)
      DMpole_xy(ind1+20,ind2+20) =   &
               DA_xy(m,i)*A_xy(m,i)*A_xy(m,j)*A_xy(n,k) +  &
               DA_xy(m,i)*A_xy(m,i)*A_xy(m,k)*A_xy(n,j) + 2.d0*  &
               DA_xy(m,i)*A_xy(m,j)*A_xy(m,k)*A_xy(n,i)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +  &
               A_xy(m,i)*DA_xy(m,i)*A_xy(m,j)*A_xy(n,k) +  &
               A_xy(m,i)*DA_xy(m,i)*A_xy(m,k)*A_xy(n,j) + 2.d0*  &
               A_xy(m,i)*DA_xy(m,j)*A_xy(m,k)*A_xy(n,i)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +  &
               A_xy(m,i)*A_xy(m,i)*DA_xy(m,j)*A_xy(n,k) +  &
               A_xy(m,i)*A_xy(m,i)*DA_xy(m,k)*A_xy(n,j) + 2.d0*  &
               A_xy(m,i)*A_xy(m,j)*DA_xy(m,k)*A_xy(n,i)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +  &
               A_xy(m,i)*A_xy(m,i)*A_xy(m,j)*DA_xy(n,k) +  &
               A_xy(m,i)*A_xy(m,i)*A_xy(m,k)*DA_xy(n,j) + 2.d0*  &
               A_xy(m,i)*A_xy(m,j)*A_xy(m,k)*DA_xy(n,i)
    enddo
  enddo

! H'_mmnn
  do ind1 = 10,12
    m = hind1(ind1)
    n = hind3(ind1)
! can put iiii & iiij together
    do ind2 = 1,9
      i = hind1(ind2)
      j = hind4(ind2)
!        Mpole_xy(ind1+20,ind2+20) = 3.d0 * (
!              A_xy(m,i)*A_xy(m,i)*A_xy(n,i)*A_xy(n,j) +
!              A_xy(m,i)*A_xy(m,j)*A_xy(n,i)*A_xy(n,i)  )
      DMpole_xy(ind1+20,ind2+20) = 3.d0 * (   &
               DA_xy(m,i)*A_xy(m,i)*A_xy(n,i)*A_xy(n,j) +   &
               DA_xy(m,i)*A_xy(m,j)*A_xy(n,i)*A_xy(n,i)  )
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +   &
          3.d0 * (                                                &
               A_xy(m,i)*DA_xy(m,i)*A_xy(n,i)*A_xy(n,j) +    &
               A_xy(m,i)*DA_xy(m,j)*A_xy(n,i)*A_xy(n,i)  )
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +   &
          3.d0 * (                                                &
               A_xy(m,i)*A_xy(m,i)*DA_xy(n,i)*A_xy(n,j) +    &
               A_xy(m,i)*A_xy(m,j)*DA_xy(n,i)*A_xy(n,i)  )
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +   &
          3.d0 * (                                                &
               A_xy(m,i)*A_xy(m,i)*A_xy(n,i)*DA_xy(n,j) +    &
               A_xy(m,i)*A_xy(m,j)*A_xy(n,i)*DA_xy(n,i)  )
    enddo
! can put iijj & iijk together
    do ind2 = 10,15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
!         Mpole_xy(ind1+20,ind2+20) =
!              A_xy(m,i)*A_xy(m,i)*A_xy(n,j)*A_xy(n,k) +
!              A_xy(m,j)*A_xy(m,k)*A_xy(n,i)*A_xy(n,i) +  2.d0*
!              A_xy(m,i)*A_xy(m,j)*A_xy(n,i)*A_xy(n,k) +  2.d0*
!              A_xy(m,i)*A_xy(m,k)*A_xy(n,i)*A_xy(n,j)
      DMpole_xy(ind1+20,ind2+20) =    &
               DA_xy(m,i)*A_xy(m,i)*A_xy(n,j)*A_xy(n,k) +    &
               DA_xy(m,j)*A_xy(m,k)*A_xy(n,i)*A_xy(n,i) +  2.d0*    &
               DA_xy(m,i)*A_xy(m,j)*A_xy(n,i)*A_xy(n,k) +  2.d0*    &
               DA_xy(m,i)*A_xy(m,k)*A_xy(n,i)*A_xy(n,j)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +    &
               A_xy(m,i)*DA_xy(m,i)*A_xy(n,j)*A_xy(n,k) +    &
               A_xy(m,j)*DA_xy(m,k)*A_xy(n,i)*A_xy(n,i) +  2.d0*    &
               A_xy(m,i)*DA_xy(m,j)*A_xy(n,i)*A_xy(n,k) +  2.d0*    &
               A_xy(m,i)*DA_xy(m,k)*A_xy(n,i)*A_xy(n,j)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +    &
               A_xy(m,i)*A_xy(m,i)*DA_xy(n,j)*A_xy(n,k) +    &
               A_xy(m,j)*A_xy(m,k)*DA_xy(n,i)*A_xy(n,i) +  2.d0*   &
               A_xy(m,i)*A_xy(m,j)*DA_xy(n,i)*A_xy(n,k) +  2.d0*   &
               A_xy(m,i)*A_xy(m,k)*DA_xy(n,i)*A_xy(n,j)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +    &
               A_xy(m,i)*A_xy(m,i)*A_xy(n,j)*DA_xy(n,k) +    &
               A_xy(m,j)*A_xy(m,k)*A_xy(n,i)*DA_xy(n,i) +  2.d0*   &
               A_xy(m,i)*A_xy(m,j)*A_xy(n,i)*DA_xy(n,k) +  2.d0*   &
               A_xy(m,i)*A_xy(m,k)*A_xy(n,i)*DA_xy(n,j)
    enddo
  enddo

! H'_mmnp
  do ind1 = 13,15
    m = hind1(ind1)
    n = hind3(ind1)
    p = hind4(ind1)
! can put iiii & iiij together
    do ind2 = 1,9
      i = hind1(ind2)
      j = hind4(ind2)
!         Mpole_xy(ind1+20,ind2+20) =
!             3.d0*A_xy(m,i)*A_xy(m,i)*A_xy(n,i)*A_xy(p,j) +
!             3.d0*A_xy(m,i)*A_xy(m,i)*A_xy(n,j)*A_xy(p,i) +
!             6.d0*A_xy(m,i)*A_xy(m,j)*A_xy(n,i)*A_xy(p,i)
      DMpole_xy(ind1+20,ind2+20) =    &
              3.d0*DA_xy(m,i)*A_xy(m,i)*A_xy(n,i)*A_xy(p,j) +    &
              3.d0*DA_xy(m,i)*A_xy(m,i)*A_xy(n,j)*A_xy(p,i) +    &
              6.d0*DA_xy(m,i)*A_xy(m,j)*A_xy(n,i)*A_xy(p,i)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +   &
              3.d0*A_xy(m,i)*DA_xy(m,i)*A_xy(n,i)*A_xy(p,j) +    &
              3.d0*A_xy(m,i)*DA_xy(m,i)*A_xy(n,j)*A_xy(p,i) +    &
              6.d0*A_xy(m,i)*DA_xy(m,j)*A_xy(n,i)*A_xy(p,i)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +   &
              3.d0*A_xy(m,i)*A_xy(m,i)*DA_xy(n,i)*A_xy(p,j) +    &
              3.d0*A_xy(m,i)*A_xy(m,i)*DA_xy(n,j)*A_xy(p,i) +    &
              6.d0*A_xy(m,i)*A_xy(m,j)*DA_xy(n,i)*A_xy(p,i)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +    &
              3.d0*A_xy(m,i)*A_xy(m,i)*A_xy(n,i)*DA_xy(p,j) +    &
              3.d0*A_xy(m,i)*A_xy(m,i)*A_xy(n,j)*DA_xy(p,i) +    &
              6.d0*A_xy(m,i)*A_xy(m,j)*A_xy(n,i)*DA_xy(p,i)
    enddo
! can put iijj & iijk together
    do ind2 = 10,15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
!         Mpole_xy(ind1+20,ind2+20) =
!              A_xy(m,i)*A_xy(m,i)*A_xy(n,j)*A_xy(p,k) +
!              A_xy(m,i)*A_xy(m,i)*A_xy(n,k)*A_xy(p,j) + 2.d0 * (
!              A_xy(m,i)*A_xy(m,j)*A_xy(n,i)*A_xy(p,k) +
!              A_xy(m,i)*A_xy(m,j)*A_xy(n,k)*A_xy(p,i) +
!              A_xy(m,i)*A_xy(m,k)*A_xy(n,i)*A_xy(p,j) +
!              A_xy(m,i)*A_xy(m,k)*A_xy(n,j)*A_xy(p,i) +
!              A_xy(m,j)*A_xy(m,k)*A_xy(n,i)*A_xy(p,i) )
      DMpole_xy(ind1+20,ind2+20) =    &
               DA_xy(m,i)*A_xy(m,i)*A_xy(n,j)*A_xy(p,k) +    &
               DA_xy(m,i)*A_xy(m,i)*A_xy(n,k)*A_xy(p,j) + 2.d0 * (   &
               DA_xy(m,i)*A_xy(m,j)*A_xy(n,i)*A_xy(p,k) +   &
               DA_xy(m,i)*A_xy(m,j)*A_xy(n,k)*A_xy(p,i) +   &
               DA_xy(m,i)*A_xy(m,k)*A_xy(n,i)*A_xy(p,j) +   &
               DA_xy(m,i)*A_xy(m,k)*A_xy(n,j)*A_xy(p,i) +   &
               DA_xy(m,j)*A_xy(m,k)*A_xy(n,i)*A_xy(p,i) )
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +   &
               A_xy(m,i)*DA_xy(m,i)*A_xy(n,j)*A_xy(p,k) +    &
               A_xy(m,i)*DA_xy(m,i)*A_xy(n,k)*A_xy(p,j) + 2.d0 * (    &
               A_xy(m,i)*DA_xy(m,j)*A_xy(n,i)*A_xy(p,k) +    &
               A_xy(m,i)*DA_xy(m,j)*A_xy(n,k)*A_xy(p,i) +    &
               A_xy(m,i)*DA_xy(m,k)*A_xy(n,i)*A_xy(p,j) +    &
               A_xy(m,i)*DA_xy(m,k)*A_xy(n,j)*A_xy(p,i) +    &
               A_xy(m,j)*DA_xy(m,k)*A_xy(n,i)*A_xy(p,i) )
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +    &
               A_xy(m,i)*A_xy(m,i)*DA_xy(n,j)*A_xy(p,k) +    &
               A_xy(m,i)*A_xy(m,i)*DA_xy(n,k)*A_xy(p,j) + 2.d0 * (    &
               A_xy(m,i)*A_xy(m,j)*DA_xy(n,i)*A_xy(p,k) +    &
               A_xy(m,i)*A_xy(m,j)*DA_xy(n,k)*A_xy(p,i) +    &
               A_xy(m,i)*A_xy(m,k)*DA_xy(n,i)*A_xy(p,j) +    &
               A_xy(m,i)*A_xy(m,k)*DA_xy(n,j)*A_xy(p,i) +    &
               A_xy(m,j)*A_xy(m,k)*DA_xy(n,i)*A_xy(p,i) )
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +    &
               A_xy(m,i)*A_xy(m,i)*A_xy(n,j)*DA_xy(p,k) +    &
               A_xy(m,i)*A_xy(m,i)*A_xy(n,k)*DA_xy(p,j) + 2.d0 * (    &
               A_xy(m,i)*A_xy(m,j)*A_xy(n,i)*DA_xy(p,k) +    &
               A_xy(m,i)*A_xy(m,j)*A_xy(n,k)*DA_xy(p,i) +    &
               A_xy(m,i)*A_xy(m,k)*A_xy(n,i)*DA_xy(p,j) +    &
               A_xy(m,i)*A_xy(m,k)*A_xy(n,j)*DA_xy(p,i) +    &
               A_xy(m,j)*A_xy(m,k)*A_xy(n,i)*DA_xy(p,i) )
    enddo
  enddo
  return
end subroutine XFORM_MPOLE_deriv_matrix
!------------------------------------------------------------
subroutine XFORM_MPOLE_field_matrix(A_xy,Field_xy,order)
! A_xy is matrix of coefficients expanding y in terms of x
! i.e. y_i = A_i1 * x_1 + A_i2 * x_2 + A_i3 * x_3
! Field_xy is resulting matrix getting expansion of field and derivs
! with basis y in terms of those with basis x
! order is order of multipoles 1 for charge, 4 for charge-dipoles
! 10 for charge,dipole,quadrupole, 20 for up to octupole and
! finally 35 for up to hexadecapole
 
! Mpole 1 is charge... Mpole 10 is x_3,x_3 quadrupolar coefficients etc.
  implicit none
  integer order
  double precision A_xy(3,3), Field_xy(order,order)
  integer i,j,k,l,m,n,p,q,ind1,ind2,jj,kk
  integer qind1(6),qind2(6),oind1(10),oind2(10),oind3(10)
  integer hind1(15),hind2(15),hind3(15),hind4(15)
  common/MxfInd/qind1,qind2,oind1,oind2,oind3,hind1,hind2,hind3,hind4
      

! CHARGE case
  Field_xy(1,1) = 1.d0
  if ( order .eq. 1 )return
  do j = 2,order
    Field_xy(1,j) = 0.d0
    Field_xy(j,1) = 0.d0
  enddo
! DIPOLES
! d_yk = A_xy(k,1)*d_x1 + A_xy(k,2)*d_x2 + A_xy(k,3)*d_x3
! D'_j
  do j = 2,4
    do k = 2,4
      Field_xy(j,k) = A_xy(j-1,k-1)
    enddo
  enddo
  if ( order .eq. 4 )return
  do j = 2,4
    do k = 5,order
      Field_xy(j,k) = 0.d0
      Field_xy(k,j) = 0.d0
    enddo
  enddo
! QUADRUPOLES
! q_ykyl = sum_i,j (A_xy(k,i)*A_xy(l,j)+A_xy(k,j)*A_xy(l,i))*(q_xixj+q_xjxi)
! Mp(5) = q_y1y1,..,Mp(7) = q_y3y3,Mp(8) = q_y1y2+q_y2y1,.,Mp(10)=q_y2y3+q_y3y2
! Q'_kk
  do ind1 = 1,3
    l = qind1(ind1)
    do ind2 = 1,3
      i = qind1(ind2)
      Field_xy(ind1+4,ind2+4) = A_xy(l,i)*A_xy(l,i)
    enddo
    do ind2 = 4,6
      i = qind1(ind2)
      j = qind2(ind2)
      Field_xy(ind1+4,ind2+4) = 2.d0*A_xy(l,i)*A_xy(l,j)
    enddo
  enddo
  do ind1 = 4,6
    l = qind1(ind1)
    m = qind2(ind1)
    do ind2 = 1,3
      i = qind1(ind2)
      Field_xy(ind1+4,ind2+4) = A_xy(l,i)*A_xy(m,i)
    enddo
    do ind2 = 4,6
      i = qind1(ind2)
      j = qind2(ind2)
      Field_xy(ind1+4,ind2+4) = (A_xy(l,i)*A_xy(m,j)+A_xy(m,i)*A_xy(l,j))
    enddo
  enddo
  if ( order .eq. 10 )return
  do j = 5,10
    do k = 11,order
      Field_xy(k,j) = 0.d0
      Field_xy(j,k) = 0.d0
    enddo
  enddo
! OCTUPOLES level field
! F'_lll
  do ind1 = 1,3
    l = oind1(ind1)
    do ind2 = 1,3
      i = oind1(ind2)
      Field_xy(ind1+10,ind2+10) = A_xy(l,i)**3
    enddo
    do ind2 = 4,9
      i = oind1(ind2)
      j = oind3(ind2)
      Field_xy(ind1+10,ind2+10) = 3.d0*A_xy(l,i)**2 * A_xy(l,j)
    enddo
    Field_xy(ind1+10,20) = 6.d0*A_xy(l,1)*A_xy(l,2)*A_xy(l,3) 
  enddo
! F'_llm
  do ind1 = 4,9
    l = oind1(ind1)
    m = oind3(ind1)
    do ind2 = 1,3
      i = oind1(ind2)
      Field_xy(ind1+10,ind2+10) = A_xy(l,i)**2 * A_xy(m,i)
    enddo
    do ind2 = 4,9
      i = oind1(ind2)
      j = oind3(ind2)
      Field_xy(ind1+10,ind2+10) = A_xy(l,i)**2 * A_xy(m,j) +   &
                                2.d0*A_xy(l,i)*A_xy(l,j)*A_xy(m,i)
    enddo
    Field_xy(ind1+10,20) = 2.d0*(A_xy(l,1)*A_xy(l,2)*A_xy(m,3)+   &
                                     A_xy(l,1)*A_xy(m,2)*A_xy(l,3)+   &
                                     A_xy(m,1)*A_xy(l,2)*A_xy(l,3) )
  enddo
! F'_123
  do ind2 = 1,3
    i = oind1(ind2)
    Field_xy(20,ind2+10) = A_xy(1,i)*A_xy(2,i)*A_xy(3,i)
  enddo
  do ind2 = 4,9
    i = oind1(ind2)
    j = oind3(ind2)
    Field_xy(20,ind2+10) = A_xy(1,i)*A_xy(2,i)*A_xy(3,j) +  &
                              A_xy(1,i)*A_xy(2,j)*A_xy(3,i) +  &
                              A_xy(1,j)*A_xy(2,i)*A_xy(3,i) 
  enddo
  Field_xy(20,20) = A_xy(1,1)*A_xy(2,2)*A_xy(3,3) +   &
                        A_xy(1,1)*A_xy(2,3)*A_xy(3,2) +   &
                        A_xy(1,2)*A_xy(2,1)*A_xy(3,3) +   &
                        A_xy(1,2)*A_xy(2,3)*A_xy(3,1) +   &
                        A_xy(1,3)*A_xy(2,1)*A_xy(3,2) +   &
                        A_xy(1,3)*A_xy(2,2)*A_xy(3,1)
  if ( order .eq. 20 )return
  do j = 11,20
    do k = 21,order
      Field_xy(k,j) = 0.d0
      Field_xy(j,k) = 0.d0
    enddo
  enddo
! HEXADECAPOLES
! F'_mmmm
  do ind1 = 1,3
    m = hind1(ind1)
    do ind2 = 1,3
      i = hind1(ind2)
      Field_xy(ind1+20,ind2+20) = A_xy(m,i)**4
    enddo
    do ind2 = 4,9
      i = hind1(ind2)
      j = hind4(ind2)
      Field_xy(ind1+20,ind2+20) = 4.d0 * A_xy(m,i)**3 * A_xy(m,j)
    enddo
    do ind2 = 10,12
      i = hind1(ind2)
      j = hind4(ind2)
      Field_xy(ind1+20,ind2+20) = 6.d0 * A_xy(m,i)**2 * A_xy(m,j)**2
    enddo
    do ind2 = 13,15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
      Field_xy(ind1+20,ind2+20) = 12.d0*A_xy(m,i)**2 * A_xy(m,j)*A_xy(m,k)
    enddo
  enddo
  do ind1 = 4,9
    m = hind1(ind1)
    n = hind4(ind1)
    do ind2 = 1,3
      i = hind1(ind2)
      Field_xy(ind1+20,ind2+20) = A_xy(m,i)**3 * A_xy(n,i)
    enddo
    do ind2 = 4,9
      i = hind1(ind2)
      j = hind4(ind2)
      Field_xy(ind1+20,ind2+20) = A_xy(m,i)**3 * A_xy(n,j) +  &
                               3.d0*A_xy(m,i)**2 * A_xy(m,j)*A_xy(n,i)
    enddo
    do ind2 = 10,12
      i = hind1(ind2)
      j = hind4(ind2)
      Field_xy(ind1+20,ind2+20) =    &
              3.d0 * A_xy(m,i)**2 * A_xy(m,j) * A_xy(n,j) +    &
              3.d0 * A_xy(m,j)**2 * A_xy(m,i) * A_xy(n,i) 
    enddo
    do ind2 = 13,15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
      Field_xy(ind1+20,ind2+20) =    &
              6.d0 * A_xy(m,i) * A_xy(m,j) * A_xy(m,k) * A_xy(n,i) +  &
              3.d0 * A_xy(m,i)**2 * A_xy(m,j) * A_xy(n,k) +  &
              3.d0 * A_xy(m,i)**2 * A_xy(m,k) * A_xy(n,j)
    enddo
  enddo
  do ind1 = 10,12
    m = hind1(ind1)
    n = hind4(ind1)
    do ind2 = 1,3
      i = hind1(ind2)
      Field_xy(ind1+20,ind2+20) = A_xy(m,i)**2 * A_xy(n,i)**2
    enddo
    do ind2 = 4,9
      i = hind1(ind2)
      j = hind4(ind2)
      Field_xy(ind1+20,ind2+20) =    &
                 2.d0 * A_xy(m,i)**2 * A_xy(n,i) * A_xy(n,j) +   &
                 2.d0 * A_xy(n,i)**2 * A_xy(m,i) * A_xy(m,j)
    enddo
    do ind2 = 10,12
      i = hind1(ind2)
      j = hind4(ind2)
      Field_xy(ind1+20,ind2+20) =    &
                       A_xy(m,i)**2 * A_xy(n,j)**2 +   &
                       A_xy(m,j)**2 * A_xy(n,i)**2 +   &
               4.d0 *  A_xy(m,i) * A_xy(m,j) *A_xy(n,i) * A_xy(n,j)
    enddo
    do ind2 = 13,15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
      Field_xy(ind1+20,ind2+20) =   &
                 2.d0 * A_xy(m,i)**2 * A_xy(n,j) * A_xy(n,k) +  &
                 2.d0 * A_xy(n,i)**2 * A_xy(m,j) * A_xy(m,k) +  &
           4.d0 * A_xy(m,i) * A_xy(m,j) * A_xy(n,i) * A_xy(n,k) +  &
           4.d0 * A_xy(m,i) * A_xy(m,k) * A_xy(n,i) * A_xy(n,j) 
    enddo
  enddo
  do ind1 = 13,15
    m = hind1(ind1)
    n = hind3(ind1)
    p = hind4(ind1)
    do ind2 = 1,3
      i = hind1(ind2)
      Field_xy(ind1+20,ind2+20) = A_xy(m,i)**2 * A_xy(n,i) * A_xy(p,i)
    enddo
    do ind2 = 4,9
      i = hind1(ind2)
      j = hind4(ind2)
      Field_xy(ind1+20,ind2+20) =    &
            2.d0 * A_xy(m,i) * A_xy(m,j) * A_xy(n,i) * A_xy(p,i) +  &
            A_xy(m,i)**2 * (A_xy(n,i)*A_xy(p,j) + A_xy(n,j)*A_xy(p,i))
    enddo
    do ind2 = 10,12
      i = hind1(ind2)
      j = hind4(ind2)
      Field_xy(ind1+20,ind2+20) =    &
               A_xy(m,i)**2 * A_xy(n,j) * A_xy(p,j) +   &
               A_xy(m,j)**2 * A_xy(n,i) * A_xy(p,i) +   &
               2.d0 * A_xy(m,i) * A_xy(m,j) *   &
                   (A_xy(n,i)*A_xy(p,j) + A_xy(n,j)*A_xy(p,i))
    enddo
    do ind2 = 13,15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
      Field_xy(ind1+20,ind2+20) =   &
           A_xy(m,i)**2 * (A_xy(n,j)*A_xy(p,k) + A_xy(n,k)*A_xy(p,j)) +   &
                      2.d0 * A_xy(m,i)*A_xy(m,j) *    &
                      (A_xy(n,i)*A_xy(p,k) + A_xy(n,k)*A_xy(p,i)) +   &
                      2.d0 * A_xy(m,i)*A_xy(m,k) *    &
                      (A_xy(n,i)*A_xy(p,j) + A_xy(n,j)*A_xy(p,i)) +   &
                   2.d0 * A_xy(m,j) * A_xy(m,k) * A_xy(n,i) * A_xy(p,i)
    enddo
  enddo

  return
end subroutine XFORM_MPOLE_field_matrix
!------------------------------------------------------------
subroutine XFORM_MPOLE(Mpole_xy,dimxy,Mpole_in,Mpole_out,order)
  implicit none
  integer order,dimxy
  double precision Mpole_xy(dimxy,dimxy)
  double precision Mpole_in(*),Mpole_out(*)

  integer i,j

  if ( order .eq. 0 )return
  Mpole_out(1) = Mpole_xy(1,1)*Mpole_in(1)
  if ( order .eq. 1 )return
! DIPOLES
  if(order == 3)then
     do i = 1,3
       Mpole_out(i) = 0.d0
       do j = 1,3
         Mpole_out(i) = Mpole_out(i) + Mpole_xy(i+1,j+1)*Mpole_in(j)
       enddo
     enddo
     return
  endif
  do i = 2,4
    Mpole_out(i) = 0.d0
    do j = 2,4
      Mpole_out(i) = Mpole_out(i) + Mpole_xy(i,j)*Mpole_in(j)
    enddo
  enddo
  if ( order .eq. 4 )return
! QUADRUPOLES
  do i = 5,10
    Mpole_out(i) = 0.d0
    do j = 5,10
      Mpole_out(i) = Mpole_out(i) + Mpole_xy(i,j)*Mpole_in(j)
    enddo
  enddo
  if ( order .eq. 10 )return
! OCTUPOLES
  do i = 11,20
    Mpole_out(i) = 0.d0
    do j = 11,20
      Mpole_out(i) = Mpole_out(i) + Mpole_xy(i,j)*Mpole_in(j)
    enddo
  enddo
  if ( order .eq. 20 )return
! HEXADECAPOLES
  do i = 21,35
    Mpole_out(i) = 0.d0
    do j = 21,35
      Mpole_out(i) = Mpole_out(i) + Mpole_xy(i,j)*Mpole_in(j)
    enddo
  enddo
  if ( order .eq. 35 )return
  return
end subroutine XFORM_MPOLE
!------------------------------------------------------------
subroutine XFORM_FIELD(Field_xy,dimxy,Field_in,Field_out,order)
  implicit none
  integer order,dimxy
  double precision Field_xy(dimxy,dimxy)
  double precision Field_in(*),Field_out(*)

  integer i,j
  !Field_out is incremented by xform of Field_in

  if ( order .eq. 0 )return
  Field_out(1) = Field_out(1) + Field_xy(1,1)*Field_in(1)
  if ( order .eq. 1 )return
! DIPOLES
  do i = 2,4
    do j = 2,4
      Field_out(i) = Field_out(i) + Field_xy(i,j)*Field_in(j)
    enddo
  enddo
  if ( order .eq. 4 )return
! QUADRUPOLES
  do i = 5,10
    do j = 5,10
      Field_out(i) = Field_out(i) + Field_xy(i,j)*Field_in(j)
    enddo
  enddo
  if ( order .eq. 10 )return
! OCTUPOLES
  do i = 11,20
    do j = 11,20
      Field_out(i) = Field_out(i) + Field_xy(i,j)*Field_in(j)
    enddo
  enddo
  if ( order .eq. 20 )return
! HEXADECAPOLES
  do i = 21,35
    do j = 21,35
      Field_out(i) = Field_out(i) + Field_xy(i,j)*Field_in(j)
    enddo
  enddo
  if ( order .eq. 35 )return
  return
end subroutine XFORM_FIELD
!------------------------------------------------------------
