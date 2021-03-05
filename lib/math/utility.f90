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
!---------------------------------------------------------------
subroutine UTIL_copy_int_array(x,y,num)
  implicit none
  integer x(*),y(*)
  integer num
  integer n
  do n = 1,num
    y(n) = x(n)
  enddo
  return
end subroutine UTIL_copy_int_array
!----------------------------------------------------------
subroutine UTIL_copy_real_array(x,y,num)
  implicit none
  double precision x(*),y(*)
  integer num
  integer n
  do n = 1,num
    y(n) = x(n)
  enddo
  return
end subroutine UTIL_copy_real_array
!----------------------------------------------------------
subroutine UTIL_add_array(x,y,z,num)
  implicit none
  double precision x(*),y(*),z(*)
  integer num
  integer n
  do n = 1,num
    z(n) = x(n) + y(n)
  enddo
  return
end subroutine UTIL_add_array
!----------------------------------------------------------
subroutine UTIL_increment_real_array(x,y,num)
  implicit none
  double precision x(*),y(*)
  integer num
  integer n
  do n = 1,num
    x(n) = x(n) + y(n)
  enddo
  return
end subroutine UTIL_increment_real_array
!----------------------------------------------------------
subroutine UTIL_zero_real_array(array,num)
  implicit none
  double precision array(*)
  integer num

  integer i
  do i = 1,num
    array(i) = 0.d0
  end do
  return
end subroutine UTIL_zero_real_array
!----------------------------------------------------------
subroutine UTIL_zero_int_array(array,num)
  implicit none
  integer array(*)
  integer num

  integer i
  do i = 1,num
    array(i) = 0
  end do
  return
end subroutine UTIL_zero_int_array
!----------------------------------------------------------
integer function UTIL_div_pow_k(n,k)
  implicit none
  integer n,k
! divide out highest power of k in n and return ratio
! e.g. if n is 24 and k is 2 result will be 3
  integer remainder,div
  remainder = 0
  div = 1
  do while(remainder == 0 )
    div = k*div
    remainder = n - div*(n/div)
  enddo
  div = div / k
  UTIL_div_pow_k = n/div
end function UTIL_div_pow_k
!----------------------------------------------------------------
subroutine UTIL_search_real_array(val,array,num,ind)
  implicit none
  integer num,ind
  double precision val,array(num)

  integer indhi,indlo
! bisection search for ind just before val (i.e. array(ind)<val<array(ind+1) )
! assume array is ordered in increasing order
  indlo = 1 
  indhi = num
  do while (indhi - indlo > 1 )
    ind = (indhi+indlo)/2
    if ( array(ind) > val )then
      indhi = ind
    else
      indlo = ind
    endif
  enddo
  ind = indlo
  return
end subroutine UTIL_search_real_array
!---------------------------------------------------------
