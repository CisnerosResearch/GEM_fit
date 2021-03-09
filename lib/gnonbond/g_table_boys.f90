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
!------------------------------------------------------------------
subroutine GNTABLE_exact(topdeg,arg,boys_table)
  implicit none
  integer topdeg
  double precision arg,boys_table(0:topdeg)
  integer n
  double precision value

  do n = 0,topdeg
    call GTABLE_F(n,arg,value)
    boys_table(n) = value
  enddo
  return
end subroutine GNTABLE_exact
!------------------------------------------------------------------
subroutine GTABLE_exact(mlist,topdeg,scratch,boys_table)
  implicit none
  integer mlist,topdeg
  double precision scratch(mlist),boys_table(mlist,0:topdeg)
  integer m,n
  double precision x,value

  do n = 0,topdeg
    do m = 1,mlist
      x = scratch(m)
      call GTABLE_F(n,x,value)
      boys_table(m,n) = value
    enddo
  enddo
  return
end subroutine GTABLE_exact
!---------------------------------------------
subroutine GTABLE_F(m,x,value)
  implicit none
  integer m
  double precision x,value
 include "math.fh"

  integer MAX
  parameter(MAX=1000)
  double precision chk,sum,summand(MAX),fac,eps,temp
  integer i,j,k,itop

  fac = 0.5d0*sqrt(M_PI)
  ! below series expansion for boys function obtained by
  ! iterating the recursion for F(j,t) in terms of  F(j+1,t)
  if ( x > 30.d0 )then
    temp = 0.5d0*sqrt(M_PI/x) 
    if ( m > 0 )then
      do j = 1,m
        temp = (2.d0*j-1.d0)*temp/(2.d0*x)
      enddo
    endif
    value = temp
    return
  endif
  eps = 1.d-20
  i = 0
  k = m + i
  fac = 2.d0*x
  chk = 1.d0 / (2.d0*k+1)
  do while ( chk .gt. eps )
    i = i + 1
    if ( i .gt. MAX )then
      write(6,*)'too many iterations in GTABLE_F: m,x,chk = ',m,x,chk
      stop
    endif
    k = m + i
    chk = chk*fac/(2.d0*k+1.d0)
    summand(i) = chk
  enddo
  itop = i
  sum = 0.d0
  ! add up in reverse order
  do k = 1,itop
    i = itop + 1 - k
    sum = sum + summand(i)
  enddo
  ! add in the initial value
  sum = sum + 1.d0 / (2.d0*m +1.d0)
  value = exp(-x)*sum

  return
end subroutine GTABLE_F
!------------------------------------------------
