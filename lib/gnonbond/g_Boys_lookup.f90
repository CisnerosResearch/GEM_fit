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
! GAC
subroutine G_boys_0to12(topdeg,arg,boys_table)
  implicit none
  integer topdeg,tabarg
  double precision arg,argstar,darg,tmp1,tmp2,boys_table(0:topdeg)
  double precision fac(0:6), facinv(0:6), signum
  integer n
  include "boys.fh"

  if (topdeg .gt. 13) then
     write(6,*)'G_boys_0to12: incorrect degree for Boys function lookup'
     stop
  endif

  if (arg .gt. 12.0d0 .or. arg .lt. 0.0d0) then
     write(6,*)'G_boys_0to12: incorrect argument for Boys function lookup',&
                & arg
     stop
  endif

  fac(0) = 1.d0
  facinv(0) = 1.d0
  do n=1,6
     fac(n) = n*fac(n-1)
     facinv(n) = 1.d0/fac(n)
  enddo

  tmp1 = 10.d0*arg
  tabarg = int(tmp1)
  !argstar = arg
  darg = 0.10d0*tabarg -arg

! Taylor expansion
  boys_table(topdeg) = B_tab(topdeg,tabarg) + &
                       darg*(B_tab(topdeg+1,tabarg)*facinv(1) + &
                       darg*(B_tab(topdeg+2,tabarg)*facinv(2) + &
                       darg*(B_tab(topdeg+3,tabarg)*facinv(3) + &
                       darg*(B_tab(topdeg+4,tabarg)*facinv(4) + &
                       darg*(B_tab(topdeg+5,tabarg)*facinv(5) + &
                       darg*(B_tab(topdeg+6,tabarg)*facinv(6)))))))
  !do n = 0,6
  !   signum = 1.0d0
  !   boys_table(topdeg) = boys_table(topdeg) + &
  !                       (signum*B_tab(topdeg+n,tabarg)*(darg)**n) / &
  !                       fac(n)
  !enddo

! downward recursion
  do n = topdeg-1, 0, -1
     boys_table(n) = (2.d0*arg*boys_table(n+1) + exp(-arg)) / &
                     (2.d0*n + 1)
  enddo

  return

end subroutine G_boys_0to12
!------------------------------------------------------------------
subroutine G_boys_12to15(topdeg,arg,boys_table)
  implicit none
  integer topdeg
  double precision arg,g,pi,boys_table(0:topdeg)
  integer n

  pi = 4.d0*(atan(1.d0))

! approximate complementary integral
  g = 0.4999489092 - 0.2473631686d0/arg + 0.321180909/(arg**2) - &
      0.3811559346/(arg**3)
  boys_table(0) = sqrt(pi)/(2.d0*sqrt(arg)) - (g*exp(-arg))/arg

! upward recursion
  do n = 1, topdeg
     boys_table(n) = (1/(2.d0*arg))*((2.d0*(n-1)+1)*boys_table(n-1)-exp(-arg))
  enddo

end subroutine G_boys_12to15
!------------------------------------------------------------------
subroutine G_boys_15to18(topdeg,arg,boys_table)
  implicit none
  integer topdeg
  double precision arg,g,pi,boys_table(0:topdeg)
  integer n

  pi = 4.d0*(atan(1.d0))

! approximate complementary integral
  g = 0.4998436875 - 0.24249438/arg + 0.24642845/(arg**2) 
  boys_table(0) = sqrt(pi)/(2.d0*sqrt(arg)) - (g*exp(-arg))/arg

! upward recursion
  do n = 1, topdeg
     boys_table(n) = (1/(2.d0*arg))*((2.d0*(n-1)+1)*boys_table(n-1)-exp(-arg))
  enddo

end subroutine G_boys_15to18
!------------------------------------------------------------------
subroutine G_boys_18to24(topdeg,arg,boys_table)
  implicit none
  integer topdeg
  double precision arg,g,pi,boys_table(0:topdeg)
  integer n

  pi = 4.d0*(atan(1.d0))

! approximate complementary integral
  g = 0.499093162 - 0.2152832/arg
  boys_table(0) = sqrt(pi)/(2.d0*sqrt(arg)) - (g*exp(-arg))/arg

! upward recursion
  do n = 1, topdeg
     boys_table(n) = (1/(2.d0*arg))*((2.d0*(n-1)+1)*boys_table(n-1)-exp(-arg))
  enddo

end subroutine G_boys_18to24
!------------------------------------------------------------------
subroutine G_boys_24to30(topdeg,arg,boys_table)
  implicit none
  integer topdeg
  double precision arg,g,pi,boys_table(0:topdeg)
  integer n

  pi = 4.d0*(atan(1.d0))

! approximate complementary integral
  g = 0.490
  boys_table(0) = sqrt(pi)/(2.d0*sqrt(arg)) - (g*exp(-arg))/arg

! upward recursion
  do n = 1, topdeg
     boys_table(n) = (1/(2.d0*arg))*((2.d0*(n-1)+1)*boys_table(n-1)-exp(-arg))
  enddo

end subroutine G_boys_24to30
!------------------------------------------------------------------
subroutine G_boys_30up(topdeg,arg,boys_table)
  implicit none
  integer topdeg
  double precision arg,g,pi,boys_table(0:topdeg)
  integer n

  pi = 4.d0*(atan(1.d0))

! approximate complementary integral
  boys_table(0) = sqrt(pi)/(2.d0*sqrt(arg))

! upward recursion
  do n = 1, topdeg
     boys_table(n) = (1/(2.d0*arg))*((2.d0*(n-1)+1)*boys_table(n-1)-exp(-arg))
  enddo

end subroutine G_boys_30up
!------------------------------------------------------------------
