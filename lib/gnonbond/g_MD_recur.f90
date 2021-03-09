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
!-------------------------------------------------------
! GAC: modified Tom's _dd_ subroutine
!subroutine GN_MD_REC(deg1,deg2,boyspar,dx,dy,dz,result,x1,y1,z1,x2,y2,z2,e1,e2)
subroutine GN_MD_REC(deg1,deg2,boyspar,dx,dy,dz,result,x1,y1,z1,x2,y2,z2)
  implicit none
  integer deg1, deg2, x1, y1, z1, x2, y2, z2, allochk! GAC
  double precision boyspar,dx,dy,dz,result

  integer topdeg,t,u,v,n,i,j,counter
  double precision dr2,pi,arg,prefac
  double precision,allocatable::R(:,:,:,:),boys_table(:)
  integer t_part(56),u_part(56),v_part(56)
  double precision signum, factor1,factor2,epsilon!,e1,e2
  parameter(epsilon = 1.0d-15)

  dr2 = dx*dx+dy*dy+dz*dz
  arg = boyspar*dr2
  if (abs(arg) .lt. epsilon) arg = 0.0d0
  topdeg = deg1+deg2
  if (topdeg .gt. 12) then
     write(6,*)'g_MD_recur: can only handle up to f functions, exiting'
     stop
  endif
  allocate(boys_table(0:topdeg), R(0:topdeg,0:topdeg,0:topdeg,0:topdeg),&
           stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'GN_MD_REC:could not allocate boys_table, R; exiting'
     stop
  endif
  boys_table(:) = 0.0d0
  if (arg .le. 12.0d0) then
     call G_boys_0to12(topdeg,arg,boys_table)
  else if (arg .le. 15.0d0) then
     call G_boys_12to15(topdeg,arg,boys_table)
  else if (arg .le. 18.0d0) then
     call G_boys_15to18(topdeg,arg,boys_table)
  else if (arg .le. 24.0d0) then
     call G_boys_18to24(topdeg,arg,boys_table)
  else if (arg .le. 30.0d0) then
     call G_boys_24to30(topdeg,arg,boys_table)
  else if (arg .gt. 30.0d0) then
     call G_boys_30up(topdeg,arg,boys_table)
  endif
  !print *,'arg = ',arg,boys_table(topdeg)
  !initialize---multiply by factor to the power of boys func degree
  ! also times the factor for electrostatics of unit charge distributions
  ! this gives the value of R^n(0,0,0) or R(0,0,0,n)
  pi = 3.14159265358979323846d0
  prefac = 2.d0*sqrt(boyspar/pi)!this is the factor for electrostatics
                                !involving unit charge distributions
  !factor1 = (pi/e1)*sqrt(pi/e1)
  !factor2 = (pi/e2)*sqrt(pi/e2)
  do n = 0,topdeg
    R(0,0,0,n) = prefac*boys_table(n)
    prefac = prefac * (-2.d0*boyspar)
  enddo
  ! next the v,z direction
  v = 1
  do n = 0,topdeg - v
    R(0,0,v,n) = dz*R(0,0,v-1,n+1)
  enddo
  do v = 2,topdeg
    do n = 0,topdeg - v 
      R(0,0,v,n) = (v-1)*R(0,0,v-2,n+1) + dz*R(0,0,v-1,n+1)
    enddo
  enddo
  ! next the u,y direction
  u = 1
  do v = 0,topdeg - u
    do n = 0,topdeg - (u+v)
      R(0,u,v,n) = dy*R(0,u-1,v,n+1)
    enddo
  enddo
  do u = 2,topdeg
    do v = 0,topdeg - u
      do n = 0,topdeg - (u+v)
        R(0,u,v,n) = (u-1)*R(0,u-2,v,n+1) + dy*R(0,u-1,v,n+1)
      enddo
    enddo
  enddo
  ! next the t,x direction
  t = 1
  do u = 0,topdeg - t
    do v = 0,topdeg - (t+u)
      do n = 0,topdeg - (t+u+v)
        R(t,u,v,n) = dx*R(t-1,u,v,n+1)
      enddo
    enddo
  enddo
  do t = 2,topdeg
    do u = 0,topdeg - t
      do v = 0,topdeg - (t+u)
        do n = 0,topdeg - (t+u+v)
          R(t,u,v,n) = (t-1)*R(t-2,u,v,n+1) + dx*R(t-1,u,v,n+1)
        enddo
      enddo
    enddo
  enddo

  if(mod((x2+y2+z2),2).eq.1)then
    signum = -1.0d0
  else
    signum = 1.0d0
  endif
  result = signum*R(x1+x2,y1+y2,z1+z2,0)
  !print *,'result = ',result
  deallocate(boys_table,R)
  return
end subroutine GN_MD_REC
!-------------------------------------------------------
! GAC: modified Tom's _dd_ subroutine
subroutine GN_MD_REC_erfc(deg1,deg2,boyspar,dx,dy,dz,result,x1,y1,z1,x2,&
                          y2,z2,e1,e2,boyspar2,beta)
  implicit none
  integer deg1, deg2, x1, y1, z1, x2, y2, z2 ! GAC
  double precision boyspar,dx,dy,dz,result,boyspar2,arg2

  integer topdeg,t,u,v,n,i,j,counter,allochk
  double precision dr2,pi,arg,prefac
  double precision,allocatable::R(:,:,:,:),R1(:,:,:,:),R2(:,:,:,:),&
                                boys_table(:),boys_table2(:)
  integer t_part(56),u_part(56),v_part(56)
  double precision signum, factor1,factor2,e1,e2,beta,prefac2,epsilon
  parameter(epsilon = 1.0d-15) 

  dr2 = dx*dx+dy*dy+dz*dz
  arg = boyspar*dr2
  if (abs(arg) .lt. epsilon) arg = 0.0d0
  arg2 = boyspar2*dr2
  if (abs(arg2) .lt. epsilon) arg2 = 0.0d0
  topdeg = deg1+deg2
  if (topdeg .gt. 12) then
     write(6,*)'g_MD_recur: can only handle up to f functions, exiting'
     stop
  endif
  allocate(boys_table(0:topdeg), R(0:topdeg,0:topdeg,0:topdeg,0:topdeg),&
           boys_table2(0:topdeg), R1(0:topdeg,0:topdeg,0:topdeg,0:topdeg),&
           R2(0:topdeg,0:topdeg,0:topdeg,0:topdeg),&
           stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'GN_MD_REC_erfc:could not allocate boys_table, R; exiting'
     stop
  endif
  boys_table(:) = 0.0d0
  if (arg .le. 12.0d0) then
     call G_boys_0to12(topdeg,arg,boys_table)
  else if (arg .le. 15.0d0) then
     call G_boys_12to15(topdeg,arg,boys_table)
  else if (arg .le. 18.0d0) then
     call G_boys_15to18(topdeg,arg,boys_table)
  else if (arg .le. 24.0d0) then
     call G_boys_18to24(topdeg,arg,boys_table)
  else if (arg .le. 30.0d0) then
     call G_boys_24to30(topdeg,arg,boys_table)
  else if (arg .gt. 30.0d0) then
     call G_boys_30up(topdeg,arg,boys_table)
  endif
  boys_table2(:) = 0.0d0
  if (arg2 .le. 12.0d0) then
     call G_boys_0to12(topdeg,arg2,boys_table2)
  else if (arg2 .le. 15.0d0) then
     call G_boys_12to15(topdeg,arg2,boys_table2)
  else if (arg2 .le. 18.0d0) then
     call G_boys_15to18(topdeg,arg2,boys_table2)
  else if (arg2 .le. 24.0d0) then
     call G_boys_18to24(topdeg,arg2,boys_table2)
  else if (arg2 .le. 30.0d0) then
     call G_boys_24to30(topdeg,arg2,boys_table2)
  else if (arg2 .gt. 30.0d0) then
     call G_boys_30up(topdeg,arg2,boys_table2)
  endif
  !print *,'arg = ',arg,boys_table(topdeg)
  !initialize---multiply by factor to the power of boys func degree
  ! also times the factor for electrostatics of unit charge distributions
  ! this gives the value of R^n(0,0,0) or R(0,0,0,n)
  pi = 3.14159265358979323846d0
  prefac = 2.d0*sqrt(boyspar/pi)!this is the factor for electrostatics
                                !involving unit charge distributions
  prefac2 = sqrt(2.0d0*boyspar2/pi)      !this is the factor for 2nd integral in erfc
  !factor1 = (pi/e1)*sqrt(pi/e1)
  !factor2 = (pi/e2)*sqrt(pi/e2)
  do n = 0,topdeg
    R1(0,0,0,n) = prefac*boys_table(n)
    R2(0,0,0,n) = prefac2*boys_table2(n)
    R(0,0,0,n) = R1(0,0,0,n) - R2(0,0,0,n)
    prefac = prefac * (-2.d0*boyspar)
    prefac2 = prefac2 * (-2.0d0*boyspar2)
  enddo
  ! next the v,z direction
  v = 1
  do n = 0,topdeg - v
    R(0,0,v,n) = dz*R(0,0,v-1,n+1)
  enddo
  do v = 2,topdeg
    do n = 0,topdeg - v 
      R(0,0,v,n) = (v-1)*R(0,0,v-2,n+1) + dz*R(0,0,v-1,n+1)
    enddo
  enddo
  ! next the u,y direction
  u = 1
  do v = 0,topdeg - u
    do n = 0,topdeg - (u+v)
      R(0,u,v,n) = dy*R(0,u-1,v,n+1)
    enddo
  enddo
  do u = 2,topdeg
    do v = 0,topdeg - u
      do n = 0,topdeg - (u+v)
        R(0,u,v,n) = (u-1)*R(0,u-2,v,n+1) + dy*R(0,u-1,v,n+1)
      enddo
    enddo
  enddo
  ! next the t,x direction
  t = 1
  do u = 0,topdeg - t
    do v = 0,topdeg - (t+u)
      do n = 0,topdeg - (t+u+v)
        R(t,u,v,n) = dx*R(t-1,u,v,n+1)
      enddo
    enddo
  enddo
  do t = 2,topdeg
    do u = 0,topdeg - t
      do v = 0,topdeg - (t+u)
        do n = 0,topdeg - (t+u+v)
          R(t,u,v,n) = (t-1)*R(t-2,u,v,n+1) + dx*R(t-1,u,v,n+1)
        enddo
      enddo
    enddo
  enddo

  if(mod((x2+y2+z2),2).eq.1)then
    signum = -1.d0
  else
    signum = 1.d0
  endif
  result = signum*R(x1+x2,y1+y2,z1+z2,0)
  !print *,'result = ',result
  deallocate(boys_table,boys_table2,R,R1,R2)
  return
end subroutine GN_MD_REC_erfc
!-------------------------------------------------------
! GAC: modified Tom's _dd_ subroutine for overlap
subroutine GN_MD_REC_OVERLAP(deg1,deg2,boyspar,dx,dy,dz,result,x1,y1,z1,&
                             x2,y2,z2,e1,e2)
  implicit none
  integer deg1, deg2, x1, y1, z1, x2, y2, z2 ! GAC
  double precision boyspar,dx,dy,dz,result

  integer topdeg,t,u,v,n,i,j,counter,allochk
  double precision dr2,pi,arg,prefac
  double precision,allocatable::R(:,:,:,:),boys_table(:)
  integer t_part(56),u_part(56),v_part(56)
  double precision signum,e1,e2,factor,factor1,factor2

  dr2 = dx*dx+dy*dy+dz*dz
  arg = boyspar*dr2
  topdeg = deg1+deg2
  if (topdeg .gt. 12) then
     write(6,*)'g_MD_recur: can only handle up to f functions, exiting'
     stop
  endif
  allocate(boys_table(0:topdeg), R(0:topdeg,0:topdeg,0:topdeg,0:topdeg),&
           stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'GN_MD_REC:could not allocate boys_table, R; exiting'
     stop
  endif
  !initialize---multiply by factor to the power of boys func degree
  ! also times the factor for electrostatics of unit charge distributions
  ! this gives the value of R^n(0,0,0) or R(0,0,0,n)
  pi = 3.14159265358979323846d0
  prefac = 1.0d0 
  factor = (boyspar/pi)*sqrt(boyspar/pi)
  !factor1 = (pi/e1)*sqrt(pi/e1)
  !factor2 = (pi/e2)*sqrt(pi/e2)
  R(0,0,0,0) = factor*exp(-arg)
  do n = 0,topdeg
    R(0,0,0,n) = prefac*R(0,0,0,0)
    prefac = prefac * (-2.d0*boyspar)
  enddo
  ! next the v,z direction
  v = 1
  do n = 0,topdeg - v
    R(0,0,v,n) = dz*R(0,0,v-1,n+1)
  enddo
  do v = 2,topdeg
    do n = 0,topdeg - v 
      R(0,0,v,n) = (v-1)*R(0,0,v-2,n+1) + dz*R(0,0,v-1,n+1)
    enddo
  enddo
  ! next the u,y direction
  u = 1
  do v = 0,topdeg - u
    do n = 0,topdeg - (u+v)
      R(0,u,v,n) = dy*R(0,u-1,v,n+1)
    enddo
  enddo
  do u = 2,topdeg
    do v = 0,topdeg - u
      do n = 0,topdeg - (u+v)
        R(0,u,v,n) = (u-1)*R(0,u-2,v,n+1) + dy*R(0,u-1,v,n+1)
      enddo
    enddo
  enddo
  ! next the t,x direction
  t = 1
  do u = 0,topdeg - t
    do v = 0,topdeg - (t+u)
      do n = 0,topdeg - (t+u+v)
        R(t,u,v,n) = dx*R(t-1,u,v,n+1)
      enddo
    enddo
  enddo
  do t = 2,topdeg
    do u = 0,topdeg - t
      do v = 0,topdeg - (t+u)
        do n = 0,topdeg - (t+u+v)
          R(t,u,v,n) = (t-1)*R(t-2,u,v,n+1) + dx*R(t-1,u,v,n+1)
        enddo
      enddo
    enddo
  enddo

  if(mod((x2+y2+z2),2).eq.1)then
    signum = -1.d0
  else
    signum = 1.d0
  endif
  result = signum*R(x1+x2,y1+y2,z1+z2,0)
  !print *,'result = ',result
  deallocate(boys_table,R)
  return
end subroutine GN_MD_REC_OVERLAP
!-------------------------------------------------------
