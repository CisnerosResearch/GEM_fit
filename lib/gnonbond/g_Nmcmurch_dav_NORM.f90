!-------------------------------------------------------
subroutine GN_MCMUR_DAV_ss_NORM(boyspar,dx,dy,dz,coeff1,coeff2,ene,vec,e1,e2)
  ! this subroutine is to test the G matrix formation to compare with
  ! cartesian, normalized hermites.
  implicit none
  double precision boyspar,dx,dy,dz,coeff1,coeff2,ene,vec(*),e1,e2,bp2
  integer site1,site2,k1,k2

  integer topdeg,i
  double precision dr2,pi,arg,prefac,boys_table,pi52
  double precision R

  dr2 = dx*dx+dy*dy+dz*dz
  arg = boyspar*dr2
  topdeg = 0
  call GNTABLE_exact(topdeg,arg,boys_table)
  !initialize---multiply by factor to the power of boys func degree
  ! also times the factor for electrostatics of unit charge distributions
  ! this gives the value of R^n(0,0,0) or R(0,0,0,n)
  pi = 3.14159265358979323846d0
  pi52=pi
  do i = 1,4
    pi52 = pi*pi52
  enddo
  pi52 = sqrt(pi52)
  !prefac = 2.d0*sqrt(boyspar/pi)!this is the factor for electrostatics
                                !involving unit charge distributions
  bp2 = (2.d0*pi52)/(e1*e2*dsqrt(e1+e2))
   
  R = prefac*boys_table
  ene = coeff1*coeff2*R
  !print *,'R,lambda,arg = ',boys_table,bp2,boyspar
  vec(1) = bp2*boys_table
  !vec(1) = R
  return
end subroutine GN_MCMUR_DAV_ss_NORM
!-------------------------------------------------------
! GAC
subroutine GN_MCMUR_DAV_sp_NORM(boyspar,dx,dy,dz,coeff1,coeff2,ene,vec,e1,e2)
  implicit none
  double precision boyspar,dx,dy,dz,coeff1,coeff2(10),ene,vec(*),e1,e2,bp2

  integer topdeg,t,u,v,n,i
  double precision dr2,pi,arg,prefac,boys_table(0:2),pi52
  double precision R(0:1,0:1,0:1,0:1),R2(0:1,0:1,0:1,0:1)

  dr2 = dx*dx+dy*dy+dz*dz
  arg = boyspar*dr2
  topdeg = 1
  call GNTABLE_exact(topdeg,arg,boys_table)
  !initialize---multiply by factor to the power of boys func degree
  ! also times the factor for electrostatics of unit charge distributions
  ! this gives the value of R^n(0,0,0) or R(0,0,0,n)
  pi = 3.14159265358979323846d0
  pi52=pi
  do i = 1,4
    pi52 = pi*pi52
  enddo
  pi52 = sqrt(pi52)
  prefac = 1
  bp2 = (2.d0*pi52)/(e1*e2*dsqrt(e1+e2))
  do n = 0,1
    R(0,0,0,n) = prefac*boys_table(n)
    prefac = prefac * (-2.d0*boyspar)
  enddo
  ! next the v,z direction
  R(0,0,1,0) = dz*R(0,0,0,1)
  ! next the u,y direction
  R(0,1,0,0) = dy*R(0,0,0,1)
  ! next the t,x direction
  R(1,0,0,0) = dx*R(0,0,0,1)
  

  !note the signs
  ene = coeff1*coeff2(1)*R(0,0,0,0) + coeff1*coeff2(2)*R(1,0,0,0) + &
        coeff1*coeff2(3)*R(0,1,0,0) + coeff1*coeff2(4)*R(0,0,1,0) 
  vec(1) = R(0,0,0,0)*bp2
  vec(2) = R(1,0,0,0)*bp2
  vec(3) = R(0,1,0,0)*bp2
  vec(4) = R(0,0,1,0)*bp2
  return
end subroutine GN_MCMUR_DAV_sp_NORM
!-------------------------------------------------------
subroutine GN_MCMUR_DAV_sd_NORM(boyspar,dx,dy,dz,coeff1,coeff2,ene,vec,e1,e2)
  implicit none
  double precision boyspar,dx,dy,dz,coeff1,coeff2(10),ene,vec(*),e1,e2,bp2

  integer topdeg,t,u,v,n,i
  double precision dr2,pi,arg,prefac,boys_table(0:2),pi52
  double precision R(0:2,0:2,0:2,0:2)

  dr2 = dx*dx+dy*dy+dz*dz
  arg = boyspar*dr2
  topdeg = 2
  call GNTABLE_exact(topdeg,arg,boys_table)
  !initialize---multiply by factor to the power of boys func degree
  ! also times the factor for electrostatics of unit charge distributions
  ! this gives the value of R^n(0,0,0) or R(0,0,0,n)
  pi = 3.14159265358979323846d0
  pi52=pi
  do i = 1,4
    pi52 = pi*pi52
  enddo
  pi52 = sqrt(pi52)
  prefac =1
  bp2 = (2.d0*pi52)/(e1*e2*dsqrt(e1+e2))

  do n = 0,2
    R(0,0,0,n) = prefac*boys_table(n)
    prefac = prefac * (-2.d0*boyspar)
  enddo
  ! next the v,z direction
  R(0,0,1,0) = dz*R(0,0,0,1)
  R(0,0,1,1) = dz*R(0,0,0,2)
  R(0,0,2,0) = R(0,0,0,1) + dz*R(0,0,1,1)
  ! next the u,y direction
  R(0,1,0,0) = dy*R(0,0,0,1)
  R(0,1,1,0) = dy*R(0,0,1,1)
  R(0,1,0,1) = dy*R(0,0,0,2)
  R(0,2,0,0) = R(0,0,0,1) + dy*R(0,1,0,1)
  ! next the t,x direction
  R(1,0,0,0) = dx*R(0,0,0,1)
  R(1,1,0,0) = dx*R(0,1,0,1)
  R(1,0,1,0) = dx*R(0,0,1,1)
  R(1,0,0,1) = dx*R(0,0,0,2)
  R(2,0,0,0) = R(0,0,0,1) + dx*R(1,0,0,1)
  
  !note the signs
  ene = coeff1*coeff2(1)*R(0,0,0,0) + coeff1*coeff2(2)*R(1,0,0,0) + &
        coeff1*coeff2(3)*R(0,1,0,0) + coeff1*coeff2(4)*R(0,0,1,0) + &
        coeff1*coeff2(5)*R(2,0,0,0) + coeff1*coeff2(6)*R(0,2,0,0) + &
        coeff1*coeff2(7)*R(0,0,2,0) + coeff1*coeff2(8)*R(1,1,0,0) + &
        coeff1*coeff2(9)*R(1,0,1,0) + coeff1*coeff2(10)*R(0,1,1,0)
  vec(1) = R(0,0,0,0)*bp2
  vec(2) = R(1,0,0,0)*bp2
  vec(3) = R(0,1,0,0)*bp2
  vec(4) = R(0,0,1,0)*bp2
  vec(5) = R(2,0,0,0)*bp2
  vec(6) = R(0,2,0,0)*bp2
  vec(7) = R(0,0,2,0)*bp2
  vec(8) = R(1,1,0,0)*bp2
  vec(9) = R(1,0,1,0)*bp2
  vec(10) = R(0,1,1,0)*bp2
  return
end subroutine GN_MCMUR_DAV_sd_NORM
!-------------------------------------------------------
! GAC (handle d only shells)
subroutine GN_MCMUR_DAV_sda_NORM(boyspar,dx,dy,dz,coeff1,coeff2,ene,vec,e1,e2)
  implicit none
  double precision boyspar,dx,dy,dz,coeff1,coeff2(10),ene,vec(*),e1,e2,bp2

  integer topdeg,t,u,v,n,i
  double precision dr2,pi,arg,prefac,boys_table(0:2),pi52
  double precision R(0:2,0:2,0:2,0:2)

  dr2 = dx*dx+dy*dy+dz*dz
  arg = boyspar*dr2
  topdeg = 2
  call GNTABLE_exact(topdeg,arg,boys_table)
  !initialize---multiply by factor to the power of boys func degree
  ! also times the factor for electrostatics of unit charge distributions
  ! this gives the value of R^n(0,0,0) or R(0,0,0,n)
  pi = 3.14159265358979323846d0
  pi52=pi
  do i = 1,4
    pi52 = pi*pi52
  enddo
  pi52 = sqrt(pi52)
  prefac = 1
  bp2 = (2.d0*pi52)/(e1*e2*dsqrt(e1+e2))

  do n = 0,2
    R(0,0,0,n) = prefac*boys_table(n)
    prefac = prefac * (-2.d0*boyspar)
  enddo
  ! next the v,z direction
  R(0,0,1,0) = dz*R(0,0,0,1)
  R(0,0,1,1) = dz*R(0,0,0,2)
  R(0,0,2,0) = R(0,0,0,1) + dz*R(0,0,1,1)
  ! next the u,y direction
  R(0,1,0,0) = dy*R(0,0,0,1)
  R(0,1,1,0) = dy*R(0,0,1,1)
  R(0,1,0,1) = dy*R(0,0,0,2)
  R(0,2,0,0) = R(0,0,0,1) + dy*R(0,1,0,1)
  ! next the t,x direction
  R(1,0,0,0) = dx*R(0,0,0,1)
  R(1,1,0,0) = dx*R(0,1,0,1)
  R(1,0,1,0) = dx*R(0,0,1,1)
  R(1,0,0,1) = dx*R(0,0,0,2)
  R(2,0,0,0) = R(0,0,0,1) + dx*R(1,0,0,1)
  
  !note the signs
  ene = coeff1*coeff2(1)*R(2,0,0,0) + coeff1*coeff2(2)*R(0,2,0,0) + &
        coeff1*coeff2(3)*R(0,0,2,0) + coeff1*coeff2(4)*R(1,1,0,0) + &
        coeff1*coeff2(5)*R(1,0,1,0) + coeff1*coeff2(6)*R(0,1,1,0)
  vec(1) = R(2,0,0,0)*bp2
  vec(2) = R(0,2,0,0)*bp2
  vec(3) = R(0,0,2,0)*bp2
  vec(4) = R(1,1,0,0)*bp2
  vec(5) = R(1,0,1,0)*bp2
  vec(6) = R(0,1,1,0)*bp2
  return
end subroutine GN_MCMUR_DAV_sda_NORM
!-------------------------------------------------------
! GAC
subroutine GN_MCMUR_DAV_ps_NORM(boyspar,dx,dy,dz,coeff1,coeff2,ene,vec,e1,e2)
  implicit none
  double precision boyspar,dx,dy,dz,coeff1(10),coeff2,ene,vec(*),e1,e2,bp2

  integer topdeg,t,u,v,n,i
  double precision dr2,pi,arg,prefac,boys_table(0:2),pi52
  double precision R(0:1,0:1,0:1,0:1)

  dr2 = dx*dx+dy*dy+dz*dz
  arg = boyspar*dr2
  topdeg = 1
  call GNTABLE_exact(topdeg,arg,boys_table)
  !initialize---multiply by factor to the power of boys func degree
  ! also times the factor for electrostatics of unit charge distributions
  ! this gives the value of R^n(0,0,0) or R(0,0,0,n)
  pi = 3.14159265358979323846d0
  pi52=pi
  do i = 1,4
    pi52 = pi*pi52
  enddo
  pi52 = sqrt(pi52)
  prefac = 1
  bp2 = (2.d0*pi52)/(e1*e2*dsqrt(e1+e2))

  do n = 0,1
    R(0,0,0,n) = prefac*boys_table(n)
    prefac = prefac * (-2.d0*boyspar)
  enddo
  ! next the v,z direction
  R(0,0,1,0) = dz*R(0,0,0,1)
  ! next the u,y direction
  R(0,1,0,0) = dy*R(0,0,0,1)
  ! next the t,x direction
  R(1,0,0,0) = dx*R(0,0,0,1)

  ene = coeff2*coeff1(1)*R(0,0,0,0) - coeff2*coeff1(2)*R(1,0,0,0) - &
        coeff2*coeff1(3)*R(0,1,0,0) - coeff2*coeff1(4)*R(0,0,1,0)
  vec(1) = R(0,0,0,0)*bp2
  vec(2) = -R(1,0,0,0)*bp2
  vec(3) = -R(0,1,0,0)*bp2
  vec(4) = -R(0,0,1,0)*bp2
  return
end subroutine GN_MCMUR_DAV_ps_NORM
!-------------------------------------------------------
subroutine GN_MCMUR_DAV_ds_NORM(boyspar,dx,dy,dz,coeff1,coeff2,ene,vec,e1,e2)
  implicit none
  double precision boyspar,dx,dy,dz,coeff1(10),coeff2,ene,vec(*),e1,e2,bp2

  integer topdeg,t,u,v,n,i
  double precision dr2,pi,arg,prefac,boys_table(0:2),pi52
  double precision R(0:2,0:2,0:2,0:2)

  dr2 = dx*dx+dy*dy+dz*dz
  arg = boyspar*dr2
  topdeg = 2
  call GNTABLE_exact(topdeg,arg,boys_table)
  !initialize---multiply by factor to the power of boys func degree
  ! also times the factor for electrostatics of unit charge distributions
  ! this gives the value of R^n(0,0,0) or R(0,0,0,n)
  pi = 3.14159265358979323846d0
  pi52=pi
  do i = 1,4
    pi52 = pi*pi52
  enddo
  pi52 = sqrt(pi52)
  prefac = 1
  bp2 = (2.d0*pi52)/(e1*e2*dsqrt(e1+e2))

  do n = 0,2
    R(0,0,0,n) = prefac*boys_table(n)
    prefac = prefac * (-2.d0*boyspar)
  enddo
  ! next the v,z direction
  R(0,0,1,0) = dz*R(0,0,0,1)
  R(0,0,1,1) = dz*R(0,0,0,2)
  R(0,0,2,0) = R(0,0,0,1) + dz*R(0,0,1,1)
  ! next the u,y direction
  R(0,1,0,0) = dy*R(0,0,0,1)
  R(0,1,1,0) = dy*R(0,0,1,1)
  R(0,1,0,1) = dy*R(0,0,0,2)
  R(0,2,0,0) = R(0,0,0,1) + dy*R(0,1,0,1)
  ! next the t,x direction
  R(1,0,0,0) = dx*R(0,0,0,1)
  R(1,1,0,0) = dx*R(0,1,0,1)
  R(1,0,1,0) = dx*R(0,0,1,1)
  R(1,0,0,1) = dx*R(0,0,0,2)
  R(2,0,0,0) = R(0,0,0,1) + dx*R(1,0,0,1)

  ene = coeff2*coeff1(1)*R(0,0,0,0) - coeff2*coeff1(2)*R(1,0,0,0) - &
        coeff2*coeff1(3)*R(0,1,0,0) - coeff2*coeff1(4)*R(0,0,1,0) + &
        coeff2*coeff1(5)*R(2,0,0,0) + coeff2*coeff1(6)*R(0,2,0,0) + &
        coeff2*coeff1(7)*R(0,0,2,0) + coeff2*coeff1(8)*R(1,1,0,0) + &
        coeff2*coeff1(9)*R(1,0,1,0) + coeff2*coeff1(10)*R(0,1,1,0)
  vec(1) = R(0,0,0,0)*bp2
  vec(2) = -R(1,0,0,0)*bp2
  vec(3) = -R(0,1,0,0)*bp2
  vec(4) = -R(0,0,1,0)*bp2
  vec(5) = R(2,0,0,0)*bp2
  vec(6) = R(0,2,0,0)*bp2
  vec(7) = R(0,0,2,0)*bp2
  vec(8) = R(1,1,0,0)*bp2
  vec(9) = R(1,0,1,0)*bp2
  vec(10) = R(0,1,1,0)*bp2
  return
end subroutine GN_MCMUR_DAV_ds_NORM
!-------------------------------------------------------
!GAC (handle only d shells)
subroutine GN_MCMUR_DAV_das_NORM(boyspar,dx,dy,dz,coeff1,coeff2,ene,vec,e1,e2)
  implicit none
  double precision boyspar,dx,dy,dz,coeff1(10),coeff2,ene,vec(*),e1,e2,bp2

  integer topdeg,t,u,v,n,i
  double precision dr2,pi,arg,prefac,boys_table(0:2),pi52
  double precision R(0:2,0:2,0:2,0:2)

  dr2 = dx*dx+dy*dy+dz*dz
  arg = boyspar*dr2
  topdeg = 2
  call GNTABLE_exact(topdeg,arg,boys_table)
  !initialize---multiply by factor to the power of boys func degree
  ! also times the factor for electrostatics of unit charge distributions
  ! this gives the value of R^n(0,0,0) or R(0,0,0,n)
  pi = 3.14159265358979323846d0
  pi52=pi
  do i = 1,4
    pi52 = pi*pi52
  enddo
  pi52 = sqrt(pi52)
  prefac = 1
  bp2 = (2.d0*pi52)/(e1*e2*dsqrt(e1+e2))

  do n = 0,2
    R(0,0,0,n) = prefac*boys_table(n)
    prefac = prefac * (-2.d0*boyspar)
  enddo
  ! next the v,z direction
  R(0,0,1,0) = dz*R(0,0,0,1)
  R(0,0,1,1) = dz*R(0,0,0,2)
  R(0,0,2,0) = R(0,0,0,1) + dz*R(0,0,1,1)
  ! next the u,y direction
  R(0,1,0,0) = dy*R(0,0,0,1)
  R(0,1,1,0) = dy*R(0,0,1,1)
  R(0,1,0,1) = dy*R(0,0,0,2)
  R(0,2,0,0) = R(0,0,0,1) + dy*R(0,1,0,1)
  ! next the t,x direction
  R(1,0,0,0) = dx*R(0,0,0,1)
  R(1,1,0,0) = dx*R(0,1,0,1)
  R(1,0,1,0) = dx*R(0,0,1,1)
  R(1,0,0,1) = dx*R(0,0,0,2)
  R(2,0,0,0) = R(0,0,0,1) + dx*R(1,0,0,1)

  ene = coeff2*coeff1(1)*R(2,0,0,0) + coeff2*coeff1(2)*R(0,2,0,0) + &
        coeff2*coeff1(3)*R(0,0,2,0) + coeff2*coeff1(4)*R(1,1,0,0) + &
        coeff2*coeff1(5)*R(1,0,1,0) + coeff2*coeff1(6)*R(0,1,1,0)
  vec(1) = R(2,0,0,0)*bp2
  vec(2) = R(0,2,0,0)*bp2
  vec(3) = R(0,0,2,0)*bp2
  vec(4) = R(1,1,0,0)*bp2
  vec(5) = R(1,0,1,0)*bp2
  vec(6) = R(0,1,1,0)*bp2
  return
end subroutine GN_MCMUR_DAV_das_NORM
!-------------------------------------------------------
subroutine GN_MCMUR_DAV_pp_NORM(boyspar,dx,dy,dz,coeff1,coeff2,ene,vec,e1,e2)
  implicit none
  double precision boyspar,dx,dy,dz,coeff1(10),coeff2(10),ene,vec(*),e1,e2

  integer topdeg,t,u,v,n,i,j,counter
  double precision dr2,pi,arg,prefac,boys_table(0:2),bp2,pi52
  double precision R(0:2,0:2,0:2,0:2)
  integer t_part(56),u_part(56),v_part(56)
  double precision signum

  call GMcMur_Dav_get_tuv_parts(t_part,u_part,v_part)
  dr2 = dx*dx+dy*dy+dz*dz
  arg = boyspar*dr2
  topdeg = 2
  call GNTABLE_exact(topdeg,arg,boys_table)
  !initialize---multiply by factor to the power of boys func degree
  ! also times the factor for electrostatics of unit charge distributions
  ! this gives the value of R^n(0,0,0) or R(0,0,0,n)
  pi = 3.14159265358979323846d0
  pi52=pi
  do i = 1,4
    pi52 = pi*pi52
  enddo
  pi52 = sqrt(pi52)
  prefac = 1
  bp2 = (2.d0*pi52)/(e1*e2*dsqrt(e1+e2))

  do n = 0,2
    R(0,0,0,n) = prefac*boys_table(n)
    prefac = prefac * (-2.d0*boyspar)
  enddo
  ! next the v,z direction
  v = 1
  do n = 0,2 - v
    R(0,0,v,n) = dz*R(0,0,v-1,n+1)
  enddo
  do v = 2,2
    do n = 0,2 - v 
      R(0,0,v,n) = (v-1)*R(0,0,v-2,n+1) + dz*R(0,0,v-1,n+1)
    enddo
  enddo
  ! next the u,y direction
  u = 1
  do v = 0,2 - u
    do n = 0,2 - (u+v)
      R(0,u,v,n) = dy*R(0,u-1,v,n+1)
    enddo
  enddo
  do u = 2,2
    do v = 0,2 - u
      do n = 0,2 - (u+v)
        R(0,u,v,n) = (u-1)*R(0,u-2,v,n+1) + dy*R(0,u-1,v,n+1)
      enddo
    enddo
  enddo
  ! next the t,x direction
  t = 1
  do u = 0,2 - t
    do v = 0,2 - (t+u)
      do n = 0,2 - (t+u+v)
        R(t,u,v,n) = dx*R(t-1,u,v,n+1)
      enddo
    enddo
  enddo
  do t = 2,2
    do u = 0,2 - t
      do v = 0,2 - (t+u)
        do n = 0,2 - (t+u+v)
          R(t,u,v,n) = (t-1)*R(t-2,u,v,n+1) + dx*R(t-1,u,v,n+1)
        enddo
      enddo
    enddo
  enddo

  counter = 0
  ene = 0.d0
  do i = 1,4
    do j = 1,4
      if ( i .eq. j .and. i .gt. 1 )then ! porque estos son los p
        signum = -1.d0
      else
        signum = 1.d0
      endif
      t = t_part(i) + t_part(j)
      u = u_part(i) + u_part(j)
      v = v_part(i) + v_part(j)
      ene = ene + signum*coeff1(i)*coeff2(j)*R(t,u,v,0)
      counter = counter + 1
      vec(counter) = signum*R(t,u,v,0)*bp2
    enddo
  enddo
  return
end subroutine GN_MCMUR_DAV_pp_NORM
!-------------------------------------------------------
subroutine GN_MCMUR_DAV_dd_NORM(boyspar,dx,dy,dz,coeff1,coeff2,ene,vec,e1,e2)
  implicit none
  double precision boyspar,dx,dy,dz,coeff1(10),coeff2(10),ene,vec(*),e1,e2

  integer topdeg,t,u,v,n,i,j,counter
  double precision dr2,pi,arg,prefac,boys_table(0:4),bp2,pi52
  double precision R(0:4,0:4,0:4,0:4)
  integer t_part(56),u_part(56),v_part(56)
  double precision signum

  call GMcMur_Dav_get_tuv_parts(t_part,u_part,v_part)
  dr2 = dx*dx+dy*dy+dz*dz
  arg = boyspar*dr2
  topdeg = 4
  call GNTABLE_exact(topdeg,arg,boys_table)
  !initialize---multiply by factor to the power of boys func degree
  ! also times the factor for electrostatics of unit charge distributions
  ! this gives the value of R^n(0,0,0) or R(0,0,0,n)
  pi = 3.14159265358979323846d0
  pi52=pi
  do i = 1,4
    pi52 = pi*pi52
  enddo
  pi52 = sqrt(pi52)
  prefac = 1
  bp2 = (2.d0*pi52)/(e1*e2*dsqrt(e1+e2))

  do n = 0,4
    R(0,0,0,n) = prefac*boys_table(n)
    prefac = prefac * (-2.d0*boyspar)
  enddo
  ! next the v,z direction
  v = 1
  do n = 0,4 - v
    R(0,0,v,n) = dz*R(0,0,v-1,n+1)
  enddo
  do v = 2,4
    do n = 0,4 - v 
      R(0,0,v,n) = (v-1)*R(0,0,v-2,n+1) + dz*R(0,0,v-1,n+1)
    enddo
  enddo
  ! next the u,y direction
  u = 1
  do v = 0,4 - u
    do n = 0,4 - (u+v)
      R(0,u,v,n) = dy*R(0,u-1,v,n+1)
    enddo
  enddo
  do u = 2,4
    do v = 0,4 - u
      do n = 0,4 - (u+v)
        R(0,u,v,n) = (u-1)*R(0,u-2,v,n+1) + dy*R(0,u-1,v,n+1)
      enddo
    enddo
  enddo
  ! next the t,x direction
  t = 1
  do u = 0,4 - t
    do v = 0,4 - (t+u)
      do n = 0,4 - (t+u+v)
        R(t,u,v,n) = dx*R(t-1,u,v,n+1)
      enddo
    enddo
  enddo
  do t = 2,4
    do u = 0,4 - t
      do v = 0,4 - (t+u)
        do n = 0,4 - (t+u+v)
          R(t,u,v,n) = (t-1)*R(t-2,u,v,n+1) + dx*R(t-1,u,v,n+1)
        enddo
      enddo
    enddo
  enddo

  counter = 0
  ene = 0.d0
  do i = 1,10
    if ( i == 2 .or. i == 3 .or. i == 4 )then ! porque estos son los p
      signum = -1.d0
    else
      signum = 1.d0
    endif
    do j = 1,10
      t = t_part(i) + t_part(j)
      u = u_part(i) + u_part(j)
      v = v_part(i) + v_part(j)
      ene = ene + signum*coeff1(i)*coeff2(j)*R(t,u,v,0)
      counter = counter + 1
      vec(counter) = signum*R(t,u,v,0)*bp2
    enddo
  enddo
  return
end subroutine GN_MCMUR_DAV_dd_NORM
!-------------------------------------------------------
subroutine GN_MCMUR_DAV_dada_NORM(boyspar,dx,dy,dz,coeff1,coeff2,ene,vec,e1,e2)
  implicit none
  double precision boyspar,dx,dy,dz,coeff1(10),coeff2(10),ene,vec(*),e1,e2

  integer topdeg,t,u,v,n,i,j,counter
  double precision dr2,pi,arg,prefac,boys_table(0:4),bp2,pi52
  double precision R(0:4,0:4,0:4,0:4)
  integer t_part(56),u_part(56),v_part(56)
  double precision signum

  call GMcMur_Dav_get_tuv_parts(t_part,u_part,v_part)
  dr2 = dx*dx+dy*dy+dz*dz
  arg = boyspar*dr2
  topdeg = 4
  call GNTABLE_exact(topdeg,arg,boys_table)
  !initialize---multiply by factor to the power of boys func degree
  ! also times the factor for electrostatics of unit charge distributions
  ! this gives the value of R^n(0,0,0) or R(0,0,0,n)
  pi = 3.14159265358979323846d0
  pi52=pi
  do i = 1,4
    pi52 = pi*pi52
  enddo
  pi52 = sqrt(pi52)
  prefac = 1
  bp2 = (2.d0*pi52)/(e1*e2*dsqrt(e1+e2))

  do n = 0,4
    R(0,0,0,n) = prefac*boys_table(n)
    prefac = prefac * (-2.d0*boyspar)
  enddo
  ! next the v,z direction
  v = 1
  do n = 0,4 - v
    R(0,0,v,n) = dz*R(0,0,v-1,n+1)
  enddo
  do v = 2,4
    do n = 0,4 - v 
      R(0,0,v,n) = (v-1)*R(0,0,v-2,n+1) + dz*R(0,0,v-1,n+1)
    enddo
  enddo
  ! next the u,y direction
  u = 1
  do v = 0,4 - u
    do n = 0,4 - (u+v)
      R(0,u,v,n) = dy*R(0,u-1,v,n+1)
    enddo
  enddo
  do u = 2,4
    do v = 0,4 - u
      do n = 0,4 - (u+v)
        R(0,u,v,n) = (u-1)*R(0,u-2,v,n+1) + dy*R(0,u-1,v,n+1)
      enddo
    enddo
  enddo
  ! next the t,x direction
  t = 1
  do u = 0,4 - t
    do v = 0,4 - (t+u)
      do n = 0,4 - (t+u+v)
        R(t,u,v,n) = dx*R(t-1,u,v,n+1)
      enddo
    enddo
  enddo
  do t = 2,4
    do u = 0,4 - t
      do v = 0,4 - (t+u)
        do n = 0,4 - (t+u+v)
          R(t,u,v,n) = (t-1)*R(t-2,u,v,n+1) + dx*R(t-1,u,v,n+1)
        enddo
      enddo
    enddo
  enddo

  counter = 0
  ene = 0.d0
  do i = 1,6
    do j = 1,6
      t = t_part(i+4) + t_part(j+4)
      u = u_part(i+4) + u_part(j+4)
      v = v_part(i+4) + v_part(j+4)
      ene = ene + signum*coeff1(i)*coeff2(j)*R(t,u,v,0)
      counter = counter + 1
      vec(counter) = signum*R(t,u,v,0)*bp2
    enddo
  enddo
  return
end subroutine GN_MCMUR_DAV_dada_NORM
!-------------------------------------------------------
subroutine GN_MCMUR_DAV_dad_NORM(boyspar,dx,dy,dz,coeff1,coeff2,ene,vec,e1,e2)
  implicit none
  double precision boyspar,dx,dy,dz,coeff1(10),coeff2(10),ene,vec(*),e1,e2

  integer topdeg,t,u,v,n,i,j,counter
  double precision dr2,pi,arg,prefac,boys_table(0:4),bp2,pi52
  double precision R(0:4,0:4,0:4,0:4)
  integer t_part(56),u_part(56),v_part(56)
  double precision signum

  call GMcMur_Dav_get_tuv_parts(t_part,u_part,v_part)
  dr2 = dx*dx+dy*dy+dz*dz
  arg = boyspar*dr2
  topdeg = 4
  call GNTABLE_exact(topdeg,arg,boys_table)
  !initialize---multiply by factor to the power of boys func degree
  ! also times the factor for electrostatics of unit charge distributions
  ! this gives the value of R^n(0,0,0) or R(0,0,0,n)
  pi = 3.14159265358979323846d0
  pi52=pi
  do i = 1,4
    pi52 = pi*pi52
  enddo
  pi52 = sqrt(pi52)
  prefac = 1
  bp2 = (2.d0*pi52)/(e1*e2*dsqrt(e1+e2))

  do n = 0,4
    R(0,0,0,n) = prefac*boys_table(n)
    prefac = prefac * (-2.d0*boyspar)
  enddo
  ! next the v,z direction
  v = 1
  do n = 0,4 - v
    R(0,0,v,n) = dz*R(0,0,v-1,n+1)
  enddo
  do v = 2,4
    do n = 0,4 - v 
      R(0,0,v,n) = (v-1)*R(0,0,v-2,n+1) + dz*R(0,0,v-1,n+1)
    enddo
  enddo
  ! next the u,y direction
  u = 1
  do v = 0,4 - u
    do n = 0,4 - (u+v)
      R(0,u,v,n) = dy*R(0,u-1,v,n+1)
    enddo
  enddo
  do u = 2,4
    do v = 0,4 - u
      do n = 0,4 - (u+v)
        R(0,u,v,n) = (u-1)*R(0,u-2,v,n+1) + dy*R(0,u-1,v,n+1)
      enddo
    enddo
  enddo
  ! next the t,x direction
  t = 1
  do u = 0,4 - t
    do v = 0,4 - (t+u)
      do n = 0,4 - (t+u+v)
        R(t,u,v,n) = dx*R(t-1,u,v,n+1)
      enddo
    enddo
  enddo
  do t = 2,4
    do u = 0,4 - t
      do v = 0,4 - (t+u)
        do n = 0,4 - (t+u+v)
          R(t,u,v,n) = (t-1)*R(t-2,u,v,n+1) + dx*R(t-1,u,v,n+1)
        enddo
      enddo
    enddo
  enddo

  counter = 0
  ene = 0.d0
  do i = 1,6
    do j = 1,10
      if ( j == 2 .or. j == 3 .or. j == 4 )then ! porque estos son los p
        signum = -1.d0
      else
        signum = 1.d0
      endif
      t = t_part(i+4) + t_part(j)
      u = u_part(i+4) + u_part(j)
      v = v_part(i+4) + v_part(j)
      ene = ene + signum*coeff1(i)*coeff2(j)*R(t,u,v,0)
      counter = counter + 1
      vec(counter) = signum*R(t,u,v,0)*bp2
    enddo
  enddo
  return
end subroutine GN_MCMUR_DAV_dad_NORM
!-------------------------------------------------------
subroutine GN_MCMUR_DAV_dda_NORM(boyspar,dx,dy,dz,coeff1,coeff2,ene,vec,e1,e2)
  implicit none
  double precision boyspar,dx,dy,dz,coeff1(10),coeff2(10),ene,vec(*),e1,e2

  integer topdeg,t,u,v,n,i,j,counter
  double precision dr2,pi,arg,prefac,boys_table(0:4),bp2,pi52
  double precision R(0:4,0:4,0:4,0:4)
  integer t_part(56),u_part(56),v_part(56)
  double precision signum

  call GMcMur_Dav_get_tuv_parts(t_part,u_part,v_part)
  dr2 = dx*dx+dy*dy+dz*dz
  arg = boyspar*dr2
  topdeg = 4
  call GNTABLE_exact(topdeg,arg,boys_table)
  !initialize---multiply by factor to the power of boys func degree
  ! also times the factor for electrostatics of unit charge distributions
  ! this gives the value of R^n(0,0,0) or R(0,0,0,n)
  pi = 3.14159265358979323846d0
  pi52=pi
  do i = 1,4
    pi52 = pi*pi52
  enddo
  pi52 = sqrt(pi52)
  prefac = 1
  bp2 = (2.d0*pi52)/(e1*e2*dsqrt(e1+e2))

  do n = 0,4
    R(0,0,0,n) = prefac*boys_table(n)
    prefac = prefac * (-2.d0*boyspar)
  enddo
  ! next the v,z direction
  v = 1
  do n = 0,4 - v
    R(0,0,v,n) = dz*R(0,0,v-1,n+1)
  enddo
  do v = 2,4
    do n = 0,4 - v 
      R(0,0,v,n) = (v-1)*R(0,0,v-2,n+1) + dz*R(0,0,v-1,n+1)
    enddo
  enddo
  ! next the u,y direction
  u = 1
  do v = 0,4 - u
    do n = 0,4 - (u+v)
      R(0,u,v,n) = dy*R(0,u-1,v,n+1)
    enddo
  enddo
  do u = 2,4
    do v = 0,4 - u
      do n = 0,4 - (u+v)
        R(0,u,v,n) = (u-1)*R(0,u-2,v,n+1) + dy*R(0,u-1,v,n+1)
      enddo
    enddo
  enddo
  ! next the t,x direction
  t = 1
  do u = 0,4 - t
    do v = 0,4 - (t+u)
      do n = 0,4 - (t+u+v)
        R(t,u,v,n) = dx*R(t-1,u,v,n+1)
      enddo
    enddo
  enddo
  do t = 2,4
    do u = 0,4 - t
      do v = 0,4 - (t+u)
        do n = 0,4 - (t+u+v)
          R(t,u,v,n) = (t-1)*R(t-2,u,v,n+1) + dx*R(t-1,u,v,n+1)
        enddo
      enddo
    enddo
  enddo

  counter = 0
  ene = 0.d0
  do i = 1,10
    if ( i == 2 .or. i == 3 .or. i == 4 )then ! porque estos son los p
      signum = -1.d0
    else
      signum = 1.d0
    endif
    do j = 1,6
      t = t_part(i) + t_part(j+4)
      u = u_part(i) + u_part(j+4)
      v = v_part(i) + v_part(j+4)
      ene = ene + signum*coeff1(i)*coeff2(j)*R(t,u,v,0)
      counter = counter + 1
      vec(counter) = signum*R(t,u,v,0)*bp2
    enddo
  enddo
  return
end subroutine GN_MCMUR_DAV_dda_NORM
!-------------------------------------------------------
! GAC
subroutine GN_MCMUR_DAV_RECp_NORM(topi,topj,boyspar,dx,dy,dz,coeff1,coeff2,&
                                  ene,vec,e1,e2)
  implicit none
  integer topi, topj ! GAC
  double precision boyspar,dx,dy,dz,coeff1(10),coeff2(10),ene,vec(*),e1,e2

  integer topdeg,t,u,v,n,i,j,counter
  double precision dr2,pi,arg,prefac,boys_table(0:4),bp2,pi52
  double precision R(0:4,0:4,0:4,0:4)
  integer t_part(56),u_part(56),v_part(56)
  double precision signum

  call GMcMur_Dav_get_tuv_parts(t_part,u_part,v_part)
  dr2 = dx*dx+dy*dy+dz*dz
  arg = boyspar*dr2
  topdeg = 3
  call GNTABLE_exact(topdeg,arg,boys_table)
  !initialize---multiply by factor to the power of boys func degree
  ! also times the factor for electrostatics of unit charge distributions
  ! this gives the value of R^n(0,0,0) or R(0,0,0,n)
  pi = 3.14159265358979323846d0
  pi52=pi
  do i = 1,4
    pi52 = pi*pi52
  enddo
  pi52 = sqrt(pi52)
  prefac = 1
  bp2 = (2.d0*pi52)/(e1*e2*dsqrt(e1+e2))

  do n = 0,3
    R(0,0,0,n) = prefac*boys_table(n)
    prefac = prefac * (-2.d0*boyspar)
  enddo
  ! next the v,z direction
  v = 1
  do n = 0,3 - v
    R(0,0,v,n) = dz*R(0,0,v-1,n+1)
  enddo
  do v = 2,3
    do n = 0,3 - v 
      R(0,0,v,n) = (v-1)*R(0,0,v-2,n+1) + dz*R(0,0,v-1,n+1)
    enddo
  enddo
  ! next the u,y direction
  u = 1
  do v = 0,3 - u
    do n = 0,3 - (u+v)
      R(0,u,v,n) = dy*R(0,u-1,v,n+1)
    enddo
  enddo
  do u = 2,3
    do v = 0,3 - u
      do n = 0,3 - (u+v)
        R(0,u,v,n) = (u-1)*R(0,u-2,v,n+1) + dy*R(0,u-1,v,n+1)
      enddo
    enddo
  enddo
  ! next the t,x direction
  t = 1
  do u = 0,3 - t
    do v = 0,3 - (t+u)
      do n = 0,3 - (t+u+v)
        R(t,u,v,n) = dx*R(t-1,u,v,n+1)
      enddo
    enddo
  enddo
  do t = 2,3
    do u = 0,3 - t
      do v = 0,3 - (t+u)
        do n = 0,3 - (t+u+v)
          R(t,u,v,n) = (t-1)*R(t-2,u,v,n+1) + dx*R(t-1,u,v,n+1)
        enddo
      enddo
    enddo
  enddo

  counter = 0
  ene = 0.d0
  do i = 1,topi
    if ( i == 2 .or. i == 3 .or. i == 4 )then ! porque estos son los p
      signum = -1.d0
    else
      signum = 1.d0
    endif
    do j = 1,topj
      t = t_part(i) + t_part(j)
      u = u_part(i) + u_part(j)
      v = v_part(i) + v_part(j)
      ene = ene + signum*coeff1(i)*coeff2(j)*R(t,u,v,0)
      counter = counter + 1
      vec(counter) = signum*R(t,u,v,0)*bp2
    enddo
  enddo
  return
end subroutine GN_MCMUR_DAV_RECp_NORM
!-------------------------------------------------------
