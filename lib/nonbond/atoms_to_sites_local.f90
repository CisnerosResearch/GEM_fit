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
!------------------------------------------------------
!------------------------------------------------------
subroutine ATMSITE_rotate_mpoles( &
                 nsites,indframe,MpoleOrder,     &                !input
                 MpoleOffset,                    &                !input
                 frame,Mpole_loc,                &
                 Mpole_xyz)                                       !output
  implicit none
  integer nsites
  integer indframe(nsites),MpoleOrder(nsites),MpoleOffset(nsites)
  double precision frame(3,3,*),Mpole_loc(*),Mpole_xyz(*)

  integer num,nf
include "mpole_sizes.fh"
  integer n,order,off,k,dimxy
  double precision Mpole_xy(MAXMP*MAXMP)
  do n = 1,nsites
    off = MpoleOffset(n)
    order = MpoleOrder(n)
    dimxy = order
    if ( order .gt. 0 )then  ! skip if order is zero
      if ( order .eq. 1 )then   !no frame needed
        Mpole_xyz(off+1) = Mpole_loc(off+1)
      else
        k = indframe(n)
        call XFORM_MPOLE_matrix(frame(1,1,k),Mpole_xy,order)
        call XFORM_MPOLE(Mpole_xy,dimxy,Mpole_loc(off+1), &
              Mpole_xyz(off+1),order)
      endif
    endif 
  enddo
  return
end subroutine ATMSITE_rotate_mpoles
!------------------------------------------------------------
subroutine ATMSITE_build_centers(numatoms,num_extra_centers,num_cenlist, &
                             cen_list,cen_wt,sitecrd)
  implicit none
  integer numatoms,num_extra_centers,num_cenlist
  integer cen_list(2,num_cenlist)
  double precision cen_wt(num_cenlist)
  double precision sitecrd(3,*)   !this array holds coords for all sites

  integer k,m,n
  double precision w

  if ( num_cenlist == 0 )return
  do n = 1,num_cenlist
    k = cen_list(1,n)  !contributing atom site
    m = cen_list(2,n)  !resulting extra center
    w = cen_wt(n)  ! weight of atom k contribution
    sitecrd(1,m) = sitecrd(1,m) + w*sitecrd(1,k)
    sitecrd(2,m) = sitecrd(2,m) + w*sitecrd(2,k)
    sitecrd(3,m) = sitecrd(3,m) + w*sitecrd(3,k)
  enddo
  return
end subroutine ATMSITE_build_centers
!----------------------------------------------------------------------
subroutine ATMSITE_build_frame_def_pts(   &
                 numfrdeflist,framedeflist,sitecrd,nframes,frdefpt,frvalid)
  implicit none
  integer numfrdeflist,nframes
  integer framedeflist(5,numfrdeflist)
  double precision sitecrd(3,*),frdefpt(3,2,nframes),epsilon
  logical frvalid(nframes)
  parameter (epsilon = 1.d-16)

  integer n,i,j,k,l,m
  double precision dx,dy,dz,wt

  do n = 1,numfrdeflist
    i = framedeflist(1,n)   ! center i gives tail of vector
    j = framedeflist(2,n)   ! center j gives head of vector
    k = framedeflist(3,n)   ! the frame number this vector helps define
    l = framedeflist(4,n)   ! the frame defpoint (1 or 2) within frame
    m = framedeflist(5,n)   ! m is number of such (unit) vectors 
                                ! averaged to get the frame defpoint k
    ! GAC : changing to match AMOEBA IN AMBER
    if (i > 0 .and. j > 0) then
       frvalid(k) = .true.
       dx = sitecrd(1,j) - sitecrd(1,i)
       dy = sitecrd(2,j) - sitecrd(2,i)
       dz = sitecrd(3,j) - sitecrd(3,i)
       wt = m*sqrt(dx*dx+dy*dy+dz*dz) ! divide by length of ij times num pairs
!    !! GAC: handle divide by 0
    if (abs(wt) .lt. epsilon) then
       frdefpt(1,l,k) = frdefpt(1,l,k) 
       frdefpt(2,l,k) = frdefpt(2,l,k) 
       frdefpt(3,l,k) = frdefpt(3,l,k)
    else
       frdefpt(1,l,k) = frdefpt(1,l,k) + dx/wt
       frdefpt(2,l,k) = frdefpt(2,l,k) + dy/wt
       frdefpt(3,l,k) = frdefpt(3,l,k) + dz/wt
    endif
    else
      frvalid(k) = .false.
      !print *,'frame not valid',k
    endif
  enddo
  return
end subroutine ATMSITE_build_frame_def_pts
!----------------------------------------------------------------------
subroutine ATMSITE_defpoints_to_frames(nframes,frame_axis,frdefpt,frame,&
                                       molframe,frvalid)
  implicit none
  integer nframes,frame_axis(3),molframe
  double precision frdefpt(3,2,nframes),frame(3,3,nframes),epsilon
  logical frvalid(nframes)
  parameter (epsilon = 1.d-16)

  integer i,n,k1,k2,k3
  double precision u(3),v(3),w(3),siz,dot

  k1 = frame_axis(1)
  k2 = frame_axis(2)
  k3 = frame_axis(3)
  do n = 1,nframes
! u is unit vector in direction of primary def pt
    if (frvalid(n)) then
       do i = 1,3
         u(i) = frdefpt(i,1,n)
       enddo
       siz = sqrt(u(1)*u(1)+u(2)*u(2)+u(3)*u(3))
       do i = 1,3
         ! GAC: handle divide by 0
         if (abs(siz) .lt. epsilon) then
            u(i) = 0.d0
         else
            u(i) = u(i) / siz
         endif
       enddo
! v is unit vector given by component of secondary pt orthog to u
       do i = 1,3
         v(i) = frdefpt(i,2,n)
       enddo
       dot = u(1)*v(1)+u(2)*v(2)+u(3)*v(3)
       do i = 1,3
         v(i) = v(i) - dot*u(i)
       enddo
       siz = sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
       do i = 1,3
         ! GAC: handle divide by 0
         if (abs(siz) .lt. epsilon) then
            v(i) = 0.d0
         else
            v(i) = v(i) / siz
         endif
       enddo
! w is u cross v
       w(1) = u(2)*v(3) - u(3)*v(2)
       w(2) = u(3)*v(1) - u(1)*v(3)
       w(3) = u(1)*v(2) - u(2)*v(1)
!  now define frame
       do i = 1,3
         frame(i,k1,n) = u(i)
         frame(i,k2,n) = v(i)
         frame(i,k3,n) = w(i)
       enddo
    endif
    !if(n == molframe)then
    !  write(12,*)'u = ',u
    !  write(12,*)'v = ',v
    !  write(12,*)'w = ',w
    !endif
  enddo

  return
end subroutine ATMSITE_defpoints_to_frames
!----------------------------------------------------------------------
subroutine ATMSITE_build_extrasites(num_centers,num_extra_sites, &
                 extra_site_list,indframe,loc_vector,frame,sitecrd)
  implicit none
  integer  num_centers,num_extra_sites
  integer extra_site_list(2,num_extra_sites),indframe(*)
  double precision loc_vector(3,num_extra_sites),frame(3,3,*),sitecrd(3,*)

  integer k,l,m,n
  double precision u,v,w
  if ( num_extra_sites == 0 )return
  do n = 1,num_extra_sites
    k = extra_site_list(1,n)   !center number for this site
    l = extra_site_list(2,n)   !frame number for this site
    m = num_centers + n  !index of the site in the coordinate list
    indframe(m) = l   ! assign this for torque->frame 
    u = loc_vector(1,n)
    v = loc_vector(2,n)
    w = loc_vector(3,n)
    sitecrd(1,m)=sitecrd(1,k) + u*frame(1,1,l)+v*frame(1,2,l)+w*frame(1,3,l)
    sitecrd(2,m)=sitecrd(2,k) + u*frame(2,1,l)+v*frame(2,2,l)+w*frame(2,3,l)
    sitecrd(3,m)=sitecrd(3,k) + u*frame(3,1,l)+v*frame(3,2,l)+w*frame(3,3,l)
  enddo
  return
end subroutine ATMSITE_build_extrasites
!------------------------------------------------------------------------------
subroutine ATMSITE_get_de_drot_mpole(nsites,Cphi,Mpole_xyz,de_drotsite, &
                                  Mpole_offset,MpOrder,coulomb_constant)
  implicit none
  integer nsites,Mpole_offset(*),MpOrder(*)
  double precision Cphi(*),Mpole_xyz(*),de_drotsite(3,*),coulomb_constant
!  get the derivative of electrostatic energy with respect to infinitesmal 
!  rotations of atomic frames about the x,y,z global axes
! Basic idea--electrostatic energy given by dot product of cartesian
! multipoles and the electrostatic potential and its cartesian derivatives
! i.e. ene = 1/2 *(q*phi + mux*dphidx + ...)
! Thus derivative obtained by rotating multipoles infinitesmally

include "mpole_sizes.fh"
  integer i,j,k,n,off,dimxy
  double precision DMP_x(MAXMP*MAXMP),DMP_y(MAXMP*MAXMP),  &
         DMP_z(MAXMP*MAXMP),A_xy(3,3),DA_xy(3,3),  &
         Tmp_x(MAXMP),Tmp_y(MAXMP),Tmp_z(MAXMP)

! to get de_drot we calculate the deriv of mpole wrt infinitesmal
! rotations about x,y and z axis
  do i = 1,3
    do j = 1,3
      A_xy(i,j) = 0.d0
    enddo
    A_xy(i,i) = 1.d0
  enddo
     
! do the maximal order
! x-axis rotation of dtheta
  do i = 1,3
    do j = 1,3
      DA_xy(i,j) = 0.d0
    enddo
  enddo
  DA_xy(3,2) = 1.d0
  DA_xy(2,3) = -1.d0
  call XFORM_MPOLE_deriv_matrix(A_xy,DA_xy,DMP_x,MAXMP)
! y-axis
  do i = 1,3
    do j = 1,3
      DA_xy(i,j) = 0.d0
    enddo
  enddo
  DA_xy(3,1) = -1.d0
  DA_xy(1,3) = 1.d0
  call XFORM_MPOLE_deriv_matrix(A_xy,DA_xy,DMP_y,MAXMP)

! z-axis
  do i = 1,3
    do j = 1,3
      DA_xy(i,j) = 0.d0
    enddo
  enddo
  DA_xy(2,1) = 1.d0
  DA_xy(1,2) = -1.d0
  call XFORM_MPOLE_deriv_matrix(A_xy,DA_xy,DMP_z,MAXMP)
 
  dimxy = MAXMP
  do n = 1,nsites
    de_drotsite(1,n) = 0.d0
    de_drotsite(2,n) = 0.d0
    de_drotsite(3,n) = 0.d0
    if ( MpOrder(n) .gt. 1 )then  !skip those with no torque
      off = Mpole_offset(n)
      call XFORM_MPOLE(DMP_x,dimxy,Mpole_xyz(off+1),Tmp_x,MpOrder(n))
      call XFORM_MPOLE(DMP_y,dimxy,Mpole_xyz(off+1),Tmp_y,MpOrder(n))
      call XFORM_MPOLE(DMP_z,dimxy,Mpole_xyz(off+1),Tmp_z,MpOrder(n))
      do k = 1,MpOrder(n)
        de_drotsite(1,n) = de_drotsite(1,n) +  &
                           coulomb_constant*Tmp_x(k)*Cphi(off+k)
        de_drotsite(2,n) = de_drotsite(2,n) +  &
                           coulomb_constant*Tmp_y(k)*Cphi(off+k)
        de_drotsite(3,n) = de_drotsite(3,n) +  &
                           coulomb_constant*Tmp_z(k)*Cphi(off+k)
      enddo
    endif
  enddo
  return
end subroutine ATMSITE_get_de_drot_mpole
!----------------------------------------------------------------------
subroutine ATMSITE_add_extra_site_contrib(  &
                num_centers,num_extra_sites,extra_site_list,sitecrd,  &
                sitefrc,de_drotsite)
  implicit none
  integer num_centers,num_extra_sites,extra_site_list(2,num_extra_sites)      
  double precision sitecrd(3,*),sitefrc(3,*),de_drotsite(3,*)
!----------------------------------------------------------------------
! add the derivs of energy due to translation of extra sites in terms
! of frame rotations. The torque is given by r x F, so the term we are looking
! for is given by minus the torque
!----------------------------------------------------------------------
  integer k,l,m,n
  double precision r(3),torque(3)
  do n = 1,num_extra_sites
    k = extra_site_list(1,n)   !center number for this site
    l = extra_site_list(1,n)   !frame number for this site
    m = num_centers + n  !index of the site in the coordinate list
    r(1) = sitecrd(1,m) - sitecrd(1,k)
    r(2) = sitecrd(2,m) - sitecrd(2,k)
    r(3) = sitecrd(3,m) - sitecrd(3,k)
    torque(1) = r(2)*sitefrc(3,m) - r(3)*sitefrc(2,m)
    torque(2) = r(3)*sitefrc(1,m) - r(1)*sitefrc(3,m)
    torque(3) = r(1)*sitefrc(2,m) - r(2)*sitefrc(1,m)
    de_drotsite(1,m) = de_drotsite(1,m) - torque(1)
    de_drotsite(2,m) = de_drotsite(2,m) - torque(2)
    de_drotsite(3,m) = de_drotsite(3,m) - torque(3)
  enddo
  return
end subroutine ATMSITE_add_extra_site_contrib
!----------------------------------------------------------------------
subroutine ATMSITE_accum_de_dframe_rot(  &
               nsites,nframes,indframe,frame_axis,    &      !input
               de_drotsite,frame,          &                       !intput
               frdefpt,                    &
               de_drotframe)                                       !output
!----------------------------------------------------------------------
  implicit none
  integer nsites,nframes,indframe(nsites),frame_axis(3)
  double precision de_drotsite(3,nsites),frame(3,3,nframes)
  double precision frdefpt(3,2,nframes)
  double precision de_drotframe(3,nframes)

  integer n,j,k,k1,k2,k3
  double precision p2unit(3),siz
  do n = 1,nframes   !clear the derivs wrt frame rotation
    do j = 1,3
      de_drotframe(j,n) = 0.d0
    enddo
  enddo
  k1 = frame_axis(1)
  k2 = frame_axis(2)
  k3 = frame_axis(3)
! deriv of energy with respect to rotation about unit vectors along
! p1 p2 and their cross product
! note that unit vector along p1 corresponds to k1st frame axis
! and unit vector in p1 x p2 direction corresponds to k3rd frame axis
! the energy derivative with respect to rotation about any unit vector
! is for each mpole given by the dot product of de_drotpole 
! (which is negative of torque due to that mpole) with the unit vector
! a frame may define several mpoles; hence the need for  pointers
  do n = 1,nsites
    k = indframe(n)
    if ( k .gt. 0 )then   !skip those with no frame pointers
      siz = sqrt(frdefpt(1,2,k)**2+frdefpt(2,2,k)**2+frdefpt(3,2,k)**2)
      do j = 1,3
        p2unit(j) = frdefpt(j,2,k) / siz
      enddo
      de_drotframe(k1,k) = de_drotframe(k1,k) +   &
                          de_drotsite(1,n)*frame(1,k1,k) +  &
                          de_drotsite(2,n)*frame(2,k1,k) +  &
                          de_drotsite(3,n)*frame(3,k1,k)
      de_drotframe(k2,k) = de_drotframe(k2,k) +   &
                          de_drotsite(1,n)*p2unit(1) +   &
                          de_drotsite(2,n)*p2unit(2) +   &
                          de_drotsite(3,n)*p2unit(3)
      de_drotframe(k3,k) = de_drotframe(k3,k) +    &
                          de_drotsite(1,n)*frame(1,k3,k) +  &
                          de_drotsite(2,n)*frame(2,k3,k) +  &
                          de_drotsite(3,n)*frame(3,k3,k)
    endif
  enddo
      
  return
end subroutine ATMSITE_accum_de_dframe_rot
!----------------------------------------------------------------------
subroutine ATMSITE_accum_de_ddefpts(   &
           nframes,frame_axis,   &                       !input
           de_drotframe,frame,defpts,  &                        !input
           de_ddefpt)                                          !output
  implicit none
  integer nframes,frame_axis(3)
  double precision de_drotframe(3,nframes),frame(3,3,nframes)
  double precision defpts(3,2,nframes)
  double precision de_ddefpt(3,2,nframes)
! get the derivs of energy with respect to movement of defpoints
! expressed in the local frame coord system
  integer n,k,j,k1,k2,k3
  double precision p1(3),p2(3),p2unit(3),p2perp1(3),p1perp2(3),  &
         u(3),v(3),w(3),dotu,dotv,dotw,  &
         sizp1perp2,sizp2perp1,dot12,dot21,  &
         sizp1,sizp2,dedrotp1,dedrotp2,dedu,dedv,dedw, &
         de_drotu,de_drotv,de_drotw

  k1 = frame_axis(1)
  k2 = frame_axis(2)
  k3 = frame_axis(3)
  do n = 1,nframes
    do j = 1,3
      p1(j) = defpts(j,1,n)
      p2(j) = defpts(j,2,n)
      u(j) = frame(j,k1,n)
      v(j) = frame(j,k2,n)
      w(j) = frame(j,k3,n)
    enddo
    de_drotu = de_drotframe(k1,n)
    de_drotv = de_drotframe(k2,n)
    de_drotw = de_drotframe(k3,n)

    sizp1 = sqrt( p1(1)**2 + p1(2)**2 + p1(3)**2 )
    sizp2 = sqrt( p2(1)**2 + p2(2)**2 + p2(3)**2 )
    do j = 1,3
      p2unit(j) = p2(j) / sizp2
!     p1unit(j) = u(j) so no need to recalculate
    enddo
    dot21 = u(1)*p2(1) + u(2)*p2(2) + u(3)*p2(3)
    dot12 = p1(1)*p2unit(1) + p1(2)*p2unit(2) + p1(3)*p2unit(3)
    do j = 1,3
      p2perp1(j) = p2(j) - dot21*u(j)
      p1perp2(j) = p1(j) - dot12*p2unit(j)
    enddo
    sizp2perp1 = sqrt(p2perp1(1)**2 + p2perp1(2)**2 + p2perp1(3)**2)
    sizp1perp2 = sqrt(p1perp2(1)**2 + p1perp2(2)**2 + p1perp2(3)**2)
! def point one is along axis one. movement du parallel to that axis does
! not rotate the frame..so deriv is zero
!       dedu = 0.d0
! movement dv in v-axis direction corresponds to rotation about local w-axis
! of dtheta = dv/sizp1; thus a change in energy of 
!    dE = dedrotw*dtheta = dedrotw*dv/sizp1
!    dE/dv = dedrotw /sizp1
    dedv = de_drotw/sizp1
! movement dw in w-axis direction corresponds to rotation about p2unit
! of dtheta = -dw/sizp1perp2 (clockwise rotation) 
    dedw = -de_drotv/sizp1perp2
! Now convert to derivs wrt x,y,z. u = p.u = x*u(1)+y*u(2)+z*u(3)
! thus dudx = u(1). similaryl dvdx = v(1)
    de_ddefpt(1,1,n) = dedv*v(1) + dedw*w(1)
    de_ddefpt(2,1,n) = dedv*v(2) + dedw*w(2)
    de_ddefpt(3,1,n) = dedv*v(3) + dedw*w(3)
 
! for point 2..any movement in the local uv plane does not affect the frame
!       dedu = 0.d0
!       dedv = 0.d0
! movement dw in w direction corresponds to rotation about local u-axis 
! of dtheta = dw/sizp2perpu
    dedw = de_drotu/sizp2perp1
    de_ddefpt(1,2,n) = dedw*w(1)
    de_ddefpt(2,2,n) = dedw*w(2)
    de_ddefpt(3,2,n) = dedw*w(3)
  enddo      
  return
end subroutine ATMSITE_accum_de_ddefpts
!----------------------------------------------------------------------
subroutine ATMSITE_trans_de_ddefpts_to_centers(  &
                numfrdeflist,framedeflist,sitecrd,sitefrc,  &
                nframes,de_ddefpt)
  implicit none
  integer numfrdeflist,nframes
  integer framedeflist(5,numfrdeflist)
  double precision sitecrd(3,*),sitefrc(3,*),de_ddefpt(3,2,nframes)

      integer n,i,j,k,l,m,p
  double precision siz,siz2,dx,dy,dz,dux_dx,dux_dy,dux_dz, &
        duy_dy,duy_dz,duz_dz,dedux,deduy,deduz,dedx,dedy,dedz,  &
        siz3inv,sizinv
  do n = 1,numfrdeflist
    i = framedeflist(1,n)   ! site i gives tail of vector
    j = framedeflist(2,n)   ! site j gives head of vector
    k = framedeflist(3,n)   ! the frame number this vector helps define
    l = framedeflist(4,n)   ! the frame defpoint (1 or 2) within frame
    m = framedeflist(5,n)   ! m is number of such (unit) vectors 
                            ! averaged to get the frame defpoint k
    dx = sitecrd(1,j) - sitecrd(1,i)
    dy = sitecrd(2,j) - sitecrd(2,i)
    dz = sitecrd(3,j) - sitecrd(3,i)
    siz2 = dx*dx+dy*dy+dz*dz
    siz = sqrt(siz2)
    siz3inv = 1.d0/(siz2*siz)
    sizinv = 1.d0/siz
! ux, uy, uz are given by dx/siz, dy/siz, and dz/siz 
    dux_dx = sizinv - dx*dx*siz3inv
    dux_dy = -dx*dy*siz3inv   ! note duy_dx = dux_dy use this below
    dux_dz = -dx*dz*siz3inv
    duy_dy = sizinv - dy*dy*siz3inv
    duy_dz = -dy*dz*siz3inv
    duz_dz = sizinv - dz*dz*siz3inv
! the derivs of E wrt coordinates of unit vector in ij direction are given
! by (1/m) times the derivs of E wrt coords of def point (l,k)
! since the def point is the simple average of m of these unit vectors
    dedux = de_ddefpt(1,l,k) / m
    deduy = de_ddefpt(2,l,k) / m
    deduz = de_ddefpt(3,l,k) / m
! now apply chain rule, using symmetry e.g. dux_dy = duy_dx
    dedx = dedux*dux_dx + deduy*dux_dy + deduz*dux_dz
    dedy = dedux*dux_dy + deduy*duy_dy + deduz*duy_dz
    dedz = dedux*dux_dz + deduy*duy_dz + deduz*duz_dz
! finally apply forces. note force is negative of deriv of energy wrt position
! also note e.g. deriv of dx wrt x position of atoms i and j is -1,+1
    sitefrc(1,i) = sitefrc(1,i) + dedx
    sitefrc(2,i) = sitefrc(2,i) + dedy
    sitefrc(3,i) = sitefrc(3,i) + dedz
    sitefrc(1,j) = sitefrc(1,j) - dedx
    sitefrc(2,j) = sitefrc(2,j) - dedy
    sitefrc(3,j) = sitefrc(3,j) - dedz
  enddo

  return
end subroutine ATMSITE_trans_de_ddefpts_to_centers
!----------------------------------------------------------------------
subroutine ATMSITE_trans_de_dcens_to_atoms(numatoms,ncenlist,  &
                                       cen_list,cen_wt,sitefrc)
  implicit none
  integer numatoms,ncenlist,cen_list(2,*)
  double precision cen_wt(*),sitefrc(3,*)

  integer n,k,m
  double precision w
  do n = 1,ncenlist
    k = cen_list(1,n)  !contributing atom site
    m = cen_list(2,n)  !resulting extra center
    w = cen_wt(n)  ! weight of atom k contribution
    sitefrc(1,k) = sitefrc(1,k) + w * sitefrc(1,m)
    sitefrc(2,k) = sitefrc(2,k) + w * sitefrc(2,m)
    sitefrc(3,k) = sitefrc(3,k) + w * sitefrc(3,m)
  enddo
  return
end subroutine ATMSITE_trans_de_dcens_to_atoms
!-----------------------------------------------------------
subroutine ATMSITE_atom_to_sitecrds(numatoms,crd,sitecrd,atomsite)
  implicit none
  integer numatoms,atomsite(numatoms)
  double precision crd(3,*),sitecrd(3,*)
  integer n,m
  do n = 1,numatoms
    m = atomsite(n)  ! site this atom is mapped to
    sitecrd(1,m) = crd(1,n)
    sitecrd(2,m) = crd(2,n)
    sitecrd(3,m) = crd(3,n)
  enddo
  return
end subroutine ATMSITE_atom_to_sitecrds
!----------------------------------------------------------------------
subroutine ATMSITE_sitefrc_to_frc(numatoms,sitefrc,frc,atomsite)
  implicit none
  integer numatoms,atomsite(numatoms)
  double precision sitefrc(3,*),frc(3,*)
  integer n,m
  do n = 1,numatoms
    m = atomsite(n)  ! site this atom is mapped to
    frc(1,n) = frc(1,n) + sitefrc(1,m)
    frc(2,n) = frc(2,n) + sitefrc(2,m)
    frc(3,n) = frc(3,n) + sitefrc(3,m)
  enddo
  return
end subroutine ATMSITE_sitefrc_to_frc
!----------------------------------------------------------------------
subroutine ATMSITE_normalize_cenwts()
end subroutine ATMSITE_normalize_cenwts
