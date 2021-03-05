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
!--------------------------------------------------------------
subroutine AHFRAME_load_deflist(site_info,parmfile,readdens,debug)
use definition
  implicit none
  type(sites_info)::site_info
  integer allochk
  character(len=*) parmfile
  include "io.fh"
  include "site_names.inc"

  integer parmunit,nsites,k,n,mdef,slo,shi,num_deflist, &
          max_num,num_pt1,num_pt2,jsite,res,tot_def,tot_def2,j
  integer AHFRAME_name_to_site
  integer,allocatable:: tmp_p2(:),tmp_p3(:)
  double precision dx,dy,dz,wt,i,vdu,r1,r2,r3,tiny,u2(3),vec(3)
  double precision, allocatable::ua(:,:),ub(:,:)
  logical debug,readdens
  character(len=4) frame_pt_sname(10)
  parameter (tiny=1.d-15)
  
  parmunit = File_Open(parmfile,"r")
  if ( parmunit < 0 )then
    write(6,*)'failed to open parmfile ',parmfile(1:Len_Trim(parmfile)),parmunit
    stop
  endif
  nsites = site_info%nsites

  max_num = 10

  !first pass get num frame def list
  num_deflist = 0
  do n = 1,nsites
    call AHFRAME_get_def_names(sitename(n),parmunit,num_pt1,num_pt2,max_num, &
                       frame_pt_sname)
    !write(6,*)'n,num_pt1,num_pt2 = ',n,num_pt1,num_pt2
    !write(6,*)'frame_pt_sname = ',(frame_pt_sname(k),k=1,num_pt1+num_pt2)
    num_deflist = num_deflist + num_pt1 + num_pt2
  enddo
  write(6,*)'num_deflist = ',num_deflist
  site_info%num_frame_deflist = num_deflist

  allocate(site_info%frame_deflist(5*num_deflist), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHFRAME_load_deflist:could not allocate frame_deflist, exiting'
     stop
  endif

  tot_def = 0
  tot_def2 = 0
  num_deflist = 0
  do n = 1,nsites
    call AHFRAME_get_def_names(sitename(n),parmunit,num_pt1,num_pt2,max_num, &
                       frame_pt_sname)
    mdef = 0
    res = site_info%residue_num(n)
    slo = site_info%residue_start(res)
    shi = site_info%residue_end(res)
    if (n == site_info%molframe) then
       tot_def = num_pt1
       tot_def2 = num_pt2
       allocate(tmp_p2(num_pt1), tmp_p3(num_pt2), ua(num_pt1,3),&
                ub(num_pt2,3),stat = allochk)
       if (allochk .gt. 0 ) then
          write(6,*)'AHFRAME_load_deflist:could not allocate tmp_p2, exiting'
          stop
       endif
    endif
    do k = 1,num_pt1
      mdef = mdef + 1
      jsite = AHFRAME_name_to_site(frame_pt_sname(mdef),slo,shi)
      if ( jsite <= 0 )then
        write(6,*)'for PT1 n,k,slo,shi = ',n,k,slo,shi
        write(6,*)'no site found for ',frame_pt_sname(mdef)
        stop
      endif
      ! fill list item
      num_deflist = num_deflist +1
      site_info%frame_deflist(5*(num_deflist-1)+1) = n !tail of vector
      site_info%frame_deflist(5*(num_deflist-1)+2) = jsite !head of vector
      site_info%frame_deflist(5*(num_deflist-1)+3) = n ! frame number
      site_info%frame_deflist(5*(num_deflist-1)+4) = 1 ! frame point 1
      site_info%frame_deflist(5*(num_deflist-1)+5) = num_pt1 !num defining point
      if (n == site_info%molframe)tmp_p2(k) = jsite
    enddo
    do k = 1,num_pt2
      mdef = mdef + 1
      jsite = AHFRAME_name_to_site(frame_pt_sname(mdef),slo,shi)
      if ( jsite <= 0 )then
        write(6,*)'for PT2 n,k,slo,shi = ',n,k,slo,shi
        write(6,*)'no site found for ',frame_pt_sname(mdef)
        stop
      endif
      ! fill list item
      num_deflist = num_deflist +1
      site_info%frame_deflist(5*(num_deflist-1)+1) = n !tail of vector
      site_info%frame_deflist(5*(num_deflist-1)+2) = jsite !head of vector
      site_info%frame_deflist(5*(num_deflist-1)+3) = n ! frame number
      site_info%frame_deflist(5*(num_deflist-1)+4) = 2 ! frame point 1
      site_info%frame_deflist(5*(num_deflist-1)+5) = num_pt2 !num defining point
      if (n == site_info%molframe)tmp_p3(k) = jsite
    enddo
  enddo
  call File_Close(parmunit)
  !do n = 1,num_deflist
  !  write(6,5006)(site_info%frame_deflist(5*(n-1)+k),k=1,5)
  !enddo
!5006    format (T2,'frame_deflist :',i4)
  ! allocate frames
  allocate(site_info%frame_def_pts(3*2*nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHFRAME_load_deflist:could not allocate frame_def_pts, exiting'
     stop
  endif
  allocate(site_info%frames(3*3*nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHFRAME_load_deflist:could not allocate frames, exiting'
     stop
  endif
  allocate(site_info%fr_valid(nsites), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'AHFRAME_load_deflist:could not allocate fr_valid, exiting'
     stop
  endif
  site_info%fr_valid(:) = .false.

  if(readdens)then

     if (tot_def == 0 .or. tot_def2==0) then
        write(6,*)'AHFRAME_load_deflist: something wrong with rotation, exiting'
        stop
     endif

     ! now calculate the vectors that define molecular frame
   
     site_info%v1(:) = 0.0d0
     site_info%v2(:) = 0.0d0
     site_info%v3(:) = 0.0d0
     u2(:) = 0.0d0
   
     do j = 1, tot_def
        jsite = tmp_p2(j)
        do k = 1, 3 
           ua(j,k) = site_info%local_crds(3*(jsite-1)+k)
        enddo
        wt = sqrt(ua(j,1)*ua(j,1)+ua(j,2)*ua(j,2)+ua(j,3)*ua(j,3))
        if (wt .gt. tiny) then
           ua(j,:) = ua(j,:)/wt
        else
           ua(j,:) = 0.d0
        endif
        do k = 1, 3
           site_info%v1(k) = site_info%v1(k) + ua(j,k)
        enddo
     enddo
     wt = sqrt(site_info%v1(1)*site_info%v1(1)+site_info%v1(2)*site_info%v1(2)+&
               site_info%v1(3)*site_info%v1(3))
     if (wt .gt. tiny) then
        site_info%v1(:) = site_info%v1(:)/wt
     else
        site_info%v1(:) = 0.d0
     endif

     do j = 1, tot_def2
        jsite = tmp_p3(j)
        do k = 1, 3 
           ub(j,k) = site_info%local_crds(3*(jsite-1)+k)
        enddo
        wt = sqrt(ub(j,1)*ub(j,1)+ub(j,2)*ub(j,2)+ub(j,3)*ub(j,3))
        if (wt .gt. tiny) then
           ub(j,:) = ub(j,:)/wt
        else
           ub(j,:) = 0.d0
        endif
        do k = 1, 3
           u2(k) = u2(k) + ub(j,k)
        enddo
     enddo
     if(debug)write(12,*),'coords v1 = ',site_info%v1
   
     vdu = dot_product(site_info%v1,u2)
     do j = 1, 3
        site_info%v2(j) = u2(j) - vdu*site_info%v1(j)
     enddo
     wt = sqrt(site_info%v2(1)*site_info%v2(1)+site_info%v2(2)*site_info%v2(2)+&
              site_info%v2(3)*site_info%v2(3))
     if (wt .gt. tiny) then
        site_info%v2(:) = site_info%v2(:)/wt
     else
        site_info%v2(:) = 0.d0
     endif
     if(debug)write(12,*),'coords v2 = ',site_info%v2
   
     site_info%v3(1) = site_info%v1(2)*site_info%v2(3) - &
                       site_info%v1(3)*site_info%v2(2)
     site_info%v3(2) = site_info%v1(3)*site_info%v2(1) - &
                       site_info%v1(1)*site_info%v2(3)
     site_info%v3(3) = site_info%v1(1)*site_info%v2(2) - &
                       site_info%v1(2)*site_info%v2(1)
     if(debug)write(12,*),'coords v3 = ',site_info%v3
   
     !if (debug) then
     !   write(12,*)'----------------------------------------------'
     !   write(12,*)'LOCAL COORDINATES BEFORE'
     !   do j = 1, nsites
     !      write (12,*)(site_info%local_crds(3*(j-1)+k)*0.529177249d0,k=1,3)
     !   enddo
     !   write(12,*)'LOCAL COORDINATES BEFORE'
     !   write(12,*)'----------------------------------------------'
     !endif
   
     do j = 1, nsites
        do k = 1, 3
           vec(k) = site_info%local_crds(3*(j-1)+k)
        enddo
        r1 = dot_product(vec,site_info%v1)
        r2 = dot_product(vec,site_info%v2)
        r3 = dot_product(vec,site_info%v3)
        do k = 1, 3
           site_info%local_crds(3*(j-1)+k) = site_info%v1(k)*r1 + &
                                             site_info%v2(k)*r2 + &
                                             site_info%v3(k)*r3
        enddo
     enddo

     if (debug) then
        write(12,*)'----------------------------------------------'
        write(12,*)'LOCAL COORDINATES '
        do j = 1, nsites
           write (12,*)(site_info%local_crds(3*(j-1)+k)*0.529177249d0,k=1,3)
        enddo
        write(12,*)'LOCAL COORDINATES '
        write(12,*)'----------------------------------------------'
     endif

  endif !readdens

  return
end subroutine AHFRAME_load_deflist
!--------------------------------------------------------
subroutine AHFRAME_build_frames(site_info)
use definition
  implicit none
  type(sites_info)::site_info

  integer n,num_deflist,nsites,molframe
  integer frame_axis(3)

  nsites = site_info%nsites

  num_deflist = site_info%num_frame_deflist
  molframe = site_info%molframe

! clear the frame def pts
  if ( nsites > 0 )then
    do n = 1, 3*2*nsites
       site_info%frame_def_pts(n) = 0.d0
    enddo
  endif
  frame_axis(1) = 3 !z then x then y
  frame_axis(2) = 1 !z then x then y
  frame_axis(3) = 2 !z then x then y
  if ( num_deflist > 0 )then
    ! first the def points
    call ATMSITE_build_frame_def_pts( &
       num_deflist,site_info%frame_deflist,site_info%site_crds,nsites,&
       site_info%frame_def_pts,site_info%fr_valid)
    ! next complete the frames
    call ATMSITE_defpoints_to_frames(nsites,frame_axis, &
               site_info%frame_def_pts,site_info%frames,molframe,&
               site_info%fr_valid)
  endif
  molframe=0
!  do n = 1, 3*nsites
!     molframe = molframe+3
!     write(6,5006)site_info%frames(molframe-2),site_info%frames(molframe-1),&
!                  site_info%frames(molframe)
!  enddo
!5006    format (T2,'frames ',f20.17,f20.17,f20.17)
  return
end subroutine AHFRAME_build_frames
!--------------------------------------------------------
subroutine AHFRAME_rotate_hermites(site_info,lhermitefile)
use definition
  implicit none
  type(sites_info)::site_info
  character(len=*) lhermitefile
  include "io.fh"
  include "site_names.inc"

  integer lhermite_unit
  integer nsites,k,m,n,scoff,dimxy,num,off,order,j
  integer ncoeff
  double precision At(3,3),Mpole_xy(20,20)

  nsites = site_info%nsites

  ncoeff = site_info%num_coefficients

  do n = 1,nsites
    num = site_info%num_primitives(n)
    off = site_info%off_primitives(n)
    call VEC3D_3x3_transpose(site_info%frames((1+(9*(n-1))):(9*n)),At)
    !order = 10 ! either 1 or 10---here use 10 to get Mpole_xy
    !call XFORM_MPOLE_matrix(At,Mpole_xy,order)
    do k = 1,num
      order = site_info%prim_order(off+k)
      if (order .gt. 20) then
         write(6,*)'Can only transform hermites up to l=3, exiting'
         stop
      endif
      if(order == 6)order=10
      call XFORM_MPOLE_matrix(At,Mpole_xy,order)
      scoff = site_info%coeff_offset(off+k)
      dimxy = order
      if ( order == 1 )then ! simply copy s-orbital coeff: local=global
        site_info%local_hermite_coeffs(scoff+1) = &
                  site_info%global_hermite_coeffs(scoff+1)
      elseif ( order > 1 .and. .not.(site_info%fr_valid(n)))then 
      ! frame is not valid!
        site_info%global_hermite_coeffs(scoff+1) = &
                  site_info%local_hermite_coeffs(scoff+1)
      elseif ( order > 1 .and. site_info%fr_valid(n))then
        call XFORM_MPOLE(Mpole_xy,dimxy,&
                   site_info%global_hermite_coeffs(scoff+1:scoff+1+order), &
                   site_info%local_hermite_coeffs(scoff+1:scoff+1+order),order)
      endif
      !print *,site_info%global_hermite_coeffs(scoff+1),site_info%local_hermite_coeffs(scoff+1)
    enddo
  enddo
  !lhermite_unit = File_Open(lhermitefile,"w")
  !if ( lhermite_unit < 0 )then
  !  write(6,*)'failed to open lhermite file ',  &
  !         lhermitefile(1:Len_Trim(lhermitefile))
  !  stop
  !endif
  lhermite_unit = 99
  open(lhermite_unit,file=lhermitefile,status="unknown")
  write(lhermite_unit,*)ncoeff
  do n = 1,ncoeff
     write(lhermite_unit,*)site_info%local_hermite_coeffs(n)
     !write(lhermite_unit,5006)site_info%local_hermite_coeffs(n)
  enddo
  close(lhermite_unit)
!5006    format (T2,f13.10)
  return
end subroutine AHFRAME_rotate_hermites
!--------------------------------------------------------
subroutine AHFRAME_rotate_Jvec(site_info,Jvec,ncoeff,debug)
use definition
  implicit none
  type(sites_info)::site_info

  integer nsites,k,m,n,scoff,dimxy,num,off,order
  integer ncoeff, counter
  double precision At(3,3),Mpole_xy(10,10)
  double precision Jvec(ncoeff)
  double precision Jvec_new(ncoeff)
  logical debug

  nsites = site_info%nsites

  Jvec_new(:) = 0.0d0
  counter = 0
  do n = 1,nsites
    num = site_info%num_primitives(n)
    off = site_info%off_primitives(n)
    call VEC3D_3x3_transpose(site_info%frames((1+(9*(n-1))):(9*n)),At)
    order = 10 ! either 1 or 10---here use 10 to get Mpole_xy
    call XFORM_MPOLE_matrix(At,Mpole_xy,order)
    do k = 1,num
      order = site_info%prim_order(off+k)
      dimxy = order
      if ( order == 1 )then ! simply copy s-orbital coeff: local=global
        counter = counter + order
        Jvec_new(counter) = Jvec(counter)
      elseif ( order > 1 )then
        counter = counter + order
        call XFORM_MPOLE(Mpole_xy,dimxy,Jvec(counter-order:counter), &
                         Jvec_new(counter-order:counter),order)
      endif
    enddo
  enddo
  Jvec(:) = Jvec_new(:)
  if (debug) then
     do n = 1,ncoeff
        write(12,1001)n,Jvec_new(n)
     enddo
  endif
1001    format ('rotated Jvec (',i4,') = ',f27.20)
  return
end subroutine AHFRAME_rotate_Jvec
!--------------------------------------------------------
subroutine AHFRAME_local_to_global(site_info)
use definition
  implicit none
  type(sites_info)::site_info
  include "io.fh"

  integer nsites, i,j
  integer k,n,num,off,order,scoff,dimxy,n_dip_sites,l
  double precision pi, expon
  double precision Mpole_xy(20,20)
  logical ind_dip,readstat

  n_dip_sites=0
  pi = 4.d0*(atan(1.0d0))

  ! check for induced dipoles in file
  inquire(file = 'ind_dip', exist = ind_dip)
  if (ind_dip) then
     open(66,file='ind_dip',status='unknown')
       read(66,'(i5)') n_dip_sites
       !if (readstat.gt.0)then
       !   print *,'something wrong reading ind_dip file, exiting'
       !   stop
       !endif
     close(66)
  endif
  nsites = site_info%nsites-n_dip_sites

  do n = 1,nsites
    num = site_info%num_primitives(n)
    off = site_info%off_primitives(n)
    ! TRANSPOSE FOR MULTIPOLES IN AMOEBA, HAVE TO TEST FOR LOCAL HERMITES!
    order = 10 ! either 1 or 10---here use 10 to get Mpole_xy
    call XFORM_MPOLE_matrix(site_info%frames((1+(9*(n-1))):(9*n)),&
                            Mpole_xy,order)
    !print *,site_info%frames(1+(9*(n-1))),site_info%frames(4+(9*(n-1))),&
    !&site_info%frames(7+(9*(n-1)))
    !print *,site_info%frames(2+(9*(n-1))),site_info%frames(5+(9*(n-1))),&
    !&site_info%frames(8+(9*(n-1)))
    !print *,site_info%frames(3+(9*(n-1))),site_info%frames(6+(9*(n-1))),&
    !&site_info%frames(9+(9*(n-1)))
    !do i = 1,order
    ! do j = 1, order
    ! print *,Mpole_xy(i,j)
    ! enddo
    !enddo
    !print *,''
    do k = 1,num
      order = site_info%prim_order(off+k)
      if(order == 6)order=10
      scoff = site_info%coeff_offset(off+k)
      !if(order .le. 10) then
         dimxy = order
      !else
      !   dimxy = 10
      !endif
      if ( order == 1 )then ! simply copy s-orbital coeff: local=global
        site_info%global_hermite_coeffs(scoff+1) = &
                  site_info%local_hermite_coeffs(scoff+1)
      elseif ( order > 1 .and. .not.(site_info%fr_valid(n)))then 
      ! frame is not valid!
        site_info%global_hermite_coeffs(scoff+1) = &
                  site_info%local_hermite_coeffs(scoff+1)
      elseif ( order > 1 .and. site_info%fr_valid(n))then 
!      elseif ( order > 1 )then
        call XFORM_MPOLE(Mpole_xy,dimxy,&
                         site_info%local_hermite_coeffs(scoff+1), &
                         site_info%global_hermite_coeffs(scoff+1),order)
      endif
      !do j = 1, 10
      !   print *,site_info%global_hermite_coeffs(scoff+1+j-1),&
      !          &site_info%local_hermite_coeffs(scoff+1+j-1)
      !enddo
      !print *,' '
    enddo
  enddo
  
  if(ind_dip)then
     do n = nsites+1,nsites+n_dip_sites
       num = site_info%num_primitives(n)
       off = site_info%off_primitives(n)
       do k = 1,num
         order = site_info%prim_order(off+k)
         scoff = site_info%coeff_offset(off+k)
         expon = site_info%prim_expo(off+k)
         do l = 1,order
         !do l = k,k+2
            site_info%global_hermite_coeffs(scoff+l) = &
                        site_info%local_hermite_coeffs(scoff+l) !&
                        !*(1.d0/(2.d0*expon))
                        !*(pi*sqrt(pi))/(2.542*2.d0*expon**2.5)
                        !*(pi/expon)*sqrt(pi/expon)*(1.d0/(2.d0*expon))&
                        !*site_info%cart_coeff_norms(scoff+l)
                        !*(expon/pi)*sqrt(expon/pi)!*(1.d0/(2.d0*expon))&
                        !*site_info%cart_coeff_norms(scoff+l)
            !print *,'dip = ',site_info%global_hermite_coeffs(scoff+l),&
            !                 site_info%local_hermite_coeffs(scoff+l)
         enddo
       enddo
     enddo
  endif
  return
end subroutine AHFRAME_local_to_global
!--------------------------------------------------------
subroutine AHFRAME_get_def_names(sname,parmunit,num_pt1,num_pt2,max_num, &
                       frame_pt_sname)
  implicit none
  character(len=4)sname
  integer parmunit,num_pt1,num_pt2,max_num
  character(len=4)frame_pt_sname(max_num)

  integer ios,ptr,j,k
  character(len=120) line,word
  integer TRDPRM_get_next_line,lenword
  rewind(parmunit)
  ios = 0
  do while (ios==0)
    line = ' '  ! to avoid junk in line after some reads
    ios = TRDPRM_get_next_line(parmunit,line)
    if ( ios == 0 )then
      ptr = 1 ! first call to str_next_token needs an initial value
      call str_next_token(line,ptr,word,lenword) !get 1st token
      if ( ptr > 0 )then !check if its 'frame'
        if ( word(1:lenword) == 'frame' )then
          call str_next_token(line,ptr,word,lenword) !get 1st token
          if ( sname(1:Len_trim(sname)) == word(1:lenword) )then
            read(line(ptr:),*)num_pt1,num_pt2
            if ( num_pt1 + num_pt2 > max_num )then
              write(6,*)'too many site names in frame def for ',sname
              stop
            endif
            ! read to start of names
            call str_next_token(line,ptr,word,lenword)
            call str_next_token(line,ptr,word,lenword)
            do k = 1,num_pt1+num_pt2
              call str_next_token(line,ptr,word,lenword)
              if ( lenword > 4 )then
                write(6,*)'name too big'
                stop
              endif
              frame_pt_sname(k)(1:lenword) = word(1:lenword)
              do j = lenword+1,4
                frame_pt_sname(k)(j:j) = ' '
              enddo
            enddo
          endif
        endif
      endif !ptr > 0
    endif!ios == 0 )
  enddo !ios==0

  return
end subroutine AHFRAME_get_def_names
!--------------------------------------------------------
integer function AHFRAME_name_to_site(sname,slo,shi)
  implicit none
  character(len=*)sname
  integer slo,shi
  include "site_names.inc"
  integer n
  do n = slo,shi
    if ( sname(1:Len_Trim(sname)) == sitename(n)(1:Len_Trim(sitename(n))) )then
      AHFRAME_name_to_site = n
      return
    endif
  enddo
  AHFRAME_name_to_site = -1
  return
end function AHFRAME_name_to_site
!--------------------------------------------------------
subroutine AHFRAME_get_molframe(site_info,parmfile)
use definition
  implicit none
  type(sites_info)::site_info
  integer allochk
  character(len=*) parmfile
  include "io.fh"
  include "site_names.inc"

  integer ios,ptr,j,k
  character(len=120) line,word
  integer TRDPRM_get_next_line,lenword
  integer parmunit,nsites,n,mdef,slo,shi,num_deflist, &
          max_num,num_pt1,num_pt2,jsite,res
  integer AHFRAME_name_to_site
  
  parmunit = File_Open(parmfile,"r")
  if ( parmunit < 0 )then
    write(6,*)'failed to open parmfile ',parmfile(1:Len_Trim(parmfile))
    stop
  endif
  nsites = site_info%nsites

  site_info%molframe = 0

  rewind(parmunit)
  ios = 0
  do while (ios==0)
    line = ' '  ! to avoid junk in line after some reads
    ios = TRDPRM_get_next_line(parmunit,line)
    if ( ios == 0 )then
      ptr = 1 ! first call to str_next_token needs an initial value
      call str_next_token(line,ptr,word,lenword) !get 1st token
      if ( ptr > 0 )then !check if its 'molframe'
        if ( word(1:lenword) == 'molframe' )then
          call str_next_token(line,ptr,word,lenword) !get 1st token
          do j = 1, nsites
            if ( sitename(j)(1:Len_trim(sitename(j))) == word(1:lenword) )then
               site_info%molframe = j
               do k = 1, 3
                  site_info%p1(k) = site_info%site_crds(3*(j-1)+k)
               enddo
            endif
          enddo
        endif
      endif !ptr > 0
    endif!ios == 0 )
  enddo !ios==0

  do j = 1, nsites
     do k = 1, 3
        site_info%local_crds(3*(j-1)+k) = site_info%site_crds(3*(j-1)+k) - &
                                          site_info%p1(k)
        !site_info%local_crds(3*(j-1)+k) = site_info%site_crds(3*(j-1)+k)
     enddo
  enddo

  call File_Close(parmunit)

  if (site_info%molframe == 0) then
     write (6,*) 'AHFRAME_get_molframe: could not find molframe, exiting'
     stop
  endif

  return
end subroutine AHFRAME_get_molframe
!--------------------------------------------------------
subroutine AHFRAME_rotate_coords(site_info, debug)
! subroutine to rotate residue to molecular frame by alining v1 along x
! and v2 along y.
use definition
  implicit none
  type(sites_info)::site_info
  integer allochk

  integer ios,ptr,i,j,k,nsite
  character(len=120) line,word
  integer parmunit,nsites,n,mdef,slo,shi,num_deflist, &
          max_num,num_pt1,num_pt2,jsite,res
  double precision :: a, b, c, d, dx, dy, dz, Ux, Uy, Uz,&
                      angle, v1bar, mag, ca, sa, tiny
  double precision, dimension(3) :: V, U, vec, ref
  double precision, dimension(3,3):: M
  logical debug
  parameter (tiny = 1.d-15)

  nsite = site_info%nsites

  if (debug) then
     write(12,*)'LOCAL COORDINATES BEFORE'
     do i = 1, nsite
        write (12,*)(site_info%local_crds(3*(i-1)+k)*0.529177249d0,k=1,3)
     enddo
     write(12,*)'LOCAL COORDINATES BEFORE'
  endif

  ! if v1 is along Z axis and v2 along X axis we're done
  if (site_info%v1(3) .eq. 1.0d0 .and. site_info%v2(1) .eq. 1.0d0 &
      .and. site_info%v3(2) .eq. 1.0d0) return 
     
  !align v1 to z

  ref(1) = 0.0d0
  ref(2) = 0.0d0
  ref(3) = 1.0d0

  V(1) = site_info%v1(2)*ref(3) - site_info%v1(3)*ref(2)
  V(2) = site_info%v1(3)*ref(1) - site_info%v1(1)*ref(3)
  V(3) = site_info%v1(1)*ref(2) - site_info%v1(2)*ref(1)

  v1bar = sqrt(site_info%v1(1)*site_info%v1(1) + &
               site_info%v1(2)*site_info%v1(2) + &
               site_info%v1(3)*site_info%v1(3))

  angle = acos(dot_product(site_info%v1,ref))

  mag = v1bar*sin(angle)

  if (mag .gt. tiny) then
     V(:) = V(:)/mag
  else
     V(:) = 0.d0
  endif

  a = V(1)
  b = V(2)
  c = V(3)
  
  ca = cos(angle)
  sa = sin(angle)
  
  M = reshape((/ca+a*a*(1-ca), (a*b)*(1-ca)+c*sa, (a*c)*(1-ca)-b*sa,&
               (a*b)*(1-ca)-c*sa, ca+b*b*(1-ca), (b*c)*(1-ca)+a*sa,&
               (a*c)*(1-ca)+b*sa, (b*c)*(1-ca)-a*sa, ca+c*c*(1-ca)/),&
               (/3,3/))

  do i = 1, nsite
     do k = 1, 3
        vec(k) = site_info%local_crds(3*(i-1)+k) 
     enddo
     vec = matmul(M,vec)
     do k = 1, 3
        site_info%local_crds(3*(i-1)+k) = vec(k)
     enddo
  enddo

  site_info%v1 = matmul(M,site_info%v1)
  site_info%v2 = matmul(M,site_info%v2)
  site_info%v3 = matmul(M,site_info%v3)

  ! if v1 is along Z axis and v2 along X axis we're done
  if (site_info%v1(3) .eq. 1.0d0 .and. site_info%v2(1) .eq. 1.0d0 &
      .and. site_info%v3(2) .eq. 1.0d0) return 
     
  !align v2 to x

  ref(1) =  1.0d0
  ref(2) =  0.0d0
  ref(3) =  0.0d0

  V(1) = site_info%v2(2)*ref(3) - site_info%v2(3)*ref(2)
  V(2) = site_info%v2(3)*ref(1) - site_info%v2(1)*ref(3)
  V(3) = site_info%v2(1)*ref(2) - site_info%v2(2)*ref(1)

  v1bar = sqrt(site_info%v2(1)*site_info%v2(1) + &
               site_info%v2(2)*site_info%v2(2) + &
               site_info%v2(3)*site_info%v2(3))

  angle = acos(dot_product(site_info%v2,ref))

  mag = v1bar*sin(angle)

  if (mag .gt. tiny) then
     V(:) = V(:)/mag
  else
     V(:) = 0.d0
  endif

  a = V(1)
  b = V(2)
  c = V(3)
  
  ca = cos(angle)
  sa = sin(angle)
  
  M = reshape((/ca+a*a*(1-ca), (a*b)*(1-ca)+c*sa, (a*c)*(1-ca)-b*sa,&
               (a*b)*(1-ca)-c*sa, ca+b*b*(1-ca), (b*c)*(1-ca)+a*sa,&
               (a*c)*(1-ca)+b*sa, (b*c)*(1-ca)-a*sa, ca+c*c*(1-ca)/),&
               (/3,3/))

  do i = 1, nsite
     do k = 1, 3
        vec(k) = site_info%local_crds(3*(i-1)+k) 
     enddo
     vec = matmul(M,vec)
     do k = 1, 3
        site_info%local_crds(3*(i-1)+k) = vec(k)
     enddo
  enddo

  site_info%v1 = matmul(M,site_info%v1)
  site_info%v2 = matmul(M,site_info%v2)
  site_info%v3 = matmul(M,site_info%v3)

  ! check to make sure it's in the positive quadrant

  !if (site_info%v1(1) .lt. 0.0d0) then
  !   do i = 1, nsite
  !      site_info%local_crds(3*(i-1)+1) = -1.d0*site_info%local_crds(3*(i-1)+1)
  !   enddo
  !else if (site_info%v2(2) .lt. 0.0d0) then
  !   do i = 1, nsite
  !      site_info%local_crds(3*(i-1)+2) = -1.d0*site_info%local_crds(3*(i-1)+2)
  !   enddo
  !else if (site_info%v3(3) .lt. 0.0d0) then
  !   do i = 1, nsite
  !      site_info%local_crds(3*(i-1)+3) = -1.d0*site_info%local_crds(3*(i-1)+3)
  !   enddo
  !endif

  if (debug) then
     write(12,*)'----------------------------------------------'
     write(12,*)'LOCAL COORDINATES AFTER'
     do i = 1, nsite
        write (12,*)(site_info%local_crds(3*(i-1)+k)*0.529177249d0,k=1,3)
     enddo
     write(12,*)'LOCAL COORDINATES AFTER'
     write(12,*)'----------------------------------------------'
  endif

  return
end subroutine AHFRAME_rotate_coords
!--------------------------------------------------------
