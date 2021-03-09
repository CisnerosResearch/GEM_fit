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
subroutine GMcMur_Dav_tabsiz(topdeg,itoptab)
  implicit none
  integer topdeg,itoptab
  integer ic,t,u,v,n
! go through stages of recursion
! first R(0,0,0,j)
  ic = 0
  do n = 0,topdeg
    ic = ic + 1
  enddo
! next R(0,0,v,j)
  do v = 1,topdeg
    do n = 0,topdeg - v
      ic = ic + 1
    enddo
  enddo
! next R(0,u,v,j)
  do u = 1,topdeg
    do v = 0,topdeg - u
      do n = 0,topdeg - (u + v)
        ic = ic + 1
      enddo
    enddo
  enddo
! next R(t,u,v,j)
  do t = 1,topdeg
    do u = 0,topdeg - t
      do v = 0,topdeg - (t + u)
        do n = 0,topdeg - (t + u + v)
          ic = ic + 1
        enddo
      enddo
    enddo
  enddo
  itoptab = ic

  return
end subroutine GMcMur_Dav_tabsiz
!-------------------------------------------------------
subroutine GMcMur_Dav_ptrs(topdeg,itoptab,  &
                          ptr1,ptr2,tuv_map,  &
                          indt,field_order1,mpole_order2, &
                          tstart,ustart,vstart)
  implicit none
  integer topdeg,itoptab
  integer ptr1(itoptab),ptr2(itoptab),tuv_map(itoptab)
  integer field_order1,mpole_order2
  integer indt(field_order1,mpole_order2)
  integer tstart,ustart,vstart

 include "mpole_sizes.fh"
  integer ic,t,u,v,n,k1,k2
  integer t_part(MAXMPG),u_part(MAXMPG),v_part(MAXMPG)
  integer scratch(0:MAXJ,0:MAXJ,0:MAXJ,0:MAXJ)

! go through stages of recursion to get index map scratch array
! first R(0,0,0,n)
  ic = 0
  do n = 0,topdeg
    ic = ic + 1
    scratch(0,0,0,n) = ic
  enddo
! next R(0,0,v,n)
  do v = 1,topdeg
    do n = 0,topdeg - v
      ic = ic + 1
      scratch(0,0,v,n) = ic
    enddo
  enddo
! next R(0,u,v,n)
  do u = 1,topdeg
    do v = 0,topdeg - u
      do n = 0,topdeg - (u + v)
        ic = ic + 1
        scratch(0,u,v,n) = ic
      enddo
    enddo
  enddo
! next R(t,u,v,n)
  do t = 1,topdeg
    do u = 0,topdeg - t
      do v = 0,topdeg - (t + u)
        do n = 0,topdeg - (t + u + v)
          ic = ic + 1
          scratch(t,u,v,n) = ic
        enddo
      enddo
    enddo
  enddo
  itoptab = ic
! now do it again to get recursion pointers
! clear the arrays
  do ic = 1,itoptab
    ptr1(ic) = 0
    ptr2(ic) = 0
    tuv_map(ic) = 0
  enddo
! first R(0,0,0,n)
  ic = 0
  do n = 0,topdeg
    ic = ic + 1
    tuv_map(ic) = 0
  enddo
! next R(0,0,v,n)
! recursion will be R(0,0,v,n) = (v-1)*R(0,0,v-2,n+1)+Z*R(0,0,v-1,n+1)
  vstart = ic + 1
  do v = 1,topdeg
    do n = 0,topdeg - v
      ic = ic + 1
      tuv_map(ic) = v
      ptr1(ic) = scratch(0,0,v-1,n+1)
      if ( v .ge. 2 )then
        ptr2(ic) = scratch(0,0,v-2,n+1)
      else 
        ptr2(ic) = 0
      endif
    enddo
  enddo
! next R(0,u,v,n)
! recursion will be R(0,u,v,n) = (u-1)*R(0,u-2,v,n+1)+Y*R(0,u-1,v,n+1)
  ustart = ic + 1
  do u = 1,topdeg
    do v = 0,topdeg - u
      do n = 0,topdeg - (u + v)
        ic = ic + 1
        tuv_map(ic) = u
        ptr1(ic) = scratch(0,u-1,v,n+1)
        if ( u .ge. 2 )then
          ptr2(ic) = scratch(0,u-2,v,n+1)
        else 
          ptr2(ic) = 0
        endif
      enddo
    enddo
  enddo
! next R(t,u,v,n)
! recursion will be R(t,u,v,n) = (t-1)*R(t-2,u,v,n+1)+X*R(t-1,u,v,n+1)
  tstart = ic + 1
  do t = 1,topdeg
    do u = 0,topdeg - t
      do v = 0,topdeg - (t + u)
        do n = 0,topdeg - (t + u + v)
          ic = ic + 1
          tuv_map(ic) = t
          ptr1(ic) = scratch(t-1,u,v,n+1)
          if ( t .ge. 2 )then
            ptr2(ic) = scratch(t-2,u,v,n+1)
          else 
            ptr2(ic) = 0
          endif
        enddo
      enddo
    enddo
  enddo
  itoptab = ic
! now get pointers for field accum
  call GMcMur_Dav_get_tuv_parts(t_part,u_part,v_part)
  do k2 = 1,mpole_order2
    do k1 = 1,field_order1
      t = t_part(k1) + t_part(k2)
      u = u_part(k1) + u_part(k2)
      v = v_part(k1) + v_part(k2)
      indt(k1,k2) = scratch(t,u,v,0)
    enddo
  enddo
  return
end subroutine GMcMur_Dav_ptrs
!----------------------------------------------------------
subroutine GMcMur_Dav_get_tuv_parts(t_part,u_part,v_part)
  implicit none
  integer t_part(*),u_part(*),v_part(*)
 include "mpole_index.fh"

! charge; no derivs
  t_part(IND_000) = 0; u_part(IND_000) = 0; v_part(IND_000) = 0
! dipoles
  t_part(IND_100) = 1; u_part(IND_100) = 0; v_part(IND_100) = 0
  t_part(IND_010) = 0; u_part(IND_010) = 1; v_part(IND_010) = 0
  t_part(IND_001) = 0; u_part(IND_001) = 0; v_part(IND_001) = 1
! quadrupoles
  t_part(IND_200) = 2; u_part(IND_200) = 0; v_part(IND_200) = 0
  t_part(IND_020) = 0; u_part(IND_020) = 2; v_part(IND_020) = 0
  t_part(IND_002) = 0; u_part(IND_002) = 0; v_part(IND_002) = 2
  t_part(IND_110) = 1; u_part(IND_110) = 1; v_part(IND_110) = 0
  t_part(IND_101) = 1; u_part(IND_101) = 0; v_part(IND_101) = 1
  t_part(IND_011) = 0; u_part(IND_011) = 1; v_part(IND_011) = 1
! octupoles
  t_part(IND_300) = 3; u_part(IND_300) = 0; v_part(IND_300) = 0
  t_part(IND_030) = 0; u_part(IND_030) = 3; v_part(IND_030) = 0
  t_part(IND_003) = 0; u_part(IND_003) = 0; v_part(IND_003) = 3
  t_part(IND_210) = 2; u_part(IND_210) = 1; v_part(IND_210) = 0
  t_part(IND_201) = 2; u_part(IND_201) = 0; v_part(IND_201) = 1
  t_part(IND_120) = 1; u_part(IND_120) = 2; v_part(IND_120) = 0
  t_part(IND_021) = 0; u_part(IND_021) = 2; v_part(IND_021) = 1
  t_part(IND_102) = 1; u_part(IND_102) = 0; v_part(IND_102) = 2
  t_part(IND_012) = 0; u_part(IND_012) = 1; v_part(IND_012) = 2
  t_part(IND_111) = 1; u_part(IND_111) = 1; v_part(IND_111) = 1
! hexadecapoles
  t_part(IND_400) = 4; u_part(IND_400) = 0; v_part(IND_400) = 0
  t_part(IND_040) = 0; u_part(IND_040) = 4; v_part(IND_040) = 0
  t_part(IND_004) = 0; u_part(IND_004) = 0; v_part(IND_004) = 4
  t_part(IND_310) = 3; u_part(IND_310) = 1; v_part(IND_310) = 0
  t_part(IND_301) = 3; u_part(IND_301) = 0; v_part(IND_301) = 1
  t_part(IND_130) = 1; u_part(IND_130) = 3; v_part(IND_130) = 0
  t_part(IND_031) = 0; u_part(IND_031) = 3; v_part(IND_031) = 1
  t_part(IND_103) = 1; u_part(IND_103) = 0; v_part(IND_103) = 3
  t_part(IND_013) = 0; u_part(IND_013) = 1; v_part(IND_013) = 3
  t_part(IND_220) = 2; u_part(IND_220) = 2; v_part(IND_220) = 0
  t_part(IND_202) = 2; u_part(IND_202) = 0; v_part(IND_202) = 2
  t_part(IND_022) = 0; u_part(IND_022) = 2; v_part(IND_022) = 2
  t_part(IND_211) = 2; u_part(IND_211) = 1; v_part(IND_211) = 1
  t_part(IND_121) = 1; u_part(IND_121) = 2; v_part(IND_121) = 1
  t_part(IND_112) = 1; u_part(IND_112) = 1; v_part(IND_112) = 2
!  5th order
  t_part(IND_500) = 5; u_part(IND_500) = 0; v_part(IND_500) = 0
  t_part(IND_050) = 0; u_part(IND_050) = 5; v_part(IND_050) = 0
  t_part(IND_005) = 0; u_part(IND_005) = 0; v_part(IND_005) = 4
  t_part(IND_410) = 4; u_part(IND_410) = 1; v_part(IND_410) = 0
  t_part(IND_401) = 4; u_part(IND_401) = 0; v_part(IND_401) = 1
  t_part(IND_140) = 1; u_part(IND_140) = 4; v_part(IND_140) = 0
  t_part(IND_041) = 0; u_part(IND_041) = 4; v_part(IND_041) = 1
  t_part(IND_104) = 1; u_part(IND_104) = 0; v_part(IND_104) = 4
  t_part(IND_014) = 0; u_part(IND_014) = 1; v_part(IND_014) = 4
  t_part(IND_320) = 3; u_part(IND_320) = 2; v_part(IND_320) = 0
  t_part(IND_302) = 3; u_part(IND_302) = 0; v_part(IND_302) = 2
  t_part(IND_230) = 2; u_part(IND_230) = 3; v_part(IND_230) = 0
  t_part(IND_032) = 0; u_part(IND_032) = 3; v_part(IND_032) = 2
  t_part(IND_203) = 2; u_part(IND_203) = 0; v_part(IND_203) = 3
  t_part(IND_023) = 0; u_part(IND_023) = 2; v_part(IND_023) = 3
  t_part(IND_311) = 3; u_part(IND_311) = 1; v_part(IND_311) = 1
  t_part(IND_131) = 1; u_part(IND_131) = 3; v_part(IND_131) = 1
  t_part(IND_113) = 1; u_part(IND_113) = 1; v_part(IND_113) = 3
  t_part(IND_221) = 2; u_part(IND_221) = 2; v_part(IND_221) = 1
  t_part(IND_212) = 2; u_part(IND_212) = 1; v_part(IND_212) = 2
  t_part(IND_122) = 1; u_part(IND_122) = 2; v_part(IND_122) = 2
   
  return
end subroutine GMcMur_Dav_get_tuv_parts
