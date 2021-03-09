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
#include "dprec.fh"
!------------------------------------------------------
subroutine ATMSITE_build_sites()
  implicit none
#include "database.fh"
#include "atoms_to_sites.inc"

  integer num_centers

! Get necessary ptrs and scalars---needed also by ATMSITE_torque_to_force 
  call ATMSITE_unpack()

! clear the sitecrd array
  call UTIL_zero_real_array(RH(p_sitecrd),3*nsites)
! clear the frame def pts
  if ( nframes .gt. 0 )then
    call UTIL_zero_real_array(RH(p_frdefpts),3*2*nframes)
  endif
! first copy the atomic coords into site coord array
  call ATMSITE_atom_to_sitecrds(numatoms,RH(p_crd),RH(p_sitecrd),  &
                                IH(p_atmsite_ind))
! next define the centers..weighted averages of atomic positions
  if ( nextra_cen > 0 )then
        call ATMSITE_build_centers(numatoms,nextra_cen,ncenlist, &
                               IH(p_cen_list),RH(p_cen_wt),RH(p_sitecrd))
      endif
! next define frames in terms of centers
  if ( nfrdeflist > 0 )then
    ! first the def points
    call ATMSITE_build_frame_def_pts( &
       nfrdeflist,IH(p_frdeflst),RH(p_sitecrd),nframes,RH(p_frdefpts))
    ! next complete the frames
    call ATMSITE_defpoints_to_frames(nframes,frame_axis, &
               RH(p_frdefpts),RH(p_frames))
  endif
! next fill out the sites of the asymmetric unit
  if ( nextra_sites > 0 )then
    num_centers = numatoms + nextra_cen
    call ATMSITE_build_extrasites(num_centers,nextra_sites, &
                              IH(p_extsitelst),IH(p_indframe), &
                              RH(p_locvec),RH(p_frames),RH(p_sitecrd))
  endif
! next fill out the unit cell using symmetry operators
! to be done...
! next rotate the local multipoles into global frame
  call ATMSITE_rotate_mpoles( &
             nsites,IH(p_indframe),IH(p_Mpole_order),  &             !input
                 IH(p_Mpole_offset),                   &             !input
                 RH(p_frames),RH(p_Mpole_Local),  &
                 RH(p_Mpole_Global))
 ! call ATM_dump_poles(  &
                 !nsites,IH(p_Mpole_order),IH(p_Mpole_offset),  &
                 !RH(p_Mpole_Local),RH(p_Mpole_Global))
  return
end subroutine ATMSITE_build_sites
!-----------------------------------------------------------------
subroutine ATMSITE_torque_to_force()
  implicit none
#include "database.fh"
#include "atoms_to_sites.inc"
  integer num_centers

  call ATMSITE_get_de_drot_mpole( &
            nsites,RH(p_field_xyz),RH(p_mpole_xyz),RH(p_de_drotsite), &
            IH(p_mpole_offset),IH(p_mpole_order),coulomb_constant)
  num_centers = numatoms + nextra_cen
  call ATMSITE_add_extra_site_contrib( &
        num_centers,nextra_sites,IH(p_extsitelst),RH(p_sitecrd), &
        RH(p_sitefrc),RH(p_de_drotsite))
  call ATMSITE_accum_de_dframe_rot( &
                 nsites,nframes,IH(p_indframe),frame_axis,  &  !input 
                 RH(p_de_drotsite),RH(p_frames),   &                 !intput
                 RH(p_frdefpts),     &
                 RH(p_de_drotframe))                              !output
  call ATMSITE_accum_de_ddefpts(  &
                 nframes,frame_axis,     &               !input
                 RH(p_de_drotframe),RH(p_frames),   &                 !intput
                 RH(p_frdefpts),                        &                !input
                 RH(p_de_ddefpt))                                  !output
  call ATMSITE_trans_de_ddefpts_to_centers(    &
           nfrdeflist,IH(p_frdeflst),RH(p_sitecrd),RH(p_sitefrc),  &
           nframes,RH(p_de_ddefpt))
  call ATMSITE_trans_de_dcens_to_atoms( &
           numatoms,ncenlist,  &
           IH(p_cen_list),RH(p_cen_wt),  &
           RH(p_sitefrc))
  call ATMSITE_sitefrc_to_frc(numatoms,RH(p_sitefrc),RH(p_frc),  &
                              IH(p_atmsite_ind))

  !call dump_atomfrc(numatoms,RH(p_frc),46)
  return
end subroutine ATMSITE_torque_to_force
!------------------------------------------------------------
subroutine ATMSITE_unpack()
  implicit none
#include "database.fh"
#include "atoms_to_sites.inc"

  integer p_topo,p_md,p_nonbond,p_atom_to_site,p_mpole_info,p_user,p_units
  integer ptr,polarizable_ff

  p_topo = DB_Get_Node(p_DBRoot,"topology")
  p_units = DB_Get_Node(p_topo,"phys_consts_and_unit_converts")
  p_md = DB_Get_Node(p_topo,"molec_dyn")
  p_user = DB_Get_Node(p_topo,"user_input")
  p_nonbond = DB_Get_Node(p_topo,"nonbond")
  p_atom_to_site = DB_Get_Node(p_nonbond,"atom_to_site")
  p_mpole_info = DB_Get_Node(p_nonbond,"multipole")

  ptr = DB_Get_Pointer(p_units,DB_real,"coulomb_law_in_kcal_per_mole",1)
  coulomb_constant = RH(ptr)
  ptr = DB_Get_pointer(p_md,DB_int,"numatoms",1)
  numatoms = IH(ptr)
  ptr = DB_Get_Pointer(p_nonbond,DB_int,"nsites",1)
  nsites = IH(ptr)
  ptr = DB_Get_Pointer(p_atom_to_site,DB_int,"nextra_centers",1)
  nextra_cen = IH(ptr)
  ptr = DB_Get_Pointer(p_atom_to_site,DB_int,"nextra_sites",1)
  nextra_sites = IH(ptr)
  ptr = DB_Get_Pointer(p_atom_to_site,DB_int,"ncenlist",1)
  ncenlist = IH(ptr)
  ptr = DB_Get_Pointer(p_atom_to_site,DB_int,"nframes",1)
  nframes = IH(ptr)
  ptr = DB_Get_Pointer(p_atom_to_site,DB_int,"frame_axis",3)
  frame_axis(1) = IH(ptr)
  frame_axis(2) = IH(ptr+1)
  frame_axis(3) = IH(ptr+2)
  ptr = DB_Get_Pointer(p_atom_to_site,DB_int,"nfrdeflist",1)
  nfrdeflist = IH(ptr)
  ptr = DB_Get_Pointer(p_mpole_info,DB_int,"mpole_top",1)
  mpole_top = IH(ptr)

  ptr = DB_Get_Pointer(p_user,DB_int,"polarizable_force_field",1)
  polarizable_ff = IH(ptr)

  p_crd = DB_Get_pointer(p_md,DB_real,"coordinates",3*numatoms)
  p_frc = DB_Get_pointer(p_md,DB_real,"forces",3*numatoms)
  p_sitecrd = DB_Get_pointer(p_nonbond,DB_real,"sitecrd",3*nsites)
  p_sitefrc = DB_Get_pointer(p_nonbond,DB_real,"sitefrc",3*nsites)
  p_atmsite_ind = DB_Get_pointer(p_atom_to_site,DB_int,  &
                                 "atom_site_index",numatoms)

  p_frdefpts = DB_Get_Pointer(p_atom_to_site,DB_real,  &
                              "frame_def_pts",3*2*nframes)
  p_frdeflst = DB_Get_Pointer(p_atom_to_site,DB_int,  &
                              "frame_def_list", 5*nfrdeflist)
  p_indframe = DB_Get_Pointer(p_atom_to_site,DB_int,"frame_index", nsites)
  p_frames = DB_Get_Pointer(p_atom_to_site,DB_real,"frames",9*nframes)
  p_cen_wt = DB_Get_Pointer(p_atom_to_site,DB_real,"center_wts",ncenlist)
  p_cen_list = DB_Get_Pointer(p_atom_to_site,DB_int,"center_list",2*ncenlist)
  if ( nextra_sites > 0 )then
    p_extsitelst = DB_Get_pointer(p_atom_to_site,DB_int,  &
                      "extra_site_list",nextra_sites)
    p_locvec = DB_Get_pointer(p_atom_to_site,DB_real,  &
                      "extra_site_loc_coords",3*3*nextra_sites)
  endif

  p_mpole_order = DB_Get_Pointer(p_mpole_info,DB_int,"mpole_order",nsites)
  p_mpole_offset = DB_Get_Pointer(p_mpole_info,DB_int,"mpole_offset",nsites)
  p_mpole_local = DB_Get_Pointer(p_mpole_info,DB_real,"local_mpole",mpole_top)
  p_mpole_global = DB_Get_Pointer(p_mpole_info,DB_real,  &
                     "global_mpole_perm",mpole_top)
!  if ( polarizable_ff == 0 )then
    p_mpole_xyz = DB_Get_Pointer(p_mpole_info,DB_real,  &
                     "global_mpole_perm",mpole_top)
    p_field_xyz = DB_Get_Pointer(p_mpole_info,DB_real,  &
                     "global_field_perm",mpole_top)
!  else
!    p_mpole_xyz = DB_Get_Pointer(p_mpole_info,DB_real,  &
!                     "global_mpole_total",mpole_top)
!    p_field_xyz = DB_Get_Pointer(p_mpole_info,DB_real,  &
!                     "global_field_total",mpole_top)
!  endif
! Stuff for converting torques back to sitefrc and thence to frc
  p_de_drotsite = DB_Get_Pointer(p_mpole_info,DB_real,  &
                     "d_ene_d_rotate_site",3*nsites)
  p_de_drotframe = DB_Get_Pointer(p_mpole_info,DB_real,  &
                     "d_ene_d_rotate_frame",3*nframes)
  p_de_ddefpt = DB_Get_Pointer(p_mpole_info,DB_real,  &
                     "d_ene_d_defpt_xyz",3*2*nframes)
  return
end subroutine ATMSITE_unpack
!---------------------------------------------------------------
subroutine ATM_dump_poles(  &
                 nsites,Mpole_order,Mpole_offset,Mpole_Local,Mpole_Global)
  implicit none
  integer nsites, Mpole_order(*),Mpole_offset(*)
  _REAL_ Mpole_Local(*),Mpole_Global(*)

  integer n,k,mpord,mpoff,i
  
  i = 0
  do n = 1,nsites
    mpoff = Mpole_offset(n)
    mpord = Mpole_order(n)
    if ( mpord > 0 )then
      i = i + 1
      write(35,*)'$$$$$$$$$$$$$$$$$$$$$$$$ i = ',i
      write(35,*)'----------------------- local ---------'
      write(35,'(5(1x,f14.10))')(Mpole_Local(mpoff+k),k=1,mpord)
      write(35,*)'----------------------- global ---------'
      write(35,'(5(1x,f14.10))')(Mpole_Global(mpoff+k),k=1,mpord)
    endif
  enddo
  return
end subroutine ATM_dump_poles
!---------------------------------------------------------
subroutine dump_atomfrc(numatoms,frc,nf)
  implicit none
  integer numatoms,nf
  _REAL_ frc(3,*)
  integer j,n
  do n = 1,numatoms
    write(nf,'(3f20.12)')(frc(j,n),j=1,3)
  enddo
  return
end subroutine dump_atomfrc
!---------------------------------------------------------
