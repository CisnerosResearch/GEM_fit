!---------------------------------------
! PRIVATE data for atoms_to_sites_main.f90
!---------------------------------------

  _REAL_                  coulomb_constant
  common/ATMSITE_scalarR/coulomb_constant
  integer                numatoms,nsites,nextra_cen,nextra_sites,  &
                         ncenlist,nframes,nfrdeflist,mpole_top, &
                         frame_axis
  common/ATMSITE_scalarI/numatoms,nsites,nextra_cen,nextra_sites,  &
                         ncenlist,nframes,nfrdeflist,mpole_top, &
                         frame_axis(3)

  integer                p_crd,p_frc,p_sitecrd,p_sitefrc,  &
                         p_atmsite_ind,p_frdefpts,p_frdeflst,p_indframe, &
                         p_frames,p_cen_wt,p_cen_list,p_extsitelst,p_locvec, &
                         p_mpole_order,p_mpole_offset,p_mpole_local,  &
                         p_mpole_global,p_mpole_xyz,p_field_xyz,  &
                         p_de_drotsite,p_de_drotframe,p_de_ddefpt
  common/ATMSITE_memptrs/p_crd,p_frc,p_sitecrd,p_sitefrc,  &
                         p_atmsite_ind,p_frdefpts,p_frdeflst,p_indframe, &
                         p_frames,p_cen_wt,p_cen_list,p_extsitelst,p_locvec, &
                         p_mpole_order,p_mpole_offset,p_mpole_local,  &
                         p_mpole_global,p_mpole_xyz,p_field_xyz,  &
                         p_de_drotsite,p_de_drotframe,p_de_ddefpt
