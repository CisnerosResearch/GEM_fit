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
module definition

implicit none

  integer,parameter::maxatom=500

  type AO_basis

       integer:: num_basis, num_primitives, num_coefficients
       integer,allocatable:: num_basis_prims(:), prim_order(:),&
                             norm_offset(:), index(:), &
                             orb_contr(:), off_basis_prims(:)
       double precision,allocatable:: expon(:), &
                                      contr(:), &
                                      norm(:)

  end type AO_basis

  type auxiliary_basis

       integer:: num_auxbasis, num_primitives, num_coefficients
       integer,allocatable:: num_basis_prims(:), off_basis_prims(:), &
                             prim_order(:), norm_offset(:)
       double precision,allocatable:: prim_expon_term(:), &
                                      prim_normalizers(:), &
                                      prim_elec(:)

  end type auxiliary_basis

  type sites_info

       integer:: nsites, num_residues, tot_num_primitives, &
                 num_coefficients, num_frame_deflist, natoms, num_stos, &
                 tot_ao_num_prim, molframe, mol_chg
       integer,allocatable:: residue_num(:), residue_start(:), &
                             residue_end(:), site_type(:), &
                             basis_index(:), num_primitives(:), &
                             off_primitives(:), prim_order(:), & 
                             coeff_offset(:), frame_deflist(:), &
                             ao_bas_ind(:), ao_num_prim(:), ao_index(:),&
                             ao_order(:), ao_off_prim(:), ao_contr_deg(:)
       double precision, dimension(3)::p1, v1, v2, v3
       double precision, allocatable:: site_crds(:), nuclear_charge(:), &
                                       nuclear_mass(:), prim_expo(:), &
                                       local_crds(:), &
                                       cartesian_coeffs(:), &
                                       global_hermite_coeffs(:), &
                                       local_hermite_coeffs(:), &
                                       cart_coeff_norms(:), &
                                       aux_elec(:), frame_def_pts(:), &
                                       frames(:), densmat(:), promol_dens(:),&
                                       ao_expo(:), ao_contr_coeff(:), &
                                       ao_norms(:)
       logical, allocatable:: fr_valid(:) 
  end type sites_info

  type aux_orbitals
 
       integer:: x, y, z
       double precision:: expo, norm, coords(3), coords_G(3)

  end type aux_orbitals

  type ao_orbitals
 
       integer:: x, y, z, deg_contr, index
       double precision:: coords(3)
       double precision,allocatable:: expo(:), contr(:), norm(:)

  end type ao_orbitals

end module                         
