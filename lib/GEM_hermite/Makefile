ifndef $(PROJECT_ROOT)
  PROJECT_ROOT=../..
endif
include $(PROJECT_ROOT)/config.mak

NAME = $(notdir $(shell pwd))

# If all three source definitions are blank, libdir.mak
# will automatically include all *.c, *.f, and *.f90 files.
F90_SOURCE= h_read_crdfile.f90 h_loadbasis.f90 h_load_site_type.f90 \
            h_load_coefficients.f90 h_frameload.f90 h_site_site.f90 \
            h_form_aux_coefs.f90 h_load_densmat.f90 \
            h_G_mat_CART.f90 h_J_vec_CART.f90 h_esp.f90 \
            h_field.f90 h_Grid_utils.f90 \
            h_exact_esp_fld.f90 h_cube_gen_exact.f90 h_esp_aux.f90 \
            h_esp_nuc.f90 h_cube_aux.f90 h_nuc_elec.f90 \
            h_density.f90 h_cube_aux_coef.f90  h_exchange.f90 \
            h_form_basis.f90 h_spherical_grid.f90 h_Grid_Gen.f90 \
            h_read_control_NEW.f90 h_read_control.f90 ch_coeff_mpoles_new.f90
F_SOURCE= 
C_SOURCE=

include ../libdir.mak
sinclude depend.mak

