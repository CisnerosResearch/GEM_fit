ifndef $(PROJECT_ROOT)
  PROJECT_ROOT=../..
endif
include $(PROJECT_ROOT)/config.mak

NAME = $(notdir $(shell pwd))

# If all three source definitions are blank, libdir.mak
# will automatically include all *.c, *.f, and *.f90 files.
!F90_SOURCE= g_load_basis.f90 g_load_hermite_coeffs.f90 g_table_boys.f90 \
!            g_mcmurch_dav.f90 g_site_site.f90 g_form_Gmat.f90 my_test.f90 \
!            g_Nmcmurch_dav.f90
F90_SOURCE= g_Nmcmurch_dav.f90 g_table_boys.f90 g_mcmurch_dav.f90 \
            g_Nmcmurch_dav_NORM.f90 g_MD_recur.f90 g_Boys_lookup.f90  \
            g_MD_recur_Jvec.f90 g_MD_recur_CART.f90 g_MD_recur_Jvec_CART.f90 \
            g_NMmD_S.f90
F_SOURCE= 
C_SOURCE=

include ../libdir.mak
sinclude depend.mak

