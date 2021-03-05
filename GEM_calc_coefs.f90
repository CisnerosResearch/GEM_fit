program GEM_calc_coefs

use definition

  implicit none

  include "io.fh"
  type(AO_basis)::atom_basis
  type(auxiliary_basis)::atom_auxis
  type(sites_info)::site_info
  type(aux_orbitals),allocatable::auxis(:)
  type(ao_orbitals),allocatable::basis(:)
  integer argc,arg,l,allochk,i,iargc,ncoeff,num_stos,inttype,cubetype
  double precision beta, scale, tot_esp, tot_fld(3), exact_esp, exact_fld(3),&
                   cutoff
  double precision,allocatable,dimension(:)::Jvec,aux_coefs
  double precision,allocatable,dimension(:,:)::Gmat
  character(len=PATH_MAX) argv
  character(len=80)parmfile,auxfile,crdfile,densfile,basfile,coeff_file,&
                   lhermitefile,mpole_file,mpole_file2,controlfile,filehead,&
                   promolfile,mpole_file3
  character(len=3)dtype
  logical constrain, readdens, readbas, debug, cartfit, cube, printcube,&
          twocentfit, sphere, promolfit,HLYfit

  controlfile = ''
  arg=1
  argc=IArgC()
  do while (arg <= argc)
    call GetArg(arg,argv)
    if(argv == '-control') then
      arg=arg+1
      call GetArg(arg,controlfile)
    endif
    arg = arg + 1
  enddo
  if (controlfile == '')then
       write(6,*)'usage: GEM_calc_coefs -control controlfile'
       stop
  endif
  cube = .false.

  call HCRD_read_control(controlfile,parmfile,auxfile,basfile,crdfile,&
                         densfile,promolfile,coeff_file,lhermitefile,&
                         mpole_file,mpole_file2,mpole_file3,&
                         constrain,readdens,readbas,debug,cartfit,&
                         inttype,beta,scale,cutoff,&
                         printcube,promolfit,twocentfit,cubetype,sphere,&
                         HLYfit)
  call AHCRD_read_file(site_info,crdfile)
  call AHBASE_load_auxnames(atom_auxis,parmfile)
  call AHBASE_load_auxbasis(atom_auxis,auxfile)
  call AHBASE_calcnorms(atom_auxis,densfile)
  if (readdens) call HLOAD_dens(site_info,densfile,dtype,.false.)
  if (readbas) then
     call AHBASE_load_basnames(atom_basis,parmfile)
     call AHBASE_load_AObasis(atom_basis,basfile)
     call AHBASE_calcAOnorms(atom_basis,densfile)
  endif
  call AHTYPE_load_site_info(site_info,atom_auxis,atom_basis,parmfile,readdens)
  call AHCOEFF_copy_norms_from_basis(site_info,atom_auxis)
  call AHFRAME_get_molframe(site_info,parmfile)
  call AHFRAME_load_deflist(site_info,parmfile,readdens,debug)
  call AHFRAME_build_frames(site_info)
  !call AHFRAME_rotate_coords(site_info,debug)

!! ALLOCATE G matrix, and Jvec
  ncoeff = site_info%num_coefficients
  num_stos = site_info%num_stos
  allocate(Gmat(ncoeff,ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'GEM_calc_coefs:could not allocate Gmat, exiting'
     stop
  endif
  allocate(Jvec(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'GEM_calc_coefs:could not allocate Jvec, exiting'
     stop
  endif
  allocate(aux_coefs(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'GEM_calc_coefs:could not allocate aux_coefs, exiting'
     stop
  endif
  allocate(auxis(ncoeff), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'GEM_calc_coefs:could not allocate auxis, exiting'
     stop
  endif
  allocate(basis(num_stos), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'GEM_calc_coefs:could not allocate basis, exiting'
     stop
  endif

     ! SUBROUTINES IN GLOBAL CARTESIAN
     call AH_form_G_mat_CART(site_info,auxis,Gmat,ncoeff,inttype,beta,&
                             debug,cube,densfile)
     call H_form_J_vec_C(site_info,auxis,basis,Jvec,ncoeff,num_stos,&
                         inttype,beta,debug,cube,.false.,twocentfit,densfile)

!! invert G matrix and form coefficients
  call AH_form_aux_coefs(site_info,Gmat,Jvec,aux_coefs,ncoeff,scale,&
                         constrain,debug)
  deallocate(aux_coefs,Jvec,Gmat)

!! NOW TRANSFORM CARTESIAN COEFS TO LOCAL HERMITE
  if (cartfit) then
     call AHCOEFF_load_cartesian(site_info,coeff_file,readdens,densfile)
  else
     !call AHCOEFF_load_cartesian(site_info,coeff_file,readdens)
     call AHCOEFF_mult_hermite(site_info,coeff_file,readdens)
  endif
  call AHFRAME_rotate_hermites(site_info,lhermitefile)
  call CH_COEFF_mpoles_NEW(atom_auxis,site_info,mpole_file,mpole_file2)
  call CH_COEFF_mpoles_GLOBAL(atom_auxis,site_info,mpole_file3)
  call H_calc_esp_CART(site_info,auxis,0.0d0,0.0d0,0.0d0,tot_esp,&
                       ncoeff,cartfit,cube)
  call H_exact_esp3(site_info,basis,0.0d0,0.0d0,0.0d0,exact_esp,num_stos,cube)
  call H_calc_fld_CART(site_info,auxis,0.0d0,0.0d0,0.0d0,tot_fld,&
                       ncoeff,cartfit,cube)
  call H_exact_fld2(site_info,basis,0.0d0,0.0d0,0.0d0,exact_fld,num_stos,cube)

!! calculate nuclear - electron intramolecular energy
  call H_nuc_elec(site_info,basis,num_stos)
  deallocate(basis)
  call H_nuc_elec_AUX(site_info,auxis,ncoeff,cartfit)

  deallocate(auxis)

  if (debug)close(12)

end program GEM_calc_coefs
