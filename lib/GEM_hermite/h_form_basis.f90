subroutine H_form_basis(site_info,basis,num_stos,densfile)
  ! GAC: subroutine to form ordered atomic orbitals for J vector
use definition
  implicit none
  integer num_stos
  type(sites_info)::site_info
  type(ao_orbitals)::basis(num_stos)

  integer site1,allochk, order1, scan_file
  integer k1,num1,off1,c_off1 !,t1,t2
  integer counter, res, c1, c2, natoms, count2
  integer slo1,shi1,ires,jres,n,i,j,k,contr
  double precision expon1(20),contr_coeff(20),norm(20)
  character(len=*) densfile
  character(len=4) tmp_name
  logical gau_dens, nw_dens

! check if it's a gaussian checkpoint
  gau_dens = .false.
  scan_file = 0
  tmp_name = 'fchk'
  scan_file = index(densfile,tmp_name)
  if (scan_file .ne. 0) gau_dens = .true.
  if(gau_dens)write(6,*)'h_form_basis: Gaussian fchk found, reordering f'

! check if it's a NWChem density file
  nw_dens = .false.
  scan_file = 0
  tmp_name = 'NW'
  scan_file = index(densfile,tmp_name)
  if (scan_file .ne. 0) nw_dens = .true.
  if(nw_dens)write(6,*)'h_form_basis: NWChem density found, reordering d,f'

  natoms = site_info%natoms
  !res = site_info%residue_num(1)
  !slo1 = site_info%residue_start(res)
  !shi1 = site_info%residue_end(res)
  counter = 0
  count2 = 0
  do site1 = 1,natoms
    num1 = site_info%ao_num_prim(site1)
    off1 = site_info%ao_off_prim(site1)
    k1 = 0
    !do k1 = 1,num1
    do while (k1 .lt. num1)
      order1 = site_info%ao_order(off1+k1+1)
      if (order1 == 1) then ! s shell
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         contr = site_info%ao_contr_deg(off1+k1+1)
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 0
         basis(counter)%y = 0
         basis(counter)%z = 0
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      elseif (order1 == 3) then ! p shell
      ! *** px orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         contr = site_info%ao_contr_deg(off1+k1+1)
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 1
         basis(counter)%y = 0
         basis(counter)%z = 0
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** py orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 0
         basis(counter)%y = 1
         basis(counter)%z = 0
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** pz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 0
         basis(counter)%y = 0
         basis(counter)%z = 1
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      elseif (order1 == 6) then ! d shell
       if (.not. nw_dens) then
      ! *** dxx orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         contr = site_info%ao_contr_deg(off1+k1+1)
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 2
         basis(counter)%y = 0
         basis(counter)%z = 0
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** dyy orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 0
         basis(counter)%y = 2
         basis(counter)%z = 0
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** dzz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 0
         basis(counter)%y = 0
         basis(counter)%z = 2
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** dxy orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 1
         basis(counter)%y = 1
         basis(counter)%z = 0
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** dxz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 1
         basis(counter)%y = 0
         basis(counter)%z = 1
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
     ! *** dyz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 0
         basis(counter)%y = 1
         basis(counter)%z = 1
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
       else ! nw_dens is true, d orbital order is different!!!
      ! *** dxx orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         contr = site_info%ao_contr_deg(off1+k1+1)
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 2
         basis(counter)%y = 0
         basis(counter)%z = 0
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** dxy orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 1
         basis(counter)%y = 1
         basis(counter)%z = 0
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** dxz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 1
         basis(counter)%y = 0
         basis(counter)%z = 1
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** dyy orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 0
         basis(counter)%y = 2
         basis(counter)%z = 0
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
     ! *** dyz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 0
         basis(counter)%y = 1
         basis(counter)%z = 1
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** dzz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 0
         basis(counter)%y = 0
         basis(counter)%z = 2
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
       endif ! nw_dens
      elseif (order1 == 10) then ! f shell
       if (gau_dens) then! it's a gaussian fchk, order of f's is different
         !print *,'fchk, reordering f'
      ! *** fxxx orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         contr = site_info%ao_contr_deg(off1+k1+1)
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 3
         basis(counter)%y = 0
         basis(counter)%z = 0
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fyyy orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 0
         basis(counter)%y = 3
         basis(counter)%z = 0
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fzzz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 0
         basis(counter)%y = 0
         basis(counter)%z = 3
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fxyy orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 1
         basis(counter)%y = 2
         basis(counter)%z = 0
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fxxy orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 2
         basis(counter)%y = 1
         basis(counter)%z = 0
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fxxz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 2
         basis(counter)%y = 0
         basis(counter)%z = 1
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fxzz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 1
         basis(counter)%y = 0
         basis(counter)%z = 2
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fyzz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 0
         basis(counter)%y = 1
         basis(counter)%z = 2
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fyyz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 0
         basis(counter)%y = 2
         basis(counter)%z = 1
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fxyz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 1
         basis(counter)%y = 1
         basis(counter)%z = 1
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
       elseif (nw_dens) then! it's a NWChem density, order of f's is different
      ! *** fxxx orbital*** 
         print *,'NWC, reordering f'
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         contr = site_info%ao_contr_deg(off1+k1+1)
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 3
         basis(counter)%y = 0
         basis(counter)%z = 0
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fxxy orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 2
         basis(counter)%y = 1
         basis(counter)%z = 0
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fxxz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 2
         basis(counter)%y = 0
         basis(counter)%z = 1
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fxyy orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 1
         basis(counter)%y = 2
         basis(counter)%z = 0
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fxyz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 1
         basis(counter)%y = 1
         basis(counter)%z = 1
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fxzz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 1
         basis(counter)%y = 0
         basis(counter)%z = 2
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fyyy orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 0
         basis(counter)%y = 3
         basis(counter)%z = 0
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fyyz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 0
         basis(counter)%y = 2
         basis(counter)%z = 1
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fyzz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 0
         basis(counter)%y = 1
         basis(counter)%z = 2
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fzzz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 0
         basis(counter)%y = 0
         basis(counter)%z = 3
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
       else
      ! *** fxxx orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         contr = site_info%ao_contr_deg(off1+k1+1)
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 3
         basis(counter)%y = 0
         basis(counter)%z = 0
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fyyy orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 0
         basis(counter)%y = 3
         basis(counter)%z = 0
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fzzz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 0
         basis(counter)%y = 0
         basis(counter)%z = 3
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fxxy orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 2
         basis(counter)%y = 1
         basis(counter)%z = 0
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fxxz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 2
         basis(counter)%y = 0
         basis(counter)%z = 1
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fxyy orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 1
         basis(counter)%y = 2
         basis(counter)%z = 0
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fyyz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 0
         basis(counter)%y = 2
         basis(counter)%z = 1
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fxzz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 1
         basis(counter)%y = 0
         basis(counter)%z = 2
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fyzz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 0
         basis(counter)%y = 1
         basis(counter)%z = 2
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
      ! *** fxyz orbital*** 
         expon1(:) = 0.0d0
         contr_coeff(:)=0.0d0
         norm(:)=0.0d0
         do i = 1, contr
            k1 = k1 + 1
            expon1(i) = site_info%ao_expo(off1+k1)
            contr_coeff(i) = site_info%ao_contr_coeff(off1+k1)
            norm(i) = site_info%ao_norms(off1+k1)
         enddo
         contr = site_info%ao_contr_deg(off1+k1)
         counter = counter + 1
         allocate(basis(counter)%expo(contr),&
                  basis(counter)%contr(contr),&
                  basis(counter)%norm(contr), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'H_form_basis:could not allocate expo, contr,&
                       & norm; exiting'
            stop
         endif
         basis(counter)%index = site1
         basis(counter)%x = 1
         basis(counter)%y = 1
         basis(counter)%z = 1
         basis(counter)%deg_contr = contr
         basis(counter)%coords(1) = site_info%site_crds(3*(site1-1)+1)
         basis(counter)%coords(2) = site_info%site_crds(3*(site1-1)+2)
         basis(counter)%coords(3) = site_info%site_crds(3*(site1-1)+3)
         basis(counter)%expo(:) = expon1(1:contr)
         basis(counter)%contr(:) = contr_coeff(1:contr)
         basis(counter)%norm(:) = norm(1:contr)
       endif ! gau_dens, nw_dens check
      endif ! select order
    enddo ! k1
    count2 = count2+k1
  enddo ! site1
  if(counter .ne. num_stos) then
    write(6,*)'AH_form_basis: something wrong while forming basis, exiting'
    write(6,*)counter,num_stos
    stop
  endif
  if(count2 .ne. site_info%tot_ao_num_prim) then
    write(6,*)'AH_form_basis: something wrong while forming basis, exiting'
    write(6,*)count2,site_info%tot_ao_num_prim
    stop
  endif
  return
end subroutine H_form_basis
!---------------------------------------------------------
subroutine H_form_exp_coef(i, j, A, B, P_x, p, E)
!
!08/04 GAC: subroutine to form expansion coefficients to transform
!           from hermite to cartesian gaussians.  See eqns 9.5.15,
!           9.5.16, 9.5.17, 9.5.20 and 9.5.21 from Helgaker.
!           This routine takes in the limits of the coefficients to
!           generate (i,j), the centers of the gaussian functions (A, B),
!           the ctr. of charge coord. (P_x), the combined exponent (p = a+b;
!           where a and b are the exponents for the corresponding gaussians)
!           and puts the resulting expansion coefficients into E(i,j,t).
!
  implicit none

  double precision A, B, P_x, p, X_pa, X_pb
  integer maxi, i, j, t, k, l, i1, j1, t1, dim
  double precision E(0:3,0:3,0:6)

! calculate geometric paramenters needed

  X_pa = P_x - A 
  X_pb = P_x - B 

! calculate seed coefficients

  E(0,0,0) = 1.d0
  E(1,0,0) = X_pa
  E(0,1,0) = X_pb
  E(0,0,1) = 0.d0
  E(1,0,1) = 1.d0/(2.d0*p)
  E(0,1,1) = 1.d0/(2.d0*p)

! next do E(i1,0,t) i1=2,...,i & t = 0,...,i1

  do i1 = 2,i
    E(i1,0,0) = X_pa*E(i1-1,0,0) + E(i1-1,0,1)
    do t = 1,i1
      E(i1,0,t) = (1.d0/(2.d0*p*t))*(i1*E(i1-1,0,t-1))
    enddo
  enddo

! next do E(0,j1,t) j1=2,...,j & t = 0,...,j1

  do j1 = 2,j
     E(0,j1,0) = X_pb*E(0,j1-1,0) + E(0,j1-1,1)
     do t = 1,j1
       E(0,j1,t) = (1.d0/(2.d0*p*t))*(j1*E(0,j1-1,t-1))
     enddo
  enddo

! finally do E(i1,j1,t) i1=1,...,i j1 = 1,...,j & t = 0,...,i1+j1

  do i1 = 1,i
    do j1 = 1,j
      E(i1,j1,0) = X_pa*E(i1-1,j1,0) + E(i1-1,j1,1)
      do t = 1,i1+j1
        E(i1,j1,t) = (1.d0/(2.d0*p*t))*(i1*E(i1-1,j1,t-1)+j1*E(i1,j1-1,t-1))
      enddo
    enddo
  enddo

  end subroutine H_form_exp_coef
!---------------------------------------------------------
subroutine H_form_exp_coef2(i, j, A, B, P_x, p, alpha, beta, mu, E)
!
!08/04 GAC: subroutine to form expansion coefficients to transform
!           from hermite to cartesian gaussians.  
!           NOTE: THIS IS DIFFERENT FROM ABOVE BECAUSE HERE THE
!           COEFFICIENTS ARE CALC. W.R.T. Ax and Bx NOT Px and X_ab
!           This routine takes in the limits of the coefficients to
!           generate (i,j), the centers of the gaussian functions (A, B),
!           the ctr. of charge coord. (P_x), the combined exponent (p = a+b;
!           where a and b are the exponents for the corresponding gaussians)
!           and puts the resulting expansion coefficients into E(i,j,t).
!
  implicit none

  double precision A, B, P_x, p, X_pa, X_pb, mu, R_ab, K_ab, pf1, pf2,&
                   alpha, beta, X_ab
  integer maxi, i, j, t, k, l, i1, j1, t1, dim, i2
  double precision E(0:3,0:3,0:6)

! calculate geometric paramenters needed

  !X_pa = P_x - A 
  !X_pb = P_x - B 
  X_ab = A - B
  R_ab = (A - B)**2
  K_ab = exp(-1.d0*mu*R_ab)
  pf1 = 0.5d0/(alpha + beta)
  pf2 = -2.d0*beta*pf1

! calculate seed coefficients

  E(0,0,0) = K_ab
  E(1,0,0) = X_ab*K_ab*pf2
  E(0,0,1) = 0.d0
  E(1,0,1) = K_ab*pf1

! next do E(i1,0,t) i1=2,...,i & t = 0,...,i1

  do i1 = 2,i
    E(i1,0,0) = pf2*X_ab*E(i1-1,0,0) + E(i1-1,0,1)
    do t = 1,i1-2
      E(i1,0,t) = pf1*E(i1-1,0,t-1) + pf2*X_ab*E(i1-1,0,t) + &
                  (t+1)*E(i1-1,0,t+1)
    enddo
    t = i1-1
    E(i1,0,t) = pf1*E(i1-1,0,t-1) + pf2*X_ab*E(i1-1,0,t)
    E(i1,0,t+1) = pf1*E(i1-1,0,t)
  enddo

! update prefactor

  pf2 = 2.d0*alpha*pf1

! next do E(0,1,0), E(0,1,1) 

  E(0,1,0) = pf2*X_ab*K_ab
  E(0,1,1) = pf1*K_ab

! finally do E(i1,j1,t) i1=1,...,i j1 = 1,...,j & t = 0,...,i1+j1

  do j1 = 1,j
    if (j1 == 1) then
       i2 = 1
    else
       i2 = 0
    endif
    do i1 = i2,i
      E(i1,j1,0) = pf2*X_ab*E(i1,j1-1,0) + E(i1,j1-1,1)
      do t = 1,i1+j1-2
        E(i1,j1,t) = pf1*E(i1,j1-1,t-1) + pf2*X_ab*E(i1,j1-1,t) + &
                     (t+1)*E(i1,j-1,t+1)
      enddo
      t = i1+j1-1
      E(i1,j1,t) = pf1*E(i1,j1-1,t-1) + pf2*X_ab*E(i1,j-1,t)
      E(i1,j1,t+1) = pf1*E(i1,j-1,t)
    enddo
  enddo

  end subroutine H_form_exp_coef2
!---------------------------------------------------------
  subroutine H_form_exp_coef_noS(i, A, p, E)
!
!08/04 GAC: subroutine to form expansion coefficients to transform
!           from hermite to cartesian gaussians.  See eqns 9.5.15,
!           9.5.16, 9.5.17, 9.5.20 and 9.5.21 from Helgaker.
!           This routine is a modification of form_exp_coef.  It deals
!           with single orbitals (no overlap). In this case the
!           resulting expansion coefficients are stored in E(i,t)
!
  implicit none

  double precision A, B, P_x, p
  integer i, t, i1, t1, dim
  double precision E(0:3,0:3)

! calculate seed coefficients

  E(0,0) = 1.d0
  E(1,0) = 0.d0
  E(0,1) = 0.d0
  E(1,1) = 1.d0/(2.d0*p)

! next do E(i1,t) i1=2,...,i & t = 0,...,i1

  do i1 = 2,i
    E(i1,0) = E(i1-1,1)
    do t = 1,i1
       E(i1,t) = (1.d0/(2*p*t))*(i1*E(i1-1,t-1))
    enddo
  enddo

  end subroutine H_form_exp_coef_noS
!---------------------------------------------------------
subroutine H_form_exp_coef_D(i, j, maxD, A, B, alpha, beta, mu, E, dE)
!
!03/07 GAC: subroutine to form expansion coefficients for derivatives of
!           the exponential prefactor w.r.t. A and B. Helgaker, eqns.
!           9.5.6 and 9.5.7 (modified to depend on X_ab instead of X_pa),
!           and puts the resulting expansion coefficients into dE(der,i,j,t).
!           NOTE: E(i,j,t) need to be precomputed and this subroutine
!           needs to be called recursively for higher derivatives.
!
  implicit none

  double precision A, B, P_x, p, p2, X_ab, tx, tx1, tx2, alpha, beta, d1,&
                   pf1, pf2, R_ab, K_ab, mu
  integer maxi, i, j, t, k, l, i1, j1, t1, dim, der, i1a, maxD, i2
  double precision E(0:3,0:3,0:6), dE(0:maxD,0:3,0:3,0:6)

! calculate geometric paramenters needed

  !X_pa = P_x - A 
  !X_pb = P_x - B 
  X_ab = A - B
  R_ab = (A - B)**2
  K_ab = exp(-1.d0*mu*R_ab)
  pf1 = 0.5d0/(alpha + beta)
  pf2 = -2.d0*beta*pf1

! calculate seed coefficients

  dE(0,0,0,0) = K_ab
  dE(0,1,0,0) = X_ab*K_ab*pf2
  dE(0,0,0,1) = 0.d0
  dE(0,1,0,1) = K_ab*pf1

! next do dE(0,i1,0,t) i1=2,...,i & t = 0,...,i1

  do i1 = 2,i
    dE(0,i1,0,0) = pf2*X_ab*dE(0,i1-1,0,0) + dE(0,i1-1,0,1)
    do t = 1,i1-2
      dE(0,i1,0,t) = pf1*dE(0,i1-1,0,t-1) + pf2*X_ab*dE(0,i1-1,0,t) + &
                  (t+1)*dE(0,i1-1,0,t+1)
    enddo
    t = i1-1
    dE(0,i1,0,t) = pf1*dE(0,i1-1,0,t-1) + pf2*X_ab*dE(0,i1-1,0,t)
    dE(0,i1,0,t+1) = pf1*dE(0,i1-1,0,t)
  enddo

! update prefactor

  pf2 = 2.d0*alpha*pf1

! next do dE(0,0,1,0), dE(0,0,1,1) 

  dE(0,0,1,0) = pf2*X_ab*K_ab
  dE(0,0,1,1) = pf1*K_ab

! finally do dE(0,i1,j1,t) i1=1,...,i j1 = 1,...,j & t = 0,...,i1+j1

  do j1 = 1,j
    if (j1 == 1) then
       i2 = 1
    else
       i2 = 0
    endif
    do i1 = i2,i
      dE(0,i1,j1,0) = pf2*X_ab*dE(0,i1,j1-1,0) + dE(0,i1,j1-1,1)
      do t = 1,i1+j1-2
        dE(0,i1,j1,t) = pf1*dE(0,i1,j1-1,t-1) + pf2*X_ab*dE(0,i1,j1-1,t) + &
                     (t+1)*dE(0,i1,j-1,t+1)
      enddo
      t = i1+j1-1
      dE(0,i1,j1,t) = pf1*dE(0,i1,j1-1,t-1) + pf2*X_ab*dE(0,i1,j-1,t)
      dE(0,i1,j1,t+1) = pf1*dE(0,i1,j-1,t)
    enddo
  enddo

! calculate needed parameters

  X_ab = A - B
  p = 0.5/(alpha + beta)
  d1 = -2.d0*(alpha*beta/(alpha+beta))

  !dE(:,:,:,:) = 0.d0
  !!dE(0,:,:,:) = E(:,:,:) ! note that E(i,j,t) need to be precomp.
  !do i1 = 0,i
  !   do j1 = 0,j
  !      do t1 = 0,i+j
  !         dE(0,i1,j1,t) = E(i1,j1,t)
  !      enddo
  !   enddo
  !enddo

! calculate seed coefficients

  if (maxD .eq. 1) then
     dE(maxD,0,0,0) = d1*X_ab*dE(0,0,0,0) 
  else 
     dE(maxD,0,0,0) = d1*(X_ab*dE(maxD-1,0,0,0) + (maxD-1)*dE(maxD-2,0,0,0))
  endif   
 
! do dE(maxD,i,0,t) i = 1, i, t = 0, t

  if (i .gt. 0) then

     p2 = -2.d0*beta*p
     if (maxD .eq. 1) then
        dE(maxD,1,0,0) = p2*(X_ab*dE(maxD,0,0,0)+dE(0,0,0,0))
     else
        dE(maxD,1,0,0) = p2*(X_ab*dE(maxD,0,0,0)+maxD*dE(maxD-1,0,0,0))
     endif
     dE(maxD,1,0,1) = p*dE(maxD,0,0,0)
 
     do i1 = 2, i

        tx1 = dE(maxD,i1-1,0,1)
        if (maxD .eq. 1) then
           tx = X_ab*dE(maxD,i1-1,0,0) + dE(0,i1-1,0,0)
        else
           tx = X_ab*dE(maxD,i1-1,0,0) + maxD*dE(maxD-1,i1-1,0,0)
        endif
        dE(maxD,i1,0,0) = p2*tx+tx1
        
        do t = 1, i1-2
           tx = dE(maxD,i1-1,0,t-1)
           tx2 = dE(maxD,i1-1,0,t+1)
           if (maxD .eq. 1) then
              tx1 = X_ab*dE(maxD,i1-1,0,t) + dE(0,i1-1,0,t)
           else
              tx1 = X_ab*dE(maxD,i1-1,0,t) + maxD*dE(maxD-1,i1-1,0,t)
           endif
           dE(maxD,i1,0,t) = p*tx + p2*tx1 + (t+1)*tx2
        enddo

        t = i1-1
        tx = dE(maxD,i1-1,0,t-1)
        if (maxD .eq. 1) then
           tx1 = X_ab*dE(maxD,i1-1,0,t) + maxD*dE(0,i1-1,0,t)
        else
           tx1 = X_ab*dE(maxD,i1-1,0,t) + maxD*dE(maxD-1,i1-1,0,t)
        endif
        dE(maxD,i1,0,t) = p*tx + p2*tx1
        tx = dE(maxD,i1-1,0,t)
        dE(maxD,i1,0,t+1) = p*tx

     enddo ! i1 = 2, i

  endif ! i > 0
  
! do dE(maxD,i,j,t) i = 1, i; j=1, j; t = 0, t

  if (j .gt. 0) then
    
     p2 = 2.d0*alpha*p
     if (maxD .eq. 1) then
        dE(maxD,0,1,0) = p2*(X_ab*dE(maxD,0,0,0)+maxD*dE(0,0,0,0))
     else
        dE(maxD,0,1,0) = p2*(X_ab*dE(maxD,0,0,0)+maxD*dE(maxD-1,0,0,0))
     endif
     dE(maxD,0,1,1) = p*dE(maxD,0,0,0)
 
     do j1 = 1, j

        if (j1 .eq. 1) then
           i1a = 1
        else
           i1a = 0
        endif

        do i1 = i1a, i

           tx1 = dE(maxD,i1,j1-1,1)
           if (maxD .eq. 1) then
              tx = X_ab*dE(maxD,i1,j1-1,0) + maxD*dE(0,i1,j1-1,0)
           else
              tx = X_ab*dE(maxD,i1-1,j1,0) + maxD*dE(maxD-1,i1-1,j1,0)
           endif
           dE(maxD,i1,j1,0) = p2*tx+tx1
        
           do t = 1, i1+j1-2
              tx = dE(maxD,i1,j1-1,t-1)
              tx2 = dE(maxD,i1,j1-1,t+1)
              if (maxD .eq. 1) then
                 tx1 = X_ab*dE(maxD,i1,j1-1,t) + maxD*dE(0,i1,j1-1,t)
              else
                 tx1 = X_ab*dE(maxD,i1,j1-1,t) + maxD*dE(maxD-1,i1,j1-1,t)
              endif
              dE(maxD,i1,j1,t) = p*tx + p2*tx1 + (t+1)*tx2
           enddo

           t = i1+j1-1
           tx = dE(maxD,i1,j1-1,t-1)
           if (maxD .eq. 1) then
              tx1 = X_ab*dE(maxD,i1,j1-1,t) + maxD*dE(0,i1,j1-1,t)
           else
              tx1 = X_ab*dE(maxD,i1,j1-1,t) + maxD*dE(maxD-1,i1,j1-1,t)
           endif
           dE(maxD,i,0,t) = p*tx + p2*tx1
           tx = dE(maxD,i1,j1-1,t)
           dE(maxD,i1,j1,t+1) = p*tx

        enddo !i1 = i1a,i
   
     enddo ! j1 = 1, j

  endif ! j > 0

  end subroutine H_form_exp_coef_D
!---------------------------------------------------------
subroutine H_basis_limits(site_info,basis,num_stos,natoms,limitvec)
  ! GAC: subroutine to form ordered atomic orbitals for J vector
use definition
  implicit none
  integer num_stos,natoms
  type(sites_info)::site_info
  type(ao_orbitals)::basis(num_stos)

  integer site1,allochk, order1 
  integer k1,num1,off1,c_off1 !,t1,t2
  integer counter, res, c1, c2, count2, count3
  integer slo1,shi1,ires,jres,n,i,j,k,contr, limitvec(natoms)
  double precision expon1(20),contr_coeff(20),norm(20)

  counter = 0
  count2 = 0
  count3 = 0
  do site1 = 1,natoms
    num1 = site_info%ao_num_prim(site1)
    off1 = site_info%ao_off_prim(site1)
    k1 = 0
    !do k1 = 1,num1
    do while (k1 .lt. num1)
      order1 = site_info%ao_order(off1+k1+1)
      if (order1 == 1) then ! s shell
         count3 = count3 + order1
         contr = site_info%ao_contr_deg(off1+k1+1)
         do i = 1, contr*order1
            k1 = k1 + 1
         enddo
      elseif (order1 == 3) then ! p shell
         count3 = count3 + order1
         contr = site_info%ao_contr_deg(off1+k1+1)
         do i = 1, contr*order1
            k1 = k1 + 1
         enddo
      elseif (order1 == 6) then ! d shell
         count3 = count3 + order1
         contr = site_info%ao_contr_deg(off1+k1+1)
         do i = 1, contr*order1
            k1 = k1 + 1
         enddo
      elseif (order1 == 10) then ! f shell
         count3 = count3 + order1
         contr = site_info%ao_contr_deg(off1+k1+1)
         do i = 1, contr*order1
            k1 = k1 + 1
         enddo
      endif ! select order
    enddo ! k1
    count2 = count2+k1
    limitvec(site1) = count3
  enddo ! site1

  if(limitvec(natoms) .ne. num_stos) then
    write(6,*)'H_basis_limits: num_stos different than limitvector, exiting'
    write(6,*)counter,num_stos
    stop
  endif
  if(count2 .ne. site_info%tot_ao_num_prim) then
    write(6,*)'H_basis_limits: total number or primitives differs, exiting'
    write(6,*)count2,site_info%tot_ao_num_prim
    stop
  endif
  return
end subroutine H_basis_limits
