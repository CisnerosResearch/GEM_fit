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
!----------------------------------------------------------
subroutine H_cube_aux(site_info,auxis,aux_cubes,ex_cubes,Gmat,Jvec,cutoff,&
                      rad_grid,wght,npts,ncoeff,debug,cubetype,sphere,&
                      promolfit,errorfit,densfile,expdens)
use definition
use cubes
  implicit none
  type(sites_info)::site_info
  type(aux_orbitals)::auxis(ncoeff)
  type(auxiliary_cubes)::aux_cubes
  type(exact_cubes)::ex_cubes

  integer i,j,k,l,m,n,i2
  integer n_cub_atoms, num_a, num_b, num_c, num_stos, allochk, istat
  integer tmp_cub_atm, tmp_num_a, ncoeff, cubetype, npts, numcoefs
  integer nsites, num, off, order, scoff, dimxy
  integer,parameter::size = 250
  double precision dens,esp,field(3),x_orig,y_orig,z_orig,x1,y1,z1,&
                   x2,y2,z2,x3,y3,z3,xin,yin,zin,tmpJ,tmpG,cutoff,&
                   dcube,ecube,fcubex,fcubey,fcubez,dens2,esp2,field2(3)
  double precision tmp_x,tmp_y,tmp_z,tmp_x1,tmp_y1,tmp_z1,dist
  double precision rhoapprox
  double precision Gmat(ncoeff,ncoeff),Jvec(ncoeff),rad_grid(npts,3),wght(npts)
  double precision,allocatable::nuc_esp(:,:,:),nuc_fld(:,:,:,:),w_r(:)
  double precision,parameter::AtoBohr = 0.529177249d0,small=1.d-6,&
                              sigma=0.42d0,logrhoref=-7.d0 ! FOR HAO_NEW
                              !sigma=0.80d0,logrhoref=-4.0d0 ! FROM PLOT
                              !sigma=0.8d0,logrhoref=-9.d0 ! ORIGINAL VALUES
  logical debug, skip, jfile, gfile, incore, update, sphere, promolfit,&
          expdens,errorfit
  character*1 int_type
  character(len=*) densfile

  jfile = .false.
  gfile = .false.
  incore = .true.
  !expdens = .false.
  !inquire(file = 'HLY_fit', exist = expdens)
  if (promolfit) expdens = .false.
  !if (sphere .and. (.not. expdens)) expdens = .true.
  if(expdens)write(6,*)'H_cube_aux: HLY_fit found, using HLY weighting &
                        &function'

! make sure only to calculate J and G on atoms
  if(.not. sphere) then
    numcoefs = ncoeff
  else
    if (promolfit) then
       numcoefs = 0
       do n = 1, site_info%natoms
          num = site_info%num_primitives(n)
          off = site_info%off_primitives(n)
          do k = 1,num
            order = site_info%prim_order(off+k)
            numcoefs = numcoefs+order
          enddo
       enddo
       print *,'numcoefs (only atoms) = ',numcoefs
    else
       numcoefs = ncoeff
    endif 
  endif

  if ((ncoeff .gt. 750) .and. (cubetype .eq. 0)) incore = .false.
  if (.not. incore) then
     write(6,*)'H_cube_aux: too many coefficients, cubes will be written to &
                &disk.'
  endif

! check cubetype is Ok
  if ((cubetype .lt. 0) .or. (cubetype .gt. 3)) then
     write(6,*)'H_cube_aux: type of cube to fit does not match valid choices;&
                & exiting',cubetype
     stop
  endif

! Create basis info
  if(.not.errorfit)call AH_form_auxis(site_info,auxis,ncoeff,densfile)

  inquire(file = 'GEM.cubeJvec', exist = jfile)
  inquire(file = 'GEM.cubeGmat', exist = gfile)
  if (jfile .and. gfile) then
     write(6,*)'H_cube_aux: Gmat and Jvec files found, skipping cube &
                &calculation'
     open (14,file='GEM.cubeJvec',status="unknown")
     do i = 1, ncoeff
        read(14,*)Jvec(i)
     enddo
     close(14)
     open (14,file='GEM.cubeGmat',status="unknown")
     do i = 1, ncoeff
        do j = 1, numcoefs
           read(14,*)Gmat(i,j)
        enddo
     enddo
     close(14)
     return
  endif

! get cube info
  if (.not. sphere) then
     n_cub_atoms = ex_cubes%n_cub_atoms
     num_a = ex_cubes%num_a
     num_b = ex_cubes%num_b
     num_c = ex_cubes%num_c
     x_orig = ex_cubes%x_orig
     y_orig = ex_cubes%y_orig
     z_orig = ex_cubes%z_orig
     x1 = ex_cubes%x1
     y1 = ex_cubes%y1
     z1 = ex_cubes%z1
     x2 = ex_cubes%x2
     y2 = ex_cubes%y2
     z2 = ex_cubes%z2
     x3 = ex_cubes%x3
     y3 = ex_cubes%y3
     z3 = ex_cubes%z3
     !!! Check if number of atoms in info file matches
     if (n_cub_atoms .ne. site_info%natoms) then
        write(6,*)'H_cube_aux: Number of atoms in cube_file &
                   &does not match, exiting'
        stop
     endif
  else
     num_a = 1
     num_b = 1
     num_c = npts
  endif

! allocate auxiliary cubes
     if (incore) then
        allocate(aux_cubes%esp_cube(ncoeff,num_c,num_b,num_a),&
                 stat=allochk)
        if (allochk .gt. 0) then
           write(6,*)'H_cube_aux: could not allocate cubes; exiting'
           stop
        endif
        aux_cubes%esp_cube(:,:,:,:) = 0.0d0
     else
       open(52,file='aux_esp_cube',status='unknown')
     endif

  if (expdens) then
     allocate(w_r(num_c),stat=allochk)
     if (allochk .gt. 0) then
        write(6,*)'H_cube_aux: could not allocate w_r; exiting'
        stop
     endif
  else
     allocate(w_r(1),stat=allochk)
     if (allochk .gt. 0) then
        write(6,*)'H_cube_aux: could not allocate w_r; exiting'
        stop
     endif
  endif
  w_r(:) = 1.d0
   
! allocate temporary cubes for nuclear contribution

     allocate(nuc_esp(num_c,num_b,num_a), stat=allochk)
     if (allochk .gt. 0) then
        write(6,*)'H_cube_aux: could not allocate nuclear cubes; exiting'
        stop
     endif

! calculate cubes
 !first do temporary cubes for nuclear contribution
  if (cubetype .ne. 1) then
     skip = .false.
     do i = 1, num_a
        do j = 1, num_b
           do k = 1, num_c
              if (.not. sphere .and. .not. expdens) then 
                 ! Calculate each point starting from origin
                 xin = x_orig + (i-1)*x1 + (j-1)*x2 + (k-1)*x3
                 yin = y_orig + (i-1)*y1 + (j-1)*y2 + (k-1)*y3
                 zin = z_orig + (i-1)*z1 + (j-1)*z2 + (k-1)*z3
                 do n = 1, n_cub_atoms
                    dist = sqrt( (site_info%site_crds(3*(n-1)+1)-xin)**2 + &
                                 (site_info%site_crds(3*(n-1)+2)-yin)**2 + &
                                 (site_info%site_crds(3*(n-1)+3)-zin)**2 ) 
                    if (promolfit) then
                       if (dist .gt. cutoff) skip = .true.
                    else
                       if (dist .lt. cutoff) skip = .true.
                    endif
                 enddo
              else
                 xin = rad_grid(k,1)
                 yin = rad_grid(k,2)
                 zin = rad_grid(k,3)
                 if (.not. expdens) then ! use regular cutoff
                    do n = 1, site_info%natoms
                       dist = sqrt( (site_info%site_crds(3*(n-1)+1)-xin)**2 + &
                                    (site_info%site_crds(3*(n-1)+2)-yin)**2 + &
                                    (site_info%site_crds(3*(n-1)+3)-zin)**2 ) 
                       !if (dist .lt. cutoff) skip = .true.
                       if (promolfit) then
                          if (dist .gt. cutoff) skip = .true.
                       else
                          if (dist .lt. cutoff) skip = .true.
                       endif
                    enddo
                 endif
              endif
              if (.not. skip) then
                    call H_esp_nuc(site_info,xin,yin,zin,esp)
                    nuc_esp(k,j,i) = esp
              else
                    nuc_esp(k,j,i) = 0.0d0
              endif
              skip = .false.
          enddo
        enddo
     enddo
  endif

 ! now do electronic for each auxiliary and calculate Jvec
  skip = .false.
  dens = 0.0d0
  esp = 0.0d0
  field(:) = 0.0d0
  do l = 1, numcoefs
     tmpJ = 0.0d0
     do i = 1, num_a
        do j = 1, num_b
           do k = 1, num_c
              if (.not. sphere) then 
                 ! Calculate each point starting from origin
                 xin = x_orig + (i-1)*x1 + (j-1)*x2 + (k-1)*x3
                 yin = y_orig + (i-1)*y1 + (j-1)*y2 + (k-1)*y3
                 zin = z_orig + (i-1)*z1 + (j-1)*z2 + (k-1)*z3
                 do n = 1, n_cub_atoms
                    dist = sqrt( (site_info%site_crds(3*(n-1)+1)-xin)**2 + &
                                 (site_info%site_crds(3*(n-1)+2)-yin)**2 + &
                                 (site_info%site_crds(3*(n-1)+3)-zin)**2 ) 
                    !if (dist .lt. cutoff) skip = .true.
                    if (promolfit) then
                       if (dist .gt. cutoff) skip = .true.
                    else
                       if (.not. expdens)then
                          if (dist .lt. cutoff) skip = .true.
                       else
                          call promol_rho(site_info,xin,yin,zin,rhoapprox)
                          w_r(k) = &
                          exp(-1.d0*sigma*(log(rhoapprox)-logrhoref)**2)
                       endif
                    endif
                 enddo
              else
                 xin = rad_grid(k,1)
                 yin = rad_grid(k,2)
                 zin = rad_grid(k,3)
                 if (.not. expdens) then ! use regular cutoff
                    do n = 1, site_info%natoms
                       dist = sqrt( (site_info%site_crds(3*(n-1)+1)-xin)**2 + &
                                    (site_info%site_crds(3*(n-1)+2)-yin)**2 + &
                                    (site_info%site_crds(3*(n-1)+3)-zin)**2 ) 
                       !if (dist .lt. cutoff) skip = .true.
                       if (promolfit) then
                          if (dist .gt. cutoff) skip = .true.
                       else
                          if (dist .lt. cutoff) skip = .true.
                       endif
                    enddo
                 else ! use Hao's weighing function
                    !w_r(k) = exp(-1.d0*sigma*&
                    !      (log(ex_cubes%dens_cube(k,j,i))-logrhoref)**2)
                    call promol_rho(site_info,xin,yin,zin,rhoapprox)
                    w_r(k) = exp(-1.d0*sigma*(log(rhoapprox)-logrhoref)**2)
                 endif
              endif
              if (.not. skip) then
                    call H_esp_aux(auxis,xin,yin,zin,esp,l)
                    if (incore) then
                       aux_cubes%esp_cube(l,k,j,i) = esp + nuc_esp(k,j,i)
                    else
                       write(52,*)esp + nuc_esp(k,j,i)
                    endif
                    if (sphere) then
                       if (expdens) then
                          tmpJ = tmpJ + (esp+nuc_esp(k,j,i))*&
                                        ex_cubes%esp_cube(k,j,i)*wght(k)*w_r(k)
                       else
                          tmpJ = tmpJ + (esp+nuc_esp(k,j,i))*&
                                        ex_cubes%esp_cube(k,j,i)*wght(k)
                       endif
                    else
                       if (expdens) then
                          tmpJ = tmpJ + (esp+nuc_esp(k,j,i))*&
                                        ex_cubes%esp_cube(k,j,i)*w_r(k)
                       else
                          tmpJ = tmpJ + (esp+nuc_esp(k,j,i))*&
                                        ex_cubes%esp_cube(k,j,i)
                       endif
                    endif
              else
                 if(.not.incore)then
                       write(52,*)0.0d0
                 endif ! not incore to write out to file
              endif
              skip = .false.
           enddo
        enddo
     enddo
     Jvec(l) = tmpJ
  enddo
  if(debug)write(12,*)'---------------- Jvec ----------------------'
  if (debug) then
     do i = 1, numcoefs
        write(12,1001)i,Jvec(i)
     enddo
  endif
1001    format ('Jvec (',i4,') = ',f27.20)
  open (14,file='GEM.cubeJvec',status="unknown")
  do i = 1, numcoefs
     write(14,*)Jvec(i)
  enddo
  close(14)

!1004 write(6,*)'H_cube_aux: reached end of J',ncoeff,num_a,num_b,num_c
!  stop

  if (.not. incore) then
        close(52)
        call system('/bin/cp aux_esp_cube aux_esp_cube2')
        open(62,file='aux_esp_cube2',status='unknown')
        open(52,file='aux_esp_cube',status='unknown')
  endif

! calculate Gmat 
  i2 = numcoefs
  update = .false.
  skip = .false.
  do m = 1, numcoefs
     update = .true.
     if (.not. incore) then
           close(62)
           open(62,file='aux_esp_cube2',status='unknown')
     endif
     if (incore) i2 = m
     do l = 1, i2
        tmpG = 0.0d0
        do i = 1, num_a
           do j = 1, num_b
              do k = 1, num_c
                 if (.not. sphere) then
                    xin = x_orig + (i-1)*x1 + (j-1)*x2 + (k-1)*x3
                    yin = y_orig + (i-1)*y1 + (j-1)*y2 + (k-1)*y3
                    zin = z_orig + (i-1)*z1 + (j-1)*z2 + (k-1)*z3
                    do n = 1, n_cub_atoms
                       dist = sqrt( (site_info%site_crds(3*(n-1)+1)-xin)**2 + &
                                    (site_info%site_crds(3*(n-1)+2)-yin)**2 + &
                                    (site_info%site_crds(3*(n-1)+3)-zin)**2 ) 
                       !if (dist .lt. cutoff) skip = .true.
                       if (promolfit) then
                          if (dist .gt. cutoff) skip = .true.
                       else
                          if (dist .lt. cutoff) skip = .true.
                       endif
                    enddo
                 else 
                    xin = rad_grid(k,1)
                    yin = rad_grid(k,2)
                    zin = rad_grid(k,3)
                    if (.not. expdens) then ! use regular cutoff
                       do n = 1, site_info%natoms
                          dist = sqrt( (site_info%site_crds(3*(n-1)+1)-xin)**2+&
                                       (site_info%site_crds(3*(n-1)+2)-yin)**2+&
                                       (site_info%site_crds(3*(n-1)+3)-zin)**2) 
                          !if (dist .lt. cutoff) skip = .true.
                          if (promolfit) then
                             if (dist .gt. cutoff) skip = .true.
                          else
                             if (dist .lt. cutoff) skip = .true.
                          endif
                       enddo
                    endif
                 endif
                 if (.not. skip) then

                       if (incore) then
                          if (sphere) then
                             if (expdens) then
                                tmpG = tmpG + aux_cubes%esp_cube(l,k,j,i)*&
                                              aux_cubes%esp_cube(m,k,j,i)*&
                                              wght(k)*w_r(k)
                             else
                                tmpG = tmpG + aux_cubes%esp_cube(l,k,j,i)*&
                                              aux_cubes%esp_cube(m,k,j,i)*&
                                              wght(k)
                             endif
                          else
                             if (expdens) then
                                tmpG = tmpG + aux_cubes%esp_cube(l,k,j,i)*&
                                              aux_cubes%esp_cube(m,k,j,i)*&
                                              w_r(k)
                             else
                                tmpG = tmpG + aux_cubes%esp_cube(l,k,j,i)*&
                                              aux_cubes%esp_cube(m,k,j,i)
                             endif
                          endif
                       else
                          if (update) then
                             read(52,*)esp
                             update = .false.
                          endif
                          read(62,*)esp2
                          tmpG = tmpG + esp*esp2
                       endif

                 else
                    if(.not.incore)then
                          read(62,*)esp2
                    endif
                 endif
                 skip = .false.
              enddo ! k
           enddo ! j
        enddo ! i
        Gmat(l,m) = tmpG
     enddo ! l
  enddo ! m
! symmetrize
  if(debug)write(12,*)'---------------- Gmat ----------------------'
  do i = 1, numcoefs
     do j = 1, numcoefs
        if(incore)Gmat(j,i) = Gmat(i,j)
        if (debug) write(12,1002)i,j,Gmat(i,j)
     enddo
  enddo
1002    format ('Gmat (',i4,',',i4,') = ',f27.20)
  open (14,file='GEM.cubeGmat',status="unknown")
  do i = 1, numcoefs
     do j = 1, numcoefs
        write(14,*)Gmat(i,j)
     enddo
  enddo
  close(14)

!1003        write(6,*)'H_cube_aux: reached end of dens2',m,l,i,j,k

! deallocate arrays
  if (incore) then
     deallocate(aux_cubes%esp_cube,nuc_esp)
     deallocate(w_r)
   else
        close(52)
        close(62)
        call system ('/bin/rm aux_esp_cube aux_esp_cube2')
  endif

end subroutine H_cube_aux
!----------------------------------------------------------
subroutine promol_rho(site_info,xin,yin,zin,rhoapprox)

use definition
  implicit none
  type(sites_info)::site_info
  integer i, j, k
  double precision dist,xin,yin,zin,dx,dy,dz,rhoapprox,rhoatom
  double precision,parameter::mone=-1.d0,tiny=1.d-12
  double precision::Hs(2,1),Cs(2,2),Ns(2,2),Os(2,2),Fs(2,2),Nes(2,2),&
                    Ps(3,2),Ss(3,2),Cls(3,2),Ars(3,2)
  ! DATA IN ANGSTROMS FROM HAO
  !data Hs/0.384137961d0,3.90762643d0/
  !data Cs/166.591448d0,3.23010126d0,29.0603279d0,5.01709331d0/
  !data Ns/256.609200d0,2.58989432d0,31.2114908d0,5.45471548d0/
  !data Os/243.630909d0,2.53736474d0,26.3836036d0,4.29335839d0/
  !data Ps/2282.83071d0,155.142338d0,1.82194667d0,73.7103367d0,15.6986998d0,&
  !        3.38628928d0/
  !data Ss/2736.19302d0,206.867393d0,2.78312612d0,78.9192252d0,17.4500522d0,&
  !        3.51974385d0/
  ! DATA IN au (EXPONENTS) FROM HAO
  !data Hs/0.384137961d0,2.0678270d0/
  !data Cs/166.591448d0,3.23010126d0,15.3780644d0,2.65493163d0/
  !data Ns/256.609200d0,2.58989432d0,16.5164108d0,2.88651133d0/
  !data Os/243.630909d0,2.53736474d0,13.9616028d0,2.27194758d0/
  !data Ps/2282.83071d0,155.142338d0,1.82194667d0,39.0058332d0,8.30739477d0,&
  !        1.79194724d0/
  !data Ss/2736.19302d0,206.867393d0,2.78312612d0,41.7622585d0,9.23417062d0,&
  !        1.86256837d0/
  ! OPTIMIZED DATA IN au
  data Hs/0.318035d0,1.999273d0/
  data Cs/127.042163d0,0.724069d0,11.982138d0,1.551190d0/
  data Ns/206.361417d0,1.617216d0,14.143747d0,1.894116d0/
  data Os/313.539616d0,3.094042d0,16.320394d0,2.231987d0/
  data Fs/452.348028d0,5.457284d0,18.519116d0,2.572569d0/
  data Nes/626.814078d0,9.026761d0,20.741032d0,2.919823d0/
  data Ps/2236.097919d0,58.843178d0,0.058126d0,32.398448d0,5.334440d0,&
          0.623309d0/
  data Ss/2732.954806d0,77.644268d0,0.133756d0,34.764797d0,5.847848d0,&
          0.774454d0/
  data Cls/3298.077636d0,100.506181d0,0.271754d0,37.143139d0,6.376868d0,&
          0.929900d0/
  data Ars/3936.019052d0,128.032322d0,0.500810d0,39.536820d0,6.923455d0,&
          1.087954d0/
  
  rhoapprox = 0.d0
  do i = 1, site_info%natoms
     dx = site_info%site_crds(3*(i-1)+1) - xin
     if (abs(dx) .lt. tiny) dx = 0.0d0
     dy = site_info%site_crds(3*(i-1)+2) - yin
     if (abs(dy) .lt. tiny) dy = 0.0d0
     dz = site_info%site_crds(3*(i-1)+3) - zin
     if (abs(dz) .lt. tiny) dz = 0.0d0
     dist = dsqrt(dx*dx+dy*dy+dz*dz)
     rhoatom = 0.d0
     if (site_info%nuclear_charge(i) .eq. 1.d0) then 
        rhoatom = rhoatom + Hs(1,1)*exp(mone*Hs(2,1)*dist)
     else if (site_info%nuclear_charge(i) .eq. 0.d0) then 
        rhoatom = 0.0d0
     else if (site_info%nuclear_charge(i) .eq. 6.d0) then 
        do j = 1, 2
           rhoatom = rhoatom + Cs(j,1)*exp(mone*Cs(j,2)*dist)
        enddo
     else if (site_info%nuclear_charge(i) .eq. 7.d0) then 
        do j = 1, 2
           rhoatom = rhoatom + Ns(j,1)*exp(mone*Ns(j,2)*dist)
        enddo
     else if (site_info%nuclear_charge(i) .eq. 8.d0) then 
        do j = 1, 2
           rhoatom = rhoatom + Os(j,1)*exp(mone*Os(j,2)*dist)
        enddo
     else if (site_info%nuclear_charge(i) .eq. 9.d0) then 
        do j = 1, 2
           rhoatom = rhoatom + Fs(j,1)*exp(mone*Fs(j,2)*dist)
        enddo
     else if (site_info%nuclear_charge(i) .eq. 10.d0) then 
        do j = 1, 2
           rhoatom = rhoatom + Nes(j,1)*exp(mone*Nes(j,2)*dist)
        enddo
     else if (site_info%nuclear_charge(i) .eq. 15.d0) then 
        do j = 1, 3
           rhoatom = rhoatom + Ps(j,1)*exp(mone*Ps(j,2)*dist)
        enddo
     else if (site_info%nuclear_charge(i) .eq. 16.d0) then 
        do j = 1, 3
           rhoatom = rhoatom + Ss(j,1)*exp(mone*Ss(j,2)*dist)
        enddo
     else if (site_info%nuclear_charge(i) .eq. 17.d0) then 
        do j = 1, 3
           rhoatom = rhoatom + Cls(j,1)*exp(mone*Cls(j,2)*dist)
        enddo
     else if (site_info%nuclear_charge(i) .eq. 18.d0) then 
        do j = 1, 3
           rhoatom = rhoatom + Ars(j,1)*exp(mone*Ars(j,2)*dist)
        enddo
     else
        write(6,*)'promol_rho: unknown atom type',site_info%nuclear_charge(i),&
                   &' exiting'
        stop
     endif
     rhoapprox = rhoapprox + rhoatom
  enddo

end subroutine promol_rho
!----------------------------------------------------------
