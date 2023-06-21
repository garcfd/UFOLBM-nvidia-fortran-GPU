C-----------------------------------------------------------------------
C     based on 2D python code by Jonas Latt Geneva University Coursera
C     code developed by Giles Richardson (GARCFD) in Cambridgeshire UK
C-----------------------------------------------------------------------
C     export NV_ACC_TIME=1
C     export ACC_NUM_CORES=8
C     export NVHPC_CUDA_HOME="/usr/local/cuda-12.1"
C
C     pgf77 -static-nvidia -acc=host      -o lbm-basic lbm-01.for
C     pgf77 -static-nvidia -acc=multicore -o lbm-mcore lbm-01.for
C     pgf77 -static-nvidia -acc -gpu=cc75 -o lbm-tesla lbm-01.for
C-----------------------------------------------------------------------
      module allocated
        INTEGER,  ALLOCATABLE :: obs(:,:,:)
        INTEGER,  ALLOCATABLE :: cellnrm(:,:,:)
        INTEGER,  ALLOCATABLE :: nrm26(:)
        INTEGER,  ALLOCATABLE :: ref(:,:)

        REAL,     ALLOCATABLE :: vel(:,:,:,:)
        REAL,     ALLOCATABLE :: fin(:,:,:,:)
        REAL,     ALLOCATABLE :: fout(:,:,:,:)
        REAL,     ALLOCATABLE :: feq(:,:,:,:)
        REAL,     ALLOCATABLE :: rho(:,:,:)
        REAL,     ALLOCATABLE :: nut(:,:,:)
        REAL,     ALLOCATABLE :: nrm(:,:)
      end module allocated
C-----------------------------------------------------------------------

      program main
      use allocated
      implicit none

CCC   integer, parameter :: ndim=2, nvec=9
      integer, parameter :: ndim=3, nvec=15, nrow=5, nvec0=26

      integer i,j,k,n,it,nits,nsav,nout,lcd,tang
      integer vec(nvec,ndim),col1(nrow),col2(nrow),col3(nrow)
      integer nexti,nextj,nextk,obsmax,nrmmax,inrm,v0,v0save
      integer imin,jmin,kmin,imax,jmax,kmax,reflect

      real vec0(nvec0,ndim)
      real Re,Lref,uLB,nuLB,omega,sum2,sum3,cu,usqr
      real cx,cy,cz,cr,xx,yy,zz,wt(nvec),vin(ndim)
      real smag,srinv,rdx,rdy,rdz,dx,dy,dz
      real dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
      real dist,distmin,nrm1,nrm2,nrm3,rr2,rr3

C-----common
      integer ni,nj,nk,sa
      common  ni,nj,nk,sa

C-----read input file
      open(20,file="lbm-inputs.txt")
      read(20,*) ni,nj,nk  ! grid dimensions
      read(20,*) nsav      ! max saves
      read(20,*) nits      ! max iterations
      read(20,*) nout      ! write every nout
      read(20,*) tang      ! tangential on/off
      read(20,*) Re        ! Re
      read(20,*) Lref      ! Lref
      read(20,*) uLB       ! velocity in lattice units
      read(20,*) smag      ! smagorinsky constant
      close(20)

C-----constants
c     ni = 1000           ! lattice nodes
c     nj = 225            ! lattice nodes
c     nk = 400            ! lattice nodes
c     nits = 100          ! max iterations
      cx = real(ni)/4.0   ! sphere x-coord
      cy = real(nj)/2.0   ! sphere y-coord
      cz = real(nk)/2.0   ! sphere z-coord
      cr = real(nj)/5.0   ! sphere radius, in terms of nj dimension

      nuLB  = uLB*Lref/Re ! viscosity in lattice units
      omega = 1.0/(3.0*nuLB + 0.5) ! relaxation parameter
      vin = [uLB,0.0,0.0] ! inlet velocity

C-----write data to screen
      write(6,'(1X,A,3I5)')'ni nj nk  = ',ni,nj,nk
      write(6,'(1X,A,I10)')'mesh size = ',ni*nj*nk
      write(6,'(1X,A,I8)') 'nsav = ',nsav
      write(6,'(1X,A,I8)') 'nits = ',nits
      write(6,'(1X,A,I8)') 'nout = ',nout
      write(6,'(1X,A,I8)') 'tang = ',tang
      write(6,*) 'Re   = ',Re
      write(6,*) 'Lref = ',Lref
      write(6,*) 'uLB  = ',uLB
      write(6,*) 'nuLB = ',nuLB
      write(6,*) 'smag = ',smag
      write(6,*) 'omega = ',omega

C=====allocate
      ALLOCATE( obs(ni,nj,nk) )       ! int
      ALLOCATE(cellnrm(ni,nj,nk) )    ! int
      ALLOCATE( vel(ni,nj,nk,ndim) )  ! real
      ALLOCATE( fin(ni,nj,nk,nvec) )  ! real
      ALLOCATE(fout(ni,nj,nk,nvec) )  ! real
      ALLOCATE( feq(ni,nj,nk,nvec) )  ! real
      ALLOCATE( rho(ni,nj,nk) )       ! real
      ALLOCATE( nut(ni,nj,nk) )       ! real

C=====initialise vec vectors
      vec(1, :) = (/ 1, 0, 0 /)
      vec(2, :) = (/ 0, 1, 0 /)
      vec(3, :) = (/ 0, 0, 1 /)
      vec(4, :) = (/ 1, 1, 1 /) 
      vec(5, :) = (/ 1, 1,-1 /)
      vec(6, :) = (/ 1,-1, 1 /) 
      vec(7, :) = (/ 1,-1,-1 /) 
      vec(8, :) = (/ 0, 0, 0 /)  
      vec(9, :) = (/-1, 1, 1 /)
      vec(10,:) = (/-1, 1,-1 /)
      vec(11,:) = (/-1,-1, 1 /)
      vec(12,:) = (/-1,-1,-1 /)
      vec(13,:) = (/ 0, 0,-1 /)
      vec(14,:) = (/ 0,-1, 0 /) 
      vec(15,:) = (/-1, 0, 0 /)

      rr2 = 0.707 ! 1/sqrt(2)
      rr3 = 0.577 ! 1/sqrt(3)

C=====initialise vec0 vectors
      vec0(1, :) = (/  1.0,  0.0,  0.0 /) ! 6 main axes
      vec0(2, :) = (/ -1.0,  0.0,  0.0 /)
      vec0(3, :) = (/  0.0,  1.0,  0.0 /)
      vec0(4, :) = (/  0.0, -1.0,  0.0 /)
      vec0(5, :) = (/  0.0,  0.0,  1.0 /) 
      vec0(6, :) = (/  0.0,  0.0, -1.0 /) 

      vec0(7, :) = (/  rr2,  rr2,  0.0 /) ! 12 in-plane diagonals
      vec0(8, :) = (/  rr2, -rr2,  0.0 /)
      vec0(9, :) = (/ -rr2,  rr2,  0.0 /)
      vec0(10,:) = (/ -rr2, -rr2,  0.0 /)
      vec0(11,:) = (/  rr2,  0.0,  rr2 /)
      vec0(12,:) = (/  rr2,  0.0, -rr2 /)
      vec0(13,:) = (/ -rr2,  0.0,  rr2 /)
      vec0(14,:) = (/ -rr2,  0.0, -rr2 /)
      vec0(15,:) = (/  0.0,  rr2,  rr2 /)
      vec0(16,:) = (/  0.0,  rr2, -rr2 /)
      vec0(17,:) = (/  0.0, -rr2,  rr2 /)
      vec0(18,:) = (/  0.0, -rr2, -rr2 /)

      vec0(19,:) = (/  rr3,  rr3,  rr3 /) ! 8 out-of-plane diagonals
      vec0(20,:) = (/ -rr3,  rr3,  rr3 /)
      vec0(21,:) = (/  rr3, -rr3,  rr3 /) 
      vec0(22,:) = (/ -rr3, -rr3,  rr3 /) 
      vec0(23,:) = (/  rr3,  rr3, -rr3 /)
      vec0(24,:) = (/ -rr3,  rr3, -rr3 /)
      vec0(25,:) = (/  rr3, -rr3, -rr3 /)
      vec0(26,:) = (/ -rr3, -rr3, -rr3 /)
 
C=====initialise weights
      wt(1)  = 1.0 /  9.0
      wt(2)  = 1.0 /  9.0
      wt(3)  = 1.0 /  9.0
      wt(4)  = 1.0 / 72.0
      wt(5)  = 1.0 / 72.0
      wt(6)  = 1.0 / 72.0
      wt(7)  = 1.0 / 72.0
      wt(8)  = 2.0 /  9.0
      wt(9)  = 1.0 / 72.0
      wt(10) = 1.0 / 72.0
      wt(11) = 1.0 / 72.0
      wt(12) = 1.0 / 72.0
      wt(13) = 1.0 /  9.0
      wt(14) = 1.0 /  9.0
      wt(15) = 1.0 /  9.0

      col1 = (/ 1,  4,  5,  6,   7 /)
      col2 = (/ 2,  3,  8,  13, 14 /)
      col3 = (/ 9,  10, 11, 12, 15 /)

      obs = 1     ! int array (default = 1)
      cellnrm = 0 ! int array (default = 0)
      fin = 0.0   ! real array
      fout= 0.0   ! real array
      feq = 0.0   ! real array
      rho = 0.0   ! real array
      vel = 0.0   ! real array
      nut = 0.0   ! real array


C-----original sphere geometry
      if (.false.) then       
      do i = 1,ni
      do j = 1,nj
      do k = 1,nk
        xx = real(i-1) + 0.5
        yy = real(j-1) + 0.5
        zz = real(k-1) + 0.5
        dx = xx-cx
        dy = yy-cy
        dz = 0.0 ! dz = zz-cz
        if ((dx*dx + dy*dy + dz*dz) .LT. (cr*cr)) then
          obs(i,j,k) = 1
        endif
      enddo
      enddo
      enddo
      endif


C-----set walls (if not periodic)
      if (.false.) then
      do i = 1,ni
        do k = 1,nk
          obs(i, 1,k) = 1 ! ymin
          obs(i,nj,k) = 1 ! ymax
        enddo
        do j = 1,nj ! side walls
          obs(i,j, 1) = 1 ! zmin
          obs(i,j,nk) = 1 ! zmax
        enddo
      enddo
      endif


C-----read obstacle.dat file from ufocfd
C-----read  normals.dat file from ufocfd
      if (.true.) then
        write(6,*) "read obstacle.dat and normals.dat"
        open(20,file="obstacle.dat")
        open(30,file="normals.dat")

        read(20,*) ! im jm km header
        read(20,*) obsmax
        read(30,*) nrmmax

        allocate( nrm(nrmmax,ndim) )
        allocate( ref(nvec0, nvec) )
        allocate( nrm26(nrmmax) )

        nrm(:,:) = 0.0
        ref(:,:) = 0
        nrm26(:) = 0

        inrm = 0

        do n = 1,obsmax

          read(20,*) lcd,i,j,k
          obs(i,j,k) = lcd ! (0 or -1)

          if (lcd.eq.0) then ! cutc
            inrm = inrm + 1
            cellnrm(i,j,k) = inrm
            read(30,*) nrm1,nrm2,nrm3
            nrm(inrm,1) = nrm1
            nrm(inrm,2) = nrm2
            nrm(inrm,3) = nrm3
          endif

        enddo

        close(20)
        close(30)

        write(6,*)"obsmax = ", obsmax
        write(6,*)"nrmmax = ", nrmmax
        write(6,*)"inrm =   ", inrm
        write(6,*)"done"
      endif


C=====convert simplify nrm(n,:) array to nrm26() array
      do n = 1,nrmmax

        nrm1 = nrm(n,1)
        nrm2 = nrm(n,2)
        nrm3 = nrm(n,3)
 
        distmin = 10.0

        do v0 = 1,nvec0

            dx = nrm1 - vec0(v0,1) 
            dy = nrm2 - vec0(v0,2) 
            dz = nrm3 - vec0(v0,3) 
            dist = sqrt(dx*dx + dy*dy + dz*dz)

            if (dist.lt.distmin) then
              v0save  = v0
              distmin = dist 
            endif

        enddo  

        nrm26(n) = v0save

      enddo




C=====specify the reflections

      ref(1, :)=(/ 15, 2, 3, 9, 10, 11, 12, 8, 4, 5, 6, 7, 13, 14, 1 /)
      ref(2, :)=(/ 15, 2, 3, 9, 10, 11, 12, 8, 4, 5, 6, 7, 13, 14, 1 /)
      ref(3, :)=(/ 1, 14, 3, 6, 7, 4, 5, 8, 11, 12, 9, 10, 13, 2, 15 /)
      ref(4, :)=(/ 1, 14, 3, 6, 7, 4, 5, 8, 11, 12, 9, 10, 13, 2, 15 /)
      ref(5, :)=(/ 1, 2, 13, 5, 4, 7, 6, 8, 10, 9, 12, 11, 3, 14, 15 /)
      ref(6, :)=(/ 1, 2, 13, 5, 4, 7, 6, 8, 10, 9, 12, 11, 3, 14, 15 /)
      ref(7, :)=(/ 14, 15, 3, 11, 12, 6, 7, 8, 9, 10, 4, 5, 13, 1, 2 /)
      ref(8, :)=(/ 2, 1, 3, 4, 5, 9, 10, 8, 6, 7, 11, 12, 13, 15, 14 /)
      ref(9, :)=(/ 2, 1, 3, 4, 5, 9, 10, 8, 6, 7, 11, 12, 13, 15, 14 /)
      ref(10,:)=(/ 14, 15, 3, 11, 12, 6, 7, 8, 9, 10, 4, 5, 13, 1, 2 /)
      ref(11,:)=(/ 13, 2, 15, 10, 5, 12, 7, 8, 9, 4, 11, 6, 1, 14, 3 /)
      ref(12,:)=(/ 3, 2, 1, 4, 9, 6, 11, 8, 5, 10, 7, 12, 15, 14, 13 /)
      ref(13,:)=(/ 3, 2, 1, 4, 9, 6, 11, 8, 5, 10, 7, 12, 15, 14, 13 /)
      ref(14,:)=(/ 13, 2, 15, 10, 5, 12, 7, 8, 9, 4, 11, 6, 1, 14, 3 /)
      ref(15,:)=(/ 1, 13, 14, 7, 5, 6, 4, 8, 12, 10, 11, 9, 2, 3, 15 /)
      ref(16,:)=(/ 1, 3, 2, 4, 6, 5, 7, 8, 9, 11, 10, 12, 14, 13, 15 /)
      ref(17,:)=(/ 1, 3, 2, 4, 6, 5, 7, 8, 9, 11, 10, 12, 14, 13, 15 /)
      ref(18,:)=(/ 1, 13, 14, 7, 5, 6, 4, 8, 12, 10, 11, 9, 2, 3, 15 /)
      ref(19,:)=(/ 7, 10, 11, 12, 13, 14, 1, 8, 15, 2, 3, 4, 5, 6, 9 /)
      ref(20,:)=(/ 4, 5, 6, 1, 2, 3, 9, 8, 7, 13, 14, 15, 10, 11, 12 /)
      ref(21,:)=(/ 5, 4, 9, 2, 1, 10, 13, 8, 3, 6, 15, 14, 7, 12, 11 /)
      ref(22,:)=(/ 6, 9, 4, 3, 11, 1, 14, 8, 2, 15, 5, 13, 12, 7, 10 /)
      ref(23,:)=(/ 6, 9, 4, 3, 11, 1, 14, 8, 2, 15, 5, 13, 12, 7, 10 /)
      ref(24,:)=(/ 5, 4, 9, 2, 1, 10, 13, 8, 3, 6, 15, 14, 7, 12, 11 /)
      ref(25,:)=(/ 4, 5, 6, 1, 2, 3, 9, 8, 7, 13, 14, 15, 10, 11, 12 /)
      ref(26,:)=(/ 7, 10, 11, 12, 13, 14, 1, 8, 15, 2, 3, 4, 5, 6, 9 /)






C=====initial velocity field
      write(6,*)"initial velocity"
      do i = 1,ni
      do j = 1,nj
      do k = 1,nk

        if (obs(i,j,k).eq.1) then ! fluid
          vel(i,j,k,1) = vin(1)
          vel(i,j,k,2) = vin(2)
          vel(i,j,k,3) = vin(3)
        endif

        if (obs(i,j,k).le.0) then ! cutcell or solid
          vel(i,j,k,1) = 0.0
          vel(i,j,k,2) = 0.0
          vel(i,j,k,3) = 0.0
        endif

      enddo
      enddo
      enddo


C=====equilibrium distribution function
C-----fin = equilibrium(1,u)
      do k = 1,nk
      do j = 1,nj
      do i = 1,ni

        usqr = vel(i,j,k,1)**2 + vel(i,j,k,2)**2 + vel(i,j,k,3)**2
        do n = 1,nvec
          cu = vec(n,1)*vel(i,j,k,1)
     &       + vec(n,2)*vel(i,j,k,2)
     &       + vec(n,3)*vel(i,j,k,3)
C-----------------------(rho=1)
          fin(i,j,k,n) = 1.0*wt(n)*(1 + 3*cu + 4.5*cu**2 - 1.5*usqr)
        enddo

      enddo
      enddo
      enddo

      write(6,*)"start iterations"


!$acc data copy(obs,vel,fin,fout,feq,rho)
!$acc data copy(vec,wt,col1,col2,col3,vin)
!$acc data copy(nut,cellnrm,nrm26,ref)

      do sa = 1, nsav
C=====main iteration loop
      do it = 1, nits

!$acc kernels loop independent
      do k = 1,nk
!$acc loop independent
      do j = 1,nj
!$acc loop independent
      do i = 1,ni

C     STEP1 - right wall outflow condition
      if (i.eq.ni) then
         do n = 1,nrow
           fin(ni,j,k,col3(n)) = fin(ni-1,j,k,col3(n))
         enddo
      endif

      rho(i,j,k) = 0.0
      vel(i,j,k,:) = 0.0

C     STEP2 - compute macroscopic variables rho and u
      do n = 1,nvec
        rho(i,j,k) = rho(i,j,k) + fin(i,j,k,n)
        vel(i,j,k,1) = vel(i,j,k,1) + vec(n,1) * fin(i,j,k,n)
        vel(i,j,k,2) = vel(i,j,k,2) + vec(n,2) * fin(i,j,k,n)
        vel(i,j,k,3) = vel(i,j,k,3) + vec(n,3) * fin(i,j,k,n)
      enddo
      
      vel(i,j,k,1) = vel(i,j,k,1) / rho(i,j,k)
      vel(i,j,k,2) = vel(i,j,k,2) / rho(i,j,k)
      vel(i,j,k,3) = vel(i,j,k,3) / rho(i,j,k)

CCCC  if (obs(i,j,k).eq.-1) vel(i,j,k,:) = 0.0

C     STEP3 - left wall (and upper/lower wall) inflow condition
      if ((i.eq.1).or.(j.eq.1).or.(j.eq.nj)) then
        vel(1,j,k,1) = vin(1)
        vel(1,j,k,2) = vin(2)
        vel(1,j,k,3) = vin(3)
        sum2 = 0.0
        sum3 = 0.0
        do n = 1,nrow
          sum2 = sum2 + fin(1,j,k,col2(n))
          sum3 = sum3 + fin(1,j,k,col3(n))
          rho(1,j,k) = (sum2 + 2.0*sum3) / (1.0-vel(1,j,k,1))
        enddo
      endif

C     STEP4 - compute equilibrium (rho, vel)
        usqr = vel(i,j,k,1)**2 + vel(i,j,k,2)**2 + vel(i,j,k,3)**2
        do n = 1,nvec
          cu = vec(n,1)*vel(i,j,k,1)
     &       + vec(n,2)*vel(i,j,k,2)
     &       + vec(n,3)*vel(i,j,k,3)
          feq(i,j,k,n) = rho(i,j,k)*wt(n)*
     &    (1.0 + 3.0*cu + 4.5*cu**2.0 - 1.5*usqr)
        enddo

C     STEP5 - calculate populations (at inlet)       
      if (i.eq.1) then
        do n = 1,nrow
          fin(1,j,k,col1(n)) = feq(1,j,k,col1(n))
     &                       + fin(1,j,k,col3(6-n))
     &                       - feq(1,j,k,col3(6-n))
        enddo
      endif

      omega = 1.0/(3.0*(nuLB+nut(i,j,k)) + 0.5) ! relaxation parameter

C     STEP6 - collision step
      do n = 1,nvec
        fout(i,j,k,n) = fin(i,j,k,n)
     &         - omega*(fin(i,j,k,n) - feq(i,j,k,n))
      enddo

C     STEP7 - obstacle wall condition
      if (obs(i,j,k).eq.0) then ! cutcell

        if (tang.eq.0) then ! no-slip bounce
          do n = 1,nvec
            fout(i,j,k,n) = fin(i,j,k,16-n)
          enddo
        endif

        if (tang.eq.1) then ! free-slip reflection
          do n = 1,nvec
            reflect = ref( nrm26( cellnrm(i,j,k) ) ,n)
            fout(i,j,k,reflect) = fin(i,j,k,n)
          enddo
        endif

      endif

      enddo
      enddo
      enddo

!$acc kernels loop independent
      do k = 1,nk
!$acc loop independent
      do j = 1,nj
!$acc loop independent
      do i = 1,ni


C     calc smagorinsky turbulence
      if (smag .gt. 0.0) then

        imin = i-1; imax = i+1; rdx = 0.5
        jmin = j-1; jmax = j+1; rdy = 0.5
        kmin = k-1; kmax = k+1; rdz = 0.5

        if (i.eq.1) then
          imin = 1; rdx  = 1.0
        endif
        if (j.eq.1) then
          jmin = 1; rdy = 1.0
        endif
        if (k.eq.1) then
          kmin = 1; rdz = 1.0
        endif
        if (i.eq.ni) then
          imax = ni; rdx = 1.0
        endif
        if (j.eq.nj) then
          jmax = nj; rdy = 1.0
        endif
        if (k.eq.nk) then
          kmax = nk; rdz = 1.0
        endif

        dudx = rdx*(vel(imax,j,k,1) - vel(imin,j,k,1))
        dudy = rdy*(vel(i,jmax,k,1) - vel(i,jmin,k,1))
        dudz = rdz*(vel(i,j,kmax,1) - vel(i,j,kmin,1))

        dvdx = rdx*(vel(imax,j,k,2) - vel(imin,j,k,2))
        dvdy = rdy*(vel(i,jmax,k,2) - vel(i,jmin,k,2))
        dvdz = rdz*(vel(i,j,kmax,2) - vel(i,j,kmin,2))

        dwdx = rdx*(vel(imax,j,k,3) - vel(imin,j,k,3))
        dwdy = rdy*(vel(i,jmax,k,3) - vel(i,jmin,k,3))
        dwdz = rdz*(vel(i,j,kmax,3) - vel(i,j,kmin,3))

        srinv = sqrt( 2.0*dudx*dudx + 2.0*dvdy*dvdy + 2.0*dwdz*dwdz +
     &        (dudz+dwdx)**2.0 + (dudy+dvdx)**2.0 + (dwdy+dvdz)**2.0 )

      else
        srinv = 0.0
      endif

      nut(i,j,k) = smag*smag*srinv


C     STEP8 - streaming step
        do n = 1,nvec

          nexti = i + vec(n,1)
          nextj = j + vec(n,2)
          nextk = k + vec(n,3)

C---------periodic boundaries
          if (.true.) then
            if (nexti.lt.1)  nexti = ni
            if (nextj.lt.1)  nextj = nj
            if (nextk.lt.1)  nextk = nk
            if (nexti.gt.ni) nexti = 1
            if (nextj.gt.nj) nextj = 1
            if (nextk.gt.nk) nextk = 1
          endif

          if ( obs(nexti,nextj,nextk).ge.0 ) then          
            fin(nexti,nextj,nextk,n) = fout(i,j,k,n)
          endif
        enddo

      enddo
      enddo
      enddo

C-----write monitor
      if (mod(it,nout).eq.0) then
        write(6,10)" sa = ",sa," it = ",it," vx = ",
     &  vel(int(1.0*ni),int(0.5*nj),int(0.5*nk),1)
      endif

      enddo ! nits

      call write_vxmax()
      call write_binary_vtk()

      enddo ! nsav
C=====end main iteration loop

!$acc end data
!$acc end data
!$acc end data

   10 format(A,I3,A,I6,A,F8.4)

      end ! main
C===============




      subroutine write_vxmax()
      use allocated
      implicit none

      integer i,j,k,ii,jj,kk
      real vxmag,vxmax
      integer ni,nj,nk
      common  ni,nj,nk

C-----calc max velocity
      vxmax = 0.0
      do k = 1,nk
      do j = 1,nj
      do i = 1,ni

        vxmag = vel(i,j,k,1)
        if (vxmag.gt.vxmax) then
          vxmax = vxmag
          ii = i
          jj = j
          kk = k
        endif

      enddo
      enddo
      enddo
      write(6,*) "vxmax = ",vxmax," at ",ii,jj,kk

      end




      subroutine write_ascii_vtk()
      use allocated
      implicit none

      integer i,j,k
      integer ni,nj,nk,sa
      common  ni,nj,nk,sa

C-----start ascii VTK file
      write(6,*)"write ascii vtk file"
      OPEN(unit=20,file='lbm.vtk')
      write(20,10)'# vtk DataFile Version 3.0'
      write(20,10)'vtk output'
      write(20,10)'ASCII'
      write(20,10)'DATASET RECTILINEAR_GRID'
      write(20,20)'DIMENSIONS ',ni+1,nj+1,nk+1
      write(20,30)'X_COORDINATES ',ni+1,' float'
      write(20,*)   (real(i-1),i=1,ni+1)
      write(20,30)'Y_COORDINATES ',nj+1,' float'
      write(20,*)   (real(j-1),j=1,nj+1)
      write(20,30)'Z_COORDINATES ',nk+1,' float'
      write(20,*)   (real(k-1),k=1,nk+1)
      write(20,40)'CELL_DATA ',ni*nj*nk
C-----obs
      write(20,10)'SCALARS obs int'
      write(20,10)'LOOKUP_TABLE default'
      write(20,*)(((obs(i,j,k),i=1,ni),j=1,nj),k=1,nk)
C-----rho
      write(20,10)'SCALARS rho float'
      write(20,10)'LOOKUP_TABLE default'
      write(20,*)(((rho(i,j,k),i=1,ni),j=1,nj),k=1,nk)
C-----nut
      write(20,10)'SCALARS nut float'
      write(20,10)'LOOKUP_TABLE default'
      write(20,*)(((nut(i,j,k),i=1,ni),j=1,nj),k=1,nk)
C-----vel
      write(20,10)'VECTORS vel float'
      write(20,*)(((vel(i,j,k,1),vel(i,j,k,2),vel(i,j,k,3),
     &                           i=1,ni),j=1,nj),k=1,nk)
      close(20)

   10 format(A)
   20 format(A,3I4)
   30 format(A,I3,A)
   40 format(A,I9)

      end



C-----write binary VTK file
      subroutine write_binary_vtk()

      use allocated
      implicit none

      integer i,j,k
      integer ni,nj,nk,sa
      common  ni,nj,nk,sa

      CHARACTER(10) outfile
      CHARACTER(3)  string
      CHARACTER(LEN=1)  :: lf
      CHARACTER(LEN=10) :: str1,str2,str3,str4
      lf = char(10)

      write(6,*)"write binary VTK"

C-----output filename
      write(unit=string, fmt='(I3.3)') sa
      outfile = 'lbm'//string//'.vtk'
      write(6,*)"outfile = ",outfile

      OPEN(unit=20, file=outfile, form='unformatted',
     &    access='stream',status='replace',convert="big_endian")
      write(20)'# vtk DataFile Version 3.0'//lf
      write(20)'vtk output'//lf
      write(20)'BINARY'//lf
      write(20)'DATASET RECTILINEAR_GRID'//lf
      write(str1(1:10),'(i10)') ni+1
      write(str2(1:10),'(i10)') nj+1
      write(str3(1:10),'(i10)') nk+1
      write(str4(1:10),'(i10)') ni*nj*nk
      write(20)'DIMENSIONS '//str1//str2//str3//lf
      write(20)'X_COORDINATES '//str1//' float'//lf
      write(20)(real(i-1),i=1,ni+1)
      write(20)'Y_COORDINATES '//str2//' float'//lf
      write(20)(real(j-1),j=1,nj+1)
      write(20)'Z_COORDINATES '//str3//' float'//lf
      write(20)(real(k-1),k=1,nk+1)
      write(20)'CELL_DATA '//str4//lf
C-----obs
      write(20)'SCALARS obs int'//lf
      write(20)'LOOKUP_TABLE default'//lf
      write(20)(((obs(i,j,k),i=1,ni),j=1,nj),k=1,nk)
C-----rho
      write(20)'SCALARS rho float'//lf
      write(20)'LOOKUP_TABLE default'//lf
      write(20)(((rho(i,j,k),i=1,ni),j=1,nj),k=1,nk)
C-----nut
      write(20)'SCALARS nut float'//lf
      write(20)'LOOKUP_TABLE default'//lf
      write(20)(((nut(i,j,k),i=1,ni),j=1,nj),k=1,nk)
C-----vel
      write(20)'VECTORS vel float'//lf
      write(20)(((vel(i,j,k,1),vel(i,j,k,2),vel(i,j,k,3),
     &             i=1,ni),j=1,nj),k=1,nk)
      close(20)
      end




