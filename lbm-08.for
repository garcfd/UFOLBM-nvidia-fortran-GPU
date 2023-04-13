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
        REAL,     ALLOCATABLE :: vel(:,:,:,:)
        REAL,     ALLOCATABLE :: fin(:,:,:,:)
        REAL,     ALLOCATABLE :: fout(:,:,:,:)
        REAL,     ALLOCATABLE :: feq(:,:,:,:)
        REAL,     ALLOCATABLE :: rho(:,:,:)
        REAL,     ALLOCATABLE :: nut(:,:,:)
      end module allocated
C-----------------------------------------------------------------------

      program main
      use allocated
      implicit none

CCC   integer, parameter :: ndim=2, nvec=9
      integer, parameter :: ndim=3, nvec=15, nrow=5

      integer n,c,it,nits,im,jm,km,nsav,nout
      integer vec(nvec,ndim),col1(nrow),col2(nrow),col3(nrow)
      integer nexti,nextj,nextk,ii,jj,kk,nn,obsmax
      integer imin,jmin,kmin,imax,jmax,kmax

      real Re,Lref,uLB,nuLB,omega,sum2,sum3,cu,usqr
      real cx,cy,cz,cr,xx,yy,zz,wt(nvec),vin(ndim)
      real vmag,vmax,smag,srinv,rdx,rdy,rdz
      real dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

      integer i,j,k,ni,nj,nk,sa
      common  i,j,k,ni,nj,nk,sa

      CHARACTER(30) outfile
      CHARACTER(6)  string

C-----read input file
      open(20,file="lbm-inputs.txt")
      read(20,*) ni,nj,nk  ! grid dimensions
      read(20,*) nsav      ! max saves
      read(20,*) nits      ! max iterations
      read(20,*) nout      ! write every nout
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
c     cx = real(ni)/4.0   ! sphere x-coord
c     cy = real(nj)/2.0   ! sphere y-coord
c     cz = real(nk)/2.0   ! sphere z-coord
c     cr = real(nj)/10.0  ! sphere radius, in terms of nj dimension

      nuLB  = uLB*Lref/Re ! viscosity in lattice units
      omega = 1.0/(3.0*nuLB + 0.5) ! relaxation parameter

      vin = [uLB,0.0,0.0] ! inlet velocity

C-----write data to screen
      write(6,'(1X,A,3I5)')'ni nj nk  = ',ni,nj,nk
      write(6,'(1X,A,I10)')'mesh size = ',ni*nj*nk
      write(6,'(1X,A,I8)') 'nsav = ',nsav
      write(6,'(1X,A,I8)') 'nits = ',nits
      write(6,'(1X,A,I8)') 'nout = ',nout
      write(6,*) 'Re   = ',Re
      write(6,*) 'Lref = ',Lref
      write(6,*) 'uLB  = ',uLB
      write(6,*) 'nuLB = ',nuLB
      write(6,*) 'smag = ',smag
      write(6,*) 'omega = ',omega

C=====allocate
      ALLOCATE( obs(ni,nj,nk) )      ! 3D
      ALLOCATE( vel(ni,nj,nk,ndim) ) ! 3D
      ALLOCATE( fin(ni,nj,nk,nvec) ) ! 3D
      ALLOCATE(fout(ni,nj,nk,nvec) ) ! 3D
      ALLOCATE( feq(ni,nj,nk,nvec) ) ! 3D
      ALLOCATE( rho(ni,nj,nk) )      ! 3D
      ALLOCATE( nut(ni,nj,nk) )      ! 3D

C=====initialise vectors
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

C=====initialise weights

      wt(1)  = 1./9.
      wt(2)  = 1./9.
      wt(3)  = 1./9.
      wt(4)  = 1./72.
      wt(5)  = 1./72.
      wt(6)  = 1./72.
      wt(7)  = 1./72.
      wt(8)  = 2./9.
      wt(9)  = 1./72.
      wt(10) = 1./72.
      wt(11) = 1./72.
      wt(12) = 1./72.
      wt(13) = 1./9.
      wt(14) = 1./9.
      wt(15) = 1./9.

      col1 = (/ 1,  4,  5,  6,   7 /)
      col2 = (/ 2,  3,  8,  13, 14 /)
      col3 = (/ 9,  10, 11, 12, 15 /)

      obs = 0    ! int  array
      fin = 0.0  ! real array
      fout= 0.0  ! real array
      feq = 0.0  ! real array
      rho = 0.0  ! real array
      vel = 0.0  ! real array
      nut = 0.0  ! real array



C-----original sphere geometry
      if (.false.) then       
      do i = 1,ni
      do j = 1,nj
      do k = 1,nk
        xx = real(i-1) + 0.5
        yy = real(j-1) + 0.5
        zz = real(k-1) + 0.5
        if (((xx-cx)**2+(yy-cy)**2+(zz-cz)**2) .LT. (cr**2)) then
          obs(i,j,k) = 0
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
c         obs(i,nj,k) = 1 ! ymax
        enddo
c        do j = 1,nj ! side walls
c          obs(i,j, 1) = 1 ! zmin
c          obs(i,j,nk) = 1 ! zmax
c        enddo
      enddo
      endif


C-----read obstacle file from ufocfd
      if (.true.) then
        write(6,*) "read obstacle"
        open(20,file="obstacle.dat")
        read(20,*) ! im jm km header
        read(20,*) obsmax
        do n = 1,obsmax
          read(20,*) i,j,k
          obs(i,j,k) = 1
        enddo
        close(20)
        write(6,*)"done"
      endif



C=====initial velocity field
      write(6,*)"initial velocity"
      do i = 1,ni
      do j = 1,nj
      do k = 1,nk

        if (obs(i,j,k).eq.1) then ! solid
          vel(i,j,k,1) = 0.0
          vel(i,j,k,2) = 0.0
          vel(i,j,k,3) = 0.0
        endif

        if (obs(i,j,k).eq.0) then ! fluid
          vel(i,j,k,1) = vin(1)
          vel(i,j,k,2) = vin(2)
          vel(i,j,k,3) = vin(3)
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


!$acc data copy(obs,vel,fin,fout,feq,rho,vec,wt,col1,col2,col3,vin,nut)

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

C     STEP3 - left wall (and upper wall) inflow condition
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

C     STEP6a - calc smagorinsky turbulence
      imin = i-1; imax = i+1; rdx = 0.5
      jmin = j-1; jmax = j+1; rdy = 0.5
      kmin = k-1; kmax = k+1; rdz = 0.5

      if (i.eq.1)  imin =  1; if (i.eq.1)  rdx = 1.0
      if (j.eq.1)  jmin =  1; if (j.eq.1)  rdy = 1.0
      if (k.eq.1)  kmin =  1; if (k.eq.1)  rdz = 1.0
      if (i.eq.ni) imax = ni; if (i.eq.ni) rdx = 1.0
      if (j.eq.nj) jmax = nj; if (j.eq.nj) rdy = 1.0
      if (k.eq.nk) kmax = nk; if (k.eq.nk) rdz = 1.0

      dudx = rdx*(vel(imax,j,k,1) - vel(imin,j,k,1))
      dudy = rdy*(vel(i,jmax,k,1) - vel(i,jmin,k,1))
      dudz = rdz*(vel(i,j,kmax,1) - vel(i,j,kmin,1))
      dvdx = rdx*(vel(imax,j,k,2) - vel(imin,j,k,2))
      dvdy = rdy*(vel(i,jmax,k,2) - vel(i,jmin,k,2))
      dvdz = rdz*(vel(i,j,kmax,2) - vel(i,j,kmin,2))
      dwdx = rdx*(vel(imax,j,k,3) - vel(imin,j,k,3))
      dwdy = rdy*(vel(i,jmax,k,3) - vel(i,jmin,k,3))
      dwdz = rdz*(vel(i,j,kmax,3) - vel(i,j,kmin,3))

      srinv = sqrt( 2*dudx*dudx + 2*dvdy*dvdy + 2*dwdz*dwdz +
     &         (dudz+dwdx)**2 + (dudy+dvdx)**2 + (dwdy+dvdz)**2 )

      nut(i,j,k) = smag*smag*srinv

      omega = 1.0/(3.0*(nuLB+nut(i,j,k)) + 0.5) ! relaxation parameter

C     STEP6b - collision step
        do n = 1,nvec
          fout(i,j,k,n) = fin(i,j,k,n)
     &           - omega*(fin(i,j,k,n)
     &                  - feq(i,j,k,n))
        enddo

C     STEP7 - obstacle bounce-back condition
        if (obs(i,j,k).eq.1) then ! solid
          do n = 1,nvec
            fout(i,j,k,n) = fin(i,j,k,16-n)
          enddo
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

          fin(nexti,nextj,nextk,n) = fout(i,j,k,n)

        enddo

      enddo
      enddo
      enddo

C-----write monitor
      if (mod(it,nout).eq.0) then
        write(6,10)" sa = ",sa," it = ",it," vx = ",
     &  vel(int(0.8*ni),int(0.5*nj),int(0.5*nk),1)
      endif

      enddo ! nits

      call write_binary_vtk()

      enddo ! nsav
C=====end main iteration loop
!$acc end data


   10 format(A,I3,A,I6,A,F8.4)


C-----calc max velocity
      vmax = 0.0
      do k = 1,nk
      do j = 1,nj
      do i = 1,ni
      do n = 1,3
        vmag = vel(i,j,k,n)
        if (vmag.gt.vmax) then
          vmax = vmag
          ii = i
          jj = j
          kk = k
          nn = n
        endif
      enddo
      enddo
      enddo
      enddo
      write(6,*) "vmax = ",vmax," at ",ii,jj,kk,"dir = ",nn

      end ! main



      subroutine write_ascii_vtk()
      use allocated
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
      write(20,10)'SCALARS obstacle int'
      write(20,10)'LOOKUP_TABLE default'
      write(20,*)(((obs(i,j,k),i=1,ni),j=1,nj),k=1,nk)
C-----rho
      write(20,10)'SCALARS density float'
      write(20,10)'LOOKUP_TABLE default'
      write(20,*)(((rho(i,j,k),i=1,ni),j=1,nj),k=1,nk)
C-----vel
      write(20,10)'VECTORS velocity float'
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

      integer i,j,k,ni,nj,nk,sa
      common  i,j,k,ni,nj,nk,sa

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
      write(20)'SCALARS obstacle int'//lf
      write(20)'LOOKUP_TABLE default'//lf
      write(20)(((obs(i,j,k),i=1,ni),j=1,nj),k=1,nk)
C-----rho
      write(20)'SCALARS density float'//lf
      write(20)'LOOKUP_TABLE default'//lf
      write(20)(((rho(i,j,k),i=1,ni),j=1,nj),k=1,nk)
C-----nut
      write(20)'SCALARS nut float'//lf
      write(20)'LOOKUP_TABLE default'//lf
      write(20)(((nut(i,j,k),i=1,ni),j=1,nj),k=1,nk)
C-----vel
      write(20)'VECTORS velocity float'//lf
      write(20)(((vel(i,j,k,1),vel(i,j,k,2),vel(i,j,k,3),
     &             i=1,ni),j=1,nj),k=1,nk)
      close(20)
      end




