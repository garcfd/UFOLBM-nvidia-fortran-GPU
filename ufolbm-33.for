C-----------------------------------------------------------------------
C     based on 2D python code by Jonas Latt Geneva University Coursera
C     code developed by Giles Richardson (GARCFD) in Cambridgeshire UK
C     code available on github under garcfd with MIT license agreement
C-----------------------------------------------------------------------
C     export NV_ACC_TIME=1
C     export ACC_NUM_CORES=8
C     export NVHPC_CUDA_HOME="/usr/local/cuda-12.1"
C
C     pgf77 -static-nvidia -acc=host      -o lbm-basic lbm-01.for
C     pgf77 -static-nvidia -acc=multicore -o lbm-mcore lbm-01.for
C     pgf77 -static-nvidia -acc -gpu=cc75 -o lbm-tesla lbm-01.for # Nvidia Tesla T4
C     pgf77 -static-nvidia -acc -gpu=cc86 -o lbm-tesla lbm-20.for # Nvidia A10
C-----------------------------------------------------------------------
      module allocated
        INTEGER,  ALLOCATABLE :: obs(:,:,:)
        REAL,     ALLOCATABLE :: vel(:,:,:,:)
        REAL,     ALLOCATABLE :: fin(:,:,:,:)
        REAL,     ALLOCATABLE :: fout(:,:,:,:)
        REAL,     ALLOCATABLE :: feq(:,:,:,:)
        REAL,     ALLOCATABLE :: rho(:,:,:)
      end module allocated
C-----------------------------------------------------------------------

      program main
      use allocated
      implicit none

      integer, parameter :: ndim=3, nvec=27, nrow=9

      integer i,j,k,n,it,nits,nsav,nout,lcd
      integer nexti,nextj,nextk,obsmax
      integer vec(nvec,ndim),col1(nrow),col2(nrow),col3(nrow)

      real wt(nvec),vin(ndim)
      real Re,Lref,uLB,nuLB,omega,sum2,sum3,cu
      real cx,cy,cz,cr,xx,yy,zz,dx,dy,dz,usqr,cs2

C-----common
      integer ni,nj,nk,sa
      common  ni,nj,nk,sa

C-----read input file
      open(10,file="ufolbm.txt")
      read(10,*) ni,nj,nk  ! grid dimensions
      read(10,*) nsav      ! max saves
      read(10,*) nits      ! max iterations
      read(10,*) nout      ! write every nout
      read(10,*) re        ! Re
      read(10,*) lref      ! Lref
      read(10,*) ulb       ! velocity in lattice units
      close(10)

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
      cs2 = 1.0/3.0       ! speed of sound squared

C-----write data to screen
      write(6,'(1X,A,3I5)')'ni nj nk  = ',ni,nj,nk
      write(6,'(1X,A,I10)')'mesh size = ',ni*nj*nk
      write(6,'(1X,A,I8)') 'nsav = ',nsav
      write(6,'(1X,A,I8)') 'nits = ',nits
      write(6,'(1X,A,I8)') 'nout = ',nout
      write(6,*) 'Re    = ',Re
      write(6,*) 'Lref  = ',Lref
      write(6,*) 'uLB   = ',uLB
      write(6,*) 'nuLB  = ',nuLB
      write(6,*) 'omega = ',omega

C=====allocate
      ALLOCATE( obs(ni,nj,nk) )       ! int
      ALLOCATE( vel(ni,nj,nk,ndim) )  ! real
      ALLOCATE( fin(ni,nj,nk,nvec) )  ! real
      ALLOCATE(fout(ni,nj,nk,nvec) )  ! real
      ALLOCATE( feq(ni,nj,nk,nvec) )  ! real
      ALLOCATE( rho(ni,nj,nk) )       ! real

C=====initialise the lbm vectors
      vec(1, :) = (/ 1,-1,-1 /) ! r3
      vec(2, :) = (/ 1,-1, 0 /)
      vec(3, :) = (/ 1,-1, 1 /) ! r3
      vec(4, :) = (/ 1, 0,-1 /) 
      vec(5, :) = (/ 1, 0, 0 /) ! x+
      vec(6, :) = (/ 1, 0, 1 /) 
      vec(7, :) = (/ 1, 1,-1 /) ! r3
      vec(8, :) = (/ 1, 1, 0 /)  
      vec(9, :) = (/ 1, 1, 1 /) ! r3

      vec(10,:) = (/ 0,-1,-1 /)
      vec(11,:) = (/ 0,-1, 0 /) ! y-
      vec(12,:) = (/ 0,-1, 1 /)
      vec(13,:) = (/ 0, 0,-1 /) ! z-
      vec(14,:) = (/ 0, 0, 0 /) 
      vec(15,:) = (/ 0, 0, 1 /) ! z+
      vec(16,:) = (/ 0, 1,-1 /)
      vec(17,:) = (/ 0, 1, 0 /) ! y+
      vec(18,:) = (/ 0, 1, 1 /)

      vec(19,:) = (/-1,-1,-1 /) ! r3
      vec(20,:) = (/-1,-1, 0 /)
      vec(21,:) = (/-1,-1, 1 /) ! r3
      vec(22,:) = (/-1, 0,-1 /) 
      vec(23,:) = (/-1, 0, 0 /) ! x-
      vec(24,:) = (/-1, 0, 1 /)
      vec(25,:) = (/-1, 1,-1 /) ! r3
      vec(26,:) = (/-1, 1, 0 /)
      vec(27,:) = (/-1, 1, 1 /) ! r3
 
C=====initialise weights
      wt(1)  = 1.0 / 216.0 ! r3
      wt(2)  = 1.0 /  54.0 ! r2
      wt(3)  = 1.0 / 216.0 ! r3
      wt(4)  = 1.0 /  54.0 ! r2
      wt(5)  = 2.0 /  27.0 ! x++
      wt(6)  = 1.0 /  54.0 ! r2
      wt(7)  = 1.0 / 216.0 ! r3
      wt(8)  = 1.0 /  54.0 ! r2
      wt(9)  = 1.0 / 216.0 ! r3

      wt(10) = 1.0 /  54.0 ! r2
      wt(11) = 2.0 /  27.0 ! y--
      wt(12) = 1.0 /  54.0 ! r2
      wt(13) = 2.0 /  27.0 ! z--
      wt(14) = 8.0 /  27.0 ! centre
      wt(15) = 2.0 /  27.0 ! z++
      wt(16) = 1.0 /  54.0 ! r2
      wt(17) = 2.0 /  27.0 ! y++
      wt(18) = 1.0 /  54.0 ! r2

      wt(19) = 1.0 / 216.0 ! r3
      wt(20) = 1.0 /  54.0 ! r2
      wt(21) = 1.0 / 216.0 ! r3
      wt(22) = 1.0 /  54.0 ! r2
      wt(23) = 2.0 /  27.0 ! x--
      wt(24) = 1.0 /  54.0 ! r2
      wt(25) = 1.0 / 216.0 ! r3
      wt(26) = 1.0 /  54.0 ! r2
      wt(27) = 1.0 / 216.0 ! r3

      col1 = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9  /)
      col2 = (/ 10,11,12,13,14,15,16,17,18 /)
      col3 = (/ 19,20,21,22,23,24,25,26,27 /)

      obs = 1     ! int array (1 = live)

      fin = 0.0   ! real array
      fout= 0.0   ! real array
      feq = 0.0   ! real array
      rho = 0.0   ! real array
      vel = 0.0   ! real array


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
        dz = zz-cz
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

        do k = 1,nk ! top/bottom walls
          obs(i, 1,k) = 0 ! ymin
          obs(i,nj,k) = 0 ! ymax
        enddo

        do j = 1,nj ! front/back walls
          obs(i,j, 1) = 0 ! zmin
          obs(i,j,nk) = 0 ! zmax
        enddo

      enddo
      endif


C-----read obstacle.dat file from ufocfd
      if (.true.) then
        write(6,*) "read obstacle.dat"
        open(20,file="obstacle.dat")

        read(20,*) ! im jm km (for user reference only)
        read(20,*) obsmax

        do n = 1,obsmax
          read(20,*) i,j,k,lcd
CCC       write(6,*) ">>>",i,j,k,lcd
          obs(i,j,k) = lcd ! (0 or -1)
        enddo

        close(20)
        close(30)

        write(6,*)"obsmax = ", obsmax
        write(6,*)"done"
      endif



C=====initial velocity field
      write(6,*)"initial velocity"
      do i = 1,ni
      do j = 1,nj
      do k = 1,nk

        if (obs(i,j,k).eq.1) then ! fluid (lcd field)
          vel(i,j,k,1) = vin(1)
          vel(i,j,k,2) = vin(2)
          vel(i,j,k,3) = vin(3)
        endif

        if (obs(i,j,k).le.0) then ! cutcell or solid (lcd field)
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


C=====saving vtk files
      do sa = 1, nsav

!$acc data copy(obs,vel,fin,fout,feq,rho)
!$acc data copy(vec,wt,col1,col2,col3,vin)

C=====main iterations
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
        vel(i,j,k,1) = vel(i,j,k,1) + vec(n,1)*fin(i,j,k,n)
        vel(i,j,k,2) = vel(i,j,k,2) + vec(n,2)*fin(i,j,k,n)
        vel(i,j,k,3) = vel(i,j,k,3) + vec(n,3)*fin(i,j,k,n)
      enddo
      
      vel(i,j,k,1) = vel(i,j,k,1) / rho(i,j,k)
      vel(i,j,k,2) = vel(i,j,k,2) / rho(i,j,k)
      vel(i,j,k,3) = vel(i,j,k,3) / rho(i,j,k)


C     STEP3 - left wall ! (and upper/lower wall) inflow condition
      if (i.eq.1) then  ! .or.(j.eq.1).or.(j.eq.nj)) then
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
     &                       + fin(1,j,k,col3(nrow+1 -n))
     &                       - feq(1,j,k,col3(nrow+1 -n))
        enddo
      endif


C     STEP6 - collision step (original LBM approach)
      do n = 1,nvec
CCC     fout(i,j,k,n) = fin(i,j,k,n) - omega*(fin(i,j,k,n) - feq(i,j,k,n))
        fout(i,j,k,n) = omega*feq(i,j,k,n) + fin(i,j,k,n)*(1.0-omega) 
      enddo


C     STEP7 - bounce-back obstacle wall condition
      if (obs(i,j,k).le.0) then ! cutcell or dead
        do n = 1,nvec
          fout(i,j,k,n) = fin(i,j,k,nvec+1-n)
        enddo
      endif  !  (obs.eq.0) cutcell



      enddo  !  i-loop
      enddo  !  j-loop
      enddo  !  k-loop


!$acc kernels loop independent
      do k = 1,nk
!$acc loop independent
      do j = 1,nj
!$acc loop independent
      do i = 1,ni

C       STEP8 - streaming step
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
     &  vel(int(0.75*ni),int(0.5*nj),int(0.5*nk),1)
      endif

      enddo ! nits

!$acc end data
!$acc end data

      call write_vxmax()
C     call write_ascii_vtk()
C     call write_binary_vtk()
      call write_binary_slice_vtk()

      enddo ! nsav
C=====end main iteration loop


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

!$acc data copy(vel)

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

!$acc end data

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

!$acc data copy(obs,vel,rho)

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
C-----vel
      write(20)'VECTORS vel float'//lf
      write(20)(((vel(i,j,k,1),vel(i,j,k,2),vel(i,j,k,3),
     &             i=1,ni),j=1,nj),k=1,nk)
      close(20)

!$acc end data

      end



C-----write binary (slice) VTK file
      subroutine write_binary_slice_vtk()

      use allocated
      implicit none

      integer i,j,k,k1,k2
      integer ni,nj,nk,sa
      common  ni,nj,nk,sa

      CHARACTER(15) outfile
      CHARACTER(3)  string
      CHARACTER(LEN=1)  :: lf
      CHARACTER(LEN=10) :: str1,str2,str3,str4
      lf = char(10)

C-----slice location
      k1 = int(nk/2)+1
      k2 = k1+1

      write(6,*)"write binary slice VTK"

C-----output filename
      write(unit=string, fmt='(I3.3)') sa
      outfile = 'lbmslice'//string//'.vtk'
      write(6,*)"outfile = ",outfile

!$acc data copy(obs,rho,vel)

      OPEN(unit=20, file=outfile, form='unformatted',
     &    access='stream',status='replace',convert="big_endian")
      write(20)'# vtk DataFile Version 3.0'//lf
      write(20)'vtk output'//lf
      write(20)'BINARY'//lf
      write(20)'DATASET RECTILINEAR_GRID'//lf
      write(str1(1:10),'(i10)') ni+1
      write(str2(1:10),'(i10)') nj+1
      write(str3(1:10),'(i10)') 2       !(slice in k)
      write(str4(1:10),'(i10)') ni*nj*1 !(slice in k)
      write(20)'DIMENSIONS '//str1//str2//str3//lf
      write(20)'X_COORDINATES '//str1//' float'//lf
      write(20)(real(i-1),i=1,ni+1)
      write(20)'Y_COORDINATES '//str2//' float'//lf
      write(20)(real(j-1),j=1,nj+1)
      write(20)'Z_COORDINATES '//str3//' float'//lf
      write(20)(real(k-1),k=k1,k2)
      write(20)'CELL_DATA '//str4//lf
C-----obs
      write(20)'SCALARS obs int'//lf
      write(20)'LOOKUP_TABLE default'//lf
      write(20)(((obs(i,j,k),i=1,ni),j=1,nj),k=k1,k1)
C-----rho
      write(20)'SCALARS rho float'//lf
      write(20)'LOOKUP_TABLE default'//lf
      write(20)(((rho(i,j,k),i=1,ni),j=1,nj),k=k1,k1)
C-----vel
      write(20)'VECTORS vel float'//lf
      write(20)(((vel(i,j,k,1),vel(i,j,k,2),
     &          vel(i,j,k,3),i=1,ni),j=1,nj),k=k1,k1)
      close(20)

!$acc end data

      end

