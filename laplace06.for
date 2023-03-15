C
C     for OpenMP
C     export OMP_NUM_THREADS=8 (takes 10sec)
C     pgf77 -fopenmp -o laplace laplace03.for
C
C     export ACC_NUM_CORES=8
C     export NVHPC_CUDA_HOME="/usr/local/cuda-12.1"
C     export NV_ACC_TIME=1
C
C pgf77 -static-nvidia -acc=host      -o laplace_basic laplace05.for
C pgf77 -static-nvidia -acc=multicore -o laplace_multi laplace05.for
C pgf77 -static-nvidia -acc -gpu=cc75 -o laplace_tesla laplace05.for

      program main

      parameter (im=200, jm=200, km=200)

      integer imm,jmm,kmm,it,nits
      real phi(im,jm,km),onesixth

      INTEGER i,j,k,iter
      CHARACTER(LEN=1)  :: lf
      CHARACTER(LEN=10) :: str1,str2,str3,str4
      CHARACTER(20) outfile
      CHARACTER(6)  string

       
      onesixth = 1.0/6.0 
       
      imm = im -1
      jmm = jm -1
      kmm = km -1

      nits = 10000
      phi  = 1.0
      phi(2:imm,2:jmm,2:kmm) = 0.0


!$acc data copy(phi(:,:,:))

      do it = 1, nits

!$acc parallel loop collapse(3) present(phi)
        do k = 2, kmm
        do j = 2, jmm
        do i = 2, imm
          phi(i,j,k) = ( phi(i-1,j,k) + phi(i+1,j,k)        
     &                 + phi(i,j-1,k) + phi(i,j+1,k)        
     &                 + phi(i,j,k-1) + phi(i,j,k+1) )*onesixth
        enddo
        enddo
        enddo

        if (mod(it,100).eq.0) then
          write(6,*) "it = ",it,"phi = ", phi(im/2,jm/2,km/2)
        endif

      enddo

!$acc end data

      write(6,*) "it = ",it,"phi = ", phi(im/2,jm/2,km/2)



      if (.false.) then
C-----write ascii VTK file
      WRITE(6,*) 'write ascii vtk'
      OPEN(UNIT=20, FILE='laplace.vtk')
      WRITE(20,10)'# vtk DataFile Version 2.0'
      WRITE(20,10)'sample rectilinear grid'
      WRITE(20,10)'ASCII'
      WRITE(20,10)'DATASET RECTILINEAR_GRID'
      WRITE(20,20)'DIMENSIONS ', im+1, jm+1, km+1
      WRITE(20,30)'X_COORDINATES ', im+1, ' float'
      WRITE(20,*)( real(i-1),   i=1,im+1)
      WRITE(20,30)'Y_COORDINATES ', jm+1, ' float'
      WRITE(20,*)( real(j-1),   j=1,jm+1)
      WRITE(20,30)'Z_COORDINATES ', km+1, ' float'
      WRITE(20,*)( real(k-1),   k=1,km+1)
      WRITE(20,40)'CELL_DATA ', im*jm*km
      WRITE(20,10)'SCALARS phi float'
      WRITE(20,10)'LOOKUP_TABLE default'
      WRITE(20,*)(((phi(i,j,k),i=1,im),j=1,jm),k=1,km)
      CLOSE(20)
   10 FORMAT(A)
   20 FORMAT(A,3I5)
   30 FORMAT(A,I5,A)
   40 FORMAT(A,I9)
      endif

      
      if (.true.) then
C-----write binary VTK file
      write(6,*) 'write binary vtk'
      lf = char(10)
      OPEN(unit=20,file='laplace.vtk',form='unformatted',
     &  access='stream',status='replace',convert="big_endian")
      write(20)'# vtk DataFile Version 3.0'//lf
      write(20)'vtk output'//lf
      write(20)'BINARY'//lf
      write(20)'DATASET RECTILINEAR_GRID'//lf
      write(str1(1:10),'(i10)') im+1
      write(str2(1:10),'(i10)') jm+1
      write(str3(1:10),'(i10)') km+1
      write(str4(1:10),'(i10)') im*jm*km
      write(20)'DIMENSIONS '//str1//str2//str3//lf
      write(20)'X_COORDINATES '//str1//' float'//lf
      write(20)(real(i-1),i=1,im+1)
      write(20)'Y_COORDINATES '//str2//' float'//lf
      write(20)(real(j-1),j=1,jm+1)
      write(20)'Z_COORDINATES '//str3//' float'//lf
      write(20)(real(k-1),k=1,km+1)
      write(20)'CELL_DATA '//str4//lf
      write(20)'SCALARS phi float'//lf
      write(20)'LOOKUP_TABLE default'//lf
      write(20)(((phi(i,j,k),i=1,im),j=1,jm),k=1,km)
      close(20)
      endif

      END


