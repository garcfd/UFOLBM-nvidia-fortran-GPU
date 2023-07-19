
      program main
      implicit none
      integer, parameter :: ndim=3, nvec=27, nnrm = 66

      integer x,y,z,cs,npts,pt,n,v,v1,v2,v2save,nfull
      integer ref(nnrm,nvec),fill(nvec),full(nnrm)

      real xx,yy,zz,xx2,yy2,zz2,r2,r3,dx,dy,dz
      real nrm1,nrm2,nrm3,unrm,vnrm,wnrm,dotp
      real vmin,vmin1,vmin2,vmin3,vmax,vmax1,vmax2,vmax3
      real n1,n2,n3,ref1,ref2,ref3,dist,distmin
      real nrm(nnrm,ndim),vec(nvec,ndim)

      if (.false.) then
C-----calc sube to sphere points
      cs = 2 ! cube size
      npts = (2*cs + 1)**3 - (2*cs - 1)**3
      write(6,*) "npts = ",npts 
      stop

C-----open ply for for paraview
      open(unit=10, file="normals218.ply")      
      write(10,10) "ply"
      write(10,10) "format ascii 1.0"
      write(10,20) "element vertex ",nnrm
      write(10,10) "property float x"
      write(10,10) "property float y"
      write(10,10) "property float z"
      write(10,10) "property uchar red"
      write(10,10) "property uchar green"
      write(10,10) "property uchar blue"
      write(10,10) "element face 0"
      write(10,10) "property list uchar int vertex_index"
      write(10,10) "end_header"

      pt = 0

C-----calc the cube to sphere mapping
      do x = -cs,cs
      do y = -cs,cs
      do z = -cs,cs

        xx = real(x)/real(cs)
        yy = real(y)/real(cs)
        zz = real(z)/real(cs)
       
        xx2 = xx*xx
        yy2 = yy*yy
        zz2 = zz*zz

        n1 = xx*sqrt(1.0 - yy2/2.0 - zz2/2.0 + yy2*zz2/3.0)
        n2 = yy*sqrt(1.0 - xx2/2.0 - zz2/2.0 + xx2*zz2/3.0) 
        n3 = zz*sqrt(1.0 - xx2/2.0 - yy2/2.0 + xx2*yy2/3.0)

        if ((x.eq. cs).or.(y.eq. cs).or.(z.eq. cs).or.
     &      (x.eq.-cs).or.(y.eq.-cs).or.(z.eq.-cs)) then        

          pt = pt + 1

          nrm(pt,:) = (/ n1,n2,n3 /)   
CCC       write(6, *) pt,nrm(pt,:)

C---------write data to ply file
          if ((x.eq.0).or.(y.eq.0).or.(z.eq.0)) then
            write(10,*) n1,n2,n3," 200 0 0 " ! red
          elseif ((abs(x).eq.cs).and.(abs(y).eq.cs)
     &       .and.(abs(z).eq.cs))then
            write(10,*) n1,n2,n3," 0 200 0 " ! green
          else
            write(10,*) n1,n2,n3," 0 0 200 " ! blue
          endif

        endif

      enddo
      enddo
      enddo

      close(10)

   10 format(A)
   20 format(A,I6)

      endif ! true/false



C-----initialise the lbm vectors
      r2 = 1.0/sqrt(2.0) ! 0.707
      r3 = 1.0/sqrt(3.0) ! 0.577

      vec(1, :) = (/ r3,-r3,-r3 /) ! r3
      vec(2, :) = (/ r2,-r2, 0. /) ! r2   
      vec(3, :) = (/ r3,-r3, r3 /) ! r3
      vec(4, :) = (/ r2, 0.,-r2 /) ! r2
      vec(5, :) = (/ 1., 0., 0. /) ! x++
      vec(6, :) = (/ r2, 0., r2 /) ! r2
      vec(7, :) = (/ r3, r3,-r3 /) ! r3
      vec(8, :) = (/ r2, r2, 0. /) ! r2
      vec(9, :) = (/ r3, r3, r3 /) ! r3

      vec(10,:) = (/ 0.,-r2,-r2 /) ! r2
      vec(11,:) = (/ 0.,-1., 0. /) ! y--
      vec(12,:) = (/ 0.,-r2, r2 /) ! r2
      vec(13,:) = (/ 0., 0.,-1. /) ! z--
      vec(14,:) = (/ 0., 0., 0. /) ! r2
      vec(15,:) = (/ 0., 0., 1. /) ! z++
      vec(16,:) = (/ 0., r2,-r2 /) ! r2
      vec(17,:) = (/ 0., 1., 0. /) ! y++
      vec(18,:) = (/ 0., r2, r2 /) ! r2

      vec(19,:) = (/-r3,-r3,-r3 /) ! r3
      vec(20,:) = (/-r2,-r2, 0. /) ! r2
      vec(21,:) = (/-r3,-r3, r3 /) ! r3
      vec(22,:) = (/-r2, 0.,-r2 /) ! r2
      vec(23,:) = (/-1., 0., 0. /) ! x--
      vec(24,:) = (/-r2, 0., r2 /) ! r2
      vec(25,:) = (/-r3, r3,-r3 /) ! r3
      vec(26,:) = (/-r2, r2, 0. /) ! r2
      vec(27,:) = (/-r3, r3, r3 /) ! r3

C-----write lbm vectors
      write(6,*)
      do v = 1,nvec
CCC     write(6,*) v,vec(v,:)
      enddo




C-----read normals from file
      if (.true.) then
      open(20, file="normals66.ply")
      do n = 1,9
        read(20,*) ! headers
      enddo
      do n = 1,nnrm
        read(20,*) nrm(n,:)
      enddo
      close(20)
      endif




C-----calc the reflections
      ref = 0 

      do n = 1,nnrm

        nrm1 = nrm(n,1)
        nrm2 = nrm(n,2)
        nrm3 = nrm(n,3)

C-------loop through vectors
        do v1 = 1,nvec ! D3Q27

          dotp = nrm1*vec(v1,1) + nrm2*vec(v1,2) + nrm3*vec(v1,3)
          unrm = nrm1*dotp ! surf normal component
          vnrm = nrm2*dotp ! surf normal component
          wnrm = nrm3*dotp ! surf normal component

          vmin1 = vec(v1,1) - 2.0*unrm
          vmin2 = vec(v1,2) - 2.0*vnrm
          vmin3 = vec(v1,3) - 2.0*wnrm
          vmin  = sqrt(vmin1*vmin1 + vmin2*vmin2 + vmin3*vmin3)

          vmax1 = vec(v1,1) + 2.0*unrm
          vmax2 = vec(v1,2) + 2.0*vnrm
          vmax3 = vec(v1,3) + 2.0*wnrm
          vmax  = sqrt(vmax1*vmax1 + vmax2*vmax2 + vmax3*vmax3)

          if (vmin.lt.vmax) then
            ref1 = vmin1
            ref2 = vmin2
            ref3 = vmin3
          else
            ref1 = vmax1
            ref2 = vmax2
            ref3 = vmax3
          endif

          distmin = 10.0

C---------loop through lbm vectors to find nearest
          do v2 = 1,nvec ! 27

            dx = ref1 - vec(v2,1) 
            dy = ref2 - vec(v2,2) 
            dz = ref3 - vec(v2,3) 
            dist = sqrt(dx*dx + dy*dy + dz*dz)

            if (dist.lt.distmin) then
              v2save = v2
              distmin = dist 
            endif

          enddo ! nvec

          ref(n,v1) = v2save

        enddo ! nvec

      enddo ! nnrm

C-----write the reflections
      do n = 1,nnrm
        write(6,'(I3,A,27I3)')n," > ",ref(n,:)
      enddo

C-----calc fullness
      nfull = 0
      
      do n = 1,nnrm

        fill = 0
        do v = 1,nvec
          fill( ref(n,v) ) = 1
        enddo

        full(n) = 0
        do v = 1,nvec
          full(n) = full(n) + fill(v)
        enddo

        if (full(n).eq.27) then
          write(6,'(I6,A,27I3,A,I3)')n," > ",fill(:)," >>> ",full(n)
          nfull = nfull+1
        endif

        if (full(n).lt.27) then
          write(6,'(I6,A,27I3,A,I3)')n," > ",fill(:),"     ",full(n)
        endif

      enddo

      write(6,*)"nnrm  = ",nnrm
      write(6,*)"nfull = ",nfull

C-----write the full normals to ply file
C-----and the reflections 

      open(unit=10, file="normals66.ply")      
      open(unit=20, file="reflect66.dat")      

      write(10,10) "ply"
      write(10,10) "format ascii 1.0"
      write(10,20) "element vertex ", nfull
      write(10,10) "property float x"
      write(10,10) "property float y"
      write(10,10) "property float z"
      write(10,10) "element face 0"
      write(10,10) "property list uchar int vertex_index"
      write(10,10) "end_header"

      nfull = 0
      do n = 1,nnrm
        if (full(n).eq.27) then
          write(10,*) nrm(n,:)
          write(20,'(27I3)') ref(n,:)
          nfull = nfull+1
        endif
      enddo
      write(6,*)"nfull check ",nfull
      close(10)
      close(20)

      end







