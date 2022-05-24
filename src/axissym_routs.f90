!
!
!
!
      subroutine get_axissym_fun_mem(nch2d,tchse,k,fcurve,np,pars, &
        rmax,npatches)
!
!  This subroutine estimates the number of patches required to discretize
!  an axissymmetric geometry where the generating curve is provided as
!  a function handle.
!  The axis of rotation is assumed to be the z
!  axis.
!
!  The mesh will with high probability be non-conforming, and 
!  each patch will have a specified max patch size
!
!  Input arguments:
!    - nch2d: integer
!        number of chunks describing the generating curve
!    - tchse: real *8 (nch2d+1)
!        starting and ending location of in parameter space
!        for each 2d chunk
!    - k: integer
!        number of points per chunk
!    - fcurve: function handle
!        function handle for corresponding to generating curve.
!        Should have calling sequence
!        fcurve(t,np,pars,r,z,drdt,dzdt,d2rdt2,d2zdt2)
!    - np: integer
!        number of parameters in fcurve
!    - pars: real *8 (np)
!        parameters of fcurve
!    - rmax: real *8
!        maximum panel dimension
!
!  Output arguments:
!    - npatches: integer
!        number of patches in the discretization
!
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: nch2d,np,k
      real *8, intent(in) :: tchse(nch2d+1),pars(np)
      real *8, intent(in) :: rmax
      integer, intent(out) :: npatches

      real *8, allocatable :: ts(:),ws(:),umat(:,:),vmat(:,:)

      external fcurve

      real *8 xs(2),rs(2)

      done = 1.0d0
      pi = atan(done)*4

      allocate(ts(k),umat(k,k),vmat(k,k),ws(k))
      
      itype = 1
      call legeexps(itype,k,ts,umat,vmat,ws)

      npatches = 0



      npatches = 0
      do ich=1,nch2d


        radmax = 0
        rlen = 0
        rs(1:2) = 0

        h = (tchse(ich+1)-tchse(ich))/2
        do j=1,k
          tt = tchse(ich) + (ts(j)+1.0d0)/2*(tchse(ich+1)-tchse(ich))
          call fcurve(tt,np,pars,r,z,drdt,dzdt,d2rdt2,d2zdt2)
          if(r.ge.radmax) r = radmax
          dsdt = sqrt(drdt**2 + dzdt**2)*h

          rlen = rlen + dsdt*ws(j)
        enddo

        call fcurve(tchse(ich),np,pars,rs(1),z, &
            drdt,dzdt,d2rdt2,d2zdt2)
        call fcurve(tchse(ich+1),np,pars,rs(2),z, &
            drdt,dzdt,d2rdt2,d2zdt2)

        if(rs(1).ge.radmax) radmax = rs(1)
        if(rs(2).ge.radmax) radmax = rs(2)

        if(rlen.ge.2*pi*radmax) then
          nt = max(ceiling(2*pi*radmax/rmax),2)
          ns = ceiling(nt*rlen/2/pi/radmax)
        else
          ns = max(ceiling(rlen/rmax),2)
          nt = ceiling(ns*2*pi*radmax/rlen)
        endif


        print *, ich,ns,nt,rlen,radmax,rmax

        npatches = npatches + 2*ns*nt
      enddo


      end subroutine get_axissym_fun_mem
!
!
!
!
!
!
!
!
!
!
!
      subroutine get_axissym_fun_geom(nch2d,tchse,k,fcurve,np,pars, &
        rmax,iort,norder,npatches,npts,norders,ixyzs, &
        iptype,srcvals,srccoefs)
!
!  This subroutine discretizes
!  an axissymmetric geometry where the generating curve is provided as
!  a function handle.
!  The axis of rotation is assumed to be the z
!  axis.
!
!  The mesh will with high probability be non-conforming, and 
!  each patch will have a specified max patch size
!
!  Input arguments:
!    - nch2d: integer
!        number of chunks describing the generating curve
!    - tchse: real *8 (nch2d+1)
!        starting and ending location of in parameter space
!        for each 2d chunk
!    - k: integer
!        number of points per chunk
!    - fcurve: function handle
!        function handle for corresponding to generating curve.
!        Should have calling sequence
!        fcurve(t,np,pars,r,z,drdt,dzdt,d2rdt2,d2zdt2)
!    - np: integer
!        number of parameters in fcurve
!    - pars: real *8 (np)
!        parameters of fcurve
!    - rmax: real *8
!        maximum panel dimension
!    - iort: integer
!        orientation of discretization
!    - npatches: integer
!        number of patches in discretization
!    - npts: integer
!        number of points in the discretization =
!          npatches*(norder+1)*(norder+2)/2
!
!  Output arguments:
!    - norders: integer(npatches)
!        discretization order of patches
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: nch2d,k
      integer, intent(in), target :: np
      integer, intent(in) :: npatches,npts,iort
      real *8, intent(in) :: tchse(nch2d+1)
      real *8, intent(in), target :: pars(np)
      real *8, intent(in) :: rmax
      integer, intent(out) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(out) :: iptype(npatches)
      real *8, intent(out) :: srcvals(12,npts),srccoefs(9,npts)
      real *8, allocatable, target :: triaskel(:,:,:)
      real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)

      real *8, allocatable :: ts(:),ws(:),umat(:,:),vmat(:,:)
      real *8, pointer :: ptr1,ptr2,ptr3,ptr4
      integer, pointer :: iptr1,iptr2,iptr3,iptr4

      external fcurve, xtri_axissym_fun_chunk

      real *8 xs(2),rs(2)

      done = 1.0d0
      pi = atan(done)*4
      allocate(ts(k),umat(k,k),vmat(k,k),ws(k))
      
      itype = 1
      call legeexps(itype,k,ts,umat,vmat,ws)



      allocate(triaskel(3,3,npatches))

      npatches0 = 0
      

      itristart = 1
      nover = 0
      do ich=1,nch2d


        radmax = 0
        rlen = 0
        rs(1:2) = 0

        h = (tchse(ich+1)-tchse(ich))/2
        do j=1,k
          tt = tchse(ich) + (ts(j)+1.0d0)/2*(tchse(ich+1)-tchse(ich))
          call fcurve(tt,np,pars,r,z,drdt,dzdt, &
            d2rdt2,d2zdt2)
          if(r.ge.radmax) r = radmax
          dsdt = sqrt(drdt**2 + dzdt**2)*h
          rlen = rlen + dsdt*ws(j)
        enddo

        call fcurve(tchse(ich),np,pars,rs(1),z, &
            drdt,dzdt,d2rdt2,d2zdt2)
        call fcurve(tchse(ich+1),np,pars,rs(2),z, &
            drdt,dzdt,d2rdt2,d2zdt2)
        if(rs(1).ge.radmax) radmax = rs(1)
        if(rs(2).ge.radmax) radmax = rs(2)

        if(rlen.ge.2*pi*radmax) then
          nt = max(ceiling(2*pi*radmax/rmax),2)
          ns = ceiling(nt*rlen/2/pi/radmax)
        else
          ns = max(ceiling(rlen/rmax),2)
          nt = ceiling(ns*2*pi*radmax/rlen)
        endif



        umin = tchse(ich)
        umax = tchse(ich+1)
        
        if(iort.eq.1) then
          vmin = 2*pi
          vmax = 0
        else
          vmin = 0
          vmax = 2*pi
        endif
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ns,nt,nover, &
          npatches,npatches0,triaskel(1,1,itristart))
        itristart = itristart + 2*ns*nt
      enddo
      ptr1 => triaskel(1,1,1)
      iptr2 => np
      ptr3 => pars(1)
      


      npols = (norder+1)*(norder+2)/2
      allocate(uvs(2,npols),wts(npols),umatr(npols,npols), &
        vmatr(npols,npols))
      call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)
      
      call getgeominfo(npatches,xtri_axissym_fun_chunk,ptr1,iptr2,ptr3, &
        fcurve,npols,uvs,umatr,srcvals,srccoefs)

      do i=1,npatches
        ixyzs(i) = (i-1)*npols + 1
        iptype(i) = 1
        norders(i) = norder
      enddo

      ixyzs(npatches+1) = npts+1

      


      end subroutine get_axissym_fun_geom
!
!
!
!
!
!
!
!
!
!
      subroutine get_axissym_mem(nch2d,npts2d,ixys2d,iptype2d,srcvals2d, &
        srccoefs2d,rmax,npatches)
!
!  This subroutine estimates the number of patches required to discretize
!  an axissymmetric geometry where the generating curve is provided in the 
!  fmm2dbie chunk format. The axis of rotation is assumed to be the z
!  axis.
!
!  The mesh will with high probability be non-conforming, and 
!  each patch will have a specified max patch size
!
!  Input arguments:
!    - nch2d: integer
!        number of chunks describing the generating curve
!    - npts2d: integer
!        number
!    - ixys2d: integer(nch2d+1)
!        starting location of points on chunk i
!    - iptype2d: integer(nch2d)
!        type of chunk
!        iptype = 1, gauss-legendre chunk with reference domain [-1,1]
!    - srcvals2d: real *8 (8,npts2d)
!        rz, drz/dt,d2rz/dt2, normals at all discretization nodes
!    - srccoefs2d: real *8 (6,npts2d)
!        Basis coefficient expansion of rz,drz/dt,d2rz/dt2
!    - rmax: real *8
!        maximum panel dimension
!
!  Output arguments:
!    - npatches: integer
!        number of patches in the discretization
!
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: nch2d,npts2d,ixys2d(nch2d+1)
      integer, intent(in) :: iptype2d(nch2d)
      real *8, intent(in) :: srcvals2d(8,npts2d),srccoefs2d(6,npts2d)
      real *8, intent(in) :: rmax
      integer, intent(out) :: npatches

      real *8, allocatable :: ts(:),ws(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: amatrint(:,:),work(:)

      real *8 xs(2),rs(2)

      done = 1.0d0
      pi = atan(done)*4
      
      itype = 1
      xs(1) = -1.0d0
      xs(2) = 1.0d0




      npatches = 0
!
!  Find max k for 2d patch discretization
!

      kmax = ixys2d(2)-ixys2d(1)
      
      do ich=1,nch2d
        k = ixys2d(ich+1)-ixys2d(ich)
        if(k.ge.kmax) kmax = k
      enddo

      allocate(ts(kmax),umat(kmax,kmax),vmat(kmax,kmax),ws(kmax))
      allocate(amatrint(2,kmax))

      lw = 4*kmax*kmax + kmax + 200
      allocate(work(lw))

      kprev = ixys2d(2)-ixys2d(1)
      k = ixys2d(2)-ixys2d(1)

      call legeexps(itype,k,ts,umat,vmat,ws)
      call lematrin(k,2,xs,amatrint,ts,work)

      npatches = 0
      do ich=1,nch2d
        k = ixys2d(ich+1)-ixys2d(ich)

        if(k.ne.kprev) then
          call legeexps(itype,k,ts,umat,vmat,ws)
          call lematrin(k,2,xs,amatrint,ts,work)
        endif

        radmax = 0
        rlen = 0
        istart = ixys2d(ich)
        rs(1:2) = 0
        do j=1,k
          rs(1:2) = rs(1:2) + amatrint(1:2,j)*srcvals2d(1,istart+j-1)
          if(srcvals2d(1,istart+j-1).ge.radmax) radmax = srcvals2d(1,istart+j-1)
          dsdt = sqrt(srcvals2d(3,istart+j-1)**2 + srcvals2d(4,istart+j-1)**2)
          rlen = rlen + dsdt*ws(j)
        enddo

        if(rs(1).ge.radmax) radmax = rs(1)
        if(rs(2).ge.radmax) radmax = rs(2)

        if(rlen.ge.2*pi*radmax) then
          nt = max(ceiling(2*pi*radmax/rmax),2)
          ns = ceiling(nt*rlen/2/pi/radmax)
        else
          ns = max(ceiling(rlen/rmax),2)
          nt = ceiling(ns*2*pi*radmax/rlen)
        endif



        npatches = npatches + 2*ns*nt
        kprev = k
      enddo


      end subroutine get_axissym_mem
!
!
!
!
!
      subroutine get_axissym_geom(nch2d,npts2d,ixys2d,iptype2d,srcvals2d, &
        srccoefs2d,rmax,iort,norder,npatches,npts,norders,ixyzs,iptype,srcvals, &
        srccoefs)
!
!  This subroutine discretizes an axissymmetric geometry where the 
!  generating curve is provided in the 
!  fmm2dbie chunk format. The axis of rotation is assumed to be the z
!  axis.
!
!  The mesh will with high probability be non-conforming, and 
!  each patch will have a specified max patch size
!
!  Input arguments:
!    - nch2d: integer
!        number of chunks describing the generating curve
!    - npts2d: integer
!        number
!    - ixys2d: integer(nch2d+1)
!        starting location of points on chunk i
!    - iptype2d: integer(nch2d)
!        type of chunk
!        iptype = 1, gauss-legendre chunk with reference domain [-1,1]
!    - srcvals2d: real *8 (8,npts2d)
!        rz, drz/dt,d2rz/dt2, normals at all discretization nodes
!    - srccoefs2d: real *8 (6,npts2d)
!        Basis coefficient expansion of rz,drz/dt,d2rz/dt2
!    - rmax: real *8
!        maximum panel dimension
!    - iort: integer
!        orientation
!    - norder: integer
!        order of discretization
!    - npatches: integer
!        number of patches in the discretization
!    - npts: integer
!        number of points in the discretization = 
!          npatches*(norder+1)*(norder+2)/2
!
!  Output arguments:
!    - norders: integer(npatches)
!        discretization order of patches
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!
!
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: nch2d,npts2d
      integer, intent(in), target :: ixys2d(nch2d+1)
      integer, intent(in) :: iptype2d(nch2d)
      real *8, intent(in) :: srcvals2d(8,npts2d)
      real *8, intent(in), target :: srccoefs2d(6,npts2d)
      real *8, intent(in) :: rmax
      integer, intent(in) :: npatches,norder,npts

      integer, intent(out) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(out) :: iptype(npatches)
      real *8, intent(out) :: srcvals(12,npts),srccoefs(9,npts)

      real *8, pointer :: ptr1,ptr2,ptr3,ptr4
      integer, pointer :: iptr1,iptr2,iptr3,iptr4

      real *8, allocatable :: ts(:),ws(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: amatrint(:,:),work(:)

      real *8 xs(2),rs(2)
      real *8, allocatable, target :: triaskel(:,:,:)
      integer, allocatable, target :: ichuse(:)
      integer, allocatable, target :: kuse(:)
      real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)
      external xtri_axissym_chunk

      done = 1.0d0
      pi = atan(done)*4
      
      itype = 1
      xs(1) = -1.0d0
      xs(2) = 1.0d0





!
!  Find max k for 2d patch discretization
!

      kmax = ixys2d(2)-ixys2d(1)
      
      do ich=1,nch2d
        k = ixys2d(ich+1)-ixys2d(ich)
        if(k.ge.kmax) kmax = k
      enddo

      allocate(ts(kmax),umat(kmax,kmax),vmat(kmax,kmax),ws(kmax))
      allocate(amatrint(2,kmax))

      lw = 4*kmax*kmax + kmax + 200
      allocate(work(lw))

      allocate(triaskel(3,3,npatches),ichuse(npatches))

      kprev = ixys2d(2)-ixys2d(1)
      k = ixys2d(2)-ixys2d(1)

      call legeexps(itype,k,ts,umat,vmat,ws)
      call lematrin(k,2,xs,amatrint,ts,work)

      npatches0 = 0
      nover = 0
      itristart = 1


      do ich=1,nch2d
        k = ixys2d(ich+1)-ixys2d(ich)

        if(k.ne.kprev) then
          call legeexps(itype,k,ts,umat,vmat,ws)
          call lematrin(k,2,xs,amatrint,ts,work)
        endif

        radmax = 0
        rlen = 0
        rs(1:2) = 0
        istart = ixys2d(ich)
        umin = -1.0d0
        umax = 1.0d0
        if(iort.eq.1) then 
          vmin = 2*pi
          vmax = 0
        else
          vmin = 0
          vmax = 2*pi
        endif

        nover = 0

        do j=1,k
          rs(1:2) = rs(1:2) + amatrint(1:2,j)*srcvals2d(1,istart+j-1)
          if(srcvals2d(1,istart+j-1).ge.radmax) radmax = srcvals2d(1,istart+j-1)
          dsdt = sqrt(srcvals2d(3,istart+j-1)**2 + srcvals2d(4,istart+j-1)**2)
          rlen = rlen + dsdt*ws(j)
        enddo
        radmax0 = radmax

        if(rs(1).ge.radmax) radmax = rs(1)
        if(rs(2).ge.radmax) radmax = rs(2)

        if(rlen.ge.2*pi*radmax) then
          nt = max(ceiling(2*pi*radmax/rmax),2)
          ns = ceiling(nt*rlen/2/pi/radmax)
        else
          ns = max(ceiling(rlen/rmax),2)
          nt = ceiling(ns*2*pi*radmax/rlen)
        endif


        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ns,nt,nover, &
          npatches,npatches0,triaskel(1,1,itristart))
        ichuse(itristart:(itristart+2*ns*nt-1)) = ich
        itristart = itristart + 2*ns*nt
        kprev = k
      enddo


      ptr1 => triaskel(1,1,1)
      ptr2 => srccoefs2d(1,1)
      iptr3 => ichuse(1)
      iptr4 => ixys2d(1)

      npols = (norder+1)*(norder+2)/2
      allocate(uvs(2,npols),wts(npols),umatr(npols,npols), &
        vmatr(npols,npols))
      call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)
      
      call getgeominfo(npatches,xtri_axissym_chunk,ptr1,ptr2,iptr3, &
        iptr4,npols,uvs,umatr,srcvals,srccoefs)

      do i=1,npatches
        ixyzs(i) = (i-1)*npols + 1
        iptype(i) = 1
        norders(i) = norder
      enddo

      ixyzs(npatches+1) = npts+1



      end subroutine get_axissym_geom
!
!
!
!
!
      subroutine xtri_axissym_chunk(itri, u, v, xyz, dxyzduv, & 
          triainfo,rzcoefs,ichuse,ixys2d)
      implicit real *8 (a-h,o-z)
      real *8 :: xyz(3), dxyzduv(3,2), triainfo(3,3,*)
      real *8 :: rzcoefs(6,*),pols(100)
      integer :: ichuse(*),ixys2d(*)


      !
      ! project the triangle itri in triainfo onto the sphere
      !
      !    Input:
      ! itri - triangle number to map
      ! u,v - local uv coordinates on triangle itri
      ! triainfo - flat skeleton triangle info
      ! rzcoefs - legendre coefficient expansions of rzvals
      ! ichuse(itri) - tells which chunk or (r,z) to use for patch itri
      !
      !    Output:
      ! xyz - point on the sphere
      ! dxyzduv - first derivative information
      !
      !


      ich = ichuse(itri)
      x0=triainfo(1,1,itri)
      y0=triainfo(2,1,itri)
      z0=triainfo(3,1,itri)

      x1=triainfo(1,2,itri)
      y1=triainfo(2,2,itri)
      z1=triainfo(3,2,itri)

      x2=triainfo(1,3,itri)
      y2=triainfo(2,3,itri)
      z2=triainfo(3,3,itri)

      !
      ! ... process the geometry, return the point location on the sphere
      ! and the derivatives with respect to u and v
      !
      s=x0+u*(x1-x0)+v*(x2-x0)
      t=y0+u*(y1-y0)+v*(y2-y0)

      dsdu = x1-x0
      dtdu = y1-y0
    
      dsdv = x2-x0
      dtdv = y2-y0
!
! s is coordinate in parametrization of r,z space
! t is coorindate in theta space
!


!
!  Compute r,z,drds,dzds
! 
!
      r = 0
      z = 0
      drds = 0
      dzds = 0
      pols = 0
      k = ixys2d(ich+1)-ixys2d(ich)
      call legepols(s,k-1,pols)
      do j=1,k
        r = r + rzcoefs(1,ixys2d(ich)+j-1)*pols(j)
        z = z + rzcoefs(2,ixys2d(ich)+j-1)*pols(j)
        drds = drds + rzcoefs(3,ixys2d(ich)+j-1)*pols(j)
        dzds = dzds + rzcoefs(4,ixys2d(ich)+j-1)*pols(j)
      enddo


      !
      !
      !
      xyz(1)= r*cos(t)
      xyz(2)= r*sin(t)
      xyz(3)= z


      ! du
      dxyzduv(1,1) = drds*cos(t)*dsdu - r*sin(t)*dtdu
      dxyzduv(2,1) = drds*sin(t)*dsdu + r*cos(t)*dtdu 
      dxyzduv(3,1) = dzds*dsdu

      ! dv
      dxyzduv(1,2) = drds*cos(t)*dsdv - r*sin(t)*dtdv
      dxyzduv(2,2) = drds*sin(t)*dsdv + r*cos(t)*dtdv 
      dxyzduv(3,2) = dzds*dsdv

      return
      end subroutine xtri_axissym_chunk


!
!
!
!
!
      subroutine xtri_axissym_fun_chunk(itri, u, v, xyz, dxyzduv, & 
          triainfo,np,pars,fcurve)
      implicit real *8 (a-h,o-z)
      real *8 :: xyz(3), dxyzduv(3,2), triainfo(3,3,*)
      real *8 :: pols(100),pars(np)
      integer :: np
      external fcurve
      



      !
      ! project the triangle itri in triainfo onto the sphere
      !
      !    Input:
      ! itri - triangle number to map
      ! u,v - local uv coordinates on triangle itri
      ! triainfo - flat skeleton triangle info
      ! rzcoefs - legendre coefficient expansions of rzvals
      ! ichuse(itri) - tells which chunk or (r,z) to use for patch itri
      !
      !    Output:
      ! xyz - point on the sphere
      ! dxyzduv - first derivative information
      !
      !


      x0=triainfo(1,1,itri)
      y0=triainfo(2,1,itri)
      z0=triainfo(3,1,itri)

      x1=triainfo(1,2,itri)
      y1=triainfo(2,2,itri)
      z1=triainfo(3,2,itri)

      x2=triainfo(1,3,itri)
      y2=triainfo(2,3,itri)
      z2=triainfo(3,3,itri)

      !
      ! ... process the geometry, return the point location on the sphere
      ! and the derivatives with respect to u and v
      !
      s=x0+u*(x1-x0)+v*(x2-x0)
      t=y0+u*(y1-y0)+v*(y2-y0)

      dsdu = x1-x0
      dtdu = y1-y0
    
      dsdv = x2-x0
      dtdv = y2-y0
!
! s is coordinate in parametrization of r,z space
! t is coorindate in theta space
!


!
!  Compute r,z,drds,dzds
! 
!
      
      call fcurve(s,np,pars,r,z,drds,dzds,tmp1,tmp1)


      !
      !
      !
      xyz(1)= r*cos(t)
      xyz(2)= r*sin(t)
      xyz(3)= z


      ! du
      dxyzduv(1,1) = drds*cos(t)*dsdu - r*sin(t)*dtdu
      dxyzduv(2,1) = drds*sin(t)*dsdu + r*cos(t)*dtdu 
      dxyzduv(3,1) = dzds*dsdu

      ! dv
      dxyzduv(1,2) = drds*cos(t)*dsdv - r*sin(t)*dtdv
      dxyzduv(2,2) = drds*sin(t)*dsdv + r*cos(t)*dtdv 
      dxyzduv(3,2) = dzds*dsdv

      return
      end subroutine xtri_axissym_fun_chunk






