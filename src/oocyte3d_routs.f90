      subroutine get_oocyte3d_mem(t,p1,p2,k,iref,nmid,nch2d, &
         rmax,npatches)
!
!  This subroutine estimates the number of patches required
!  for parametrizing an oocyte like geometry parameterized via
!  the parameters t,p1,p2
!
!  r(z) = t*2^(p1+p2-1)*z^(p1)*(1-z)^p2, z \in [0,1]
!
!  The generating curve is discretized via the parameters,
!  iref, nmid, leading to a total number of chunks
!  2*(iref+1) + nmid. The interval [0,1] is first uniformly
!  discretized using nmid+2 patches, and the first and last
!  patch are then dyadically refined for iref levels
!
!  Input arguments:
!    - t: real *8
!         parameter governing the scaling of r(z)
!    - p1: real *8
!         power of the singularity at z=0
!    - p2: real *8
!         power of the singularity at z=1 
!    - k: integer
!        order of discretization for chunking the generating curve
!    - iref: integer
!        number of dyadic refinements at either end points
!    - nmid: integer
!        number of panels in the middle section
!    - nch2d: integer
!        number of chunks on the generating curve = 2*(iref+1) + nmid
!    - rmax: real *8
!        max size of patches in the 3d discretization
!  
!  Output arguments:
!    - npatches: integer
!        number of patches in the 3d discretization
!
!
!
      implicit real *8 (a-h,o-z)
      real *8 t,p1,p2,rmax
      integer nch2d,npatches,k
      real *8, allocatable :: srcvals2d(:,:),srccoefs2d(:,:)
      integer, allocatable :: ixys2d(:),iptype2d(:)

      npts2d = k*nch2d
      allocate(srcvals2d(8,npts2d),srccoefs2d(6,npts2d))
      allocate(ixys2d(nch2d+1),iptype2d(nch2d))

      call get_oocyte3d_chunks(t,p1,p2,k,iref,nmid,nch2d,npts2d, &
        iptype2d,ixys2d,srccoefs2d,srcvals2d)

      call get_axissym_mem(nch2d,npts2d,ixys2d,iptype2d,srcvals2d, &
        srccoefs2d,rmax,npatches)
      
      return
      end subroutine get_oocyte3d_mem
!
!
!
!
      subroutine get_oocyte3d_geom(t,p1,p2,k,iref,nmid,nch2d, &
         rmax,norder,npatches,npts,norders,ixyzs,iptype,srcvals, &
         srccoefs)
!
!  This subroutine discretizes an oocyte like geometry parameterized via
!  the parameters t,p1,p2
!
!  r(z) = t*2^(p1+p2-1)*z^(p1)*(1-z)^p2, z \in [0,1]
!
!  The generating curve is discretized via the parameters,
!  iref, nmid, leading to a total number of chunks
!  2*(iref+1) + nmid. The interval [0,1] is first uniformly
!  discretized using nmid+2 patches, and the first and last
!  patch are then dyadically refined for iref levels
!
!  Input arguments:
!    - t: real *8
!         parameter governing the scaling of r(z)
!    - p1: real *8
!         power of the singularity at z=0
!    - p2: real *8
!         power of the singularity at z=1 
!    - k: integer
!        order of discretization for chunking the generating curve
!    - iref: integer
!        number of dyadic refinements at either end points
!    - nmid: integer
!        number of panels in the middle section
!    - nch2d: integer
!        number of chunks on the generating curve = 2*(iref+1) + nmid
!    - rmax: real *8
!        max size of patches in the 3d discretization
!  
!  Output arguments:
!    - npatches: integer
!        number of patches in the 3d discretization
!
!
!
      implicit real *8 (a-h,o-z)
      real *8 t,p1,p2,rmax
      integer nch2d,npatches,k
      integer npts,norder
      integer norders(npatches),ixyzs(npatches+1),iptype(npatches)
      real *8 srcvals(12,npts),srccoefs(9,npts)
      real *8, allocatable :: srcvals2d(:,:),srccoefs2d(:,:)
      integer, allocatable :: ixys2d(:),iptype2d(:)

      npts2d = k*nch2d
      allocate(srcvals2d(8,npts2d),srccoefs2d(6,npts2d))
      allocate(ixys2d(nch2d+1),iptype2d(nch2d))

      call get_oocyte3d_chunks(t,p1,p2,k,iref,nmid,nch2d,npts2d, &
        iptype2d,ixys2d,srccoefs2d,srcvals2d)

      call get_axissym_geom(nch2d,npts2d,ixys2d,iptype2d,srcvals2d, &
        srccoefs2d,rmax,norder,npatches,npts,norders,ixyzs,iptype, &
        srcvals,srccoefs)
      
      return
      end subroutine get_oocyte3d_geom
!
!
!
!
!
!
      subroutine get_oocyte3d_chunks(t,p1,p2,k,iref,nmid,nch2d,npts2d, &
        iptype2d,ixys2d,srccoefs2d,srcvals2d)
      implicit real *8 (a-h,o-z)
      real *8 t,p1,p2
      integer k,nch2d,npts2d
      integer iptype2d(nch2d),ixys2d(nch2d+1)
      real *8 srccoefs2d(6,npts2d),srcvals2d(8,npts2d)
      real *8, allocatable :: ts(:),ws(:),umat(:,:),vmat(:,:)
      real *8 pars(3)
      real *8, allocatable :: tchse(:)
      
      allocate(tchse(nch2d+1))
      
      tchse(1) = 0.0d0
      hpan = 1.0d0/(nmid+2)
      rpan = hpan*(0.5d0)**(iref)
      tchse(2) = tchse(1) + rpan
      do i=2,iref+1
        tchse(i+1) = tchse(i) + rpan
        rpan = rpan*2
      enddo

      do i=iref+2,iref+1+nmid
        tchse(i+1) = tchse(i) + hpan
      enddo

      rpan = hpan/2
      do i=iref+2+nmid,nch2d-1
        tchse(i+1) = tchse(i) + rpan
        rpan = rpan/2
      enddo
      tchse(nch2d+1) = 1.0d0

      itype = 2
      allocate(ts(k),ws(k),umat(k,k),vmat(k,k))
      call legeexps(itype,k,ts,umat,vmat,ws)

      do i=1,nch2d
        ixys2d(i) = (i-1)*k + 1
        iptype2d(i) = 1
      enddo

      ixys2d(nch2d+1) = npts2d+1
      pars(1) = t
      pars(2) = p1
      pars(3) = p2
      
      h = 1.0d0/nch2d
      do ich=1,nch2d
        tchstart = tchse(ich) 
        tchend = tchse(ich+1) 

        hh = (tchend-tchstart)/2.0d0

        do j=1,k
          ipt = (ich-1)*k + j
          tt = tchstart + (ts(j)+1.0d0)/2.0d0*(tchend-tchstart)

          call funcurve_oocyte(tt,pars,x,y,dx,dy,dx2,dy2)
          srcvals2d(1,ipt) = x
          srcvals2d(2,ipt) = y
          srcvals2d(3,ipt) = dx*hh
          srcvals2d(4,ipt) = dy*hh
          srcvals2d(5,ipt) = dx2*hh**2
          srcvals2d(6,ipt) = dy2*hh**2
          ds = sqrt(dx**2 + dy**2)
          srcvals2d(7,ipt) = dy/ds
          srcvals2d(8,ipt) = -dx/ds
        enddo
        do j=1,k
          jpt = (ich-1)*k + j
          srccoefs2d(1:6,jpt) = 0
          do l=1,k
            lpt = (ich-1)*k + l
            srccoefs2d(1:6,jpt) = srccoefs2d(1:6,jpt) + &
               umat(j,l)*srcvals2d(1:6,lpt)
          enddo
        enddo
      enddo


      return
      end subroutine get_oocyte3d_chunks
!
!
!
!
!

      subroutine funcurve_oocyte(t,pars,x,y,dxdt,dydt,d2xdt2,d2ydt2)
      implicit real *8 (a-h,o-z)
      real *8 t, pars(3),x,y,dxdt,dydt,d2xdt2,d2ydt2
      
      y = t
      dydt = 1
      d2ydt2 = 0
      
      tt = pars(1)
      p1 = pars(2)
      p2 = pars(3)
      
      c = tt*2.0d0**(p1+p2)/2.0d0
      x1 = t**p1
      dx1dt = p1*t**(p1-1.0d0)
      dx1dt2 = p1*(p1-1.0d0)*t**(p1-2.0d0)
      
      x2 = (1.0d0-t)**p2
      dx2dt = -p2*(1.0d0-t)**(p2-1.0d0)
      dx2dt2 = p2*(p2-1.0d0)*(1.0d0-t)**(p2-2.0d0)

      x = c*x1*x2
      dxdt = c*(dx1dt*x2 + x1*dx2dt)
      d2xdt2 = c*(dx1dt2*x2 + 2*dx1dt*dx2dt + x1*dx2dt2)
      
      return
      end subroutine funcurve_oocyte
