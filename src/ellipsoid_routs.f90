
      subroutine get_ellipsoid_mem(a,b,c,rmax,ifc,npatches)
!
!  This subroutine estimates the number of patches required
!  discretizing an ellipsoid with half axis lengths (a,b,c).
!
!  The mesh can either be conforming or non-conforming with
!  a specified maximum patch size
!
!
!  Input arguments:
!    - a: real *8
!        half axis length in the x direction
!    - b: real *8
!        half axis length in the y direction
!    - c: real *8
!        half axis length in the z direction
!    - rmax: real *8
!        maximum patch size
!    - ifc: integer
!        flag for determining if the mesh should be conforming
!        or not
!        ifc = 1, then mesh is conforming
!        ifc = 0, then mesh is non-conforming
!
!  Output arguments:
!    - npatches: integer
!        number of patches in the discretization
!
      implicit real *8 (a-h,o-z)
      real *8, intent(in) :: a,b,c,rmax
      integer, intent(in) :: ifc
      integer, intent(out) :: npatches

      done = 1.0d0
      pi = atan(done)*4

      npatches = 0

      nthet = ceiling(2*c/rmax)
      if(ifc.eq.1) then

        alpha = a
        beta = b
        hh = (alpha-beta)**2/(alpha+beta)**2 

        ellip_p = pi*(alpha + beta)* &
           (1.0d0 + 3*hh/(10.0d0 + sqrt(4.0d0-3*hh)))
        
        nphi = ceiling(ellip_p/rmax)
        npatches = 2*nthet*nphi
      else
        hthet = pi/(nthet+0.0d0)
        do ithet=1,nthet
          t0 = (ithet-1)*hthet
          t1 = (ithet)*hthet

          tuse = t0
          if(abs(t0-pi/2).ge.abs(t1-pi/2)) tuse = t1
          if((t0-pi/2)*(t1-pi/2).le.0) tuse = pi/2


          alpha = a*sin(tuse)
          beta = b*sin(tuse)
          hh = (alpha-beta)**2/(alpha+beta)**2 

          ellip_p = pi*(alpha + beta)* &
             (1.0d0 + 3*hh/(10.0d0 + sqrt(4.0d0-3*hh)))
        
          nphi = ceiling(ellip_p/rmax)
          npatches = npatches + 2*nphi
        enddo
      endif

      return
      end subroutine get_ellipsoid_mem




      subroutine get_ellipsoid_geom(a,b,c,rmax,ifc,norder,npatches, & 
        npts,norders,ixyzs,iptype,srcvals,srcoefs)
!
!  This subroutine discretizes 
!  an ellipsoid with half axis lengths (a,b,c).
!
!
!  Input arguments:
!    - a: real *8
!        half axis length in the x direction
!    - b: real *8
!        half axis length in the y direction
!    - c: real *8
!        half axis length in the z direction
!    - rmax: real *8
!        maximum patch size
!    - ifc: integer
!        flag for determining if the mesh should be conforming
!        or not
!        ifc = 1, then mesh is conforming
!        ifc = 0, then mesh is non-conforming
!    - norder:
!        order of discretization
!    - npatches: integer
!        number of patches in the discretization can be computed
!        by a call to get_ellipsoid_mem
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
      real *8, intent(in) :: a,b,c,rmax
      integer, intent(in) :: ifc,npatches,norder,npts
      integer, intent(out) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(out) :: iptype(npatches)
      real *8, intent(out) :: srcvals(12,npts),srccoefs(9,npts)


      real *8 v1(3),v2(3),v3(3),v4(3),xyz0(1:3)
      real *8, allocatable, target :: triaskel(:,:,:)
      real *8, pointer :: ptr1,ptr2,ptr3,ptr4
      real *8, allocatable :: uvs(:,:),wts(:),umatr(:,:),vmatr(:,:)
      real *8, target :: p1(10),p2(10),p3(10),p4(10)
      external xtri_ellipsoid_eval

      npols = (norder+1)*(norder+2)/2
      do i=1,npatches
        norders(i) = norder
        ixyzs(i) = (i-1)*npols + 1
        iptype(i) = 1
      enddo
      ixyzs(npatches+1) = npts+1

      allocate(triaskel(3,3,npatches))

      done = 1.0d0
      pi = atan(done)*4


      nthet = ceiling(2*c/rmax)
      vmin = 0
      vmax = 2*pi
      nover = 0

      if(ifc.eq.1) then

        n0 = 0
        alpha = a
        beta = b
        hh = (alpha-beta)**2/(alpha+beta)**2 

        ellip_p = pi*(alpha + beta)* &
           (1.0d0 + 3*hh/(10.0d0 + sqrt(4.0d0-3*hh)))
        
        umin = 0
        umax = pi
        nphi = ceiling(ellip_p/rmax)
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,nthet,nphi,nover, &
          npatches,n0,triaskel)
      else
        hthet = pi/(nthet+0.0d0)
        istart = 1
        nthet0 = 1
        do ithet=1,nthet
          n0 = 0
          umin = (ithet-1)*hthet
          umax = (ithet)*hthet

          tuse = umin
          if(abs(t0-pi/2).ge.abs(t1-pi/2)) tuse = umax
          if((t0-pi/2)*(t1-pi/2).le.0) tuse = pi/2


          alpha = a*sin(tuse)
          beta = b*sin(tuse)
          hh = (alpha-beta)**2/(alpha+beta)**2 

          ellip_p = pi*(alpha + beta)* &
             (1.0d0 + 3*hh/(10.0d0 + sqrt(4.0d0-3*hh)))
        
          nphi = ceiling(ellip_p/rmax)
          call xtri_rectmesh_ani(umin,umax,vmin,vmax,nthet0,nphi,nover, &
            npatches,n0,triaskel(1,1,istart))
          istart = istart + 2*nphi
        enddo
      endif
      xyz0(1:3) = 0

      p2(1) = a
      p2(2) = b
      p2(3) = c

      p3(1:3) = xyz0(1:3)
      
      ptr1 => triaskel(1,1,1)
      ptr2 => p2(1)
      ptr3 => p3(1)
      ptr4 => p4(1)

      npols = (norder+1)*(norder+2)/2
      allocate(uvs(2,npols),wts(npols),umatr(npols,npols), &
        vmatr(npols,npols))
      call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)
      call getgeominfo(npatches,xtri_ellipsoid_eval,ptr1,ptr2,ptr3, &
        ptr4,npols,uvs,umatr,srcvals,srccoefs)


      return
      end subroutine get_ellipsoid_geom
!
!
!
!
!



      subroutine xtri_ellipsoid_eval(itri, u, v, xyz, dxyzduv, & 
          triainfo,p2, p3, p4)
      implicit real *8 (a-h,o-z)
      real *8 :: xyz(3), dxyzduv(3,2), triainfo(3,3,*),p2(3),p3(3)


      !
      ! Evaluate chart of ellipsoid in polar coordinates, u
      ! is elevation, and phi is azimuth
      !
      !    Input:
      ! itri - triangle number to map
      ! u,v - local uv coordinates on triangle itri
      ! triainfo - flat skeleton triangle info
      ! p2,p3,p4 - dummy parameters
      !
      !    Output:
      ! xyz - point on the ellipsoid
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
      x=x0+u*(x1-x0)+v*(x2-x0)
      y=y0+u*(y1-y0)+v*(y2-y0)
      z=z0+u*(z1-z0)+v*(z2-z0)

      dxdu = x1-x0
      dydu = y1-y0
      dzdu = z1-z0
    
      dxdv = x2-x0
      dydv = y2-y0
      dzdv = z2-z0

      !
      ! second derivatives are zero...
      !

      !
      ! evaluate ellipsoid chart 
      !
      xyz(1)=p2(1)*sin(x)*cos(y) + p3(1)
      xyz(2)=p2(2)*sin(x)*sin(y) + p3(2)
      xyz(3)=p2(3)*cos(x) + p3(3)


      ! du
      dxyzduv(1,1) = p2(1)*(cos(x)*cos(y)*dxdu - sin(x)*sin(y)*dydu)
      dxyzduv(2,1) = p2(2)*(cos(x)*sin(y)*dxdu + sin(x)*cos(y)*dydu)
      dxyzduv(3,1) = p2(3)*(-sin(x)*dxdu)

      ! dv
      dxyzduv(1,2) = p2(1)*(cos(x)*cos(y)*dxdv - sin(x)*sin(y)*dydv)
      dxyzduv(2,2) = p2(2)*(cos(x)*sin(y)*dxdv + sin(x)*cos(y)*dydv)
      dxyzduv(3,2) = p2(3)*(-sin(x)*dxdv)

      return
      end subroutine xtri_ellipsoid_eval








