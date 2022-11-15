      subroutine get_torus_mem(rmajor,rminor,rmax,npatches)
!
!  This subroutine estimates the number of patches required in
!  the discretization of a standard torus defined via the major and
!  minor axes given by rmajor=r1, and rminor=r2. The parametrization
!  of the surface is given by
!
!  X(u,v) = (r1 + r2*cos(v))*cos(u)
!  Y(u,v) = (r1 + r2*cos(v))*sin(u)
!  Z(u,v) = r2*sin(v)
!
!  Input arguments:
!    - rmajor: real *8
!        major axis of the torus
!    - rminor: real *8
!        minor axis of the torus
!    - rmax: real *8
!        max size of patch (this is just an estimate, rmax 
!        is used as a heuristic to divide the u axis in 2*pi*(r1+r2)/rmax
!        chunks, and the v axis into 2*pi*r2/rmax chunks
!
!  Output arguments:
!    - npatches: integer
!        number of patches used in the discretization
!
!
      implicit none
      real *8, intent(in) :: rmajor,rminor,rmax
      integer, intent(out) :: npatches

      real *8 pi,done
      integer nu,nv


      done = 1.0d0
      pi = atan(done)*4.0d0

      nu = max(ceiling(2*pi*(rmajor+rminor)/rmax),3)
      nv = max(ceiling(2*pi*rminor/rmax),3)
      npatches = 2*nu*nv

      return
      end



      subroutine get_torus_geom(rmajor,rminor,rmax,norder,npatches, &
        npts,norders,ixyzs,iptype,srcvals,srccoefs)
!
!  This subroutine estimates the number of patches required in
!  the discretization of a standard torus defined via the major and
!  minor axes given by rmajor=r1, and rminor=r2. 
!  The parametrization of the surface is given by
!
!  X(u,v) = (r1 + r2*cos(v))*cos(u)
!  Y(u,v) = (r1 + r2*cos(v))*sin(u)
!  Z(u,v) = r2*sin(v)
!
!  Input arguments:
!    - rmajor: real *8
!        major axis of the torus
!    - rminor: real *8
!        minor axis of the torus
!    - rmax: real *8
!        max size of patch (this is just an estimate, rmax 
!        is used as a heuristic to divide the u axis in 2*pi*(r1+r2)/rmax
!        chunks, and the v axis into 2*pi*r2/rmax chunks
!    - norder:
!        order of discretization
!    - npatches: integer
!        number of patches used in the discretization
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
!
      implicit real *8 (a-h,o-z)
      real *8, intent(in) :: rmajor,rminor,rmax
      integer, intent(in) :: npatches,npts,norder
      integer, intent(out) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(out) :: iptype(npatches)
      real *8, intent(out) :: srcvals(12,npts),srccoefs(9,npts)
      
      real *8 pi,done
      integer nu,nv
      real *8 umin,umax,vmin,vmax
      real *8, allocatable, target :: triaskel(:,:,:)
      integer nover
      real *8, target :: p1(3),p2(3),p3(1),p4(1)
      real *8, pointer :: ptr1,ptr2,ptr3,ptr4

      real *8, allocatable :: uvs(:,:),wts(:),umatr(:,:),vmatr(:,:)

      procedure (), pointer :: xtri_geometry

      external xtri_wtorus_eval
      


      done = 1.0d0
      pi = atan(done)*4.0d0

      nu = max(ceiling(2*pi*(rmajor+rminor)/rmax),3)
      nv = max(ceiling(2*pi*rminor/rmax),3)
      npols = (norder+1)*(norder+2)/2
      allocate(uvs(2,npols),wts(npols),umatr(npols,npols))
      allocate(vmatr(npols,npols))

      call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)

      umin = 0
      umax = 2*pi
      vmin = 0
      vmax = 2*pi

      allocate(triaskel(3,3,npatches))
      nover = 0
      npatches0 = 0 
      call xtri_rectmesh_ani(umin,umax,vmin,vmax,nu,nv,nover,npatches, &
        npatches0,triaskel)
      
      p1(1) = rminor
      p1(2) = rmajor
      p1(3) = 0.0d0
      
      p2(1) = 1.0d0
      p2(2) = 1.0d0
      p2(3) = 1.0d0

      p4(1) = 0.0d0
      ptr1 => triaskel(1,1,1)
      ptr2 => p1(1)
      ptr3 => p2(1)
      ptr4 => p4(1)

      xtri_geometry => xtri_wtorus_eval
      call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4, &
        npols,uvs,umatr,srcvals,srccoefs)
      
      do i=1,npatches
        norders(i) = norder
        ixyzs(i) = (i-1)*npols + 1
        iptype(i) = 1
      enddo
      ixyzs(npatches+1) = npts+1


      return
      end

!
!
!
!
!
      subroutine torus_interp(nd,rmajor,rminor,rmax,ntarg,targs,npatches, &
        norders,ixyzs,iptype,npts,u,uinterp)
!
!  This function interpolated a given function on surface 
!  to a collection of targets on a triangulated torus with major
!  and minor axis rmajor,rminor, maximum patch size rmax, at a collection
!  of user specified targets on the ellipsoid provided in their cartesian 
!  coordinates.
!
!  The parametrization of the surface is given by
!
!  X(u,v) = (r1 + r2*cos(v))*cos(u)
!  Y(u,v) = (r1 + r2*cos(v))*sin(u)
!  Z(u,v) = r2*sin(v)
!
!
!
!  Input arguments:
!    - nd: integer
!       number of functions
!    - rmajor: real *8
!       major axis of torus 
!    - rminor: real *8
!       minor axis of torus 
!    - rmax: real *8
!        maximum patch size
!    - ntarg: integer
!        number of targets
!    - targs: real *8 (3,ntarg)
!        xyz coordinates of the targets
!    - npatches: integer
!        number of patches in the discretization can be computed
!        by a call to get_ellipsoid_mem
!    - norders: integer(npatches)
!        discretization order of patches
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        number of points in the discretization = 
!          npatches*(norder+1)*(norder+2)/2
!    - u: real *8 (nd,npts)
!        function values on surface
!
!
!  Output arguments:
!    - uinterp: real *8 (nd,ntarg)
!        interpolated function values on surface
!

      implicit real *8 (a-h,o-z)
      real *8, intent(in) :: rmajor,rminor,rmax
      integer, intent(in) :: ntarg,nd,npts
      real *8, intent(in) :: targs(3,ntarg)
      integer, intent(in) :: npatches,norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) :: u(nd,npts)
      real *8, intent(out) :: uinterp(nd,ntarg)

      integer, allocatable :: ipatchtarg(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8, allocatable :: ucoefs(:,:),pols(:)

      real *8 pi,done,alpha,beta
      integer nu,nv

      allocate(ipatchtarg(ntarg),uvs_targ(2,ntarg))
      call torus_local_coord_targ(rmajor,rminor,rmax,ntarg,targs, &
        npatches,norders,ixyzs,iptype,npts,ipatchtarg,uvs_targ)
      allocate(ucoefs(nd,npts))

      call surf_vals_to_coefs(nd,npatches,norders,ixyzs,iptype,npts,&
        u,ucoefs)
      npols = (norders(1)+1)*(norders(2)+2)/2  

      nordermax = maxval(norders(1:npatches))
      npmax = (nordermax+1)*(nordermax+2)/2
      allocate(pols(npmax))

      incx = 1
      incy = 1
      alpha = 1.0d0
      beta = 0.0d0


!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,norder,npols,ii,pols)
      do i=1,ntarg
        norder = norders(ipatchtarg(i))
        npols = (norder+1)*(norder+2)/2
        call koorn_pols(uvs_targ(1,i),norder,npols,pols)
        ii = ixyzs(ipatchtarg(i))
        
        call dgemv_guru('n',nd,npols,alpha,ucoefs(1,ii),nd,pols,incx,beta, &
          uinterp(1,i),incy)

      enddo
!$OMP END PARALLEL DO      

      return
      end

!
!
!
!
!

      subroutine torus_local_coord_targ(rmajor,rminor,rmax,ntarg,targs, &
        npatches,norders,ixyzs,iptype,npts,ipatchtarg,uvs_targ)
      
!
!  This function extracts the local patch id and uv coordinates
!  for a collection of targets on a triangulated torus with major
!  and minor axis rmajor,rminor, maximum patch size rmax, at a collection
!  of user specified targets on the ellipsoid provided in their cartesian 
!  coordinates.
!
!  The parametrization of the surface is given by
!
!  X(u,v) = (r1 + r2*cos(v))*cos(u)
!  Y(u,v) = (r1 + r2*cos(v))*sin(u)
!  Z(u,v) = r2*sin(v)
!
!
!
!  Input arguments:
!    - rmajor: real *8
!       major axis of torus 
!    - rminor: real *8
!       minor axis of torus 
!    - rmax: real *8
!        maximum patch size
!    - ntarg: integer
!        number of targets
!    - targs: real *8 (3,ntarg)
!        xyz coordinates of the targets
!    - npatches: integer
!        number of patches in the discretization can be computed
!        by a call to get_ellipsoid_mem
!    - norders: integer(npatches)
!        discretization order of patches
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        number of points in the discretization = 
!          npatches*(norder+1)*(norder+2)/2
!
!  Output arguments:
!    - ipatchtarg: integer(ntarg)
!        ipatchtarg(i) is the patch on which target i is located
!    - uvs_targ: real *8 (2,ntarg)
!        local uv coordinates on patch
!

      implicit real *8 (a-h,o-z)
      real *8, intent(in) :: rmajor,rminor,rmax
      integer, intent(in) :: ntarg,npts
      real *8, intent(in) :: targs(3,ntarg)
      integer, intent(in) :: npatches,norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      integer, intent(out) :: ipatchtarg(ntarg)
      real *8, intent(out) :: uvs_targ(2,ntarg)

      real *8 pi,done
      integer nu,nv


      done = 1.0d0
      pi = atan(done)*4.0d0

      nu = max(ceiling(2*pi*(rmajor+rminor)/rmax),3)
      nv = max(ceiling(2*pi*rminor/rmax),3)
      
      hu = 2*pi/(nu+0.0d0)
      hv = 2*pi/(nv+0.0d0)


!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,cv,sv,v,cu,su,u) &
!$OMP PRIVATE(iu,iv,u0,v0,uuse,vuse)
      do i=1,ntarg
        cv = sqrt(targs(1,i)**2 + targs(2,i)**2) - rmajor
        sv = targs(3,i)
        v = atan2(sv,cv)
        if(v.lt.0) v = v + 2*pi
        cu = targs(1,i)
        su = targs(2,i)
        u = atan2(su,cu)
        if(u.lt.0) u = u + 2*pi

        iu = ceiling(u/hu)
        iv = ceiling(v/hv)
        if(iu.eq.0) iu = 1
        if(iv.eq.0) iv = 1

        u0 = (iu-1)*hu
        v0 = (iv-1)*hv
        uuse = (u-u0)/hu
        vuse = (v-v0)/hv


        if(uuse+vuse.le.1) then
          ipatchtarg(i) = (iu-1)*nv*2 + 2*(iv-1) + 1
          uvs_targ(1,i) = uuse
          uvs_targ(2,i) = vuse
        else
          ipatchtarg(i) = (iu-1)*nv*2 + 2*(iv-1) + 2
          uvs_targ(1,i) = 1.0d0-uuse
          uvs_targ(2,i) = 1.0d0-vuse
        endif
      enddo
!$OMP END PARALLEL DO      

      return
      end
