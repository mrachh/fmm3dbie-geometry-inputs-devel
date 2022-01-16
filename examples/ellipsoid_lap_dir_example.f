      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:)
      character *100 fname
      integer ipars(2)

      real *8, allocatable :: targs(:,:),uvs_targ(:,:)
      integer, allocatable :: ipatch_id(:)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3)
      real *8, allocatable :: sigma(:),rhs(:)
      real *8 did
      real *8, allocatable :: errs(:)
      real *8 dpars(2)
      integer numit,niter

      real *8 pot,potex
      complex *16 ztmp,ima

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4

      a = 1.1d0
      b = 2.3d0
      c = 1.3d0

      rmax = 1.0d0
      ifc = 0
      call get_ellipsoid_mem(a,b,c,rmax,ifc,npatches)
      call prinf('npatches=*',npatches,1)



      norder = 4 
      npols = (norder+1)*(norder+2)/2

      npts = npatches*npols
      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))

      call get_ellipsoid_geom(a,b,c,rmax,ifc,norder,npatches,
     1   npts,norders,ixyzs,iptype,srcvals,srccoefs)
      call prinf('npatches=*',npatches,1)
      call prinf('npts=*',npts,1)


      ifinout = 1

      r1 = 1.5d0 + 0.5d0*hkrand(0)
      r2 = 0.5d0*hkrand(0)
      if(ifinout.eq.1) then
        rin = r1
        rout = r2
      else
        rin = r2
        rout = r1
      endif

      thet = hkrand(0)*pi
      phi = hkrand(0)*2*pi
      xyz_in(1) = a*rin*sin(thet)*cos(phi)
      xyz_in(2) = b*rin*sin(thet)*sin(phi)
      xyz_in(3) = c*rin*cos(thet)

      thet = hkrand(0)*pi
      phi = hkrand(0)*2*pi
      xyz_out(1) = a*rout*sin(thet)*cos(phi)
      xyz_out(2) = b*rout*sin(thet)*sin(phi)
      xyz_out(3) = c*rout*cos(thet)

      dpars(1) = 1.0d0
      dpars(2) = 1.0d0

      allocate(wts(npts))
      allocate(targs(3,npts),ipatch_id(npts),uvs_targ(2,npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)


      ndtarg = 3
     
      do i=1,npts
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
        targs(3,i) = srcvals(3,i)
      enddo

      do i=1,npts
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo

      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts, 
     1         ipatch_id,uvs_targ)


      allocate(sigma(npts),rhs(npts))

      do i=1,npts
        rhs(i) = 0
      enddo


      do i=1,npts
        call l3d_slp(xyz_out,3,srcvals(1,i),0,dpars,0,zpars,0,ipars,
     1     rhs(i))
      enddo

      

      eps = 0.51d-6
      numit = 50
      allocate(errs(numit+1))
      

      niter = 0
      rres = 0
      call lap_comb_dir_solver(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,eps,dpars,numit,ifinout,rhs,eps,
     2  niter,errs,rres,sigma)



c
c       test solution at a point
c
      call l3d_slp(xyz_out,3,xyz_in,0,dpars,0,zpars,0,ipars,potex)
      pot = 0
      do i=1,npts
        call l3d_comb(srcvals(1,i),3,xyz_in,2,dpars,0,zpars,0,ipars,
     1     ztmp)
        pot = pot + sigma(i)*wts(i)*ztmp
      enddo

      erra = abs(pot-potex)/abs(potex)
      call prin2('relative error in solve=*',erra,1)

      stop
      end




      subroutine setup_geom(igeomtype,norder,npatches,ipars, 
     1    srcvals,srccoefs,ifplot,fname)
      implicit real *8 (a-h,o-z)
      integer igeomtype,norder,npatches,ipars(*),ifplot
      character (len=*) fname
      real *8 srcvals(12,*), srccoefs(9,*)
      real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)

      real *8, pointer :: ptr1,ptr2,ptr3,ptr4
      integer, pointer :: iptr1,iptr2,iptr3,iptr4
      real *8, target :: p1(10),p2(10),p3(10),p4(10)
      real *8, allocatable, target :: triaskel(:,:,:)
      real *8, allocatable, target :: deltas(:,:)
      integer, allocatable :: isides(:)
      integer, target :: nmax,mmax

      procedure (), pointer :: xtri_geometry


      external xtri_stell_eval,xtri_sphere_eval
      
      npols = (norder+1)*(norder+2)/2
      allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
      allocate(wts(npols))

      call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)

      if(igeomtype.eq.1) then
        itype = 2
        allocate(triaskel(3,3,npatches))
        allocate(isides(npatches))
        npmax = npatches
        ntri = 0
        call xtri_platonic(itype, ipars(1), npmax, ntri, 
     1      triaskel, isides)

        xtri_geometry => xtri_sphere_eval
        ptr1 => triaskel(1,1,1)
        ptr2 => p2(1)
        ptr3 => p3(1)
        ptr4 => p4(1)


        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1         ptr3,ptr4, norder,'Triangulated surface of the sphere')
        endif


        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)
      endif

      if(igeomtype.eq.2) then
        done = 1
        pi = atan(done)*4
        umin = 0
        umax = 2*pi
        vmin = 2*pi
        vmax = 0
        allocate(triaskel(3,3,npatches))
        nover = 0
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),
     1     nover,npatches,npatches,triaskel)

        mmax = 2
        nmax = 1
        xtri_geometry => xtri_stell_eval

        allocate(deltas(-1:mmax,-1:nmax))
        deltas(-1,-1) = 0.17d0
        deltas(0,-1) = 0
        deltas(1,-1) = 0
        deltas(2,-1) = 0

        deltas(-1,0) = 0.11d0
        deltas(0,0) = 1
        deltas(1,0) = 4.5d0
        deltas(2,0) = -0.25d0

        deltas(-1,1) = 0
        deltas(0,1) = 0.07d0
        deltas(1,1) = 0
        deltas(2,1) = -0.45d0

        ptr1 => triaskel(1,1,1)
        ptr2 => deltas(-1,-1)
        iptr3 => mmax
        iptr4 => nmax

        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1         iptr3,iptr4, norder,
     2         'Triangulated surface of the stellarator')
        endif

        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,iptr3,iptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)
      endif
      
      return  
      end

