      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:)
      character *100 fname
      integer ipars(2)

      real *8, allocatable :: targs(:,:)
      integer, allocatable :: ipatch_id(:)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3)
      real *8, allocatable :: rhs(:,:),rhs_interp(:,:),rhs_ex(:,:)
      real *8 did
      real *8, allocatable :: errs(:)
      real *8 dpars(2)
      integer numit,niter

      real *8, allocatable :: uvs_targ(:,:)
      integer, allocatable :: ipatchtarg(:),ixmattarg(:)
      real *8, allocatable :: xmattarg(:)

      real *8, allocatable :: xyztarg(:,:),ftarg(:)
      integer, allocatable :: itarg(:)

      real *8 pot,potex
      complex *16 ztmp,ima

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4

      a = 1.1d0
      b = 2.3d0
      c = 1.3d0

      rmax = 1.0d0/2
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


      rout = 1.5d0 + 0.5d0*hkrand(0)
      rin = 0.5d0*hkrand(0)

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

      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)
      nd = 2
!
!
!  Create two test functions due to point source inside and outside
!
      allocate(rhs(2,npts))

      ra = 0

      do i=1,npts
        call l3d_slp(xyz_in,3,srcvals(1,i),0,dpars,0,zpars,0,ipars,
     1     rhs(1,i))
        call l3d_slp(xyz_out,3,srcvals(1,i),0,dpars,0,zpars,0,ipars,
     1     rhs(2,i))
        ra = ra + (rhs(1,i)**2 + rhs(2,i)**2)*wts(i)
      enddo
      ra = sqrt(ra)

      call prin2('ra=*',ra,1)



      ntarg = 100
      allocate(targs(3,ntarg),rhs_ex(nd,ntarg),rhs_interp(nd,ntarg))
      do i=1,ntarg
        thet = hkrand(0)*pi
        phi = hkrand(0)*2*pi
        targs(1,i) = a*sin(thet)*cos(phi)
        targs(2,i) = b*sin(thet)*sin(phi)
        targs(3,i) = c*cos(thet)
        call l3d_slp(xyz_in,3,targs(1,i),0,dpars,0,zpars,0,ipars,
     1     rhs_ex(1,i))
        call l3d_slp(xyz_out,3,targs(1,i),0,dpars,0,zpars,0,ipars,
     1     rhs_ex(2,i))
      enddo

      call ellipsoid_interp(nd,a,b,c,rmax,ifc,ntarg,targs,npatches,
     1  norders,ixyzs,iptype,npts,rhs,rhs_interp)
      erra = 0
      do i=1,ntarg
        do idim=1,nd
          erra = erra + (rhs_interp(idim,i) - rhs_ex(idim,i))**2
        enddo
      enddo

      erra = sqrt(erra/ntarg)
      erra = erra
      call prin2('relative l2 error in interpolated value at targets=*',
     1   erra,1)

      allocate(ipatchtarg(ntarg),uvs_targ(2,ntarg))
      call ellipsoid_local_coord_targ(a,b,c,rmax,ifc,ntarg,targs,
     1   npatches,norders,ixyzs,iptype,npts,ipatchtarg,uvs_targ)
      

      lmem = 0

      call get_surf_interp_mat_targ_mem(npatches,ixyzs,ntarg,
     1  ipatchtarg,lmem)
      call prinf('lmem =*',lmem,1)
      call prinf('expected lmem=*',ntarg*npols,1)

      allocate(ixmattarg(ntarg+1),xmattarg(lmem))
      call get_surf_interp_mat_targ(npatches,norders,ixyzs,iptype,npts,
     1 ntarg,ipatchtarg,uvs_targ,lmem,xmattarg,ixmattarg)

c
c  check interpolation error using xmattarg
c
c
      rhs_interp = 0

      do i=1,ntarg
        ipatch = ipatchtarg(i)
        istart = ixyzs(ipatch)
        istart2 = ixmattarg(i)
        npols_use = ixyzs(ipatch+1) - ixyzs(ipatch) 
        do l=1,npols
          do idim = 1,nd
            rhs_interp(idim,i) = rhs_interp(idim,i) +
     1         xmattarg(istart2+l-1)*rhs(idim,istart+l-1)
          enddo
        enddo
      enddo


      erra = 0
      do i=1,ntarg
        do idim=1,nd
          erra = erra + (rhs_interp(idim,i) - rhs_ex(idim,i))**2
        enddo
      enddo

      erra = sqrt(erra/ntarg)
      erra = erra
      call prin2('relative l2 error in interpolated value at targets=*',
     1   erra,1)

      
      ntarg = 1000
      allocate(xyztarg(3,ntarg))

      allocate(itarg(ntarg),ftarg(ntarg))

      xmax = 2*a
      xmin = -2*a
      ymax = 2*b
      ymin = -2*b
      zmax = 2*c
      zmin = -2*c

      do i=1,ntarg
        xyztarg(1,i) = hkrand(0)*(xmax-xmin) + xmin
        xyztarg(2,i) = hkrand(0)*(ymax-ymin) + ymin
        xyztarg(3,i) = hkrand(0)*(zmax-zmin) + zmin
        call check_ellipsoid_interior(a,b,c,xyztarg(1,i),itarg(i))
        ftarg(i) = itarg(i) + 0.0d0
      enddo

      call surf_vtk_plot(npatches,norders,ixyzs,iptype,npts,srccoefs,
     1   srcvals,'ellipsoid.vtk','a')
      call vtk_scatter_plot_scalar(ntarg,xyztarg,
     1   ftarg,'ellip_pts.vtk','a')
      


      

      stop
      end


