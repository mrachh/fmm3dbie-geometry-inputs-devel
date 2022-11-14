      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:)
      character *100 fname
      integer ipars(2)

      real *8, allocatable :: targs(:,:)
      integer, allocatable :: ipatch_id(:)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3)
      real *8 pars(3)
      
      real *8, allocatable :: rhs(:,:),rhs_interp(:,:),rhs_ex(:,:)
      real *8 did
      real *8, allocatable :: errs(:)
      real *8, allocatable :: tchse(:)
      integer, allocatable :: ipatchtarg(:),ixmattarg(:)
      real *8, allocatable :: xmattarg(:),uvs_targ(:,:)
      real *8 dpars(2)
      integer numit,niter

      real *8 pot,potex
      complex *16 ztmp,ima
      procedure (), pointer :: funcurve

      external funcurve_oocyte,funcurve_oocyte_riemann

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4

      ifriemann = 0

      if(ifriemann.eq.1) then

        a = 0.25d0
        b = 0.1d0

        np = 2

        pars(1) = a
        pars(2) = b
        funcurve => funcurve_oocyte_riemann
      else
        tt = 0.72d0
        p1 = 0.4d0
        p2 = 0.2d0

        np = 3
        pars(1) = tt
        pars(2) = p1
        pars(3) = p2
        funcurve => funcurve_oocyte


      endif

      k = 16
      iref = 0
      nmid = 1
      nch2d = 2*(iref+1) + nmid
      npts2d = nch2d*k

      allocate(tchse(nch2d+1))
      call get_oocyte3d_tchse(iref,nmid,nch2d,tchse)

      rmax = 1.0d0

      call get_axissym_fun_mem(nch2d,tchse,k,funcurve,
     1   np,pars,rmax,npatches)
      call prinf('npatches=*',npatches,1)
      norder = 4
      npols = (norder+1)*(norder+2)/2

      npts = npatches*npols
      print *, "npts=",npts
      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))

      iort = 1
      call get_axissym_fun_geom(nch2d,tchse,k,funcurve,
     1  np,pars,rmax,iort,norder,npatches,npts,norders,ixyzs,iptype,
     2  srcvals,srccoefs)
      call prinf('npatches=*',npatches,1)



      ifinout = 1

      r1 = 2.5d0 + 0.5d0*hkrand(0)
      r2 = 0.1d0*hkrand(0)
      if(ifinout.eq.1) then
        rin = r2
        rout = r1
      else
        rin = r1
        rout = r2
      endif

      thet = hkrand(0)*pi
      phi = hkrand(0)*2*pi
      xyz_in(1) = rin*sin(thet)*cos(phi)
      xyz_in(2) = rin*sin(thet)*sin(phi)
      xyz_in(3) = rin*cos(thet)

      thet = hkrand(0)*pi
      phi = hkrand(0)*2*pi
      xyz_out(1) = rout*sin(thet)*cos(phi)
      xyz_out(2) = rout*sin(thet)*sin(phi)
      xyz_out(3) = rout*cos(thet)


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

      call prin2('rhs=*',rhs,48)



      ntarg = 10000
      allocate(targs(3,ntarg),rhs_ex(nd,ntarg),rhs_interp(nd,ntarg))
      do i=1,ntarg
        tval = hkrand(0) 
        thet = hkrand(0)*2*pi
        call funcurve(tval,np,pars,rr,zz,tmp,tmp,tmp,
     1    tmp)
        targs(1,i) = rr*cos(thet)
        targs(2,i) = rr*sin(thet)
        targs(3,i) = zz
        call l3d_slp(xyz_in,3,targs(1,i),0,dpars,0,zpars,0,ipars,
     1     rhs_ex(1,i))
        call l3d_slp(xyz_out,3,targs(1,i),0,dpars,0,zpars,0,ipars,
     1     rhs_ex(2,i))
      enddo


      call axissym_fun_interp(nd,nch2d,tchse,k,funcurve,
     1  np,pars,rmax,iort,ntarg,targs,npatches,norders,ixyzs,iptype,
     2  npts,rhs,rhs_interp)

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
      call axissym_fun_local_coord_targ(nch2d,tchse,k,
     1  funcurve,
     1  np,pars,rmax,iort,ntarg,targs,npatches,norders,ixyzs,iptype,
     2  npts,ipatchtarg,uvs_targ)
      
cc      call prinf('ipatchtarg=*',ipatchtarg,ntarg)
cc      call prin2('uvs_targ=*',uvs_targ,2*ntarg)
      

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

      stop
      end


