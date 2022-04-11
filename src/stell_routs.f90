!
!
!
!
!


subroutine setup_stell_geom(norder,ns,nt,npatches,nfp,deltas,m,norders, &
  ixyzs,iptype,npts,srccoefs,srcvals)
!-----------------------
!  This subroutine discretizes stellarator like geometries which
!  are maps from [0,2\pi]^2 of the form
!
!  R(s,t) = \sum_{i=-m}^{m} \sum_{j=-m,m} deltas(i,j,1) cos(j*s - nfp*i*t)
!  Z(s,t) = \sum_{i=-m}^{m} \sum_{j=-m,m} deltas(i,j,2) cos(j*s - nfp*i*t)
!
!  x(s,t) = R(s,t)*cos(t)
!  y(s,t) = R(s,t)*sin(t)
!  z(s,t) = Z(s,t)
!  
!
!  Input arguments:
!    - norder: integer
!         order of discretization for the patches
!    - ns: integer
!         number of patches to be used in the s variable above 
!    - nt: integer
!         number of patches to be used in the t variable above
!    - npatches: integer
!         must be equal to 2*ns*nt
!    - nfp: integer
!         nfp above
!    - deltas: real *8 (-m:m,-m:m,2)
!         the coeffs deltas above
!    - m: integer
!         number of coeffs in deltas above
!    - npts: integer
!         must be equal to npatches*(norder+1)*(norder+2)/2
!
!
!  Output arguments:
!    - norders: integer(npatches)
!        discretization order of patches
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!        
!----------------------------
 
  implicit none
  integer, intent(in) :: norder,ns,nt,npatches,npts
  real *8, intent(in), target :: deltas(-m:m,-m:m,2)
  integer, intent(in), target :: m,nfp
  integer, intent(out) :: norders(npatches),ixyzs(npatches+1)
  integer, intent(out) :: iptype(npatches)
  real *8, intent(out) :: srccoefs(9,npts),srcvals(12,npts)


  real *8, allocatable, target :: triaskel(:,:,:)
  integer, pointer :: iptr1,iptr2,iptr3,iptr4
  real *8, pointer :: ptr1,ptr2,ptr3,ptr4
  real *8 smin,smax,tmin,tmax,done,pi
  integer nover,npols
  real *8, allocatable :: uvs(:,:),wts(:),umatr(:,:),vmatr(:,:)

  external xtri_biest_stell_eval

  done = 1.0d0
  pi = atan(done)*4
  
  smin = 0
  smax = 2*pi
  tmin = 0
  tmax = 2*pi

  allocate(triaskel(3,3,npatches))
  nover = 0
  call xtri_rectmesh_ani(umin,umax,vmin,vmax,ns,nt,nover,npatches, &
    npatches,triaskel)

  ptr1 => triaskel(1,1,1)
  ptr2 => deltas(-m,-m,1)
  iptr3 => m
  iptr4=>nfp

  npols = (norder+1)*(norder+2)/2
  allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
  allocate(wts(npols))

  call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)


  call getgeominfo(npatches,xtri_biest_stell_eval,ptr1,ptr2,iptr3, &
    iptr4,npols,uvs,umatr,srcvals,srccoefs) 

  do i=1,npatches
    iptype(i) = 1
    norders(i) = norder
    ixyzs(i) = (i-1)*npols + 1
  enddo
  ixyzs(npatches+1) = npts+1


return
end subroutine setup_stell_geom
!
!
!
!
!
subroutine setup_stell_geom_adap(norder,ns0,nt0,npatchesmax,nptsmax,nfp, &
  deltas,m,npatches,norders,ixyzs,iptype,npts,srccoefs,srcvals)
!-----------------------
!  This subroutine discretizes stellarator like geometries which
!  are maps from [0,2\pi]^2 of the form
!
!  R(s,t) = \sum_{i=-m}^{m} \sum_{j=-m,m} deltas(i,j,1) cos(j*s - nfp*i*t)
!  Z(s,t) = \sum_{i=-m}^{m} \sum_{j=-m,m} deltas(i,j,2) cos(j*s - nfp*i*t)
!
!  x(s,t) = R(s,t)*cos(t)
!  y(s,t) = R(s,t)*sin(t)
!  z(s,t) = Z(s,t)
!
!  in an adaptive manner
!  
!
!  Input arguments:
!    - norder: integer
!         order of discretization for the patches
!    - ns0: integer
!         initial number of patches to be used in the s variable above 
!    - nt0: integer
!         initial number of patches to be used in the t variable above
!    - npatchesmax: integer
!         max number of patches
!    - nptsmax: integer
!         max number of points
!    - nfp: integer
!         nfp above
!    - deltas: real *8 (-m:m,-m:m,2)
!         the coeffs deltas above
!    - m: integer
!         number of coeffs in deltas above
!
!
!  Output arguments:
!    - npatches: integer
!        numbers of patches
!    - norders: integer(npatches)
!        discretization order of patches
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!        
!----------------------------

  implicit none
  integer, intent(in) :: norder,ns0,nt0,npatchesmax,nptsmax
  real *8, intent(in), target :: deltas(-m:m,-m:m,2)
  integer, intent(in), target :: m,nfp
  integer, intent(out) :: norders(npatchesmax),ixyzs(npatchesmax+1)
  integer, intent(out) :: iptype(npatchesmax)
  real *8, intent(out) :: srccoefs(9,nptsmax),srcvals(12,nptsmax)


  real *8, allocatable, target :: triaskel(:,:,:)
  integer, pointer :: iptr1,iptr2,iptr3,iptr4
  real *8, pointer :: ptr1,ptr2,ptr3,ptr4
  real *8 smin,smax,tmin,tmax,done,pi
  integer nover,npols
  real *8, allocatable :: uvs(:,:),wts(:),umatr(:,:),vmatr(:,:)

  external xtri_biest_stell_eval

  done = 1.0d0
  pi = atan(done)*4
  
  smin = 0
  smax = 2*pi
  tmin = 0
  tmax = 2*pi

  allocate(triaskel(3,3,npatchesmax))
  nover = 0
  call xtri_rectmesh_ani(umin,umax,vmin,vmax,ns0,nt0,nover,npatches, &
    npatches,triaskel)

  ptr1 => triaskel(1,1,1)
  ptr2 => deltas(-m,-m,1)
  iptr3 => m
  iptr4=>nfp

  npols = (norder+1)*(norder+2)/2
  allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
  allocate(wts(npols))

  call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)

!
!  get initial collection of panels
! 
  call getgeominfo(npatches,xtri_biest_stell_eval,ptr1,ptr2,iptr3, &
    iptr4,npols,uvs,umatr,srcvals,srccoefs) 

  do i=1,npatches
    iptype(i) = 1
    norders(i) = norder
    ixyzs(i) = (i-1)*npols + 1
  enddo
  ixyzs(npatches+1) = npts+1


  nrefmax = 52

  do ilev=1,nrefmax
!
!
!  routine to identify whether a patch needs refinement
!  relative error to surface area -> surf_routs
!
!  routine to refine skeleton mesh and identify triangle numbers -> xtri_routs
!
!  routine to rewrite srcvals and other arrays -> Here
!
  enddo

return
end subroutine setup_stell_geom_adap


!
!
!
!
!
subroutine xtri_biest_stell_eval(itri, u, v, xyz, dxyzduv, triainfo, &
    deltas, m, nfp)
  implicit real *8 (a-h,o-z)
  integer nfp,m
  real *8 :: xyz(3), dxyzduv(3,2), triainfo(3,3,*), deltas(-m:m,-m:m,2)
  real *8 :: dxyzds(3),dxyzdt(3)

  !
  ! project the triangle itri in triainfo onto a stellarator
  !
  !    Input:
  ! itri - triangle number to map
  ! u,v - local uv coordinates on triangle itri
  ! triainfo - flat skeleton triangle info
  !
  !    Output:
  ! xyz - point on the sphere
  ! dxyzduv - first and second derivative information
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
  ! ... process the geometry, return the point location on the almond
  ! and the derivatives with respect to u and v
  !
  s = x0+u*(x1-x0)+v*(x2-x0)
  t = y0+u*(y1-y0)+v*(y2-y0)

  xyz(1) = 0
  xyz(2) = 0
  xyz(3) = 0

  dxyzds(1) = 0
  dxyzds(2) = 0
  dxyzds(3) = 0

  dxyzdt(1) = 0
  dxyzdt(2) = 0
  dxyzdt(3) = 0
  !call prin2('deltas = *', deltas, (m+2)*(n+2))
  !stop
  !
  ! In BIEST notation, t is theta and s is phi

  ct = cos(t)
  st = sin(t) 
  do i = -m,m
    do j = -m,m

      cst = cos(-nfp*i*t + j*s)
      sst = sin(-nfp*i*t + j*s)
      xyz(1) = xyz(1) + ct*deltas(i,j,1)*cst
      xyz(2) = xyz(2) + st*deltas(i,j,1)*cst
      xyz(3) = xyz(3) + deltas(i,j,2)*sst


      dxyzds(1) = dxyzds(1) - (j+0.0d0)*ct*deltas(i,j,1)*sst
      dxyzds(2) = dxyzds(2) - (j+0.0d0)*st*deltas(i,j,1)*sst
      dxyzds(3) = dxyzds(3) + (j+0.0d0)*deltas(i,j,2)*cst

      dxyzdt(1) = dxyzdt(1) + deltas(i,j)*(-st*cst -sst*ct*nfp*i)
      dxyzdt(2) = dxyzdt(2) + deltas(i,j)*(ct*cst - sst*st*nfp*i)
      dxyzdt(3) = dxyzdt(3) + deltas(i,j)*cst*nfp*i

    end do
  end do


  dsdu = (x1-x0)
  dsdv = (x2-x0)
  dtdu = (y1-y0)
  dtdv = (y2-y0)

  dxyzduv(1,1) = dxyzds(1)*dsdu + dxyzdt(1)*dtdu
  dxyzduv(2,1) = dxyzds(2)*dsdu + dxyzdt(2)*dtdu
  dxyzduv(3,1) = dxyzds(3)*dsdu + dxyzdt(3)*dtdu


  dxyzduv(1,2) = (dxyzds(1)*dsdv + dxyzdt(1)*dtdv)
  dxyzduv(2,2) = (dxyzds(2)*dsdv + dxyzdt(2)*dtdv)
  dxyzduv(3,2) = (dxyzds(3)*dsdv + dxyzdt(3)*dtdv)

  return
end subroutine xtri_biest_stell_eval










