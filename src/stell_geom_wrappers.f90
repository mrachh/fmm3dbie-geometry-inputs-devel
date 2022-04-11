

subroutine setup_stell_simsopt_geom(istelltype,norder,ns,nt,npatches, &
  norders,ixyzs,iptype,npts,srccoefs,srcvals)
!-----------------------
!  This subroutine discretizes one of the two stellarator types
!  w7x or ncsx in the fmm3dbie discretization format. The fourier coeffs
!  of the r and z components of the torus are preloaded, if you wish
!  to use a different stellarator design, use the guru version of the
!  discretizer setup_stell_geom_guru
!
!
!  For all of these geometries, the surface is defined via a parametrization
!  of the form
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
!
!    - istelltype: integer
!         type of stellrator
!         istelltype = 1, w7x
!         istelltype = 2, ncsx
!    - norder: integer
!         order of discretization for the patches
!    - ns: integer
!         number of patches to be used in the s variable above 
!    - nt: integer
!         number of patches to be used in the t variable above
!    - npatches: integer
!         must be equal to 2*ns*nt
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
  integer, intent(in) :: istelltype,norder,ns,nt,npatches,npts
  integer, intent(out) :: norders(npatches),ixyzs(npatches+1)
  integer, intent(out) :: iptype(npatches)
  real *8, intent(out) :: srccoefs(9,npts),srcvals(12,npts)

  real *8, allocatable :: deltas(:,:,:)
  integer m,nfp

  if(istelltype.eq.1) then
    m = 10
    nfp = 5
    allocate(delta(-m:m,-m:m,2))
    call load_w7x_coeffs(deltas)
  endif

  if(istelltype.eq.2) then
    nfp = 3
    m = 10
    allocate(delta(-m:m,-m:m,2))
    call load_ncsx_coeffs(deltas)
  endif


  call setup_stell_geom(norder,ns,nt,npatches,nfp,deltas,m,norders, &
    ixyzs,iptype,npts,srccoefs,srcvals)


return
end subroutine setup_stell_simsopt_geom
!
!
!
!
!
