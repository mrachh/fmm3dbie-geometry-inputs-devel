      implicit real *8 (a-h,o-z)
      real *8, allocatable :: tchse(:)
      real *8, allocatable :: stargex(:),starg(:),xyztarg(:,:)
      real *8, allocatable :: restarg(:)
      real *8 pars(2)
      external funcurve_oocyte_riemann

      call prini(6,13)
      

      k = 16
      iref = 0
      nmid = 10
      nch2d = 2*(iref+1)+nmid
      allocate(tchse(nch2d+1))
      call get_oocyte3d_tchse(iref,nmid,nch2d,tchse)

      np = 2
      pars(1) = 0.25d0
      pars(2) = 0.1d0

      ntarg = 100
      allocate(stargex(ntarg),starg(ntarg),xyztarg(3,ntarg))
      allocate(restarg(ntarg))

      do i=1,ntarg
        stargex(i) = hkrand(0)
        call funcurve_oocyte_riemann(stargex(i),np,pars,r,z,tmp,tmp,
     1     tmp,tmp)
        thet = hkrand(0)
        xyztarg(1,i) = r*cos(thet)
        xyztarg(2,i) = r*sin(thet)
        xyztarg(3,i) = z
      enddo

      
      call get_param_vals(nch2d,tchse,k,funcurve_oocyte_riemann,np,
     1   pars,ntarg,xyztarg,starg,restarg)

      nprin = min(ntarg,12)
      call prin2('starg=*',starg,nprin)
      call prin2('stargex=*',stargex,nprin)
      call prin2('restarg=*',restarg,nprin)

      erra = 0
      ra = 0
      do i=1,ntarg
        ra = ra + abs(stargex(i))**2
        erra = erra + abs(starg(i)-stargex(i))**2
      enddo

      erra = sqrt(erra/ra)
      call prin2('err in svals =*',erra,1)

      stop
      end



