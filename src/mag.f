c
c     this is C interface 
c
      subroutine magFunc(mass1, sep, xsc, ysc, radius, 
     $  Gam, amp, eps, errval)
      include 'global.h'
      integer errval
      integer i
      real*8 separation
      integer flag
      real*8 A, Acenter, A2rho2, A4rho4
      real*8 eps
      real*8 xsc, ysc, radius, Gam
      real*8 sep, mass1
      real*8 amp
      integer nphi

      logical firstTime
      data firstTime /.true./

      if (errval .ne. 0) then
         firstTime = .true.
      endif

      errval = 0

c     pass the information on global variables
      m1 = mass1
      separation = sep
      Gamma = Gam
      xsCenter = xsc
      ysCenter = ysc
      rs = radius


c
      m2 = 1.0d0 - m1
      z1 = -m1*separation
      z2 =  m2*separation
      debug = .false.

c     initialise the image positions
      if (firstTime) then
        nimage = 0
        do i=1, 5
           images(i) = dcmplx(1.0d-20, 1.0d-20)
        enddo
        firstTime = .false.
      endif

      nphi = 0
      call GouldApproximation(A, Acenter, A2rho2, A4rho4,errval)

c     if the approximation is bad, e.g., we are close to caustics, then use the image track method
      if (dabs(A4rho4) .gt. 1.0d-3) then
         A = Acenter 
         call uniformFunc(A, flag, nphi, eps, errval)
      endif

      amp = A


      if (Gamma .gt. 1.0d-10) then
         errval = 3001
c        stop 'oops, limb darkening profile not yet implemented (mag.f)'
         return
      endif

      return
      end
