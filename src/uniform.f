c
c     we first try to find the magnification for a uniform sphere with Stoke's theorem.
c
c     Sometimes it fails when the limb of the source is very close to two converging fold caustics
c     or a cusp. If this happens, we make the source first slightly larger and then smaller, obtain 
c     their magnifications, and then do a linear interpolation
c
      subroutine uniformFunc(mag, flag, nphi, eps,errval)
      integer errval
      real*8 mag
      integer nphi
      integer flag
      real*8 eps

      errval = 0

      call uniform(mag, flag, nphi, eps, errval)


c     failed using Stokes's theorem, use interpolation
      if (flag .lt. 1) then
         print *, 'using interpolation'
         call uniformInterpolate(mag, flag, nphi, eps,errval)
         if (flag .lt. 1) then 
            errval=3202
c           stop 'failed in uniformInterpolate'
            return
         endif
      endif

      return
      end

c
c     calcualte the magnification for a uniform extended source
c
      subroutine uniform(mag, flag, nphi, eps, errval)
      include 'global.h'

      integer errval
      real*8 mag
      integer nphi
      integer flag
      real*8 eps
      integer i, j, k

      complex*16 ztrack(nmax, 5)
c      complex*16 complexTrack(nmax, 5)
      integer npTrack(5)

      real*8 areaFunc, areaFuncPolar
      external areaFunc, areaFuncPolar
      real*8 sum, mag0
      integer ntrack

      real*8 area(5)
      integer sign(5)           !parity for image tracks
      integer inside, ccw
      external inside, ccw
      integer index(5)

      nphi = 2**6 
      mag0 = mag

      errval = 0

c     this is a point source
      if (rs .lt. 1.0d-10) return

      flag = 0
      do i=1, bisectMax-6
         call findTrack(nphi, ntrack, ztrack, npTrack, flag, errval)
         if (errval .gt. 0) then
c     stop
            return
         endif
         

c     we failed, double number of points, try again
         if (flag .lt. 0) then
            nphi = nphi*2 - 1
            cycle
         endif

c     calculate areas
         sum = 0.0d0
         do k=1, ntrack
            area(k) = dabs(areaFunc(zTrack(1, k), npTrack(k)))
            sum = sum + area(k)
            sign(k) = +1
         enddo


c     the following seeming confusing code deals with nested image tracks
         if (ntrack .eq. 5) then     !five tracks, they cannot be nested
            mag = sum/(pi*rs*rs)
            
         else
c     obtain the index table in the rising magnification order
            call indexx(ntrack, area, index)
c     check whether the smaller tracks are inside the larger tracks
            sign(index(ntrack)) = +1
            do k=ntrack, 2, -1
               do j=1, k-1
                  if (inside(ztrack(1, index(j)), 
     $                 ztrack(1, index(k)), npTrack(index(k)))
     $                 .eq. 1) then
                     sign(index(j)) = -sign(index(j))
                     if (debug) then
                        print *, 'track ', j, ' is inside track ',k
                     endif
                  endif
               enddo
            enddo

            sum = 0.0d0
            do k=1, ntrack
               sum = sum + sign(index(k)) * area(index(k))
            enddo
         endif

         mag = sum/(pi*rs*rs)

         if (dabs(mag-mag0) .lt. eps*mag) then !we have reached the accuracy we desired, return
            if (debug) print *, 'done mag0, mag=', mag0, mag
            flag = 1

            return
         else                   !not yet done, double points and try again
            if (debug) print *, 'not yet done -- mag0, mag=', mag0, mag
            mag0 = mag          
            nphi = 2*nphi - 1
         endif

         if (nphi .gt. nmax) then
            print *, 'trying to use too many points in uniform'
            flag = -1

            return
         endif

      enddo

c     if we reach here, we have failed in our tricks
      flag = -1

      return
      end

c
c     The source may have encountered a narrow caustic crossing, we apply linear interpolation
c
      subroutine uniformInterpolate(mag, flag, nphi, eps, errval)

      include 'global.h'
      integer errval
      real*8 mag
      integer nphi
      integer flag
      real*8 eps

      real*8 rs0, A1, A2
      real*8 fraction

      fraction = 0.02d0
      rs0 = rs

      errval = 0
c     decrease the source radius by fraction
      rs = (1.0d0-fraction) * rs0
      A1 = mag
      call uniform(A1, flag, nphi, eps, errval)

      if (flag .lt. 1)then
         errval=3201
c        stop 'linear interpolation failed'
         return
      endif

 

c     increase the source radius by fraction
      rs = (1.0d0+fraction) * rs0
      A2 = mag
      call uniform(A2, flag, nphi, eps, errval)

      if (flag .lt. 1) then
         errval=3201
c        stop 'linear interpolation failed'
         return
      endif

    

c     reset the source radius
      rs = rs0

c     take the magnification as the average
      mag = (A1+A2)*0.5d0
      
      return
      end
