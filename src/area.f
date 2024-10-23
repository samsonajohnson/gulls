c
c     calcualte the area contained within a track
c
      function areaFunc(zTrack, n)
      implicit none

      integer n
      complex*16 zTrack(n)
      real*8 areaFunc

      integer i
      real*8 x1, x2, y1, y2

      
c     determine center of light

c     calculate the area enclosed by the image track
      areaFunc = 0.0d0
      x1 = dreal(ztrack(1))
      y1 = dimag(ztrack(1))

      do i=2, n
         x2 = dreal(zTrack(i))
         y2 = dimag(zTrack(i))

         areaFunc = areaFunc +
     $        0.25d0*( (y2+y1)*(x2-x1)-(y2-y1)*(x2+x1))

         x1 = x2
         y1 = y2
      enddo

      return
      end



c
c     calcualte the area contained within a track using polar coordinates
c     centroiding may be incorrect due to non-even sampling on the circles
c
      function areaFuncPolar(zTrack, n)
      implicit none

      integer i
      integer n
      complex*16 zTrack(n)
      real*8 x1, x2, y1, y2
      real*8 r1, r2, phi1, phi2

      real*8 pi2
      parameter (pi2=2.0d0*3.14159265358979323846264338327950288d0)
      real*8 areaFuncPolar
      real*8 x0, y0

c     determine centre of light
      x0 = 0.0d0
      y0 = 0.0d0
      do i=1, n
         x0 = x0 + dreal(zTrack(i))
         y0 = y0 + dimag(zTrack(i))
      enddo
      x0 = x0/n
      y0 = y0/n

c     calculate the area enclosed by the image track
      areaFuncPolar = 0.0d0
      x1 = dreal(ztrack(1)) - x0
      y1 = dimag(ztrack(1)) - y0
      do i=2, n
         x2 = dreal(zTrack(i)) - x0
         y2 = dimag(zTrack(i)) - y0

         r1 = dsqrt(x1**2+y1**2)
         r2 = dsqrt(x2**2+y2**2)

         phi1 = datan2(y1, x1)
         phi2 = datan2(y2, x2)

         if (phi1 .lt. 0.0d0) phi1 = phi1 + pi2
         if (phi2 .lt. 0.0d0) phi2 = phi2 + pi2

         areaFuncPolar = areaFuncPolar + 0.125d0*(r1+r2)**2*(phi2-phi1)
c         areaFunc = areaFunc +
c     $        0.25*( (y2+y1)*(x2-x1)-(y2-y1)*(x2+x1))
         x1 = x2
         y1 = y2
      enddo

      print *, 'x0, y0, n = ', x0, y0, n, areaFuncPolar

      return
      end
