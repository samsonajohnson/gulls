c
c     routine to check whether a point is inside a polygon, given
c     by a complex array z(n). We use the algorithm of Sedgewick:
c     Algorithm in C, Addison Wesley (1990), chapter 24
c
      integer function inside(z0, z, n)
      implicit none

      integer intercept
      external intercept

      complex*16 z0
      integer n
      complex*16 z(n)
      complex*16 z1, z2, zmax

      integer i, count

c     zmax is an arbitrary very large complex number
      zmax = dcmplx(100.0d0, 1.0d0)

      count = 0

      do i=1, n-1
         z1 = z(i)
         z2 = z(i+1)

         if (intercept(z1, z2, z0, zmax) .eq. 1) then
           count = count + 1
         endif
      enddo

c      print *, 'count=', count

      if (mod(count, 2) .eq. 1) then
         inside = + 1
      else
         inside = -1
      endif

      return
      end

c
c     does two lines intercept? the first line is given by z1 and z2,
c     the second line is given by z3 and z4
c
      integer function intercept(z1, z2, z3, z4)
      implicit none
      complex*16 z1, z2, z3, z4

      integer ccw
      external ccw

      if (ccw(z1, z2, z3)*ccw(z1, z2, z4) .le. 0 .and.
     $     ccw(z3, z4, z1)*ccw(z3, z4, z2) .le. 0) then
         intercept = 1
      else
         intercept = -1
      endif

      return
      end

c
c     check whether three points form a counter-clock-wise triangle
c
      integer function ccw(z0, z1, z2)
      implicit none

      complex*16 z0, z1, z2
      real*8 dx1, dx2, dy1, dy2

      dx1=dreal(z1-z0)
      dy1=dimag(z1-z0)
      dx2=dreal(z2-z0)
      dy2=dimag(z2-z0)
c     print *, 'z0,z1,z2', z0,z1,z2


      if (dx1*dy2 .gt. dy1*dx2) then
         ccw = +1
         return
      else
         ccw = -1
         return
      endif

      if (dx1*dx2 .lt. 0.0 .or. dy1*dy2 .lt. 0.0) then
         ccw = -1
         return
      endif

      if (dx1*dx1+dy1*dy1 .lt. dx2*dx2+dy2*dy2) then
         ccw=+1
         return
      endif

c     this should never happen
      ccw = 0
      print *, 'this should never happen in ccw unless linear'

      return
      end
