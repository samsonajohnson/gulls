c
c     calculate the magnification for a limb-darkened star using qgaus.f from Numerical Recipies
c
      subroutine limb(mag, flag, nphi, eps)
      include 'global.h'

      integer j
      real*8 mag
      integer nphi
      integer flag
      real*8 eps

      REAL*8 dx,xm,xr,w(5),x(5)
      SAVE w,x
      real*8 sum
      real*8 left, right
      real*8 rs0
      real*8 u1, u2

      DATA w/.2955242247,.2692667193,.2190863625,.1494513491,
     *.0666713443/
      DATA x/.1488743389,.4333953941,.6794095682,.8650633666,
     *.9739065285/
      
      rs0 = rs

      xm=0.5 * (pi/2.0d0-0.0d0)
      xr=0.5 * (pi/2.0d0+0.0d0)

      sum = 0.0d0
      do j=1, 5
        dx=xr*x(j)

        rs = rs0 * dsin(xm - dx)
        left = mag
        call uniform(left, flag, nphi, eps)

        rs = rs0 * dsin(xm + dx)
        right = mag
        call uniform(right, flag, nphi, eps)

        sum =sum + w(j)*(left * dsin(xm-dx)**3
     $       + right * dsin(xm+dx)**3)
      enddo

      sum = xr*sum

      u1 = -Gamma
      u2 = 1.5d0*Gamma
      mag = ((1.0+u1)*mag + u2*sum )/(1.0d0+u1+2.0d0/3.0d0*u2)
      
      print *, 'sum = ', sum

      return
      end
