c
c     this subroutine calculates the coeffecients of the 5-th order
c     complex lens equation of a binary lens
c     input: 
c	source position (xs, ys), passed in global.h via a common block
c     output:
c     	coefficients c(6)
c
c     see mag.f for coordinate systems and conventions
c
      subroutine coefficients(c)
      include 'global.h'

c     source position in complex*16 notation
      complex*16 zs, zsc
      complex*16 zs1, zs2, zsm, zsa

c     local variables
      real*8 a, b
      real*8 sq(5)
      complex*16 tmp(5)
      complex*16 c(6)

c     source position in complex notation
      zs = dcmplx(xs, ys)
      zsc = dcmplx(xs, -ys)

c     coefficients for (z-z1)(z-z2)=z^2+a z + b
      a = -(z1+z2)
      b = z1 * z2

c     coefficients for sq == (z-z1)^2(z-z2)^2
      sq(5) = 1.0d0
      sq(4) = a+a
      sq(3) = b+b+a*a
      sq(2) = 2.0d0*a*b
      sq(1) = b*b

      zs1 = zsc-z1
      zs2 = zsc-z2
      zsa = zs1+zs2
      zsm = zs1*zs2

c     tmp = (zsc-z2)(zsc-z1) sq + (2zsc-z1-z2)(z-z1)(z-z2) z +z^2
      tmp(5) = zsm*sq(5)
      tmp(4) = zsm*sq(4) + zsa

      tmp(3) = zsm*sq(3) + zsa*a + 1.0d0
      tmp(2) = zsm*sq(2) + zsa*b
      tmp(1) = zsm*sq(1)

c     (z-zs) tmp - zsc(z-z1)^2(z-z2)^2 - z(z-z1)(z-z2)
      c(6) = tmp(5)
      c(5) = tmp(4) - zs*tmp(5) - zsc*sq(5)
      c(4) = tmp(3) - zs*tmp(4) - zsc*sq(4) - 1.0d0
      c(3) = tmp(2) - zs*tmp(3) - zsc*sq(3) - a
      c(2) = tmp(1) - zs*tmp(2) - zsc*sq(2) - b
      c(1) = -zs*tmp(1) - zsc*sq(1)

      return
      end
