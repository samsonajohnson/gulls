c
c     Gould approximation
c
      subroutine GouldApproximation(mag, A0, A2rho2, A4rho4,errval)
      include 'global.h'

      integer errval
      real*8 mag

      integer i
c      real*8 phi, dphi
      real*8 muRho(8)
      real*8 muRhoHalf(4)
      real*8 AhalfPlus, APlus, Across, A0
      real*8 A2rho2, A4rho4
      real*8 cosphi(8)
      data cosphi/1.0d0, 0.7071067811865476d0, 0.0d0, 
     $     -0.7071067811865475d0, -1.0d0, -0.7071067811865477d0, 
     $     0.0d0, 0.7071067811865474d0/
      real*8 sinphi(8)
      data sinphi/0.0d0, 0.7071067811865475d0, 1.0d0, 
     $     0.7071067811865476d0, 0.0d0, -0.7071067811865475d0, -1.0d0, 
     $     -0.7071067811865477d0/
      real*8 hcosphi(4)
      data hcosphi/0.5d0, 0.0d0, -0.5d0, 0.0d0/
      real*8 hsinphi(4)
      data hsinphi/0.0d0, 0.5d0, 0.0d0, -0.5d0/

      


c     central images
      xs = xsCenter
      ys = ysCenter

      call findImages(errval)
      A0 = muTotal

c      dphi=pi*0.25d0            !pi/4
      do i=1, 8
c         phi = dphi*dble(i-1)
c         xs = xsCenter + rs*dcos(phi)
         xs = xsCenter + rs*cosphi(i)
c         ys = ysCenter + rs*dsin(phi)
         ys = ysCenter + rs*sinphi(i)
         call findImages(errval)
         
         muRho(i) = muTotal
      enddo

c     four points at half radius
c      dphi=pi*0.5d0             !pi/2
      do i=1, 4
c         phi = dphi*dble(i-1)
c         xs = xsCenter + 0.5d0*rs*dcos(phi)
         xs = xsCenter + rs*hcosphi(i)
c         ys = ysCenter + 0.5d0*rs*dsin(phi)
         ys = ysCenter + rs*hsinphi(i)
         call findImages(errval)
         
         muRhoHalf(i) = muTotal
      enddo

      AhalfPlus = (muRhoHalf(1)+muRhoHalf(2)+
     $     muRhoHalf(3)+muRhoHalf(4))*0.25d0 - A0

      Aplus  = (muRho(1)+muRho(3)+muRho(5)+muRho(7))*0.25d0-A0
      Across = (muRho(2)+muRho(4)+muRho(6)+muRho(8))*0.25d0-A0

      A2rho2 = (16.0d0*AhalfPlus-Aplus)/3.0d0
      A4rho4 = (Aplus + Across)*0.5d0 - A2rho2

c      mag = A0 + A2rho2*0.5d0 * (1.0d0-0.2d0*Gamma) 
c     $     + A4rho4/3.0d0 * (1.0d0-11.0d0/35.0d0*Gamma)
c     No limb darkening version
      mag = A0 + A2rho2*0.5d0 + A4rho4/3.0d0

      return
      end
