c
c     find the images for a binary lens using a combination of
c     Newton-Raphson method and polynomial solver of the 5-th order 
c     complex polynomial
c
c     the lens parameters passed through common blocks in global.h
c
      subroutine findImages(errval)
      include 'global.h'

      integer i, k
      integer nimage0

      complex*16 c(6)
      complex*16 roots(5)
      complex*16 rem, swap
      complex*16 Delta

      complex*16 dz, z
      real*8 jacobian
      real*8 uu,ths,rrs

      logical imageFound

      logical firstTime
      data firstTime /.true./
      integer errval

c     Deal with the single lens case first
      if (dabs(m2) .lt. 1.0d-8 .or. dabs(z2-z1) .lt. 1.0d-5) then
         nimage=2
         uu=dsqrt(xs**2 + ys**2)
         ths=datan2(ys,xs)
         rrs=0.5d0*(uu+dsqrt(uu**2+4.0d0))
         images(1)=dcmplx(rrs*dcos(ths),rrs*dsin(ths))
         rrs=0.5d0*(uu-dsqrt(uu**2+4.0d0))
         images(2)=dcmplx(-rrs*dcos(ths),-rrs*dsin(ths))
         muTotal=(uu**2+2.0d0)/(uu*(uu**2+4.0d0))
         return
      endif


c     check the coordinate system
      if (dabs(m1*z2+m2*z1) .gt. 1.0d-10) then
c         print *, 'coordinate system incorrect, m1 z2+m2 z1 !=0 (findImages.f)'
         errval=3111
c        stop
         return
      endif

c     first time, we are entering the findImage routine
      if (firstTime) then
         nimage = 0
         muTotal = 0.0d0
         do i=1, 5
            images(i) = dcmplx(0.0d0, 0.0d0)
         enddo

         call findCloseImages(errval)

         call polynomialSolver(errval)

         if (debug) then
            print *, 'entering findImage the first time',xs,ys
            print *, 'first time, ', nimage, ' images found'
            print *, 'muTotal = ', muTotal
         endif

         nimage0 = nimage
         firstTime = .false.

         call checkImageProperties(errval)

         return
      else
         nimage0 = nimage
      endif

      muTotal = 0.0d0           !total magnification is zero
      nimage = 0                    !reset the current image counter as 0

c     use previous image positions as initial guess in the
c     Newton-Raphson iteration
      do i=1, nimage0
         z = images(i)
         jacobian = mu(i)

         call newton(z, imageFound, jacobian)
         if (imageFound) call newImage(z, jacobian,errval)
      enddo

c     no image found!
      if (nimage .eq. 0) then
c         print *,
c     $   'Warning: no images found with initial guesses'
      endif
      
c     few than three images found, try to find the image close to
c     the source and two close to the lenses
      if (nimage .lt. 3) then
         call findCloseImages(errval)
      endif

      if (nimage .eq. 5) return     !all images found, we are done

c     either the true number of images is fewer than 5
c     or we are missing some real images -- deflate the polynomial
c     and find the solutions to the deflated polynomials and check
c     whether they are real images

      call coefficients(c)      !calculate the complex 5-th order
                                !complex polynomial coefficients

      do k=1, nimage                !deflate the polynomial, which will be 
                                !of order 5-ni
         rem = c(6-k+1)
         do i=6-k, 1, -1
            swap = c(i)
            c(i) = rem
            rem = swap + rem * images(k)
         enddo
      enddo
      
      if (nimage .eq. 3) then       !3 images found, solve the quadratic equation
         Delta = cdsqrt(c(2)**2-4.0d0*c(1)*c(3))
         roots(1) = (-c(2)+Delta)/(2.0d0*c(3))
         roots(2) = (-c(2)-Delta)/(2.0d0*c(3))
      else if (nimage .eq. 4) then  !4 images found, the deflated polynomial is linear
         roots(1) = -c(1)/c(2)
      else                      !missing 3 or more images, fall back to the polynomial solver
         do i=1, 5-nimage       !initialise possible solutions
            roots(i) = dcmplx(0.0d0, 0.0d0)
	 enddo

         call zroots(c, 5-nimage, roots, .true.)
      endif

c
c     check whether the newly found solutions to (the deflated polynomial) 
c     is a true solution to the lens equation
c
      do k=1, 5-nimage
         z = roots(k)
c     pre-copy the array
c         images(nimage+k) = roots(k)
         call realLensEquation(z, jacobian, dz)

         if (cdabs(dz) .lt. rootEPS) call newImage(z, jacobian,errval)
      enddo

c     we are desperate, try again using polynomial solver
      if (nimage .ne. 3 .and. nimage .ne. 5) then
         call polynomialSolver(errval)
      endif

c     check whether the image properties are OK
      call checkImageProperties(errval)

      return
      end

c     
c     use Newton-Ralphson iteration to find an image close to an initial guess
c     
      subroutine newton(z, imageFound, jacobian)
      include 'global.h'

      integer i
      logical imageFound
      complex*16 z, dz
      real*8 jacobian
      real*8 jacobian0
      complex*16 z0

      jacobian0 = jacobian
      
      z0 = z
      do i=1, 100       !maximum iteration is 30
         call realLensEquation(z, jacobian, dz)
         
         z = z + dz

         if (cdabs(dz) .lt. rootEPS) then
            imageFound = .true.

            return
         endif
c         continue
      enddo

c     no image found if reaching here
      imageFound = .false.

c     restore original guess
      z = z0
      jacobian = jacobian0

      return
      end
      
c     
c     check whether the found image is a new image. If yes
c     update the global arrays that hold the image information
c     such as images(5), mu(5) and nimage counter
c
      subroutine newImage(z, jacobian, errval)
      include 'global.h'

      integer errval
      complex*16 z
      real*8 jacobian

      integer i
      errval =0 

c     number of images should between 1 to 5
      if (nimage .lt. 0 .or. nimage .gt. 5) then
         print *, 'nimages = ', nimage
         errval=3113
c        stop 'image number incorrect in newImage'
         return
      endif

      do i=1, nimage
         if (cdabs(z-images(i)) .lt. 2.0d0*rootEPS) return
      enddo

      nimage = nimage + 1       !update the image counter
      images(nimage) = z
      mu(nimage) = jacobian
      muTotal = muTotal + dabs(jacobian)

      return
      end


c
c     check the global properties of all images found so far
c
      subroutine checkImageProperties(errval)
      include 'global.h'

c      integer k
      real*8 muTotal2
      integer errval

      errval = 0

      if (nimage .ne. 3 .and. nimage .ne. 5) then
c         print *, 'nimages incorrect, must be 3 or 5, but is ', nimage
c         print *, 'xs, ys=', xs, ys
c         print *, 'm1, m2=', m1, m2
c         print *, 'z1, z2=', z1, z2

c         do k=1, nimage
c            print *, k, images(k), mu(k)
c         enddo
         errval=3114
c         stop 'in checkImageProperties 888'
         return
      else                      !image number OK
         if (muTotal .lt. 1.0d0) then
c            print *, 'total magnification < 1, impossible'
            errval=3115
c            stop
            return
         else if (nimage .eq. 5) then !signed total must be one for 5-image systems
            muTotal2 = mu(1)+mu(2)+mu(3)+mu(4)+mu(5)
            if (dabs(muTotal2-1.0d0) .gt. 1.0d-4*muTotal) then
               print *, 'Warning: magnification relation inexact'
               print *, 'signed mu_tot & mu_tot:', muTotal2, muTotal
               print *, 'images: ', images(1)
               print *, 'images: ', images(2)
               print *, 'images: ', images(3)
               print *, 'images: ', images(4)
               print *, 'images: ', images(5)
               print *, 'mu: ', mu(1),mu(2),mu(3),mu(4),mu(5)
               print *, 'lens position and mass', z1, z2, m1, m2
               print *, 'source position', xs, ys
            endif
         endif
      endif
      
      return
      end


c
c     find default images: once close to the source and two close to
c     the lenses. These must be the images when the source is far away
c
      subroutine findCloseImages(errval)

      include 'global.h'
      integer i
      complex*16 zs, zsc
      complex*16 r(3), z
      logical imageFound
      real*8 jacobian
      integer errval

      zs = dcmplx(xs, ys)
      zsc = dcmplx(xs, -ys)
      errval = 0

      if (nimage .gt. 2) then
c        print *, 'more than 2 images already found'
c        print *, 'should not enter findCloseImages'
         errval=3112
c        stop
         return
      endif

      r(1) = zs - zsc/((zsc-z1+1.0d-20)*(zsc-z2+1.0d-20))
      r(2) = z1 + m1/(z1-zsc-m2/(z1-z2 ))
      r(3) = z2 + m2/(z2-zsc-m1/(z2-z1 ))

      do i=1, 3
         z = r(i)
         jacobian = 1.0d0       !assign an arbitrary parity

         call newton(z, imageFound, jacobian)
         if (imageFound) call newImage(z, jacobian,errval)
      enddo

      if (debug) then
         print *, 'entered findCloseImages'
      endif
         

      return
      end

c     
c     solve the lens equation using real notations
c     
      subroutine realLensEquation(z,jacobian,dz)
      include 'global.h'

      complex*16 z, dz
      real*8 jacobian
      real*8 x, y
      real*8 dx1, dx2
      real*8 dx1_2, dx2_2
      real*8 dr1_2, dr2_2, m1_dr1_4, m2_dr2_4
      real*8 dxs, dys
      real*8 dx, dy
      real*8 y2
      real*8 sum2, sq

      x = dreal(z)
      y = dimag(z)

      dx1 = x - z1
      dx2 = x - z2

      dx1_2 = dx1*dx1
      dx2_2 = dx2*dx2

      y2 = y*y

      dr1_2 = 1.0d0/(dx1_2 + y2)
      dr2_2 = 1.0d0/(dx2_2 + y2)

      m1_dr1_4 = dr1_2 * dr1_2 * m1
      m2_dr2_4 = dr2_2 * dr2_2 * m2

      dxs = xs - (x- m1*dx1*dr1_2 - m2*dx2*dr2_2)
      dys = ys - y*(1.0d0- m1*dr1_2 - m2*dr2_2)

      phixx=1.0d0-(y2-dx1_2)*m1_dr1_4 -(y2-dx2_2)*m2_dr2_4
c     there is no convergence, so phixx+phiyy = 2
      phiyy = 2.0d0- phixx

      phixy = 2.0d0*y*(dx1*m1_dr1_4 + dx2*m2_dr2_4)

      jacobian = 1.0d0/(phixx * phiyy - phixy*phixy)

      dx = jacobian*( phiyy * dxs - phixy * dys)
      dy = jacobian*(-phixy * dxs + phixx * dys)

c     figure out eigen values
c      sum2=(phixx+phiyy)/2.0d0
c      sq=dsqrt(sum2*sum2-1.0d0/jacobian)

c      lambda1 = sum2 + sq
c      lambda2 = sum2 - sq

c     rotation angle:
c      theta = 0.5d0 * datan(2.0d0*phixy/(phiyy-phixx)+1.d-50) !avoid zero

c      print *, 'theta =', x, y, theta*57.3, lambda1, lambda2
      dz = dcmplx(dx, dy)
c      stop
      return
      end


      subroutine polynomialSolver(errval)
      include 'global.h'

      integer errval
      integer i
      complex*16 c(6)
      complex*16 roots(5)
      complex*16 z, dz
      real*8 jacobian
      real*8 sum

      do i=1, nimage
         roots(i) = images(i)
      enddo
      
      do i=nimage+1, 5
         roots(i) = dcmplx(0.0d0, 0.0d0) !initialise root array to zero
      enddo
      
      call coefficients(c)      !calculate the complex 5-th order
      call zroots(c, 5, roots, .true.)
      
      do i=1, 5
         z = roots(i)
         call realLensEquation(z, jacobian, dz)
         
         if (cdabs(dz) .lt. rootEPS) call newImage(z, jacobian,errval)
      enddo

      if (nimage .eq. 5) then
         sum = mu(1)+mu(2)+mu(3)+mu(4)+mu(5)
c         print *, 'image number = ', nimage, sum-1.0
      endif

      return
      end
