c
c     find image tracks
c
      subroutine findTrack(nphi, ntrack, ztrack,
     $           npTrack, flag, errval)
      include 'global.h'
      integer errval

      integer i, k
      integer nphi
      integer ntrack            !number of complex tracks, <=5
      complex*16 ztrack(nmax, 5) !this is the track for all 5 complex solutions
      integer npTrack(5)        !number of points in tracks
      integer flag
      complex*16 complexTrack(nmax, 5)
      complex*16 segments(nmax, nsegmentsMax) !
      integer npSegment(nsegmentsMax)
      integer nSegments


      flag = 1
      errval =0
    
c     first find the complex tracks adaptively 
      call adaptiveFindComplexTracks(nphi, complexTrack, flag, errval)
      if (flag .lt. 0) return

      if (errval .gt. 0) then
c      stop
      return
      endif

c     find image segments
      call findSegments(nphi, complexTrack, flag,
     $     segments, npSegment, nSegments)

      if (flag .lt. 0) return

c     connect segments
      call connectSegments(segments, npSegment, nSegments, flag, errval)
      if (flag .lt. 0) return

c     collect the connected tracks into the array zTrack
      ntrack = 0
      do i=1, nSegments
         if (npSegment(i) .gt. 0) then

            ntrack = ntrack + 1

            do k=1, npSegment(i)
               ztrack(k, ntrack) = segments(k, i)
            enddo
            npTrack(ntrack) = npSegment(i)
         endif
      enddo

c     if debugging, output the image tracks
      if (debug) then
         print *, 'ntrack = ', ntrack
         call outputTrack(ztrack, npTrack, ntrack)
      endif

      return
      end

c     
c     find the image tracks for all complex solutions 
c     to the fifth order polynomial
c     notice not all the solutions are real images. This subroutine uses 
c     fix step size
c
      subroutine findComplexTracks(nphi, ztrack, flag, errval)
      include 'global.h'

      integer errval
      integer nphi
      complex*16 ztrack(nmax, 5) !this is the track for all 5 complex solutions
c      real*8 jacob(nmax, 5)
      integer flag
      real*8 Jacobian(5)

      integer i, k
      real*8 dphi, phi

      complex*16 roots(5)

      dphi = 8.0d0*datan(1.0d0)/dble(nphi-1)

c     start at positive x
      xs = xsCenter + rs
      ys = ysCenter

      nimage = 0
      do i=1, 5
         roots(i) = dcmplx(0.0d0, 0.0d0)
      enddo

      errval=0
      call complexSolve (roots, Jacobian, errval)
      if(errval .gt. 0) return

      do i=1, 5
         ztrack(1, i) = roots(i)
c         jacob(1, i) = Jacobian(i)
      enddo

      flag = 1
      if (debug) open(33, file='source.dat', status="old", 
     $     position="append", action="write")
      if (debug) open(34, file='images.dat', status="old", 
     $     position="append", action="write")

      do k=2, nphi
	 phi = dble(k-1)*dphi
         xs = xsCenter + rs*dcos(phi)
         ys = ysCenter + rs*dsin(phi)

         if (debug) write(33, *) "findtrack ", xs, ys

         call complexSolve(roots, jacobian, errval)
         if(errval .gt. 0) return

         do i=1, 5
            ztrack(k, i) = roots(i)
            if (debug) write(34, *) "findtrack ", xs, ys
         enddo
      enddo
      if (debug) close(33)

      return
      end


c     
c     find true image segments from the complex tracks found above
c
      subroutine findSegments(nphi, ztrack, flag,
     $     segments, npSegment, nSegments)

      include 'global.h'

      integer nphi
      complex*16 ztrack(nmax, 5) !this is the track for all 5 complex solutions
c      real*8 jacob(nmax, 5)
      integer flag

      integer i, k

      complex*16 zs
      real*8 jacobian
      real*8 radius2
      
      integer np
      complex*16 segments(nmax, nsegmentsMax) !
      logical newSegment
      integer npSegment(nsegmentsMax)
      integer nSegments

      flag = 1
      nsegments = 0
      do i=1, 5
         newSegment = .true.
         np = 0

         do k=1, nphi
            call complexLensEquation(ztrack(k, i), zs, jacobian)
            radius2 = (dreal(zs)-xsCenter)**2+(dimag(zs)-ysCenter)**2

            if (dabs(dsqrt(radius2)/rs-1.0d0) .lt. 1.0d-8) then
               np = np + 1
               if (newSegment) nsegments = nsegments + 1

               if (nsegments .gt. nsegmentsMax) then
                  print *, 'too many segments'

                  flag = -1
                  return
               endif

               if (np .gt. nmax) then
                  print *, 'too many points in segment'

                  flag = -1
                  return
               endif

               segments(np, nsegments) = ztrack(k, i)
               newSegment = .false. !this segment is now old


               npSegment(nsegments) = np
            else
               newSegment = .true.
               np = 0
            endif
         enddo
      enddo

      if (debug) print *, 'number of segments: ', nsegments

      return
      end


c     complex lens equation. For a given src position
c     find the corresponding source position and Jacobian

      subroutine complexLensEquation(z, zs, jacobian)
      include 'global.h'
      
      complex*16 z, zs
      real*8 jacobian

      zs = z - m1/dconjg(z-z1) - m2/dconjg(z-z2)
      jacobian = 1.0d0-cdabs(m1/(z-z1)**2+m2/(z-z2)**2)**2

      return
      end


c     
c     output segments
c     
      subroutine outputSegments(segments, npSegment, nsegments)
      include 'global.h'

      integer i, k
      integer nsegments
      complex*16 segments(nmax, nsegmentsMax)
      integer npSegment(nsegmentsMax)
      
      open(11, file='segments.dat')
      write(11, *) nsegments
      do i=1, nsegments
         write(11, *) npSegment(i)
      enddo

      do i=1, nsegments
         do k=1, npSegment(i)
            write(11, *) dreal(segments(k, i)), dimag(segments(k,i))
         enddo
      enddo
      close(11)

      do i=1, nsegments
         do k=1, npSegment(i)
            write(30+i, *) dreal(segments(k, i)), dimag(segments(k, i))
         enddo
      enddo

      return
      end

c
c     connect different segments
c
      subroutine connectSegments(segments, npSegment, 
     $            nSegments, flag, errval)
      include 'global.h'

      integer i, j
      complex*16 segments(nmax, nsegmentsMax)
      integer npSegment(nsegmentsMax)
      integer nSegments
      integer ntrack
      logical selfClosed(nsegmentsMax)
      integer flag
      integer errval
      flag = 1

      errval = 0

c     find self-closed segments
      do i=1, nSegments
         selfClosed(i) = .false.

c     we do not close single point track which may happen the source
c     crosses a very caustic
         if (npSegment(i) .eq. 1) then
            if (debug)print *, 'single point track ', i, segments(1, i)
            cycle
         endif

         if (cdabs(segments(1, i)-segments(npSegment(i), i)) .lt.
     $        rootEPS) then

            if (debug) print *, 'segment ', i, ' is self closed'

            selfClosed(i) = .true.
         endif
      enddo


      do i=nSegments, 1, -1
         if (selfClosed(i)) cycle !already closed, do nothing

         do j=1, nSegments
            if (selfClosed(j)) cycle !already closed do nothing

            if (i .ne. j .and. npSegment(j) .gt. 0) then
               call connectEnds(segments, npSegment,
     $              nSegments, i, j, errval)
            endif
         enddo
      enddo

      if (debug) then
         call outputSegments(segments, npSegment, nsegments)
      endif


c
c     now connect ends that require us to jump over caustics
c
      do i=nSegments, 1, -1
         if (selfClosed(i)) cycle

         do j=1, nSegments
            if (selfClosed(j)) cycle

            if (i .ne. j .and. npSegment(j) .gt. 0) then
               call jumpConnects(segments, npSegment,
     $              nSegments, i, j, errval)
            endif
         enddo
      enddo


c     close the tracks
      ntrack = 0
      do i=1, nSegments
         if (npSegment(i) .ne. 0) then
            ntrack = ntrack + 1
            if (ntrack .gt. 5) then
               print *, 'too many tracks', ntrack

               flag = -1
               return
            endif

            if (cdabs(segments(1, i)-segments(npSegment(i), i))
     $           .gt. rootEPS) then
               if (npSegment(i)+1 .gt. nmax) then
                 errval=3216
c                 stop 'too many points+1'
                 return
               endif

               segments(npSegment(i)+1, i) = segments(1, i)
               npSegment(i) = npSegment(i) + 1

            endif
         else
            !do nothing, already closed
         endif
      enddo

      return
      end


c
c     connect segments that have mutual ends
c
      subroutine connectEnds(segments, npSegment, nSegments, 
     $           i, j, errval)

      include 'global.h'

      complex*16 segments(nmax, nsegmentsMax)
      integer npSegment(nsegmentsMax)
      integer nSegments
      integer i, j, k
      integer errval
      complex*16 zi1, zi2
      complex*16 zj1, zj2      
      errval = 0

      k=nSegments !do nothing line to stop a warning

c     no segments to connect, return
      if (npSegment(i) .eq. 0 .or. npSegment(j) .eq. 0) return

c
      zi1 = segments(1, i)      
      zi2 = segments(npSegment(i), i)      

      zj1 = segments(1, j)
      zj2 = segments(npSegment(j), j)

c     ends of i meet with beginning of j
      if (cdabs(zi2 - zj1) .lt. rootEPS) then
         if (npSegment(i)+npSegment(j)-1 .gt. nmax) then
	    print *, npSegment(i)+npSegment(j)-1
           errval=3211
c           stop 'too many points in i & j'
           return
         endif

c     notice the first, identical element in j is not counted!
         do k=2, npSegment(j) 
            segments(npSegment(i)+k-1, i) = 
     $           segments(k, j)
         enddo

         if (debug) print *, 'connecting ends of ', i, j,
     $        npSegment(i), npSegment(j)

         npSegment(i) = npSegment(i)+npSegment(j)-1
         npSegment(j) = 0      !set segment j to zero

c     ends of i meet with end of j, last element of segment in j is not counted
      else if (cdabs(zi2 - zj2) .lt. rootEPS ) then

         if (npSegment(i)+npSegment(j)-1 .gt. nmax) then
c	    print *, npSegment(i)+npSegment(j)-1
            errval=3211
c     stop 'too many points in i & j'
            return
         endif

         do k=1, npSegment(j)-1
            segments(npSegment(i)+k, i) = 
     $           segments(npSegment(j)-k, j)
         enddo

         if (debug) print *, 'connecting ends of ', i, j,
     $        npSegment(i), npSegment(j)

         npSegment(i) = npSegment(i)+npSegment(j)-1
         npSegment(j) = 0      !set segment j to zero
      endif

      return
      end


c     
c     find the image tracks for all complex solutions to the fifth order polynomial
c     notice not all the solutions are real images. This subroutine uses adaptive step
c
      subroutine adaptiveFindcomplextracks(nphi, ztrack, flag, errval)
      include 'global.h'

      integer nphi
      integer nsub 
      parameter (nsub = 10)
      complex*16 ztrack(nmax, 5) !this is the track for all 5 complex solutions
c      real*8 jacob(nmax, 5)

      integer flag
      integer i, j, k
      real*8 dphi, phi
      integer errval
      complex*16 roots(5)
      real*8 mu0
      integer count
      real*8 Jacobian(5)
      real*8 muTotal0

c     initialise
      nimage = 0
      do i=1, 5
         roots(i) = dcmplx(0.0d0, 0.0d0)
      enddo

c     start at positive x
      xs = xsCenter + rs
      ys = ysCenter

      errval = 0
      call complexSolve(roots, Jacobian, errval)
      if(errval .gt. 0) return

      count = 1
      do i=1, 5
         ztrack(1, i) = roots(i)
c         jacob(1, i) = Jacobian(i)
      enddo

      mu0 = muTotal
c
      flag = 1
      if (debug) open(33, file='source.dat')

      dphi = 8.0d0*datan(1.0d0)/dble(nphi-1)
      do k=2, nphi
         count = count + 1

	 phi = dble(k-1)*dphi
         xs = xsCenter + rs*dcos(phi)
         ys = ysCenter + rs*dsin(phi)

         if (debug) write(33, *) xs, ys

c     in this step we also return muTotal
         call complexSolve(roots, Jacobian, errval)
         if (errval .gt. 0) return

         muTotal0 = muTotal

         do i=1, 5
            ztrack(count, i) = roots(i)
c            jacob(count, i) = Jacobian(i)
         enddo


c     magnification gradient is too large, further divide into nsub steps
         if (muTotal/mu0.gt.1.1d0 .or. mu0/muTotal.gt.1.1d0) then

c     take one step back, and subdivide
            phi = phi - dphi
            count = count - 1

            do j=1, nsub
               count = count + 1

               phi = phi + dphi/dble(nsub)
               xs = xsCenter + rs*dcos(phi)
               ys = ysCenter + rs*dsin(phi)            

               call complexSolve(roots, Jacobian, errval)
               if(errval .gt. 0) return

               if (count .gt. nmax) then
c     stop 'too many points in adaptive'
                  errval=3212
                  return
               endif

               do i=1, 5
                  ztrack(count, i) = roots(i)
c                  jacob(count, i) = Jacobian(i)
               enddo
            enddo
         endif

         mu0 = muTotal
      enddo
      nphi = count

      if (debug) close(33)

      return
      end

c
c     jump over caustics and connect different segments
c
      subroutine jumpConnects(segments, npSegment,
     $           nSegments, i, j, errval)
      include 'global.h'

      complex*16 segments(nmax, nsegmentsMax)
      integer npSegment(nsegmentsMax)
      integer nSegments
      integer i, j, k
      integer errval

      complex*16 zs_i2, zi1, zi2
      complex*16 zs_j1, zs_j2, zj1, zj2      
      real*8 J_i2, J_j1, J_j2
      
      k=nSegments !do nothing line to stop a warning
      errval=0

c     empty segments, do nothing & return
      if (npSegment(i) .eq. 0 .or. npSegment(j) .eq. 0) return

      if (debug) print *, 'attempting to jump connect ', i, j

      zi1 = segments(1, i)      
      zi2 = segments(npSegment(i), i)      

      zj1 = segments(1, j)
      zj2 = segments(npSegment(j), j)

      call complexLensEquation(zi2, zs_i2, J_i2)
      call complexLensEquation(zj1, zs_j1, J_j1)
      call complexLensEquation(zj2, zs_j2, J_j2)

c     ends of segment i jumps to meet beginning of j
      if ( cdabs(zs_i2 - zs_j1) .lt. rootEPS) then

         if (npSegment(i)+npSegment(j) .gt. nmax) then
c	    print *, npSegment(i)+npSegment(j)-1
            errval=3213
c           stop 'too many points in i & j'
            return
         endif

         do k=1, npSegment(j)
            segments(npSegment(i)+k, i) = segments(k, j)
         enddo

         if (debug) print *, 'jump connected ', i, j, J_i2, J_j1,
     $        npSegment(i), npSegment(j)

         npSegment(i) = npSegment(i)+npSegment(j)
         npSegment(j) = 0      !segment j now contains no points

c     ends of i meet with end of j
      else if (cdabs(zs_i2 - zs_j2) .lt. rootEPS) then
         if (npSegment(i)+npSegment(j) .gt. nmax) then
c	    print *, npSegment(i)+npSegment(j)-1
            errval=3213
c           stop 'too many points in i & j'
            return
         endif
         
         do k=1, npSegment(j)
            segments(npSegment(i)+k, i) = 
     $           segments(npSegment(j)-k+1, j)
            
         enddo

         if (debug) print *, 'jump connected ', i, j, J_i2, J_j1,
     $        npSegment(i), npSegment(j)

         npSegment(i) = npSegment(i)+npSegment(j)
         npSegment(j) = 0      !segment j now contains no points
      endif

      return
      end


c
c     output imageTracks 
c
      subroutine outputTrack(ztrack, npTrack, ntrack)
      include 'global.h'

      complex*16 ztrack(nmax, 5)
c     number of points in a track
      integer npTrack(5)
c     total number of tracks
      integer ntrack
c     
      integer i, k


c     write out the tracks
      open(11, file='track.dat')
      write(11, *) ntrack
      do i=1, ntrack
         write(11, *) npTrack(i)
      enddo
      do i=1, ntrack
         do k=1, npTrack(i)
            write(11, *) dreal(ztrack(k, i)), dimag(ztrack(k,i))
         enddo
      enddo
      close(11)

      return
      end

c
c     solve the complex lens equation using zroots from Numerical Recipies, we keep all
c     the complex roots, whether they are true solutions to the lens equation or not
c     
      subroutine complexSolve(roots, jacob, errval)
      include 'global.h'

      complex*16 roots(5)
      real*8 jacob(5)

      integer i
      complex*16 c(6)
      complex*16 zs
      real*8 radius2

      integer errval

      integer zrootFlag
      common /zrootsPass/ zrootFlag

      call coefficients(c)
      call zroots(c, 5, roots, .true.)

      errval = 0

      if (zrootFlag .lt. 1) then
        errval=3214
c       stop 'complex root finding failed'
        return
      endif

      nimage = 0
      muTotal = 0.0d0

      do i=1, 5
         call complexLensEquation(roots(i), zs, jacob(i))
         radius2 = (dreal(zs)-xsCenter)**2+(dimag(zs)-ysCenter)**2
            
         if (dabs(dsqrt(radius2)/rs-1.0d0) .lt. 1.0d-8) then
            nimage = nimage + 1
            images(nimage) = roots(i)
            mu(nimage) = 1.0d0/jacob(i)
            muTotal = muTotal + dabs(mu(nimage))
         endif
      enddo


  

c     check image numbers
      if (nimage .ne. 3 .and. nimage .ne. 5) then
c         print *, 'wrong image number (complexSolve)', nimage
         errval=3215
c         stop
         return

      endif

      return
      end
