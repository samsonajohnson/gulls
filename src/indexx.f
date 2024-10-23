c
c     indexing a small array
c
      SUBROUTINE indexx(n,arr,indx)
      implicit none
      INTEGER n,indx(n)
      REAL*8 arr(n)
      INTEGER i,indxt,j,l
      REAL*8 a

      do j=1,n
        indx(j)=j
      enddo

      l=1

      do j=l+1, n
         indxt=indx(j)
         a=arr(indxt)
         do i=j-1,l,-1
            if (arr(indx(i)).le.a) goto 2
            indx(i+1)=indx(i)
         enddo
         i=l-1
 2       indx(i+1)=indxt
      enddo
      
      return
      end

