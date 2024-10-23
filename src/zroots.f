      SUBROUTINE zroots(a,m,roots,polish)
      integer zrootFlag
      common /zrootsPass/ zrootFlag
      INTEGER m,MAXM
      REAL*8 EPS
      COMPLEX*16 a(m+1),roots(m)
      LOGICAL polish
      PARAMETER (EPS=1.d-15,MAXM=101)
CU    USES laguer
      INTEGER j,jj,its
      COMPLEX*16 ad(MAXM),x,b,c
      do 11 j=1,m+1
        ad(j)=a(j)
11    continue
      do 13 j=m,1,-1
c         x = dcmplx(0.0d0, 0.0d0)
        x=roots(j)
        call laguer(ad,j,x,its)
        if(dabs(dimag(x)).le.2.0d0*EPS**2*dabs(dreal(x))) then
           x=dcmplx(dreal(x),0.d0)
        endif
        roots(j)=x
        b=ad(j+1)
        do 12 jj=j,1,-1
          c=ad(jj)
          ad(jj)=b
          b=x*b+c
12      continue
13    continue
      if (polish) then
        do 14 j=1,m
          call laguer(a,m,roots(j),its)
14      continue
      endif

      return
      END



      SUBROUTINE laguer(a,m,x,its)
      integer zrootFlag
      common /zrootsPass/ zrootFlag
      INTEGER m,its,MAXIT,MR,MT
      REAL*8 EPSS
      COMPLEX*16 a(m+1),x
      PARAMETER (EPSS=1.0d-15,MR=8,MT=10,MAXIT=MT*MR)
      INTEGER iter,j
      REAL*8 abx,abp,abm,err,frac(MR)
      COMPLEX*16 dx,x1,b,d,f,g,h,sq,gp,gm,g2
      SAVE frac
      DATA frac /.5,.25,.75,.13,.38,.62,.88,1./
      zrootFlag = 1
      do 12 iter=1,MAXIT
        its=iter
        b=a(m+1)
        err=cdabs(b)
        d=dcmplx(0.d0,0.d0)
        f=dcmplx(0.d0,0.d0)
        abx=cdabs(x)
        do 11 j=m,1,-1
          f=x*f+d
          d=x*d+b
          b=x*b+a(j)
          err=cdabs(b)+abx*err
11      continue
        err=EPSS*err
        if(cdabs(b).le.err) then
          return
        else
          g=d/b
          g2=g*g
          h=g2-2.*f/b
          sq=cdsqrt((m-1)*(m*h-g2))
          gp=g+sq
          gm=g-sq
          abp=cdabs(gp)
          abm=cdabs(gm)
          if(abp.lt.abm) gp=gm
          if (max(abp,abm).gt.0.d0) then
            dx=m/gp
          else
            dx=exp(dcmplx(dlog(1.+abx),dble(iter)))
          endif
        endif
        x1=x-dx
        if(x.eq.x1)return
        if (mod(iter,MT).ne.0) then
          x=x1
        else
          x=x-dx*frac(iter/MT)
        endif
12    continue

      zrootFlag = 0
c      pause 'too many iterations in laguer'
      return
      END
