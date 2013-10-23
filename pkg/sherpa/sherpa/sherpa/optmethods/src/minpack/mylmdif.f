      subroutine mylmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
     *     diag,mode,factor,nprint,info,nfev,fjac,ldfjac,
     *     ipvt,qtf,wa1,wa2,wa3,wa4,lb,ub,fmin,multicore,myfdjac,
     *     lowtri,ifault,covarerr)
c     **********
c
c     subroutine mylmdif
      implicit none
      integer m,n,maxfev,mode,nprint,info,nfev,ldfjac,multicore
      integer ifault,ipvt(n)
      double precision ftol,xtol,gtol,epsfcn,factor
      double precision x(n),fvec(m),diag(n),fjac(ldfjac,n),qtf(n),
     *                 wa1(n),wa2(n),wa3(n),wa4(m)
      double precision lb(n),ub(n),lowtri(n*(n+1)/2),covarerr(n)
      double precision fmin,enorm
      integer iflag
      external fcn,myfdjac
      iflag = 1

      call lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
     *     diag,mode,factor,nprint,info,nfev,fjac,ldfjac,
     *     ipvt,qtf,wa1,wa2,wa3,wa4,lb,ub,multicore,myfdjac,lowtri)
      fmin = enorm( m, fvec )**2.0
      call calccovar(n,wa1,ifault,lowtri,covarerr)
      return
c
c     last card of subroutine mylmdif.
c
      end
      subroutine calccovar(n,wa,ifault,lowtri,covarerr)
c     **********
c
c     subroutine calccovar
      implicit none
      integer n,ifault
      double precision wa(n),lowtri(n*(n+1)/2)
      double precision covarerr(n)
      external symmatmult, syminv
      integer ii,nullty
      double precision rmax
      
      call syminv(lowtri,n,lowtri,wa,nullty,ifault,rmax)

      if ( ifault .eq. 0 ) then
         do ii = 1, n
            if ( lowtri(ii*(ii+1)/2) .gt. 0.0d0 ) then
               covarerr(ii) = dsqrt( lowtri(ii*(ii+1)/2) )
            end if
         end do
      endif

      return
c
c     last card of subroutine calccovar.
c
      end
      subroutine symmatmult(m,n,a,b,c)
c     **********
c
c     subroutine symmatmult
      implicit none
      integer m,n
      double precision a(m,n),b(m,n),c(n*(n+1)/2)
      integer rr,cc,kk,ll
      double precision tmp
      ll = 1
      do rr = 1, n
         do cc = 1, rr
            tmp = 0.0d0
            do kk = 1, m
               tmp = tmp + a(kk,rr)*b(kk,cc)
            end do
            c( ll ) = tmp
            ll = ll + 1
         end do
      end do
      return
c
c     last card of subroutine symmatmult.
c
      end
