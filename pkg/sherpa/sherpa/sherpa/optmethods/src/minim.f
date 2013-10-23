c
      subroutine minim(p,step,nop,func,max,iprint,stopcr,nloop,iquad,
     1  simp,var,functn,ifault,neval,
     2     lb,ub,g,h,pbar,pstar,pstst,aval,bmat,pmin,vc,temp)
      implicit none
c     
c     a program for function minimization using the simplex method.
c     the minimum found will often be a local, not a global, minimum.
c
c     for details, see nelder & mead, the computer journal, january 1965
c
c     programmed by d.e.shaw,
c     csiro, division of mathematics & statistics
c     p.o. box 218, lindfield, n.s.w. 2070
c
c     with amendments by r.w.m.wedderburn
c     rothamsted experimental station
c     harpenden, hertfordshire, england
c
c     further amended by alan miller,
c     csiro, division of mathematics & statistics
c     private bag 10, clayton, vic. 3168
c
c     arguments:-
c     p()     = input, starting values of parameters
c               output, final values of parameters
c     step()  = input, initial step sizes
c     nop     = input, no. of parameters, incl. any to be held fixed
c     func    = output, the function value corresponding to the final
c               parameter values
c     max     = input, the maximum no. of function evaluations allowed
c     iprint  = input, print control parameter
c                     < 0 no printing
c                     = 0 printing of parameter values and the function
c                         value after initial evidence of convergence.
c                     > 0 as for iprint = 0 plus progress reports after
c                         every iprint evaluations, plus printing for the
c                         initial simplex.
c     stopcr  = input, stopping criterion
c     nloop   = input, the stopping rule is applied after every nloop
c               function evaluations.
c     iquad   = input, = 1 if the fitting of a quadratic surface is required
c                      = 0 if not
c     simp    = input, criterion for expanding the simplex to overcome
c               rounding errors before fitting the quadratic surface.
c     var()   = output, contains the diagonal elements of the inverse of
c               the information matrix.
c     functn  = input, name of the user's subroutine - arguments (p,func)
c               which returns the function value for a given set of
c               parameter values in array p.
c****   functn must be declared external in the calling program.
c       ifault  = output, = 0 for successful termination
c                         = 1 if maximum no. of function evaluations exceeded
c                         = 2 if information matrix is not +ve semi-definite
c                         = 3 if nop < 1
c                         = 4 if nloop < 1
c
c       advice on usage:
c       if the function minimized can be expected to be smooth in the vicinity
c       of the minimum, users are strongly urged to use the quadratic-surface
c       fitting option.   this is the only satisfactory way of testing that the
c       minimum has been found.   the value of simp should be set to at least
c       1000 times the rounding error in calculating the fitted function.
c       e.g. in double precision on a micro- or mini-computer with about 16
c       decimal digit representation of floating-point numbers, the rounding
c       errors in calculating the objective function may be of the order of
c       1.e-12 say in a particular case.   a suitable value for simp would then
c       be 1.e-08.   however, if numerical integration is required in the
c       calculation of the objective function, it may only be accurate to say
c       1.e-05 and an appropriate value for simp would be about 0.1.
c       if the fitted quadratic surface is not +ve definite (and the function
c       should be smooth in the vicinity of the minimum), it probably means
c       that the search terminated prematurely and you have not found the
c       minimum.
c
c       n.b. p, step and var (if iquad = 1) must have dimension at least nop
c            in the calling program.
c       the dimensions below are for a maximum of 20 parameters.
c      the dimension of bmat should be at least nop*(nop+1)/2.
c
c****      n.b. this version is in double precision throughout
c
c       latest revision - 11 august 1991
c
c*****************************************************************************
c
      integer*4 nop,max,iprint,nloop,iquad,ifault
      real*8 p(nop),step(nop),func,stopcr,simp,var(nop)
      real*8 g(nop+1,nop),h(nop+1),pbar(nop),pstar(nop),
     1  pstst(nop),aval(nop),bmat(nop*(nop+1)/2),pmin(nop),
     1  vc(nop*(nop+1)/2),temp(nop),lb(nop),ub(nop)
      real*8 functn
      external functn

      external ordersimplex
c      integer*4  NVMAX
c      parameter (NVMAX=128)
c
c      real*8 g(NVMAX+1,NVMAX),h(NVMAX+1),pbar(NVMAX),pstar(NVMAX),
c     1  pstst(NVMAX),aval(NVMAX),bmat(NVMAX*(NVMAX+1)/2),pmin(NVMAX),
c     1  vc(NVMAX*(NVMAX+1)/2),temp(NVMAX)
      integer*4 i,iflag,imax,imin,j
      integer*4 loop,nap,neval,np1
      real*8 fnp1,hmax,hmean,hmin,hstar,hstd,hstst,savemn
c
c     a = reflection coefficient, b = contraction coefficient, and
c     c = expansion coefficient.
c
      integer*4 lout
      real*8 a,b,c
      real*8 zero,half
      data zero/0.d0/, half/0.5d0/
      data a,b,c/1.d0, 0.5d0, 2.d0/
c
c     set lout = logical unit no. for output
c
      data lout/6/
      savemn = 0.0d0
c
c     if progress reports have been requested, print heading
c
      if(iprint.gt.0) write(lout,1000) iprint
 1000 format(' progress report every',i4,' function evaluations'/,
     1  ' eval.  func.',15x,'parameter values')
c
c     check input arguments
c
      ifault=0
      if(nop.le.0) ifault=3
      if(nloop.le.0) ifault=4
      if(ifault.ne.0) return
c
c     set nap = no. of parameters to be varied, i.e. with step.ne.0
c
      nap=0
      loop=0
      iflag=0
      do 10 i=1,nop
        if(step(i).ne.zero) nap=nap+1
   10 continue
c
c     if nap = 0 evaluate function at the starting point and return
c
      if(nap.gt.0) go to 30
      func=functn(p,nop)
      return
c
c     set up the initial simplex
c
 30   call initsimplex(nop,p,step,g,h,lb,ub,
     & nap,np1,lout,iprint,neval,functn)
c
c     start of main cycle.
c
c     find max. & min. values for current simplex (hmax & hmin).
c
  100 loop=loop+1
c
c     find max. & min. values for current simplex (hmax & hmin).
c
      call ordersimplex(np1,h,imax,imin,hmax,hmin)
c
c     find the centroid of the vertices other than p(imax)
c
      call findcentroid(nap,nop,np1,imax,g,pbar)
c
c     reflect maximum through pbar to pstar,
c     hstar = function value at pstar.
c
      do 170 i=1,nop
  170 pstar(i)=a*(pbar(i)-g(imax,i))+pbar(i)
      hstar=functn(pstar,nop)
      neval=neval+1
      if(iprint.le.0) go to 180
      if(mod(neval,iprint).eq.0) write(lout,1010) neval,hstar,
     1  (pstar(j),j=1,nop)
c
c     if hstar < hmin, reflect pbar through pstar,
c     hstst = function value at pstst.
c
  180 if(hstar.ge.hmin) go to 220
      do 190 i=1,nop
  190 pstst(i)=c*(pstar(i)-pbar(i))+pbar(i)
      hstst=functn(pstst,nop)
      neval=neval+1
      if(iprint.le.0) go to 200
      if(mod(neval,iprint).eq.0) write(lout,1010) neval,hstst,
     1  (pstst(j),j=1,nop)
c
c     if hstst < hmin replace current maximum point by pstst and
c     hmax by hstst, then test for convergence.
c
  200 if(hstst.ge.hmin) go to 320
      do 210 i=1,nop
        if(step(i).ne.zero) g(imax,i)=pstst(i)
  210 continue
      h(imax)=hstst
      go to 340
c
c     hstar is not < hmin.
c     test whether it is < function value at some point other than
c     p(imax).   if it is replace p(imax) by pstar & hmax by hstar.
c
  220 do 230 i=1,np1
        if(i.eq.imax) go to 230
        if(hstar.lt.h(i)) go to 320
  230 continue
c
c     hstar > all function values except possibly hmax.
c     if hstar <= hmax, replace p(imax) by pstar & hmax by hstar.
c
      if(hstar.gt.hmax) go to 260
      do 250 i=1,nop
        if(step(i).ne.zero) g(imax,i)=pstar(i)
  250 continue
      hmax=hstar
      h(imax)=hstar
c
c     contracted step to the point pstst,
c     hstst = function value at pstst.
c
  260 do 270 i=1,nop
  270 pstst(i)=b*g(imax,i) + (1.d0-b)*pbar(i)
      hstst=functn(pstst,nop)
      neval=neval+1
      if(iprint.le.0) go to 280
      if(mod(neval,iprint).eq.0) write(lout,1010) neval,hstst,
     1  (pstst(j),j=1,nop)
c
c     if hstst < hmax replace p(imax) by pstst & hmax by hstst.
c
  280 if(hstst.gt.hmax) go to 300
      do 290 i=1,nop
        if(step(i).ne.zero) g(imax,i)=pstst(i)
  290 continue
      h(imax)=hstst
      go to 340
c
c     hstst > hmax.
c     shrink the simplex by replacing each point, other than the current
c     minimum, by a point mid-way between its current position and the
c     minimum.
c
  300 do 315 i=1,np1
        if(i.eq.imin) go to 315
        do 310 j=1,nop
          if(step(j).ne.zero) g(i,j)=(g(i,j)+g(imin,j))*half
          p(j)=g(i,j)
  310   continue
       h(i)=functn(p,nop)
        neval=neval+1
        if(iprint.le.0) go to 315
        if(mod(neval,iprint).eq.0) write(lout,1010) neval,h(i),
     1              (p(j),j=1,nop)
  315 continue
      go to 340
c
c     replace maximum point by pstar & h(imax) by hstar.
c
  320 do 330 i=1,nop
        if(step(i).ne.zero) g(imax,i)=pstar(i)
  330 continue
      h(imax)=hstar
c
c     if loop = nloop test for convergence, otherwise repeat main cycle.
c
  340 if(loop.lt.nloop) go to 100
c
c     calculate mean & standard deviation of function values for the
c     current simplex.
c
      call meanstddev(np1,h,hmean,hstd)
c
c     if the rms > stopcr, set iflag & loop to zero and go to the
c     start of the main cycle again.
c
      if(hstd.le.stopcr.or.neval.gt.max) go to 410
      iflag=0
      loop=0
      go to 100
c
c     find the centroid of the current simplex and the function value there.
c
  410 do 380 i=1,nop
        if(step(i).eq.zero) go to 380
        p(i)=zero
        do 370 j=1,np1
  370   p(i)=p(i)+g(j,i)
        fnp1 = np1
        p(i)=p(i)/fnp1
  380 continue
      func=functn(p,nop)
      neval=neval+1
      if(iprint.le.0) go to 390
      if(mod(neval,iprint).eq.0) write(lout,1010) neval,func,
     1  (p(j),j=1,nop)
c
c     test whether the no. of function values allowed, max, has been
c     overrun; if so, exit with ifault = 1.
c
  390 if(neval.le.max) go to 420
      ifault=1
      if(iprint.lt.0) return
      write(lout,1020) max
 1020 format(' no. of function evaluations exceeds',i5)
      write(lout,1030) hstd
 1030 format(' rms of function values of last simplex =',g14.6)
      write(lout,1040)(p(i),i=1,nop)
 1040 format(' centroid of last simplex =',4(/1x,6g13.5))
      write(lout,1050) func
 1050 format(' function value at centroid =',g14.6)
      return
c
c     convergence criterion satisfied.
c     if iflag = 0, set iflag & save hmean.
c     if iflag = 1 & change in hmean <= stopcr then search is complete.
c
  420 if(iprint.lt.0) go to 430
      write(lout,1060)
 1060 format(/' evidence of convergence')
      write(lout,1040)(p(i),i=1,nop)
      write(lout,1050) func
  430 if(iflag.gt.0) go to 450
      iflag=1
  440 savemn=hmean
      loop=0
      go to 100
  450 if(abs(savemn-hmean).ge.stopcr) go to 440
      if(iprint.lt.0) go to 460
      write(lout,1070) neval
 1070 format(//' minimum found after',i5,' function evaluations')
      write(lout,1080)(p(i),i=1,nop)
 1080 format(' minimum at',4(/1x,6g13.6))
      write(lout,1090) func
 1090 format(' function value at minimum =',g14.6)
  460 if(iquad.le.0) return
      call quadsurf(lout,p,step,nop,func,iprint,
     1     simp,var,functn,ifault,neval,
     2     g,h,pbar,pstar,pstst,aval,bmat,pmin,vc,temp,max)
 1010      format(/i4, 2x, g12.5, 2x, 5g12.5, 3(/20x, 5g12.5))
      return
      end

c-------------------------------------------------------------------
c
c     quadratic surface fitting
c
c-------------------------------------------------------------------

      subroutine quadsurf(lout,p,step,nop,func,iprint,
     1     simp,var,functn,ifault,neval,
     2     g,h,pbar,pstar,pstst,aval,bmat,pmin,vc,temp,maxeval)
      implicit none
      integer*4 lout,nop,iprint,ifault
      real*8 p(nop),step(nop),func,simp,var(nop)
      real*8 g(nop+1,nop),h(nop+1),pbar(nop),pstar(nop),
     1  pstst(nop),aval(nop),bmat(nop*(nop+1)/2),pmin(nop),
     1  vc(nop*(nop+1)/2),temp(nop),maxeval
      real*8 functn
      external functn

c
c     quadratic surface fitting
c
      integer*4 i,i1,i2,ii,ij,ijk,irank,j,j1,jj
      integer*4 k,l,nap,neval,np1,nullty
      real*8 a0,hstar,hstst,rmax
      real*8 test,ymin
c
c     a = reflection coefficient, b = contraction coefficient, and
c     c = expansion coefficient.
c
      real*8 zero,one,two,three,half
      data zero/0.d0/, one/1.d0/, two/2.d0/, three/3.d0/, half/0.5d0/

      nap = nop
      np1 = nop + 1
  410 do 380 i=1,nop
        if(step(i).eq.zero) go to 380
        p(i)=zero
        do 370 j=1,np1
  370   p(i)=p(i)+g(j,i)
c        fnp1 = np1
        p(i)=p(i)/np1
  380 continue
      func=functn(p,nop)
      neval=neval+1
      if ( neval .ge. maxeval ) return
      if(iprint.ge.0) write(lout,1110)
 1110 format(/' quadratic surface fitting about supposed minimum'/)
c
c     expand the final simplex, if necessary, to overcome rounding
c     errors.
c
c      neval=0
      do 490 i=1,np1
  470   test=abs(h(i)-func)
        if(test.ge.simp) go to 490
        do 480 j=1,nop
          if(step(j).ne.zero) g(i,j)=(g(i,j)-p(j))+g(i,j)
          pstst(j)=g(i,j)
  480   continue
        h(i)=functn(pstst,nop)
        neval=neval+1
        if ( neval .ge. maxeval ) return
        go to 470
  490 continue
c
c     function values are calculated at an additional nap points.
c
      do 510 i=1,nap
        i1=i+1
        do 500 j=1,nop
  500   pstar(j)=(g(1,j)+g(i1,j))*half
        aval(i)=functn(pstar,nop)
        neval=neval+1
        if ( neval .ge. maxeval ) return
  510 continue
c
c     the matrix of estimated second derivatives is calculated and its
c     lower triangle stored in bmat.
c
      a0=h(1)
      do 540 i=1,nap
        i1=i-1
        i2=i+1
        if(i1.lt.1) go to 540
        do 530 j=1,i1
          j1=j+1
          do 520 k=1,nop
  520     pstst(k)=(g(i2,k)+g(j1,k))*half
          hstst=functn(pstst,nop)
          neval=neval+1
          if ( neval .ge. maxeval ) return
          l=i*(i-1)/2+j
          bmat(l)=two*(hstst+a0-aval(i)-aval(j))
  530   continue
  540 continue
      l=0
      do 550 i=1,nap
        i1=i+1
        l=l+i
        bmat(l)=two*(h(i1)+a0-two*aval(i))
  550 continue
c
c     the vector of estimated first derivatives is calculated and
c     stored in aval.
c
      do 560 i=1,nap
        i1=i+1
        aval(i)=two*aval(i)-(h(i1)+three*a0)*half
  560 continue
c
c     the matrix q of nelder & mead is calculated and stored in g.
c
      do 570 i=1,nop
  570 pmin(i)=g(1,i)
      do 580 i=1,nap
        i1=i+1
        do 580 j=1,nop
        g(i1,j)=g(i1,j)-g(1,j)
  580 continue
      do 590 i=1,nap
        i1=i+1
        do 590 j=1,nop
          g(i,j)=g(i1,j)
  590 continue
c
c     invert bmat
c
      call syminv(bmat,nap,bmat,temp,nullty,ifault,rmax)
      if(ifault.ne.0) go to 600
      irank=nap-nullty
      go to 610
  600 if(iprint.ge.0) write(lout,1120)
 1120 format(/' matrix of estimated second derivatives not +ve defn.'/
     1  ' minimum probably not found'/)
      ifault=2
      return
c
c     bmat*a/2 is calculated and stored in h.
c
  610 do 650 i=1,nap
        h(i)=zero
        do 640 j=1,nap
          if(j.gt.i) go to 620
          l=i*(i-1)/2+j
          go to 630
  620     l=j*(j-1)/2+i
  630     h(i)=h(i)+bmat(l)*aval(j)
  640   continue
  650 continue
c
c     find the position, pmin, & value, ymin, of the minimum of the
c     quadratic.
c
      ymin=zero
      do 660 i=1,nap
  660 ymin=ymin+h(i)*aval(i)
      ymin=a0-ymin
      do 670 i=1,nop
        pstst(i)=zero
        do 670 j=1,nap
  670 pstst(i)=pstst(i)+h(j)*g(j,i)
      do 680 i=1,nop
  680 pmin(i)=pmin(i)-pstst(i)
      if(iprint.lt.0) go to 682
      write(lout,1130) ymin,(pmin(i),i=1,nop)
 1130 format(' minimum of quadratic surface =',g14.6,' at',
     1  4(/1x,6g13.5))
      write(lout,1150)
 1150 format(' if this differs by much from the minimum estimated',
     1  1x,'from the minimization,'/
     2  ' the minimum may be false &/or the information matrix may be',
     3  1x,'inaccurate'/)
c
c     calculate true function value at the minimum of the quadratic.
c
  682 neval = neval + 1
      hstar=functn(pmin,nop)
      if ( neval .ge. maxeval ) return
c
c     if hstar < func, replace search minimum with quadratic minimum.
c
      if (hstar .ge. func) go to 690
      func = hstar
      do 684 i = 1, nop
  684 p(i) = pmin(i)
c      write(lout, 1140) func
 1140 format(' true func. value at minimum of quadratic = ', g14.6/)
c
c     q*bmat*q'/2 is calculated & its lower triangle stored in vc
c
  690 do 760 i=1,nop
        do 730 j=1,nap
          h(j)=zero
          do 720 k=1,nap
            if(k.gt.j) go to 700
            l=j*(j-1)/2+k
            go to 710
  700       l=k*(k-1)/2+j
  710       h(j)=h(j)+bmat(l)*g(k,i)*half
  720     continue
  730   continue
        do 750 j=i,nop
          l=j*(j-1)/2+i
          vc(l)=zero
          do 740 k=1,nap
  740     vc(l)=vc(l)+h(k)*g(k,j)
  750   continue
  760 continue
c
c     the diagonal elements of vc are copied into var.
c
      j=0
      do 770 i=1,nop
        j=j+i
        var(i)=vc(j)
  770    continue
      if(iprint.lt.0) return
      write(lout,1160) irank
 1160 format(' rank of information matrix =',i3/
     1  ' generalized inverse of information matrix:-')
      ijk=1
      go to 880
  790 continue
      write(lout,1170)
 1170 format(/' if the function minimized was -log(likelihood),'/
     1  ' this is the covariance matrix of the parameters'/
     2  ' if the function was a sum of squares of residuals'/
     3  ' this matrix must be multiplied by twice the estimated',
     4  1x,'residual variance'/' to obtain the covariance matrix.'/)
      call syminv(vc,nap,bmat,temp,nullty,ifault,rmax)
c
c     bmat now contains the information matrix
c
      write(lout,1190)
 1190 format(' information matrix:-'/)
      ijk=3
      go to 880
c
c     calculate correlations of parameter estimates, put into vc.
c
  800 ijk=2
      ii=0
      ij=0
      do 840 i=1,nop
        ii=ii+i
        if(vc(ii).gt.zero) then
          vc(ii)=one/sqrt(vc(ii))
        else 
          vc(ii)=zero
	end if
        jj=0
        do 830 j=1,i-1
          jj=jj+j
          ij=ij+1
          vc(ij)=vc(ij)*vc(ii)*vc(jj)
  830   continue
        ij=ij+1
  840 continue
      write(lout,1200)
 1200 format(/' correlation matrix:-')
      ii=0
      do 850 i=1,nop
        ii=ii+i
        if(vc(ii).ne.zero) vc(ii)=one
  850 continue
      go to 880
  860 write(lout,1210) neval
 1210 format(/' a further',i4,' function evaluations have been used'/)
      return
c
c     pseudo-subroutine to print vc if ijk = 1 or 2, or
c     bmat if ijk = 3.
c
  880 l=1
  890 if(l.gt.nop) go to (790,860,800),ijk
      ii=l*(l-1)/2
      do 910 i=l,nop
        i1=ii+l
        ii=ii+i
        i2=min(ii,i1+5)
        if(ijk.eq.3) go to 900
        write(lout,1230)(vc(j),j=i1,i2)
        go to 910
  900   write(lout,1230)(bmat(j),j=i1,i2)
  910 continue
 1230 format(1x,6g13.5)
      write(lout,1240)
 1240 format(/)
      l=l+6
      go to 890
      return
      end

      subroutine findcentroid(nap,nop,np1,imax,g,pbar)
      implicit none
      integer nap,nop,np1,imax
      real*8 g(nop+1,nop),pbar(nop)

      integer*4 i,j
      real*8 fnap
      
      do 130 i=1,nop
  130 pbar(i)=0.0d0
      do 150 i=1,np1
        if(i.eq.imax) go to 150
        do 140 j=1,nop
  140   pbar(j)=pbar(j)+g(i,j)
  150 continue
      do 160 j=1,nop
      fnap = nap
  160 pbar(j)=pbar(j)/fnap
      return
      end

      subroutine initsimplex(nop,p,step,g,h,lb,ub,
     &     nap,np1,lout,iprint,neval,functn)
      implicit none
      integer*4 nop,nap,np1,lout,iprint,neval
      real*8 p(nop),step(nop),g(nop+1,nop),h(nop+1),ub(nop),lb(nop)
      real*8 functn
      external functn

      integer*4 i,j,irow

   30 do 40 i=1,nop
   40 g(1,i)=p(i)
      irow=2
      do 60 i=1,nop
        if(step(i).eq.0.0d0) go to 60
        do 50 j=1,nop
   50   g(irow,j)=p(j)
c        g(irow,i)=p(i)+step(i)
c dtn
        g(irow,i) = dmax1(lb(i),dmin1(p(i)+step(i),ub(i)))
c dtn
        irow=irow+1
   60 continue
      np1=nap+1
      neval=0
      do 90 i=1,np1
        do 70 j=1,nop
 70        p(j)=g(i,j)
           h(i)=functn(p,nop)
           neval=neval+1
           if(iprint.le.0) go to 90
           write(lout,1010) neval,h(i),(p(j),j=1,nop)
 1010      format(/i4, 2x, g12.5, 2x, 5g12.5, 3(/20x, 5g12.5))
 90   continue

      return
      end

      subroutine meanstddev(np1,h,hmean,hstd)
      implicit none

      integer np1
      real*8 h(np1),hmean,hstd

      integer i
      real*8 fnp1

      hstd=0.0d0
      hmean=0.0d0
      do 350 i=1,np1
  350 hmean=hmean+h(i)
      fnp1 = np1
      hmean=hmean/fnp1
      do 360 i=1,np1
  360 hstd=hstd+(h(i)-hmean)**2
      hstd=sqrt(hstd/fnp1)
      return
      end

      subroutine ordersimplex(np1,h,imax,imin,hmax,hmin)
      implicit none
      integer*4 np1,imax,imin
      real*8 h(np1),hmax,hmin

      integer*4 i

      imax=1
      imin=1
      hmax=h(1)
      hmin=h(1)
      do 120 i=2,np1
        if(h(i).le.hmax) go to 110
        imax=i
        hmax=h(i)
        go to 120
  110   if(h(i).ge.hmin) go to 120
        imin=i
        hmin=h(i)
  120 continue

      return
      end
