c
      subroutine kalman(f,dta,y,p,it,w,v,x,xt,imax,m,atime,stn)
c
c*********************************************************************
c
c     routine to apply kalman filter to a set of obs for qc purposes.
c     
c     original: john mcginley, noaa/fsl  spring 1998
c     changes:
c       21 aug 1998  peter stamus, noaa/fsl
c          make code dynamic, housekeeping changes, for use in laps.
c       09 dec 1999  john mcginley and peter stamus, noaa/fsl
c          new version; additional housekeeping changes too.
c
c*********************************************************************
c
      real k(m,m),p(m,m),f(m,m),ii(m,m)
      real h(m,m),ht(m,m),ft(m,m),g(m)
      real x(m),y(m),xt(m)
      real pt(m,m),zz(m,m)
      real a(m,m),b(m),c(m),z(m),d(m,m),e(m,m)
      real w(m,m),v(m,m)
      real dta(m),uu(m,m),vv(m,m)
      integer sca(2,2),scb(2,2),scf(2,2),on,off
      character atime*(*),stn(m)*5
c
c initialize arrays
c
      on=1
      off=0
      iiii=20998
c
c     fill initial matrix values
c
      call zero(ii, m,m)
      call zero( h, m,m)
      call zero(zz, m,m)
      call zero(pt, m,m)
      call zero(a, m,m)
      call zero(e, m,m)
      call zero(ht,m,m)
      
c
      do i=1,2
      do j=1,2
         sca(i,j) = 0
         scf(i,j) = 0 
         scb(i,j) = 0
      enddo !j
      enddo !i
c
c writeout parameter settings
c
      do i=1,imax
         z(i) = 0.
         ii(i,i) = 1.
         h(i,i) = 1.
      enddo !i
c
c     set obs 
c first guess - initial
c
      do i=1,imax
         x(i) = dta(i)  
      enddo !i
c
c xt=fx
c
      call mvmult(f,x,xt,imax,imax,1,m)
      call trans(f,ft,imax,imax,m)
c     call writev(f,imax,imax,m,'   f        ',atime,on ,1.0)
c     call writev(x,imax,1,m,'   x       ',atime,on,10000.)
c     call writev(xt,imax,1,m,'   xt       ',atime,on,10000.)
c
c pt=fpft+t
c
      call mvmult(f,p,a,imax,imax,imax,m)
      call mvmult(a,ft,pt,imax,imax,imax,m)
      call addmv(pt,w,pt,imax,imax,m)
      call writev(pt,imax,imax,m,'   pt       ',atime,off,0.)
c
c k=pth/(hptht+v)
c
      call mvmult(h,pt,a,imax,imax,imax,m)
      call trans(h,ht,imax,imax,m)
      call mvmult(a,ht,e,imax,imax,imax,m)
      call addmv(e,v,zz,imax,imax,m)
      call writev(zz,imax,imax,m,'hptht+v 2inv',atime,off ,0.)
      idiag=0
      call matrixanal(zz,imax,imax,m,idiag, ' hptht+v  ')
      if(idiag.eq.1) then 
       call fastinv(zz,imax,imax,m)
       call mvmult(pt,h,a,imax,imax,imax,m)
       call mvmult(a,zz,k,imax,imax,imax,m)
       go to 34
      endif
      call trans(zz,a,imax,imax,m)
      call replace(a,uu,imax,imax,m,m)
      call svdcmp(uu,imax,imax,m,m,b,vv,m)
c     call writev(b ,imax,1,m,'diag wj     ',atime,on,0.)
c     call writev(uu,imax,imax,m,'  uu svdcmp ',atime,off,0.)
c     call writev(vv,imax,imax,m,'  vv svdcmp ',atime,off,0.)
      wmax=0.
      do j=1,imax
         if(b(j) .gt. wmax) wmax = b(j)
         g(j) = b(j) 
      enddo !j
      wmin=wmax*1.e-6
      do j=1,imax
         if(b(j) .lt. wmin) b(j) = 0.  
      enddo !j
      call mvmult(pt,h,a,imax,imax,imax,m)
      call trans(a,d ,imax,imax,m)

      do j=1,imax
         do i=1,imax
            c(i) = d(i,j)
         enddo !i 
         call svbksb(uu,b,vv,imax,imax,m,m,c,z,m)
         do i=1,imax
            zz(i,j) = z(i)
         enddo !i
      enddo!on j
      call trans(zz,k,imax,imax,m)
c     call invert(a,imax,m,e,m)  
c     
c.....  this is single value decomposition solution for a
c
      call trans(uu,zz,imax,imax,m)
      call zero(e,m,m)     
      call mvmult(e,zz,uu,imax,imax,imax,m)
      call mvmult(vv,uu,zz,imax,imax,imax,m) 
      call trans(zz,a ,imax,imax,m)
c     call writev(d,imax,imax,m,'pt trans    ',atime,off,0.)
c     call writev(a,imax,imax,m,'at inverted ',atime,off,0.)
 34   call writev(k,imax,imax,m,'kalman gain ',atime,off,0.)
c
c.....  estimate obs loop
c
c x=xt+k(y-hxt)      
c
      call mvmult(h,xt,c,imax,imax,1,m)
      call submv(y,c,b,imax,1,m)
      call mvmult(k,b,c,imax,imax,1,m)
      call addmv(xt,c,x,imax,1,m)
c
c p=(i-k)pt
c
      call mvmult(k,h,e,imax,imax,imax,m)
      call submv(ii,e,a,imax,imax,m)
      call mvmult(a,pt,p,imax,imax,imax,m)
c
      sum = 0.
      write(6,2000) 
 2000 format(1x,' stn indx',' kalman x ',' forecast ',' observatn'
     &,' kalmgn','     w    ','     v    ')
      do i=1,imax
         write(6,1098) stn(i),i,x(i)-10000.,xt(i)-10000.,
     &              y(i)-10000.,k(i,i),w(i,i),v(i,i)
 1098 format(1x,a5,i4,3f10.3,f7.4,2f10.3)
         sum=sum+k(i,i)
      enddo !i
      print*, 'mean kalman ',sum/float(imax)
      write(6,*) 'mean kalman ',sum/float(imax)
      call writev(p,imax,imax,m,'anal cov err',atime,off,0.)
c
      return
      end
c
c
      subroutine kalmod(f,yta,byta,dta,ta,wmt,wot,wbt,offset,
     &                  imax,mwt,m)
c
c*********************************************************************
c
c     kalman filter tool.
c     
c     original: john mcginley, noaa/fsl  december 1999
c     changes:
c
c       09 dec 1999  peter stamus, noaa/fsl
c          housekeeping changes.
c
c*********************************************************************
c
      real yta(m),byta(m),dta(m),ta(m),wmt(m),wot(m),wbt(m)
      real mwt(m,m),f(m,m),a,b,c
c
      do i=1,imax
         sum=wmt(i)+wot(i)+wbt(i)
         a=0.5*(wmt(i)+wbt(i))/sum
         b=0.5*(wot(i)+wmt(i))/sum
         c=0.5*(wbt(i)+wot(i))/sum
         sum=0.
         sum1=0.
         if(mwt(i,i).eq.1.) then
             print*,'station ',i,' is isolated: set buddy trend to 0'
             sum=0.
          else
            do j = 1,imax
             if(i.eq.j) go to 1
             sum=mwt(i,j)/(1.-mwt(i,i))*yta(j)+sum
             f(i,j)=0.
 1           continue
            enddo
         endif
         byta(i)=sum
         f(i,i)=a*(1.+yta(i)/(ta(i)+offset)) + 
     &          c*(1.+dta(i)/(ta(i)+offset)) +
     &          b*(1+byta(i)/(ta(i)+offset))
      enddo
c
      return
      end

