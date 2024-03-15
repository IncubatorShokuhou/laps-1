

      subroutine posint(t,planet,r1,r2,r3)
      implicit real*8 (a,b,c,d,e,f,g,h,o,p,q,r,s,t,u,v,w,x,y,z)
      include '../../include/astparms.for'
      integer*4 aw,pl,dg
      parameter (aw = 59)
      parameter (pl = 13)
      parameter (dg = 13)
      integer s,arycen,planet,apoint,aryhw,degree,degm1,topof,apntp1
     .,apntp2,apntp3,apntm1,apntm2
      real*8 lma,ksq,lat,lon,lst,mx,my,mz
      dimension a1(dg),ap1(dg),b1(dg),g1(dg),blc(dg),er(40),ev(39),s(14)
      common x(pl,aw),y(pl,aw),z(pl,aw),xdd(pl,aw),ydd(pl,aw),zdd(pl,aw)
     .,xd(pl),yd(pl),zd(pl),ksq,a(dg),ap(dg),b(dg),g(dg),tcent,h
     .,arg3                ,iegree,npinms
     .,init,maxpnt,minpnt,arycen,aryhw,number,nrecen,numb1,iarg1,iarg3

      character*14 c_filename,c_fileref
      character*7  c_file1 
      data c_filename /'posint.2440800'/
      data c_file1    /'posout.'/

      data init1/0/,nplan/13/,nmass/12/,moon/1/,topof/0/,nrel/9/
      data emratp1/82.3007d0/
      data ap1/4621155471343.d0,-13232841914856.d0,47013743726958.d0
     .,-114321700672600.d0,202271967611865.d0,-266609549584656.d0
     .,264429021895332.d0,-197106808276656.d0,108982933333425.d0
     .,-43427592828040.d0,11807143978638.d0,-1962777574776.d0
     .,150653570023.d0/
      data b1/9072652009253.d0,-39726106418680.d0,140544566352762.d0
     .,-344579280210120.d0,613137294629235.d0,-811345852376496.d0
     .,807012356281740.d0,-602852367932304.d0,333888089374395.d0
     .,-133228219027160.d0,36262456774618.d0,-6033724094760.d0
     .,463483373517.d0/
      data a1/150653570023.d0,2662659061044.d0,-1481863453062.d0
     .,3926822700380.d0,-6604398106155.d0,8380822992264.d0
     .,-8088023425188.d0,5907495735864.d0,-3215663657055.d0
     .,1265630766980.d0,-340671801462.d0,56165516844.d0,-4281164477.d0/
      data g1/8153167962181.d0,-30282470168220.d0,101786970400854.d0
     .,-244754720525140.d0,430901143447815.d0,-566360031923352.d0
     .,560707147204260.d0,-417423287061288.d0,230582338000515.d0
     .,-91816356051340.d0,24948866958246.d0,-4145485036740.d0
     .,318065528209.d0/
      data blc/13.,-78.,286.,-715.,1287.,-1716.,1716.,-1287.,715.,-286.,
     .78.,-13.,1./
      data er/2440800.5d0
     .,0.6388226331871791d+0,-.7234663691390341d+0,-.3137200168890611d+0
     .,-.3472826956196578d+0,-.2538198621886492d+0,-.1002526349881503d+0
     .,-.2729258630496380d+0,-.6191874090265185d+0,-.2617889156852868d+0
     .,-.1043274450714903d+1,0.1149609567068845d+1,0.5556917256839074d+0
     .,-.4241139166447770d+1,-.3149106786784197d+1,-.1247304720843395d+1
     .,0.6459610758313189d+1,0.6101130124894796d+1,0.2243382247325508d+1
     .,-.1814311868095884d+2,-.2460531337350525d+1,-.8217265400205982d+0
     .,-.1530833404391927d+2,-.2435427381300139d+2,-.9591188954641517d+1
     .,-.3031628535964815d+2,-.1792426075360332d+1,0.8620253359830489d+1
     .,0.2789435760038078d+1,0.8315178079235009d+0,-.1770597331965483d+0
     .,0.2838518810188359d+1,-.1731274670350053d+1,0.1349237040557463d+0
     .,-.2229853716668692d+1,-.3544645600066997d+0,0.1505957250879900d+0
     .,0.2055647541706483d+1,0.4646808549219674d+0,0.571194552728441d-2/
      data ev/
     . 0.1308702020817010d-1,0.9880582391240620d-2,0.4284464891607220d-2
     .,0.1164351062650545d-1,-.1796891116169349d-1,-.1081911400871376d-1
     .,0.1860072818427010d-1,-.6590509148918250d-2,-.4144660354448460d-2
     .,-.1028575666404713d-1,-.7076375256663700d-2,-.2972120784045060d-2
     .,0.4626311466814210d-2,-.5059626431851830d-2,-.2283653988263900d-2
     .,-.4252234340822560d-2,0.3559265916616560d-2,0.1655565829382040d-2
     .,0.5261924464533600d-3,-.3735782236086040d-2,-.1644317466285200d-2
     .,0.2693302696796940d-2,-.1429477504289330d-2,-.6534280574085100d-3
     .,0.3986115842506400d-3,-.3146654475399760d-2,-.1114235221665060d-2
     .,-.2763866219486554d-2,0.8274414179803879d-2,0.4463057993782906d-2
     .,0.3424650232882695d-2,0.7498263250144949d-2,-.1733648725951915d-2
     .,0.2299023705728994d-2,-.1052177497695365d-1,-.4500463020367574d-2
     .,-.4717858754284870d-2,0.1190825966047914d-1,0.242366396205390d-2/
      data clty,cltz/-.3978811766d0,.9174369566d0/
!     data k/.01720209895d0/,difmax/0.d0/,lat/40.72d0/,lon/-73.87d0/
      data difmax/0.d0/,lat/40.72d0/,lon/-73.87d0/

      if(init.ne.0)go to 1000

!     read in initial conditions
!     read(11,11)c_filename
      c_filename = c_fileref
11    format(a14)
      if(c_filename .ne. c_fileref)then
          write(6,*)' opening planet file for initial conditions: '
     1                ,c_filename
          open(12,file='/home/fab/albers/ast/planets/'
     1                 //c_filename,status='old')
          read(12,12)er,ev
          write(6,12)er,ev
12        format(//1x,d25.10/26(/1x,3d26.19))
          close(12)
      else
          write(6,*)' using default initial conditions'
      endif ! we are reading from an external file (not using default data)
!     read(11,*)topof,lat,lon
      write(6,*)' topof,lat,lon',topof,lat,lon

      h=1.d0
      arycen=aw/2 + 1
      degree=13
      degm1=degree-1
      hsq=h*h
      aryhw=arycen-1
      iarg2=degree/2
      nrecen=arycen+iarg2
      iarg1=iarg2+1
      iarg3=arycen-iarg1
      arg3=float(iarg3)*h
      number=arycen
      numb2=number+1
      numb3=number+degm1
      numb1=number-1
      ksq=k*k
      npinms=min0(nplan-1,nmass)
      phi=lat*rpd

      do i=1,13
          a(i)=a1(i)/2615348736000.d0*hsq
          ap(i)=ap1(i)/2615348736000.d0*hsq
          b(i)=b1(i)/1743565824000.d0*hsq
          g(i)=g1(i)/5230697472000.d0*hsq
      enddo

      idir=0
c
c initialize the integration
      tcent=er(1)
      iscr=1

      do j=1,nplan
          x(j,number)=er(iscr+1)
          y(j,number)=er(iscr+2)
          z(j,number)=er(iscr+3)

          xd(j)=ev(iscr)
          yd(j)=ev(iscr+1)
          zd(j)=ev(iscr+2)

          iscr=iscr+3

      enddo ! j

      call accel(number,nplan,nrel,nmass,npinms,tcent)
      go to 100
c
c extrapolate the accelerations backwards
25    idir=1
      difmx1=difmax
      difmax=0.
      do np=1,degm1
          n=number-np
          n1=n+1

          do j=1,nplan
              xddold=xdd(j,n)
              yddold=ydd(j,n)
              zddold=zdd(j,n)
              xdd(j,n)=blc(1)*xdd(j,n1)
              ydd(j,n)=blc(1)*ydd(j,n1)
              zdd(j,n)=blc(1)*zdd(j,n1)

              do i=2,degree
                  xdd(j,n)=xdd(j,n)+blc(i)*xdd(j,n+i)
                  ydd(j,n)=ydd(j,n)+blc(i)*ydd(j,n+i)
                  zdd(j,n)=zdd(j,n)+blc(i)*zdd(j,n+i)

              enddo ! i

              if(init1.ge.5)then
                  difmax=
     1     dmax1(difmax,abs(xdd(j,n)-xddold),abs(ydd(j,n)-yddold),
     .             abs(zdd(j,n)-zddold))
              endif

          enddo ! j

      enddo ! np

      write(6,9101)difmax,difmx1
      if(init1.gt.25)then
          write(6,*)' maximum iterations exceeded'
          stop
      endif
      if(                     difmax.lt.1d-9.and.init1.gt.5)go to 990
!     if(difmax.gt.difmx1.and.difmax.lt.1d-9.and.init1.gt.5)go to 990
100   init1=init1+1
c
c calculate x(-1), y(-1), z(-1)
      do j=1,nplan
          x(j,numb1)=x(j,number)-h*ev(j*3-2)
          y(j,numb1)=y(j,number)-h*ev(j*3-1)
          z(j,numb1)=z(j,number)-h*ev(j*3)

          do n=1,degree
              iarg=number+(n-1)*idir
              x(j,numb1)=x(j,numb1)+g(n)*xdd(j,iarg)
              y(j,numb1)=y(j,numb1)+g(n)*ydd(j,iarg)
              z(j,numb1)=z(j,numb1)+g(n)*zdd(j,iarg)

          enddo ! n
      enddo ! j
c
c integrate forward 12 steps
      do n=numb2,numb3
          if(init1.ne.1)then

              do m=1,degree
                  s(m)=n-m+1
              enddo ! m

          else

              do m=1,degree
                  s(m)=number
              enddo ! m

          endif ! init1 .ne. 1


          nm1=n-1
          nm2=n-2

          do j=1,nplan
              x(j,n)=x(j,nm1)+x(j,nm1)-x(j,nm2)
              y(j,n)=y(j,nm1)+y(j,nm1)-y(j,nm2)
              z(j,n)=z(j,nm1)+z(j,nm1)-z(j,nm2)

              do m=1,degree
                  x(j,n)=x(j,n)+a(m)*xdd(j,s(m))
                  y(j,n)=y(j,n)+a(m)*ydd(j,s(m))
                  z(j,n)=z(j,n)+a(m)*zdd(j,s(m))

              enddo ! m
          enddo ! j
      enddo ! n

      do n=numb2,numb3

          if(nrel.ne.0)then
              do j=1,nrel
                  xd(j)=x(j,n-1)-x(j,n-2)
                  yd(j)=y(j,n-1)-y(j,n-2)
                  zd(j)=z(j,n-1)-z(j,n-2)

                  do m=1,degree
                      nmm=n-m
                      if(init1.eq.1)nmm=number
                      xd(j)=xd(j)+b(m)*xdd(j,nmm)
                      yd(j)=yd(j)+b(m)*ydd(j,nmm)
                      zd(j)=zd(j)+b(m)*zdd(j,nmm)

                  enddo ! m

                  xd(j)=xd(j)/h
                  yd(j)=yd(j)/h
                  zd(j)=zd(j)/h

              enddo ! j
          endif ! nrel .ne. 0

          t1=er(1)+float(n-arycen)*h
610       call accel(n,nplan,nrel,nmass,npinms,t1)

      enddo ! n

      go to 25

990   init=1
      maxpnt=numb3-arycen
      minpnt=number-arycen
c
c determine whether to integrate or interpolate
1000  rsteps=(t-tcent)/h
      isteps=int(rsteps)

      if(isteps.eq.0)then
          if(rsteps.ne.0.d0)then
              idpl=int(rsteps/abs(rsteps))

          else
              idpl=-1

          endif
      else
          idpl=isteps/iabs(isteps)

      endif ! isteps .eq. 0

      p=abs(rsteps-float(isteps))

      if(planet.eq.0)then
          maxint=isteps
          if(idpl.eq.-1)maxint=maxint+degree
          minint=maxint-degree
          idr1=idpl

      else
          maxint=isteps+(5+idpl)/2
          minint=isteps-(5-idpl)/2

      endif ! planet .eq. 0


1050  if(maxint.gt.maxpnt)then
          maxpnt=maxpnt+1
          ipoint=maxpnt
          idir=1

      else
          if(minint.lt.minpnt)then
              minpnt=minpnt-1
              ipoint=minpnt
              idir=-1
          else
              goto4000
          endif

      endif

      if(iabs(ipoint).gt.aryhw)then
c
c rearrange arrays
          n1=arycen*(1-idir)
          n2=arycen-iarg1*idir

          do n=1,nrecen
              n1=n1+idir
              n2=n2+idir
              do j=1,nplan
                  x(j,n1)=x(j,n2)
                  y(j,n1)=y(j,n2)
                  z(j,n1)=z(j,n2)
                  xdd(j,n1)=xdd(j,n2)
                  ydd(j,n1)=ydd(j,n2)
                  zdd(j,n1)=zdd(j,n2)
              enddo
          enddo

          if(idir.le.0)then
              ipoint=-iarg1
              tcent=tcent-arg3
              isteps=isteps+iarg3
              maxint=maxint+iarg3
              minint=minint+iarg3
              maxpnt=aryhw
              minpnt=ipoint

          else
              ipoint=iarg1
              tcent=tcent+arg3
              isteps=isteps-iarg3
              maxint=maxint-iarg3
              minint=minint-iarg3
              maxpnt=ipoint
              minpnt=-aryhw

          endif ! idir .le. 0

      endif ! rearrange arrays

c
c integrate forward or backward one step
      s(1)=ipoint+arycen
      s(2)=s(1)-idir
      s(3)=s(2)-idir
      s(4)=s(3)-idir
      s(5)=s(4)-idir
      s(6)=s(5)-idir
      s(7)=s(6)-idir
      s(8)=s(7)-idir
      s(9)=s(8)-idir
      s(10)=s(9)-idir
      s(11)=s(10)-idir
      s(12)=s(11)-idir
      s(13)=s(12)-idir
      s(14)=s(13)-idir
      do 2500 i=1,nplan
      x(i,s(1))=x(i,s(2))+x(i,s(2))-x(i,s(3))+ap(1)*xdd(i,s(2))
     .+ap(2)*xdd(i,s(3))+ap(3)*xdd(i,s(4))+ap(4)*xdd(i,s(5))
     .+ap(5)*xdd(i,s(6))+ap(6)*xdd(i,s(7))+ap(7)*xdd(i,s(8))
     .+ap(8)*xdd(i,s(9))+ap(9)*xdd(i,s(10))+ap(10)*xdd(i,s(11))
     .+ap(11)*xdd(i,s(12))+ap(12)*xdd(i,s(13))+ap(13)*xdd(i,s(14))
      y(i,s(1))=y(i,s(2))+y(i,s(2))-y(i,s(3))+ap(1)*ydd(i,s(2))
     .+ap(2)*ydd(i,s(3))+ap(3)*ydd(i,s(4))+ap(4)*ydd(i,s(5))
     .+ap(5)*ydd(i,s(6))+ap(6)*ydd(i,s(7))+ap(7)*ydd(i,s(8))
     .+ap(8)*ydd(i,s(9))+ap(9)*ydd(i,s(10))+ap(10)*ydd(i,s(11))
     .+ap(11)*ydd(i,s(12))+ap(12)*ydd(i,s(13))+ap(13)*ydd(i,s(14))
      z(i,s(1))=z(i,s(2))+z(i,s(2))-z(i,s(3))+ap(1)*zdd(i,s(2))
     .+ap(2)*zdd(i,s(3))+ap(3)*zdd(i,s(4))+ap(4)*zdd(i,s(5))
     .+ap(5)*zdd(i,s(6))+ap(6)*zdd(i,s(7))+ap(7)*zdd(i,s(8))
     .+ap(8)*zdd(i,s(9))+ap(9)*zdd(i,s(10))+ap(10)*zdd(i,s(11))
     .+ap(11)*zdd(i,s(12))+ap(12)*zdd(i,s(13))+ap(13)*zdd(i,s(14))
      if(nrel.lt.i)go to 2500
      xd(i)=(x(i,s(2))-x(i,s(3))+b(1)*xdd(i,s(2))+b(2)*xdd(i,s(3))
     .+b(3)*xdd(i,s(4))+b(4)*xdd(i,s(5))+b(5)*xdd(i,s(6))
     .+b(6)*xdd(i,s(7))+b(7)*xdd(i,s(8))+b(8)*xdd(i,s(9))
     .+b(9)*xdd(i,s(10))+b(10)*xdd(i,s(11))+b(11)*xdd(i,s(12))
     .+b(12)*xdd(i,s(13))+b(13)*xdd(i,s(14)))/h
      yd(i)=(y(i,s(2))-y(i,s(3))+b(1)*ydd(i,s(2))+b(2)*ydd(i,s(3))
     .+b(3)*ydd(i,s(4))+b(4)*ydd(i,s(5))+b(5)*ydd(i,s(6))
     .+b(6)*ydd(i,s(7))+b(7)*ydd(i,s(8))+b(8)*ydd(i,s(9))
     .+b(9)*ydd(i,s(10))+b(10)*ydd(i,s(11))+b(11)*ydd(i,s(12))
     .+b(12)*ydd(i,s(13))+b(13)*ydd(i,s(14)))/h
      zd(i)=(z(i,s(2))-z(i,s(3))+b(1)*zdd(i,s(2))+b(2)*zdd(i,s(3))
     .+b(3)*zdd(i,s(4))+b(4)*zdd(i,s(5))+b(5)*zdd(i,s(6))
     .+b(6)*zdd(i,s(7))+b(7)*zdd(i,s(8))+b(8)*zdd(i,s(9))
     .+b(9)*zdd(i,s(10))+b(10)*zdd(i,s(11))+b(11)*zdd(i,s(12))
     .+b(12)*zdd(i,s(13))+b(13)*zdd(i,s(14)))/h
c     print*,i,s(1),x(i,s(1))
2500  continue

      t1=tcent+h*float(ipoint)
c     print*,nplan,s(1),x(1,s(1)),x(13,s(1))
      call accel(s(1),nplan,nrel,nmass,npinms,t1)
      go to 1050

4000  apoint=isteps+arycen
      if(planet.eq.0)go to 9000
      psq=p*p
      pcubd=psq*p
      apntm1=apoint-idpl
      apntp1=apoint+idpl
      apntp2=apntp1+idpl
c
c interpolate to find position of planet at time t
      b2=.25*(psq-p)
      b3=(.5*p-1.5*psq+pcubd)/6.d0
      dx01=x(planet,apntp1)-x(planet,apoint)
      dy01=y(planet,apntp1)-y(planet,apoint)
      dz01=z(planet,apntp1)-z(planet,apoint)
      dx02=x(planet,apntp1)-2.*x(planet,apoint)+x(planet,apntm1)
      dy02=y(planet,apntp1)-2.*y(planet,apoint)+y(planet,apntm1)
      dz02=z(planet,apntp1)-2.*z(planet,apoint)+z(planet,apntm1)
      dx12=x(planet,apntp2)-2.*x(planet,apntp1)+x(planet,apoint)
      dy12=y(planet,apntp2)-2.*y(planet,apntp1)+y(planet,apoint)
      dz12=z(planet,apntp2)-2.*z(planet,apntp1)+z(planet,apoint)
      dx03=dx12-dx02
      dy03=dy12-dy02
      dz03=dz12-dz02
      r1=x(planet,apoint)+p*dx01+b2*(dx02+dx12)+b3*dx03
      r2=y(planet,apoint)+p*dy01+b2*(dy02+dy12)+b3*dy03
      r3=z(planet,apoint)+p*dz01+b2*(dz02+dz12)+b3*dz03
      if(planet.ne.2)go to 4400
      apntp3=apntp2+idpl
      apntm2=apntm1-idpl
      b4=p*(pcubd-psq-psq-p+2.d0)/48.d0
      ddx4=x(planet,apntp3)-2.*x(planet,apntp2)+x(planet,apntp1)-dx12
     .    +x(planet,apoint)-2.*x(planet,apntm1)+x(planet,apntm2)-dx02
      ddy4=y(planet,apntp3)-2.*y(planet,apntp2)+y(planet,apntp1)-dy12
     .    +y(planet,apoint)-2.*y(planet,apntm1)+y(planet,apntm2)-dy02
      ddz4=z(planet,apntp3)-2.*z(planet,apntp2)+z(planet,apntp1)-dz12
     .    +z(planet,apoint)-2.*z(planet,apntm1)+z(planet,apntm2)-dz02
      r1=r1+b4*ddx4
      r2=r2+b4*ddy4
      r3=r3+b4*ddz4
4400  continue


      if(planet.eq.1.and.moon.eq.1)then

!         add lunar perturbation to earth's position

!         get earth moon vector for equinox of date
          call moon_brwn(t,mx,my,mz)

!         convert vector to 1950 coordinates
          call preces(t,t1950,mx,my,mz,1)

          bx2 = -mx/emratp1
          by2 = -my/emratp1
          bz2 = -mz/emratp1

          r1=r1+bx2
          r2=r2+by2
          r3=r3+bz2


          if(topof.eq.1)then
!             convert to topocentric position
              ut1 = t - 48./86400.
              call topo(phi,lon,ut1,tx,ty,tz)
!             convert topocentric position to 1950 coordinates
              call preces(t,t1950,tx,ty,tz,1)
              r1=r1+tx
              r2=r2+ty
              r3=r3+tz
          endif

      endif ! planet .eq. 1 .and. moon .eq. 1


!     call rotate(r1,r2,r3,-23.45)

      return


c
c write out positions and velocities for a time near time t
9000  tt=tcent+float(isteps)*h
      write(c_filename,9001)c_file1,int(tt - 0.5)
9001  format(a7,i7)

      open(2,file='/home/fab/albers/ast/planets/'
     1            //c_filename//'_new',status='new')

      write(6,9100)tt
      write(2,9100)tt
9100  format(//1x,f25.10/)
9101  format(1x,3d26.19)

      do j=1,nplan
          write(2,9101)x(j,apoint),y(j,apoint),z(j,apoint)
          write(6,9101)x(j,apoint),y(j,apoint),z(j,apoint)

      enddo ! j

      do j=1,nplan
          xvel=x(j,apoint-idr1)-x(j,apoint-idr1-idr1)
          yvel=y(j,apoint-idr1)-y(j,apoint-idr1-idr1)
          zvel=z(j,apoint-idr1)-z(j,apoint-idr1-idr1)

          do i=1,degree
              xvel=xvel+b(i)*xdd(j,apoint-i*idr1)
              yvel=yvel+b(i)*ydd(j,apoint-i*idr1)
              zvel=zvel+b(i)*zdd(j,apoint-i*idr1)
          enddo ! i


          xvel=xvel/h*float(idr1)
          yvel=yvel/h*float(idr1)
          zvel=zvel/h*float(idr1)

          write(2,9101)xvel,yvel,zvel
          write(6,9101)xvel,yvel,zvel
      enddo ! j

      do i=2,nplan
          call magnitude(i,0,x(1,apoint),y(1,apoint),z(1,apoint)
     .    ,x(i,apoint),y(i,apoint),z(i,apoint),value,dum)
          write(6,9600)i,value
9600      format(1x,i4,f15.8)
      enddo ! i

      close(2)


9999  return
      end

