
      subroutine cjymd(t,iyear,month,date)
      implicit real*8(a,b,c,d,e,f,g,h,o,p,q,r,s,t,u,v,w,x,y,z)
      dimension k(49)
      data k/0,31,60,91,121,152,182,213,244,274,305,335,366,397,425,456
     .,486,517,547,578,609,639,670,700,731,762,790,821,851,882,912,943
     .,974,1004,1035,1065,1096,1127,1155,1186,1216,1247,1277,1308,1339
     .,1369,1400,1430,1461/

      time=t+0.5d0
      i=idint(time)

      if(i .gt. 2299160)then
          i=i+10
          lp1=idint((dfloat(i)-2305517.d0)/36524.25d0)
          i=i+lp1-lp1/4
      endif


      lp=i/1461
      j=i-lp*1461
      month=j/30+1

      if(k(month) .gt. j)then
          month=month-1
      endif

      date=j-k(month)+1+time-int(time)
      iyear=-4712+lp*4+(month-1)/12
      month=month-12*((month-1)/12)

      if(iyear.le.0)then
          iyear=iyear-1
      endif

11    return
      end
