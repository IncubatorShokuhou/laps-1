
      subroutine clock(h,c5_string)
      real*8 h,rhour
      character*5 c5_string
      integer hour,minute

!     hour is the input time in radians

      rhour=dmod(h*3.819718634d0+48.d0,24.d0)
      hour=idint(rhour)

      minute=nint((rhour-hour)*60.d0)

      if(minute .eq. 60)then
          minute = 0
          hour = hour + 1
      endif

      min1 = minute/10
      min2 = minute - min1 * 10

      write(c5_string,1)hour,min1,min2
1     format(i2,':',i1,i1)

      return
      end
