      subroutine bcdtoint (input1, input2, output)
c
c *****************************************************************************
c     subroutine to convert goes/gvar "bcd time" into time value which noaa's
c     "gimloc" software requires.
c  
c     original:	bruce h. thomas, the aerospace corporation
c		defense meteorological satellite program (dmsp)
c		environmental applications center, march 1997
c
c     note:     the second half of the "bcd" date/time (oadata(13)) is empty!!!
c		goes/gvar spec indicates hours in 1's, min, sec, milli sec
c		should be present.  work later if required!!! 
c *****************************************************************************
c
      implicit none
c
      logical print_flag
      integer iy1000,iy100,iy10,iy1,idoy100,idoy10,idoy1,ihr10
      integer input1,input2,ibyte1,ibyte2,ibyte3,ibyte4,ihb1,ihb2
      integer iyear,iday,ihour,ihr1
      real*8 output
c
      data print_flag / .false. /
c
      ibyte1 = ibits(input1,00,8)
      ihr10 = ibits(ibyte1,0,4)
      idoy1 = ibits(ibyte1,4,4)
c
      ibyte2 = ibits(input1,08,8)
      idoy10 = ibits(ibyte2,0,4)
      idoy100 = ibits(ibyte2,4,4)
c
      ibyte3 = ibits(input1,16,8)
      iy1 = ibits(ibyte3,0,4)
      iy10 = ibits(ibyte3,4,4)
c
      ibyte4 = ibits(input1,24,8)
      iy100 = ibits(ibyte4,0,4)
      iy1000 = ibits(ibyte4,4,4)
c
      iyear = iy1000 * 1000 + iy100 * 100 + iy10 * 10 + iy1
      iday = idoy100 * 100 + idoy10 * 10 + idoy1
      ihour = ihr10 * 10 + ihr1
c
      if (print_flag) then
        print *,'bcd_to_int: year = ',iyear
        print *,'bcd_to_int: day of year = ',iday
        print *,'bcd_to_int: hour of day = ',ihour
      endif
c
      ibyte1 = ibits(input2,00,8)
      ihb1 = ibits(ibyte1,0,4)
      ihb2 = ibits(ibyte1,4,4)
c
      if (print_flag) then
        print *,'bcd_to_int: oadata(13) information.... '
        print *,'bcd_to_int: input2 byte 1 = ',ibyte1
        print *,'bcd_to_int: half byte 1 = ',ihb1
        print *,'bcd_to_int: half byte 2 = ',ihb2
      endif
c
      ibyte2 = ibits(input2,08,8)
      ihb1 = ibits(ibyte2,0,4)
      ihb2 = ibits(ibyte2,4,4)
c
      if (print_flag) then
        print *,'bcd_to_int: input2 byte 2 = ',ibyte2
        print *,'bcd_to_int: half byte 1 = ',ihb1
        print *,'bcd_to_int: half byte 2 = ',ihb2
      endif
c
      ibyte3 = ibits(input2,16,8)
      ihb1 = ibits(ibyte3,0,4)
      ihb2 = ibits(ibyte3,4,4)
c
      if (print_flag) then
        print *,'bcd_to_int: half byte 1 = ',ihb1
        print *,'bcd_to_int: half byte 2 = ',ihb2
        print *,'bcd_to_int: input2 byte 3 = ',ibyte3
      endif
c
      ibyte4 = ibits(input2,24,8)
      ihb1 = ibits(ibyte4,0,4)
      ihb2 = ibits(ibyte4,4,4)
c
      if (print_flag) then
        print *,'bcd_to_int: input2 byte 4 = ',ibyte4
        print *,'bcd_to_int: half byte 1 = ',ihb1
        print *,'bcd_to_int: half byte 2 = ',ihb2
      endif
c
      return
      end
