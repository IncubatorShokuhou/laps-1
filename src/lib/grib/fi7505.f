      subroutine fi7505 (n,nbits)
c$$$  subprogram documentation block
c subprogram:    fi7505      determine number of bits to contain value
c   prgmmr: cavanaugh        org: w/nmc42    date: 93-06-23
c
c abstract: calculate number of bits to contain value n, with a
c            maximum of 32 bits.
c
c program history log:
c   93-06-23  cavanaugh
c   95-10-31  iredell     removed saves and prints
c
c usage:    call fi7505 (n,nbits)
c   input argument list:
c     n        - integer value
c
c   output argument list:      (including work arrays)
c     nbits    - number of bits to contain n
c
c remarks: subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: ibm vs fortran 77, cray cft77 fortran
c   machine:  hds, cray c916/256, y-mp8/64, y-mp el92/256
c
c$$$
      integer        n,nbits
      integer        ibits(31)
c
      data           ibits/1,3,7,15,31,63,127,255,511,1023,2047,
     *               4095,8191,16383,32767,65535,131071,262143,
     *               524287,1048575,2097151,4194303,8388607,
     *               16777215,33554431,67108863,134217727,268435455,
     *               536870911,1073741823,2147483647/
c  ----------------------------------------------------------------
c
      do 1000 nbits = 1, 31
          if (n.le.ibits(nbits)) then
              return
          end if
 1000 continue
      return
      end
