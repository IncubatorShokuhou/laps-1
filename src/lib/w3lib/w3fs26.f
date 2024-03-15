       subroutine w3fs26(jldayn,iyear,month,iday,idaywk,idayyr)
c$$$   subprogram  documentation  block
c
c subprogram: w3fs26         year, month, day from julian day number
c   author: jones,r.e.       org: w342       date: 87-03-29
c
c abstract: computes year (4 digits), month, day, day of week, day
c   of year from julian day number. this subroutine will work
c   from 1583 a.d. to 3300 a.d.
c
c program history log:
c   87-03-29  r.e.jones
c   89-10-25  r.e.jones   convert to cray cft77 fortran
c
c usage:  call w3fs26(jldayn,iyear,month,iday,idaywk,idayyr)
c
c   input variables:
c     names  interface description of variables and types
c     ------ --------- -----------------------------------------------
c     jldayn arg list  integer   julian day number
c
c   output variables:
c     names  interface description of variables and types
c     ------ --------- -----------------------------------------------
c     iyear  arg list  integer   year  (4 digits)
c     month  arg list  integer   month
c     iday   arg list  integer   day
c     idaywk arg list  integer   day of week (1 is sunday, 7 is sat)
c     idayyr arg list  integer   day of year (1 to 366)
c
c   remarks: a julian day number can be computed by using one of the
c     following statement functions. a day of week can be computed
c     from the julian day number. a day of year can be computed from
c     a julian day number and year.
c
c      iyear (4 digits)
c
c      jdn(iyear,month,iday) = iday - 32075
c    &            + 1461 * (iyear + 4800 + (month - 14) / 12) / 4
c    &            + 367 * (month - 2 - (month -14) / 12 * 12) / 12
c    &            - 3 * ((iyear + 4900 + (month - 14) / 12) / 100) / 4
c
c      iyr (4 digits) , idyr(1-366) day of year
c
c      julian(iyr,idyr) = -31739 + 1461 * (iyr + 4799) / 4
c    &                    -3 * ((iyr + 4899) / 100) / 4 + idyr
c
c      day of week from julian day number, 1 is sunday, 7 is saturday.
c
c      jdaywk(jldayn) = mod((jldayn + 1),7) + 1
c
c      day of year from julian day number and 4 digit year.
c
c      jdayyr(jldayn,iyear) = jldayn -
c     &  (-31739+1461*(iyear+4799)/4-3*((iyear+4899)/100)/4)
c
c      the first function was in a letter to the editor communications
c      of the acm  volume 11 / number 10 / october, 1968. the 2nd
c      function was derived from the first. this subroutine was also
c      included in the same letter. julian day number 1 is
c      jan 1,4713 b.c. a julian day number can be used to replace a
c      day of century, this will take care of the date problem in
c      the year 2000, or reduce program changes to one line change
c      of 1900 to 2000. julian day numbers can be used for finding
c      record numbers in an archive or day of week, or day of year.
c
c attributes:
c   language: cray cft77 fortran
c   machine:  cray y-mp8/864
c
c$$$
c
       l      = jldayn + 68569
       n      = 4 * l / 146097
       l      = l - (146097 * n + 3) / 4
       i      = 4000 * (l + 1) / 1461001
       l      = l - 1461 * i / 4 + 31
       j      = 80 * l / 2447
       iday   = l - 2447 * j / 80
       l      = j / 11
       month  = j + 2 - 12 * l
       iyear  = 100 * (n - 49) + i + l
       idaywk = mod((jldayn + 1),7) + 1
       idayyr = jldayn -
     &  (-31739 +1461 * (iyear+4799) / 4 - 3 * ((iyear+4899)/100)/4)
       return
       end
