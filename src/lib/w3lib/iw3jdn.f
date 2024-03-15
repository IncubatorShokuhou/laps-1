       function iw3jdn(iyear,month,iday)
c$$$   subprogram  documentation  block
c
c subprogram: iw3jdn         compute julian day number
c   author: jones,r.e.       org: w342       date: 87-03-29
c
c abstract: computes julian day number from year (4 digits), month,
c   and day. iw3jdn is valid for years 1583 a.d. to 3300 a.d.
c   julian day number can be used to compute day of week, day of
c   year, record numbers in an archive, replace day of century,
c   find the number of days between two dates.
c
c program history log:
c   87-03-29  r.e.jones
c   89-10-25  r.e.jones   convert to cray cft77 fortran
c
c usage:   ii = iw3jdn(iyear,month,iday)
c
c   input variables:
c     names  interface description of variables and types
c     ------ --------- -----------------------------------------------
c     iyear  arg list  integer   year           ( 4 digits)
c     month  arg list  integer   month of year   (1 - 12)
c     iday   arg list  integer   day of month    (1 - 31)
c
c   output variables:
c     names  interface description of variables and types
c     ------ --------- -----------------------------------------------
c     iw3jdn funtion   integer   julian day number
c                      jan. 1,1960 is julian day number 2436935
c                      jan. 1,1987 is julian day number 2446797
c
c   remarks: julian period was devised by joseph scaliger in 1582.
c     julian day number #1 started on jan. 1,4713 b.c. three major
c     chronological cycles begin on the same day. a 28-year solar
c     cycle, a 19-year luner cycle, a 15-year indiction cycle, used
c     in ancient rome to regulate taxes. it will take 7980 years
c     to complete the period, the product of 28, 19, and 15.
c     scaliger named the period, date, and number after his father
c     julius (not after the julian calendar). this seems to have
c     caused a lot of confusion in text books. scaliger name is
c     spelled three different ways. julian date and julian day
c     number are interchanged. a julian date is used by astronomers
c     to compute accurate time, it has a fraction. when truncated to
c     an integer it is called an julian day number. this function
c     was in a letter to the editor of the communications of the acm
c     volume 11 / number 10 / october 1968. the julian day number
c     can be converted to a year, month, day, day of week, day of
c     year by calling subroutine w3fs26.
c
c attributes:
c   language: cray cft77 fortran
c   machine:  cray y-mp8/864, cray y-mp el2/256
c
c$$$
c
       iw3jdn  =    iday - 32075
     &            + 1461 * (iyear + 4800 + (month - 14) / 12) / 4
     &            + 367 * (month - 2 - (month -14) / 12 * 12) / 12
     &            - 3 * ((iyear + 4900 + (month - 14) / 12) / 100) / 4
       return
       end
