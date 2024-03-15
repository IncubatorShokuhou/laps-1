      subroutine w3fs21(idate, nmin)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:   w3fs21       number of minutes since jan 1, 1978
c   prgmmr: rejones          org: nmc421     date: 89-07-17
c
c abstract: calculates the number of minutes since 0000,
c   1 january 1978.
c
c program history log:
c   84-06-21  a. desmarais
c   89-07-14  r.e.jones    convert to cyber 205 fortran 200,
c                          change logic so it will work in
c                          21 century.
c   89-11-02  r.e.jones    convert to cray cft77 fortran
c
c usage:    call w3fs21 (idate, nmin)
c   input argument list:
c     idate    - integer  size 5 array containing year of century,
c                month, day, hour and minute.  idate(1) may be
c                a two digit year or 4. if 2 digits and ge than 78
c                1900 is added to it. if lt 78 then 2000 is added
c                to it. if 4 digits the subroutine will work
c                correctly to the year 3300 a.d.
c
c   output argument list:
c     nmin     - integer number of minutes since 1 january 1978
c
c   subprograms called:
c     library:
c       w3lib    - iw3jdn
c
c attributes:
c   language: cray cft77 fortran
c   machine:  cray y-mp8/832
c
c$$$
c
      integer  idate(5)
      integer  nmin
      integer  jdn78
c
      data  jdn78 / 2443510 /
c
c***   idate(1)       year of century
c***   idate(2)       month of year
c***   idate(3)       day of month
c***   idate(4)       hour of day
c***   idate(5)       minute of hour
c
      nmin  = 0
c
      iyear = idate(1)
c
      if (iyear.le.99) then
        if (iyear.lt.78) then
          iyear = iyear + 2000
        else
          iyear = iyear + 1900
        endif
      endif
c
c     compute julian day number from year, month, day
c
      ijdn  = iw3jdn(iyear,idate(2),idate(3))
c
c     subtract julian day number of jan 1,1978 to get the
c     number of days between dates
c
      ndays = ijdn - jdn78
c
c***  number of minutes
c
      nmin = ndays * 1440 + idate(4) * 60 + idate(5)
c
      return
      end
