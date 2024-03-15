!-----------------------------------------------------------------------
      subroutine w3doxdat(idat,jdow,jdoy,jday)
!$$$   subprogram  documentation  block
!
! subprogram: w3doxdat       return week day, year day, and julian day
!   author: mark iredell     org: wp23       date: 98-01-05
!
! abstract: this subprogram returns the integer day of week, the day
!   of year, and julian day given an ncep absolute date and time.
!
! program history log:
!   98-01-05  mark iredell
!
! usage:  call w3doxdat(idat,jdow,jdoy,jday)
!
!   input variables:
!     idat       integer (8) ncep absolute date and time
!                (year, month, day, time zone,
!                 hour, minute, second, millisecond)
!
!   output variables:
!     jdow       integer day of week (1-7, where 1 is sunday)
!     jdoy       integer day of year (1-366, where 1 is january 1)
!     jday       integer julian day (day number from jan. 1,4713 b.c.)
!
! subprograms called:
!     iw3jdn         compute julian day number     
!     w3fs26         year, month, day from julian day number
!
! attributes:
!   language: fortran 90
!
!$$$
      integer idat(8)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  get julian day and then get day of week and day of year
      jday=iw3jdn(idat(1),idat(2),idat(3))
      call w3fs26(jday,jy,jm,jd,jdow,jdoy)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end
