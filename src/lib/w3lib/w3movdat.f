!-----------------------------------------------------------------------
      subroutine w3movdat(rinc,idat,jdat)
!$$$   subprogram  documentation  block
!
! subprogram: w3movdat       return a date from a time interval and date
!   author: mark iredell     org: wp23       date: 98-01-05
!
! abstract: this subprogram returns the date and time that is a given
!   ncep relative time interval from an ncep absolute date and time.
!   the output is in the ncep absolute date and time data structure.
!
! program history log:
!   98-01-05  mark iredell
!
! usage:  call w3movdat(rinc,idat,jdat)
!
!   input variables:
!     rinc       real (5) ncep relative time interval
!                (days, hours, minutes, seconds, milliseconds)
!     idat       integer (8) ncep absolute date and time
!                (year, month, day, time zone,
!                 hour, minute, second, millisecond)
!
!   output variables:
!     jdat       integer (8) ncep absolute date and time
!                (year, month, day, time zone,
!                 hour, minute, second, millisecond)
!                (jdat is later than idat if time interval is positive.)
!
! subprograms called:
!     iw3jdn         compute julian day number     
!     w3fs26         year, month, day from julian day number
!     w3reddat       reduce a time interval to a canonical form
!
! attributes:
!   language: fortran 90
!
!$$$
      real rinc(5)
      integer idat(8),jdat(8)
      real rinc1(5),rinc2(5)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  add the interval to the input time of day and put into reduced form
!  and then compute new date using julian day arithmetic.
      rinc1(1)=rinc(1)
      rinc1(2:5)=rinc(2:5)+idat(5:8)
      call w3reddat(-1,rinc1,rinc2)
      jldayn=iw3jdn(idat(1),idat(2),idat(3))+nint(rinc2(1))
      call w3fs26(jldayn,jdat(1),jdat(2),jdat(3),jdow,jdoy)
      jdat(4)=idat(4)
      jdat(5:8)=nint(rinc2(2:5))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end
