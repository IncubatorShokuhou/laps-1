!-----------------------------------------------------------------------
      subroutine w3difdat(jdat,idat,it,rinc)
!$$$   subprogram  documentation  block
!
! subprogram: w3difdat       return a time interval between two dates
!   author: mark iredell     org: wp23       date: 98-01-05
!
! abstract: this subprogram returns the elapsed time interval from
!   an ncep absolute date and time given in the second argument until
!   an ncep absolute date and time given in the first argument.
!   the output time interval is in one of seven canonical forms
!   of the ncep relative time interval data structure.
!
! program history log:
!   98-01-05  mark iredell
!
! usage:  call w3difdat(jdat,idat,it,rinc)
!
!   input variables:
!     jdat       integer (8) ncep absolute date and time
!                (year, month, day, time zone,
!                 hour, minute, second, millisecond)
!     idat       integer (8) ncep absolute date and time
!                (year, month, day, time zone,
!                 hour, minute, second, millisecond)
!     it         integer relative time interval format type
!                (-1 for first reduced type (hours always positive),
!                 0 for second reduced type (hours can be negative),
!                 1 for days only, 2 for hours only, 3 for minutes only,
!                 4 for seconds only, 5 for milliseconds only)
!
!   output variables:
!     rinc       real (5) ncep relative time interval
!                (days, hours, minutes, seconds, milliseconds)
!                (time interval is positive if jdat is later than idat.)
!
! subprograms called:
!     iw3jdn         compute julian day number     
!     w3reddat       reduce a time interval to a canonical form
!
! attributes:
!   language: fortran 90
!
!$$$
      integer jdat(8),idat(8)
      real rinc(5)
      real rinc1(5)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  difference the days and time and put into canonical form
      rinc1(1)=iw3jdn(jdat(1),jdat(2),jdat(3))-
     &         iw3jdn(idat(1),idat(2),idat(3))
      rinc1(2:5)=jdat(5:8)-idat(5:8)
      call w3reddat(it,rinc1,rinc)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end
