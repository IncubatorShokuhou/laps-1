      subroutine w3reddat(it,rinc,dinc)
!$$$   subprogram  documentation  block
!
! subprogram: w3reddat       reduce a time interval to a canonical form
!   author: mark iredell     org: wp23       date: 98-01-05
!
! abstract: this subprogram reduces an ncep relative time interval
!   into one of seven canonical forms, depending on the input it value.
!
!   first reduced format type (it=-1):
!        rinc(1) is an arbitrary integer.
!        rinc(2) is an integer between 00 and 23, inclusive.
!        rinc(3) is an integer between 00 and 59, inclusive.
!        rinc(4) is an integer between 00 and 59, inclusive.
!        rinc(5) is an integer between 000 and 999, inclusive.
!      if rinc(1) is negative, then the time interval is negative.
!    
!   second reduced format type (it=0):
!      if the time interval is not negative, then the format is:
!        rinc(1) is zero or a positive integer. 
!        rinc(2) is an integer between 00 and 23, inclusive.
!        rinc(3) is an integer between 00 and 59, inclusive.
!        rinc(4) is an integer between 00 and 59, inclusive.
!        rinc(5) is an integer between 000 and 999, inclusive.
!      otherwise if the time interval is negative, then the format is:
!        rinc(1) is zero or a negative integer. 
!        rinc(2) is an integer between 00 and -23, inclusive.
!        rinc(3) is an integer between 00 and -59, inclusive.
!        rinc(4) is an integer between 00 and -59, inclusive.
!        rinc(5) is an integer between 000 and -999, inclusive.
!    
!   days format type (it=1):
!        rinc(1) is arbitrary.
!        rinc(2) is zero.
!        rinc(3) is zero.
!        rinc(4) is zero.
!        rinc(5) is zero.
!    
!   hours format type (it=2):
!        rinc(1) is zero.
!        rinc(2) is arbitrary.
!        rinc(3) is zero.
!        rinc(4) is zero.
!        rinc(5) is zero.
!      (this format should not express time intervals longer than 300 years.)
!    
!   minutes format type (it=3):
!        rinc(1) is zero.
!        rinc(2) is zero.
!        rinc(3) is arbitrary.
!        rinc(4) is zero.
!        rinc(5) is zero.
!      (this format should not express time intervals longer than five years.)
!    
!   seconds format type (it=4):
!        rinc(1) is zero.
!        rinc(2) is zero.
!        rinc(3) is zero.
!        rinc(4) is arbitrary.
!        rinc(5) is zero.
!      (this format should not express time intervals longer than one month.)
!    
!   milliseconds format type (it=5):
!        rinc(1) is zero.
!        rinc(2) is zero.
!        rinc(3) is zero.
!        rinc(4) is zero.
!        rinc(5) is arbitrary.
!     (this format should not express time intervals longer than one hour.)
!
! program history log:
!   98-01-05  mark iredell
!
! usage:  call w3reddat(it,rinc,dinc)
!
!   input variables:
!     it         integer relative time interval format type
!                (-1 for first reduced type (hours always positive),
!                 0 for second reduced type (hours can be negative),
!                 1 for days only, 2 for hours only, 3 for minutes only,
!                 4 for seconds only, 5 for milliseconds only)
!     rinc       real (5) ncep relative time interval
!                (days, hours, minutes, seconds, milliseconds)
!
!   output variables:
!     dinc       real (5) ncep relative time interval
!                (days, hours, minutes, seconds, milliseconds)
!
! subprograms called:
!
! attributes:
!   language: fortran 90
!
!$$$
      real rinc(5),dinc(5)
!  parameters for number of units in a day
!  and number of milliseconds in a unit
!  and number of next smaller units in a unit, respectively
      integer,dimension(5),parameter:: itd=(/1,24,1440,86400,86400000/),
     &                                 itm=itd(5)/itd
      integer,dimension(4),parameter:: itn=itd(2:5)/itd(1:4)
      integer,parameter:: np=16
      integer iinc(4),jinc(5),kinc(5)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  first reduce to the first reduced form
      iinc=floor(rinc(1:4))
!  convert all positive fractional parts to milliseconds
!  and determine canonical milliseconds
      jinc(5)=nint(dot_product(rinc(1:4)-iinc,real(itm(1:4)))+rinc(5))
      kinc(5)=modulo(jinc(5),itn(4))
!  convert remainder to seconds and determine canonical seconds
      jinc(4)=iinc(4)+(jinc(5)-kinc(5))/itn(4)
      kinc(4)=modulo(jinc(4),itn(3))
!  convert remainder to minutes and determine canonical minutes
      jinc(3)=iinc(3)+(jinc(4)-kinc(4))/itn(3)
      kinc(3)=modulo(jinc(3),itn(2))
!  convert remainder to hours and determine canonical hours
      jinc(2)=iinc(2)+(jinc(3)-kinc(3))/itn(2)
      kinc(2)=modulo(jinc(2),itn(1))
!  convert remainder to days and compute milliseconds of the day
      kinc(1)=iinc(1)+(jinc(2)-kinc(2))/itn(1)
      ms=dot_product(kinc(2:5),itm(2:5))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  next reduce to either single value canonical form
!  or to one of the two reduced forms
      if(it.ge.1.and.it.le.5) then
!  ensure that exact multiples of 1./np are expressed exactly
!  (other fractions may have precision errors)
        rp=(np*ms)/itm(it)+mod(np*ms,itm(it))/real(itm(it))
        dinc=0
        dinc(it)=real(kinc(1))*itd(it)+rp/np
      else
!  the reduced form is done except the second reduced form is modified
!  for negative time intervals with fractional days
        dinc=kinc
        if(it.eq.0.and.kinc(1).lt.0.and.ms.gt.0) then
          dinc(1)=dinc(1)+1
          dinc(2:5)=mod(ms-itm(1),itm(1:4))/itm(2:5)
        endif
      endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end
