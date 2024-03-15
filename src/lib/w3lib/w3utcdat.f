!-----------------------------------------------------------------------
      subroutine w3utcdat(idat)
!$$$   subprogram  documentation  block
!
! subprogram: w3utcdat       return the utc date and time
!   author: mark iredell     org: wp23       date: 98-01-05
!
! abstract: this subprogram returns the utc (greenwich) date and time
!   in the ncep absolute date and time data structure.
!
! program history log:
!   98-01-05  mark iredell
! 1999-04-28  gilbert         - added a patch to check for the proper
!                               utc offset.  needed until the ibm bug
!                               in date_and_time is fixed.  the patch
!                               can then be removed.  see comments in
!                               the section blocked with "&&&&&&&&&&&".
! 1999-08-12  gilbert         - changed so that czone variable is saved
!                               and the system call is only done for
!                               first invocation of this routine.
!
! usage:  call w3utcdat(idat)
!
!   output variables:
!     idat       integer (8) ncep absolute date and time
!                (year, month, day, time zone,
!                 hour, minute, second, millisecond)
!
! subprograms called:
!     date_and_time  fortran 90 system date intrinsic
!     iw3jdn         compute julian day number     
!     w3fs26         year, month, day from julian day number
!
! attributes:
!   language: fortran 90
!
!$$$
      integer idat(8)
      character cdate*8,ctime*10,czone*5
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  get local date and time but use the character time zone
      call date_and_time(cdate,ctime,czone,idat)
      read(czone,'(i5)') idat(4)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  convert to hours and minutes to utc time
!  and possibly adjust the date as well
      idat(6)=idat(6)-mod(idat(4),100)
      idat(5)=idat(5)-idat(4)/100
      idat(4)=0
      if(idat(6).lt.00) then
        idat(6)=idat(6)+60
        idat(5)=idat(5)-1
      elseif(idat(6).ge.60) then
        idat(6)=idat(6)-60
        idat(5)=idat(5)+1
      endif
      if(idat(5).lt.00) then
        idat(5)=idat(5)+24
        jldayn=iw3jdn(idat(1),idat(2),idat(3))-1
        call w3fs26(jldayn,idat(1),idat(2),idat(3),idaywk,idayyr)
      elseif(idat(5).ge.24) then
        idat(5)=idat(5)-24
        jldayn=iw3jdn(idat(1),idat(2),idat(3))+1
        call w3fs26(jldayn,idat(1),idat(2),idat(3),idaywk,idayyr)
      endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end
