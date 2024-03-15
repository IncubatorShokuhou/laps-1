!-----------------------------------------------------------------------
      subroutine w3locdat(idat)
!$$$   subprogram  documentation  block
!
! subprogram: w3locdat       return the local date and time
!   author: mark iredell     org: wp23       date: 98-01-05
!
! abstract: this subprogram returns the local date and time
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
! usage:  call w3locdat(idat)
!
!   output variables:
!     idat       integer (8) ncep absolute date and time
!                (year, month, day, time zone,
!                 hour, minute, second, millisecond)
!
! subprograms called:
!     date_and_time  fortran 90 system date intrinsic
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
      end
