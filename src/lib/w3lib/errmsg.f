c-----------------------------------------------------------------------
      subroutine errmsg(cmsg)
c$$$  subprogram documentation block
c
c subprogram: errmsg         write a message to stderr
c   prgmmr: iredell          org: w/nmc23     date: 95-10-31
c
c abstract: write a message to stderr.
c
c program history log:
c   95-10-31  iredell
c
c usage:    call errmsg(cmsg)
c   input arguments:
c     cmsg         character*(*) message to write
c
c remarks: this is a machine-dependent subprogram.
c
c attributes:
c   language: fortran
c   machine:  cray
c
c$$$
      character*(*) cmsg
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      write(0,'(a)') cmsg
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
