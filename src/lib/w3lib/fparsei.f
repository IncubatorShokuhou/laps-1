c-----------------------------------------------------------------------
      subroutine fparsei(carg,marg,karg)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:  fparser       parse integers from a character string
c   prgmmr: iredell          org: np23        date:1998-09-03
c
c abstract: this subprogram extracts integers from a free-format
c   character string.  it is useful for parsing command arguments.
c
c program history log:
c 1998-09-03  iredell  
c
c usage:  call fparsei(carg,marg,karg)
c
c   input argument list:
c     carg     - character*(*) string of ascii digits to parse.
c                integers may be separated by a comma or by blanks.
c     marg     - integer maximum number of integers to parse.
c
c   output argument list:
c     karg     - integer (marg) numbers parsed.
c                (from 0 to marg values may be returned.)
c
c remarks:
c   to determine the actual number of integers found in the string,
c   karg should be set to fill values before the call to fparsei and
c   the number of non-fill values should be counted after the call.
c
c attributes:
c   language: fortran 90
c
c$$$
      character*(*) carg
      integer karg(marg)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      read(carg,*,iostat=ios) karg
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end
