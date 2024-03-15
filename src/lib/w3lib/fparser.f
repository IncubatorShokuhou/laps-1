c-----------------------------------------------------------------------
      subroutine fparser(carg,marg,rarg)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:  fparser       parse real numbers from a character string
c   prgmmr: iredell          org: np23        date:1998-09-03
c
c abstract: this subprogram extracts real numbers from a free-format
c   character string.  it is useful for parsing command arguments.
c
c program history log:
c 1998-09-03  iredell  
c
c usage:  call fparser(carg,marg,rarg)
c
c   input argument list:
c     carg     - character*(*) string of ascii digits to parse.
c                real numbers may be separated by a comma or by blanks.
c     marg     - integer maximum number of real numbers to parse.
c
c   output argument list:
c     rarg     - real (marg) numbers parsed.
c                (from 0 to marg values may be returned.)
c
c remarks:
c   to determine the actual number of real numbers found in the string,
c   rarg should be set to fill values before the call to fparser and
c   the number of non-fill values should be counted after the call.
c
c attributes:
c   language: fortran 90
c
c$$$
      character*(*) carg
      real rarg(marg)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      read(carg,*,iostat=ios) rarg
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end
