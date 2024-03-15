      subroutine w3fi01(lw)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    w3fi01      determines machine word length in bytes
c   prgmmr: keyser           org: w/nmc22    date: 06-29-92
c
c abstract: determines the number of bytes in a full word for the
c   particular machine (ibm or cray).
c
c program history log:
c   92-01-10  r. kistler (w/nmc23)
c   92-05-22  d. a. keyser -- docblocked/commented
c   95-10-31  iredell     removed saves and prints
c
c usage:    call w3fi01(lw)
c   output argument list:      (including work arrays)
c     lw       - machine word length in bytes
c
c remarks: subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: fortran 77
c   machine:  cray, workstations
c
c$$$
c
      character*8  ctest1,ctest2
      character*4  cprint(2)
c
      integer      itest1,itest2
c
      equivalence  (ctest1,itest1),(ctest2,itest2)
c
      data  ctest1/'12345678'/
c
      itest2 = itest1
      if (ctest1 .eq. ctest2) then
        lw = 8
      else
        lw = 4
      end if
      return
      end
