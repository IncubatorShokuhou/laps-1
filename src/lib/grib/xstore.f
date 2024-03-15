      subroutine xstore(cout,con,mwords)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    xstore      stores a constant value into an array
c   prgmmr: keyser           org: w/nmc22    date: 07-02-92
c
c abstract: stores an 4-byte (fullword) value through consecutive
c   storage locations.  (moving is accomplished with a do loop.)
c
c program history log:
c   92-07-02  d. a. keyser (w/nmc22)
c   93-03-29  r.e.jones   add save statement
c   93-04-15  r.e.jones   changes for microsoft fortran 5.0
c
c usage:    call xstore(cout,con,mwords)
c   input argument list:
c     con      - constant to be stored into "mwords" consecutive
c                fullwords beginning with "cout" array
c     mwords   - number of fullwords in "cout" array to store "con";
c                must be .gt. zero (not checked for this)
c
c   output argument list:      (including work arrays)
c     cout     - starting address for array of "mwords" fullwords
c                set to the contents of the value "con"
c
c remarks: the version of this subroutine on the hds common library
c   is nas-specific subr. written in assembly lang. to allow fast
c   computation time.  subr. placed in cray w3lib to allow codes to
c   compile on both the hds and cray machines.
c
c attributes:    
c   language: microsoft fortran 5.0 optimizing compiler
c   machine:  ibm pc, at, ps/2, 386, 486, 586, clones.
c
c$$$
c
      dimension  cout(*)
c
      save
c
      do 1000  i = 1,mwords
        cout(i) = con
1000  continue   
c
      return     
      end

