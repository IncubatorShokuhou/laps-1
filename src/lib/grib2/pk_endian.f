      logical function pk_endian( )
c
c        march  2001    mattison/lawrence  gsc/mdl   original coding
c        october 2002   lawrence           whfs/ohd  
c                       fixed a bug that was preventing the
c                       determination of the correct operating 
c                       system.
c
c        purpose
c           returns a value of "true" if the machine that 
c           this routine is being run on is using 
c           big endian memory architecture.  this routine
c           will return a value of "false" if the machine
c           that this routine is being run on is using
c           little endian memory architecture.
c
c           on big endian systems (such as the hp-9000 workstations),
c           the last byte of an integer*4 value is the least 
c           significant byte.  on little endian systems (such as
c           intel pentium systems) the last byte of an integer*4
c           value is the most significant byte.
c
c           this routine functions as follows:
c
c           the integer*4 variable, "i", is set to "0".  this
c           zeros out each of the four bytes in "i".  this also
c           sets each of the characters in "letter" to "0" since
c           this character string is equivalenced to "i".  then,
c           the fourth character of "letter" is set to the character
c           with an ascii value of "1".  because of the equivalence,
c           this sets the fourth byte of "i" to contain a value
c           of "1".
c
c           on a big endian system, with i's last (least significant)
c           byte being "1", "i" is seen to have the value "1".  
c           however, on a little endian system with i's last
c           (most significant) byte being "1", "i" is seen to have 
c           the value of 2**24 (16777216).  thus, by testing on the
c           value of "i", this routine can determine which type
c           of hardware this routine is being run on.
c
c        data set use
c           none
c        variables
c           none
c
c             local variables
c                   i = contains the value of "1" to be used in
c                       determining the "endianess" of the 
c                       hardware that this routine is being
c                       run on.  (integer*4)
c              letter = equivalenced to "i".  this four-byte
c                       character string allows this routine
c                       to manipulate the individual byte
c                       values in "i".
c
c        non system subroutines called
c           none

      integer*4 i
      character*4 letter
      equivalence (letter, i)
c
      i = 0
      letter(4:4) = char(1)
c
      if (i.eq.1) then
         pk_endian=.true.
      else
         pk_endian=.false.
      endif
c
      return
      end
