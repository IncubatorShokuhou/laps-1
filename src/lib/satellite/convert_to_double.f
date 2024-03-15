      subroutine convert_to_double(integrin1,integrin2,realout)

c***********************************************************************
c  purpose: properly interpret the bit string that represents a vax 
c  double precision number to obtain the correct number on the target 
c  system.   
c
c  method: the bit string that represents the vax double precision 
c  number is passed into the routine as two integers.  note the bytes 
c  have not been swapped.  the data is assumed to be in big endian 
c  ieee format.  the 4 bytes of the input integers are swaped, using 
c  mvbits, so the result is a little endian representation stored in an
c  integer array of dimension two.  this bit string is then transfered 
c  to a  variable of type double precision (real*8) using the routine 
c  transfer.  once the bit string is transfered it is properly 
c  interpreted as a double precision number. 
c
c    
c  references:
c  1.  final interface specification for the satellite data handling 
c  system communications network (sdhs-comnet), 01 february 1988
c  2.  afgwc/dons pv-wave program auto_convert.pro 
c  3.  vax fortran reference manual from the sdhs programmer library  
c***********************************************************************

      integer   integrin1,integrin2 !input integers that contain the 
                                    !the bit strings
      real*8 realout     !the output double precision number
      integer,dimension(2) ::  integrtemp !a working integer array 
                                          !that holds the double
                                          !precision bit string after
                                          !swaping of the words and 
                                          
c***********************************************************************
c  on big endian systems the first word contains the least significant
c  bits.  on a little endian system the second (high address 32 bit word)
c  contains the least significant bits.  so in addtion to swaping bytes
c  within a word the word themselves must be swapped.  both actions are
c  accomplished with the use of mvbits
c***********************************************************************                                          !bytes

c     call mvbits(integrin1,24,8,integrtemp(2),0)
c     call mvbits(integrin1,16,8,integrtemp(2),8)
c     call mvbits(integrin1,8,8,integrtemp(2),16)
c     call mvbits(integrin1,0,8,integrtemp(2),24)
c     call mvbits(integrin2,24,8,integrtemp(1),0)
c     call mvbits(integrin2,16,8,integrtemp(1),8)
c     call mvbits(integrin2,8,8,integrtemp(1),16)
c     call mvbits(integrin2,0,8,integrtemp(1),24)

      integrtemp(2)=integrin1
      integrtemp(1)=integrin2
      realout=transfer(integrtemp,realout)
 
      return
      end

