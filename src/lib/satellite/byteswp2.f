      function byteswp2(vaxint)

c***********************************************************************
c  purpose:  perform the functional equivalent of a byte swap of a vax 
c  two byte integer
c
c  method:  the bit strings contained in the bytes of the input integer
c  vaxint are parsed into the integer variables byte1,  and byte2 
c  these bits strings are used to compute the integer that 
c  they represent.  this method is used rather than just swapping bytes
c  to make the program system independent.  no knowledge of how the
c  target system stores its integers is needed.  one only needs to know 
c  how the vax system stores integers.  for a vax system the least
c  significant bits are stored in byte 1 and the most significant in byte
c  2--the reverse is true for most unix system. 
c  
c  references:
c  1.  final interface specification for the satellite data handling 
c  system communications network (sdhs-comnet), 01 february 1988
c  2.  afgwc/dons pv-wave program auto_convert.pro 
c  3.  vax fortran reference manual from the sdhs programmer library  
c***********************************************************************
      implicit none
      integer byteswp2
      integer vaxint
      integer byte1
      integer byte2

c***********************************************************************
c  pick off 8-bit bytes, using the standard fortran routine ibits, from
c  the input variable and store the information in an integer variables 
c  the most significant bits are stored in byte2, and the
c  least significant in byte1.  the input vax bit string has the least
c  significant bits in the first byte and the most significant in
c  the second byte.  since most unix sytem expect the bytes to be stored
c  in reverse order, a call to ibits to get bits 0 to 8 on a unix machine
c  means the routine goes to the 2th word to read the the bit string.  
c  becuase the the 2th byte on a vax machine has the most significant 
c  bits a call to ibits with the argument 0, 8 gets the most significant
c  bits of the integer.  sound complicated? it is. the bit string 
c  diagram given below may help.
c   
c        byte1   byte2   
c   vax  07-00   15-08    
c   unix 15-08   07-00 

c***********************************************************************

cc      byte2 = ibits(vaxint,0,8)
cc      byte1 = ibits(vaxint,8,8)
      byte1=0
      byte2=0
      byte2 = ibits(vaxint,16,8)
      byte1 = ibits(vaxint,24,8)
     

***********************************************************************
c  recontruct the integer by multipling the integer representation of an
c  8-bit byte by the appropriate power of 2 and adding the two bytes.
c  for example:  byte2 represents the most significant string of
c  8 bits and must be multiplied by 2**8 to obtain the correct integer
c  representation.   the  sum of byte1 and byte2 represents the 
c  integer.  note: if negative values are stored in 2's or 1's 
c  complement form then the integer computed will not be a correct
c  representation, but the bit string will be equal to the input vax
c  string with the bytes swapped.
c*********************************************************************** 
      if(byte2.lt.128) then
         byteswp2 = byte1 + byte2*256 
      else
         write(6,*) 'warning: byteswp2 incountered value < 0'
     +             ,' the value of this translation is doubtful'
         byteswp2 = -byte1 + (128-byte2)*256
      endif

      return
      end
