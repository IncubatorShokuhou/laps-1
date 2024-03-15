      function byteswp4(vaxint)
c***********************************************************************
c  purpose:  perform the functional equivalent of a byte swap of a vax 
c  four byte integer
c
c  method:  the bit strings contained in the bytes of the input integer
c  vaxint are parsed into the integer variables byte1, byte2, byte3, 
c  and byte4.  these bits strings are used to compute the integer that 
c  they represent.  this method is used rather than just swapping bytes
c  to make the program system independent.  no knowledge of how the
c  target system stores its integers is needed.  one only needs to know 
c  how the vax system stores integers.  for a vax system the least
c  significant bits are stored in byte 1 and the most significant in byte
c  4--the reverse is true for most unix system.      
c 
c  references:
c  1.  final interface specification for the satellite data handling 
c  system communications network (sdhs-comnet), 01 february 1988
c  2.  afgwc/dons pv-wave program auto_convert.pro 
c  3.  vax fortran reference manual from the sdhs programmer library  
c***********************************************************************

      implicit none
      integer byteswp4
      integer vaxint
      integer byte1
      integer byte2
      integer byte3
      integer byte4

c***********************************************************************
c  pick off 8-bit bytes, using the standard fortran routine ibits, from
c  the input variable and store the information in integer variables 
c  the most significant bits are stored in byte4, the next most 
c  significant in byte3, the next in byte2, and the least significant in
c  in byte1.  the input vax bit string has the most significant bits in
c  the fourth byte, the next most significant in the third byte, the 
c  the next in the second byte and the least significant in the first
c  byte.  since most unix sytem expect the bytes to be stored in reverse
c  order, a call to ibits to get bits 0 to 8 on a unix machine means the
c  routine goes to the 4th word to read the the bit string.  becuase the 
c  the 4th byte on a vax machine has the most significant bits a call to 
c  ibits with the argument 0, 8 gets the most significant bits of the 
c  byte.  sound complicated? it is. the bit string diagram given below
c  may help.
c   
c        byte1   byte2   byte3   byte4
c   vax  07-00   15-08   23-16   31-24 
c   unix 31-24   23-16   15-08   07-00
c     
c***********************************************************************

      byte1=0
      byte2=0
      byte3=0
      byte4=0

      byte4 = ibits(vaxint,0,8)  
      byte3 = ibits(vaxint,8,8)
      byte2 = ibits(vaxint,16,8)
      byte1 = ibits(vaxint,24,8)

c***********************************************************************
c  recontruct the integer by multipling the integer representation of an
c  8-bit byte by the appropriate power of 2 and adding the four bytes.
c  for example:  byte2 represents the second least significant string of
c  8 bits and must be multiplied by 2**8 to obtain the correct integer
c  representation. byte3 must be multiplied by 2**16 and byte4 by 2**24
c  to obtain the proper integer representation.  their sum represents 
c  the integer.  note: if negative values are stored in 2's or 1's 
c  complement form then the integer computed will not be a correct
c  representation, but the bit string will be equal to the input vax
c  string with the bytes swapped.
c***********************************************************************
      
      byteswp4 = byte1 + byte2*256 + byte3*65536 + byte4*16777216

      return
      end

