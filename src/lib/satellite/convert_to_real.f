      subroutine convert_to_real(integrin,realout)

c***********************************************************************
c  purpose:  properly interpret the bit string for a input
c  floating point number to obtain the correct floating point number on 
c  the target system.   
c
c  method: the bit string that represents the vax floating point number is
c  passed into the routine as an integer.  note the bytes have not been 
c  swapped.  the data is assumed to be in big endian ieee format.  the 
c  4 bytes of the input integer are swaped, using mvbits, so the result 
c  is a little endian representation.  this bit string is then transfered
c  to a  variable of type real using the function transfer.  once the bit 
c  string is transfered it is properly interpreted as a real number. 
c
c  references:
c  1.  final interface specification for the satellite data handling 
c  system communications network (sdhs-comnet), 01 february 1988
c  2.  afgwc/dons pv-wave program auto_convert.pro 
c  3.  vax fortran reference manual from the sdhs programmer library  
c***********************************************************************

      integer   integrin   !the input integer that contains the bit string
      real      realout    !the output floating point number
c     integer*4 integrtemp  !a temporary integer that contains the swapped
                            !bytes of the input interger integrin 

c***********************************************************************
c  call the intrinsic routine mvbits to swap bytes 4 to 1, 1 to 4, 2 to 3
c  and 3 to 2.  call the intrinsic function to transfer the bit string
c  from the integer where the bytes have been swaped to a variable of 
c  type real.  
c***********************************************************************

c     call mvbits(integrin,24,8,integrtemp,0)
c     call mvbits(integrin,16,8,integrtemp,8)
c     call mvbits(integrin,8,8,integrtemp,16) 
c     call mvbits(integrin,0,8,integrtemp,24)
c     realout=transfer(integrtemp,realout)
      realout=transfer(integrin,realout)

      return
      end

