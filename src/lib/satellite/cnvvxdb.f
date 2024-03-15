      subroutine cnvvxdb(integrin1,integrin2,realout)

c***********************************************************************
c  purpose: properly interpret the bit string that represents a vax 
c  double precision number to obtain the correct number on the target 
c  system.   
c
c  method: the bit string that represents the vax double precision number
c  is passed into the routine as two integers.  note the bytes have not
c  been swapped.  the sign bit, the bits that represent the exponent or 
c  characteristic, and the bits that represent the fraction are picked off
c  using the standard fortran function ibits.  the following table show how
c  the bits and what they represent (sign, exponent etc) map to most unix 
c  systems.
c    
c                        first integer   
c          byte1       byte2       byte3     byte4
c   vax    07 06-00    15 14-8     23-16     31-24 
c          e1 f7-f1     s e8-e2    f15-f8    f23-f16
c   unix   31 30-24    23 22-16    15-08     07-00
c
c                        second integer   
c          byte1      byte2        byte3     byte4
c   vax    07-00      15-08        23-16     31-24
c          f31-24     f39-f32      f47-40    f55-48  
c   unix   31-24      23-16        15-08     07-00
c
c  here f1 stands for the first bit that is used to represent the 
c  fraction, f2 for the second bit used to represent the fraction, etc., 
c  e1 represents the first exponent bit, e2 the second, and s represents 
c  the sign bit.  a vax floating point number is represented by :
c
c  floating point number = -1**signbit * fraction * 2 **(exponent-128)
c    
c  !!!!!!note: to save space the vax system assumes the most significant
c              fraction bit is set.  this fact is taken into account when
c              computing the fraction.  also, as appears to be standard
c              practice, the exponent is biased by a fix number(128 in 
c              this case). 
c    
c  references:
c  1.  final interface specification for the satellite data handling 
c  system communications network (sdhs-comnet), 01 february 1988
c  2.  afgwc/dons pv-wave program auto_convert.pro 
c  3.  vax fortran reference manual from the sdhs programmer library  
c***********************************************************************

      integer integrin1,integrin2 !input integers that contain the 
                                    !the bit strings
      real*8 realout     !the output double precision number
      integer signbit  !the sign bit
      integer characts !the integer representation of the exponent
                         !before it is biased


      real*8 m1,m2,m3,m4 ! components of the fraction;m1 is most 
                         !significant

      real*8 m5,m6,m7    ! components of the fraction;m1 is most 
                         !significant

      real*8 fraction    !the fraction portion of the double precision
                         !number 

      real*8 exponent    !the exponent portion of the double precision 
                         !number

c***********************************************************************
c  pick off the sign bit and convert it to an integer (a zero or a one)
c***********************************************************************

      signbit = ibits(integrin1,23,1)
      

c***********************************************************************
c  pick off the exponent bits (7 bits from one byte, 1 bit from another)
c  convert the two bit strings to their integer representation and add
c  to get the unbiased exponent value  
c***********************************************************************

      characts = ibits(integrin1,16,7)*2 + ibits(integrin1,31,1)

c***********************************************************************
c  bias the exponent
c***********************************************************************
     
      exponent = dfloat (characts - 128)

c***********************************************************************
c pick off the bits strings from the three difference bytes that contain
c fraction bits.  convert the bit strings to the appropriate fraction.
c to account for the fact that the most significant fraction bit is
c assumed to be set, 2**7 is added to the most significant fractional 
c component m1. 
c***********************************************************************
     
      m1 = dfloat( (ibits( integrin1,24,7 ) + 2**7 ) )/(2.0d0)**8
      m2 = dfloat(ibits(integrin1,8,8))/(2.0d0)**24 
      m3 = dfloat(ibits(integrin1,0,8))/(2.0d0)**16
      m4 = dfloat( ibits(integrin2,16,8) )/(2.0d0)**32
      m5 = dfloat(ibits(integrin2,24,8))/(2.0d0)**40
      m6 = dfloat(ibits(integrin2,0,8))/(2.0d0)**48
      m7 = dfloat(ibits(integrin2,8,8))/(2.0d0)**56

c***********************************************************************
c  compute the fraction portion of the double precision number by adding
c  the fractional components
c***********************************************************************

      fraction = m1+m2+m3+m4+m5+m6+m7

c***********************************************************************
c  compute the double precision number according to the defintion of a 
c  vax double precision number.
c***********************************************************************
     
      realout= dfloat(1-2*signbit)*fraction*2.0d0**exponent

c***********************************************************************
c  since one must assume that the most significant fraction bit is set.
c  zero in the above computation will be computed as .5 i.e. all fraction
c  bits were zero but the implied bit resulted in the fraction being .5.
c  to ensure zero is properly represented the following code checks for 
c  when the fraction is .5 and set the floating point number to zero
c***********************************************************************

      if(fraction .eq. .5)then
        realout=0.0
      endif
 
      return
      end

