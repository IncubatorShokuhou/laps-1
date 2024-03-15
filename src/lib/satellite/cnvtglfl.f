      subroutine cnvtglfl(integrin,realout)

c***********************************************************************
c  purpose: properly interpret the bit string that represents a gould 
c  floating point number to obtain the correct floating point number on 
c  the target system.   
c
c  method: the bit string that represents the gould floating point number
c  is passed into the routine as an integer. the sign bit, the bits that
c  represent the exponent or characteristic, and the bits that represent 
c  the fraction are picked off using the standard fortran function ibits.
c  the following table show the bits and what they represent (sign, 
c  exponent etc).
c    
c          byte1       byte2    byte3     byte4
c          31 30-24    23-16    15-08     07-00
c          s  e7-e1    f8-f1    f16-f9    f24-f17
c
c    
c
c  here f1 stands for the first bit that is used to represent the 
c  fraction, f2 for the second bit used to represent the fraction, etc., 
c  e1 represents the first exponent bit, e2 the second, and s represents 
c  the sign bit.  a gould floating point number is represented by :
c
c  floating point number = -1**signbit * fraction * 2 **((exponent-64)*4)
c    
c  !!!!!!note: for negative numbers the bit string is stored in 2's
c              complement form.  this form must be inverted before the
c              bit string can be properly interpreted
c    
c  reference: operations ground equipment interface specification drl 
c             504-02-1 part1, document no. e007020;space systems/loral
c   
c***********************************************************************

      integer integrin !the input integer that contains the bit string
      real*8 realout !the output floating point number
      integer signbit !the sign bit
      integer characts !the integer representation of the exponent 
                         !before it is biased

      real*8 m1,m2,m3 !components of the fraction;m1 is most significant
      real*8 fraction !the fraction portion of the floating point number
      real*8 exponent !the exponent portion of the floating point number

c***********************************************************************
c  pick off the sign bit and convert it to an integer (a zero or a one)
c***********************************************************************
      
      signbit = ibits(integrin,31,1)

c***********************************************************************
c  check to see if the number is negative, if so, invert the 2's 
c  complement notation.  the inversion is done by clearing the sign bit,
c  subtracting one and flipping the bits
c***********************************************************************

      if (signbit .eq. 1) then
        integrin = ibclr(integrin,31)
        integrin = not((integrin-1))
      endif

c***********************************************************************
c  pick off the exponent bits (7 bits) convert the two bit strings to
c  their integer representation and add to get the unbiased exponent
c  value  
c***********************************************************************

      characts = ibits(integrin,24,7)

c***********************************************************************
c  bias the exponent
c***********************************************************************

      exponent = dfloat ((characts - 64)*4 )

c***********************************************************************
c pick off the bits strings from the three difference bytes that contain
c fraction bits.  convert the bit strings to the appropriate fraction.
c***********************************************************************

      m1 = dfloat( (ibits( integrin,16,8 ) ) )/(2.0d0)**8
      m2 = dfloat(ibits(integrin,8,8))/(2.0d0)**16 
      m3 = dfloat(ibits(integrin,0,8))/(2.0d0)**24

c***********************************************************************
c  compute the fraction portion of the floating point number by adding
c  the three fractional components
c***********************************************************************
 
      fraction = m1 +m2 + m3

c***********************************************************************
c  compute the floating point number according to the defintion of a 
c  gould floating point number.
c***********************************************************************

      realout= float(1-2*signbit)*fraction*2.0**exponent


      return
      end
