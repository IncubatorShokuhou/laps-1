      subroutine cnvtint(intin)

c***********************************************************************
c  purpose: properly interpret the bit string that represents a negative 
c  integer in 2's complement form.   
c
c  method: if the sign bit of the input integer is set, then the number
c  is negative and is assumed to be in 2's complement form.  the 
c  appropriate representation of the integer on the target system is 
c  obtained by inverting the 2's complement form using the standard 
c  fortran functions ibclr (clear a bit) and not (a bitwise complement)
c    
c  reference: operations ground equipment interface specification drl 
c             504-02-1 part1, document no. e007020;space systems/loral
c   
c***********************************************************************

      integer intin
      integer signbit

c***********************************************************************
c  pick off the sign bit and convert it to an integer (a zero or a one)
c***********************************************************************

      signbit = ibits(intin,31,1)

c***********************************************************************
c  check to see if the number is negative, if so, invert the 2's 
c  complement notation.  the inversion is done by clearing the sign bit,
c  subtracting one , flipping the bits, and reclearing the sign bit. 
c  the result is the magnitude of the negative number
c***********************************************************************

      if(signbit .eq. 1) then
        intin = ibclr(intin,31)
        intin = not((intin-1))
        intin = ibclr(intin,31)
      endif

c***********************************************************************
c  compute the integer based on the sign bit and the magnitude of the
c  input number. 
c***********************************************************************

      intin = (1-2*signbit)*intin

      return
      end

