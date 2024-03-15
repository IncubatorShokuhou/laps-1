      subroutine pk_smple(kfildo,ia,nval,ipack,nd5,locn,ipos,ibit,
     1                    l3264b,ier,*)
c
c        march    2000   glahn   tdl   hp   for grib2 
c        january  2001   glahn   changed algorithm for computing ibit;
c                                comments
c        november 2001   glahn   put * in front of 900 in last call
c                                to pkbg
c        january  2002   glahn   added error ier = 706
c
c        purpose
c            packs data at "units" resolution provided in
c            ia( ) using the 'simple' packing method detailed in
c            the grib2 wmo standards.
c
c            the following equation is used to pack the data:
c               x1 = [(y - r) * (2 ** -e) * (10 ** -d)]
c                    x1 = the packed value
c                     y = the value we are packing
c                     r = the reference value (first order minima)
c                     e = the binary scale factor
c                     d = the decimal scale factor
c            r has already been removed upon entry.
c
c
c        data set use
c           kfildo - unit number for output (print) file. (output) 
c
c        variables 
c              kfildo = unit number for output (print) file.  (input) 
c               ia(k) = data to pack (k=1,nval). (input)
c                nval = number of values in ia( ).  (input)
c            ipack(j) = the array to hold the actual packed message
c                       (j=1,max of nd5).  (input/output)
c                 nd5 = dimension of ipack( ).  (input)
c                locn = the word position to place the next value.
c                       (input/output)
c                ipos = the bit position in locn to start placing
c                       the next value. (input/output)
c                ibit = the number of bits required to pack each
c                       value in ia( ).  (output) 
c              l3264b = contains the number of bits in a word
c                       implemented on this particular platform.
c                       (input).
c                 ier = status error return. (output)
c                         0 = good return.
c                         1 = packing would overflow ipack( ).
c                         2 = ipos not in range 1 to l3264b.
c                         3 = nbit not in range 0 to 32.
c                         4 = nbit eq 0, but nvalue ne 0.
c                       706 = value will not pack in 30 bits.
c                   * = alternate return when ier ne 0.
c
c         local variables
c               ifill = number of bits to pad section 7 to a whole
c                       octet.
c               izero = contains 0.
c
c        non system subroutines called
c           pkbg, pk_s7
c
      dimension ia(nval)
      dimension ipack(nd5)
c
      data izero/0/
c
c        determine ibit, the number of bits required to pack ia( ).
c        the initial loop is to see whether there are any non-zero
c        values.  if there aren't, then only the reference is
c        needed, and the nval points are packed zero bits each.
c
      ier=0
c
      do 110 k=1,nval
         if(ia(k).gt.0)go to 112
 110  continue
c
c        drop through here means all values are zero.
c 
      ibit=0
      go to 130 
c
 112  ibit=1
      ibxx2=2
c
 114  do 120 l=k,nval
         if(ia(l).lt.ibxx2)go to 120
         ibit=ibit+1
c 
         if(ibit.ge.31)then
            ier=706
c           write(kfildo,115)ia(l),ier
c115        format(' ****error in pk_smple.  value ='i12,
c    1             ' will not pack in 30 bits.  error code =',i5)    
            go to 900
         endif
c                 
         ibxx2=ibxx2*2
         go to 114
c
 120  continue      
c
c        pack the values in ia( ).
c
 130  call pk_s7(kfildo,ipack,nd5,locn,ipos,ia,nval,ibit,
     1           l3264b,ier,*900)
c
c        pad with zeros to fill out an octet, if necessary.
c
      ifill=mod(33-ipos,8)
c
      if(ifill.ne.0)then
         call pkbg(kfildo,ipack,nd5,locn,ipos,izero,ifill,l3264b,
     1             ier,*900)
      endif
c
 900  if(ier.ne.0)return1
c
      return
      end
