      subroutine pkbg(kfildo,ipack,nd5,loc,ipos,nvalue,nbit,l3264b,
     1                ier,*)
c
c        december 1994   glahn   tdl   mos-2000
c        may      1997   glahn   modified to use mvbits rather than
c                                shifting and oring, and eliminated
c                                use of mod function
c 
c        purpose 
c            packs nbit bits in the positive integer nvalue into array
c            ipack(nd5) starting in word loc, bit ipos.  the word
c            pointer loc and bit position pointer ipos are updated
c            as necessary.  packing will not occur if ipack( ) would
c            be overflowed.  in that case, return is with ier=1
c            rather than for the good return ier=0.  when nbit eq 0
c            and nvalue eq 0, no packing is done.  this routine acts
c            as "insertion" rather than "addition."  that is, the
c            bits, if any, to the right of the packed value
c            are retained.  this means that the ipack( ) array should
c            be zeroed out before using.  also, any insertion must be
c            be into an area that has all zero bits.  the integer word
c            length of the machine being used is l32b4b.  this packing
c            routine will work on either a 32- or 64-bit machine.
c            a maximum of 32 bits can be packed on a single call.
c
c        data set use 
c           kfildo - unit number for output (print) file. (output) 
c
c        variables 
c              kfildo = unit number for output (print) file.  (input)
c            ipack(j) = array to pack into (j=1,nd5).  (input-output)
c                 nd5 = dimension of ipack( ).  (input)
c                 loc = word in ipack( ) to start packing.  updated
c                       as necessary after packing is completed.
c                       (input-output)
c                ipos = bit position (counting leftmost bit in word
c                       as 1) to start packing.  must be ge 1 and
c                       le l3264b.  updated as necessary
c                       after packing is completed.  (input-output)
c              nvalue = the rightmost nbit bits in nvalue will
c                       be packed.  (input)
c                nbit = see nvalue.  must be ge 0 and le 32.  (input)   
c              l3264b = integer word length of machine being used.
c                       (input)
c                 ier = status return:
c                       0 = good return.
c                       1 = packing would overflow ipack( ).
c                       2 = ipos not in range 1 to l3264b.
c                       3 = nbit not in range 0 to 32.
c                       4 = nbit eq 0, but nvalue ne 0.
c                   * = alternate return when ier ne 0.
c
c        non system subroutines called
c            none
c
      dimension ipack(nd5)
c
c        check correctness of input and set status return.
c
      ier=0
      if(nbit.eq.0.and.nvalue.eq.0)go to 150
      if(nbit.ne.0)go to 111
      ier=4
c        when nbit=0, nvalue must be also.
 111  if(loc.lt.1.or.loc.gt.nd5)ier=1
c        packing would overflow ipack( ).
      if(ipos.le.0.or.ipos.gt.l3264b)ier=2
      if(nbit.lt.0.or.nbit.gt.32)ier=3
      if(ier.ne.0)return 1      
c
      newipos=ipos+nbit
c
c        when newipos le l3264+1, then only one word is packed into.
c        else two words are involved.
c
      if(newipos.le.l3264b+1)then
         call mvbits(nvalue,0,nbit,ipack(loc),l3264b+1-newipos)
c
         if(newipos.le.l3264b)then
            ipos=newipos
         else
            ipos=1
            loc=loc+1
         endif
c
      else
         nbit1=l3264b+1-ipos
         nbit2=nbit-nbit1
         call mvbits(nvalue,nbit2,nbit1,ipack(loc),0)
         loc=loc+1
c
         if(loc.le.nd5)go to 130
         ier=1
         return 1
c
 130     call mvbits(nvalue,0,nbit2,ipack(loc),l3264b-nbit2)
         ipos=nbit2+1
      endif
c
 150  return
      end
