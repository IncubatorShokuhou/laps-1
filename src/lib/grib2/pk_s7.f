      subroutine pk_s7(kfildo,ipack,nd5,locn,ipos,ia,nval,ibit,
     1                 l3264b,ier,*)
c
c        may     1997   glahn   tdl   hp
c        march   2000   glahn   changed name from pks4lx;
c                               ic( ) to ia( ); nxy to nval;
c                               * = return1
c        january 2001   glahn   comments; ier = 1 changed to 705;
c                               standardized return
c
c        purpose
c            packs nval values into ipack( ).  the packed values
c            are taken from ia( ).  pk_s7 eliminates the calling of
c            pkbg, and rather incorporates it into the loop.
c            since this is a highly used routine, all reasonable
c            attempts at efficiency must be pursued.  this is for
c            simple packing, the counterpart of pk_c7 for complex
c            packing.  pk_s7 accommodates ibit = 0
c
c        data set use 
c           kfildo - unit number for output (print) file. (output) 
c
c        variables 
c              kfildo = unit number for output (print) file.  (input) 
c            ipack(j) = the array holding the actual packed message
c                       (j=1,max of nd5).  (input/output)
c                 nd5 = dimension of ipack( ).  (input)
c                locn = the word position to place the next value.
c                       (input/output)
c                ipos = the bit position in locn to start placing
c               ia(k) = data to pack (k=1,nval).  (input)
c                nval = dimension of ia( ).  the number of values 
c                       to be packed.  (input)
c                ibit = the number of bits used to pack each value.
c                       (input)
c              l3264b = integer word length of machine being used.
c                       (input)
c                 ier = error return.
c                          2 = ipos not in the range 1-l3264b.
c                          3 = ibit not in the range 0-32.
c                        705 = nd5 is not large enough to accommodate the
c                              bits necessary to pack nval values 
c                              starting at the values locn and ipos.
c                       (output)
c                   * = alternate return when ier ne 0.
c
c        local variables
c         ibit1,ibit2 = used in packing the data using mvbits. they
c                       keep track of temporary bit positions.
c             newipos = used to keep track of the bit position to
c                       start packing at.
c
c        non system subroutines called
c           none
c
      dimension ipack(nd5)
      dimension ia(nval)
c
      ier=0
c 
      if(ibit.eq.0)go to 900
c        when ibit = 0, no values are packed.
c
c        check legitimate value of ipos.
c
      if(ipos.le.0.or.ipos.gt.l3264b)then
         ier=2
c        write(kfildo,101)ipos,ier
c101     format(/' ipos = 'i6,' not in the range 1 to l3264b',
c    1           ' in pk_s7.  return from pk_s7 with ier = 'i4)
         go to 900 
      endif
c
c        check legitimate value of ibit.
c
      if(ibit.lt.0.or.ibit.gt.32)then
         ier=3
c        write(kfildo,102)ibit,ier
c102     format(/' ibit = 'i6,' not in the range 0 to 32 in pk_s7.',
c    1           ' return from pk_s7 with ier = 'i4)
         go to 900
      endif
c
c        check whether nd5 is sufficient for all nval values.
c
      if(ibit*nval.gt.(l3264b+1-ipos)+(nd5-locn)*l3264b)then
         ier=705
c        write(kfildo,103)nval,ibit,locn,ipos,nd5,ier
c103     format(/' nval = 'i9,' and ibit = 'i6,' require more bits',
c    1           ' than are available in ipack( ),',
c    2           ' with locn ='i8,', ipos ='i4,', and nd5 ='i8,'.'/
c    3           ' return from pk_s7 with ier ='i4)
         go to 900
      endif
c     
      do 300 k=1,nval 
c
      newipos=ipos+ibit
c
      if(newipos.le.l3264b+1)then
         call mvbits(ia(k),0,ibit,ipack(locn),l3264b+1-newipos)
c
         if(newipos.le.l3264b)then
            ipos=newipos
         else
            ipos=1
            locn=locn+1
         endif
c
      else
         ibit1=l3264b+1-ipos
         ibit2=ibit-ibit1
         call mvbits(ia(k),ibit2,ibit1,ipack(locn),0)
         locn=locn+1
         call mvbits(ia(k),0,ibit2,ipack(locn),l3264b-ibit2)
         ipos=ibit2+1
      endif
c
 300  continue
c
 900  if(ier.ne.0)return 1
c
      return
      end
