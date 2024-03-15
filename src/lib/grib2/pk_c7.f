      subroutine pk_c7(kfildo,ipack,nd5,locn,ipos,ia,nxy,nov,lbit,lx,
     1                 l3264b,ncount,ier,*)
c
c        may     1997   glahn    tdl   hp
c        july    1999   lawrence updated this routine so that it returns a 
c                                count of the values packed into ipack.
c                                this is needed by pk_cmplx and it is
c                                sharing this routine with pk52. 
c        august  1999   lawrence updated the documentation in this routine 
c                                in keeping with tdl standards.
c        march   2000   glahn    updated for grib2; name changed from;
c                                pkc4lx; changed loc to locn
c        january 2001   glahn    comments; ier = 1 changed to ier = 705;
c                                added return1; changed ier from 932 to 2
c                                etc.
c
c        purpose 
c            packs up to nxy values into ipack( ).  the packed values
c            are taken from ia( ) with no reference value or 
c            scaling considered.  pk_c7 eliminates the calling of
c            pkbg, and rather incorporates it into the loop.
c            since this is a highly used routine, all reasonable
c            attempts at efficiency must be pursued.  this is
c            for complex packing, the counterpart of pk_s7 for
c            simple packing.
c
c        data set use 
c           kfildo - unit number for output (print) file. (output) 
c
c        variables 
c              kfildo = unit number for output (print) file.  (input) 
c            ipack(j) = the array holding the actual packed message
c                       (j=1,max of nd5).  (input/output)
c                 nd5 = dimension of ipack( ).  (input)
c                locn = holds word position in ipack of next value to
c                       pack.  (input/output)
c                ipos = holds bit position in ipack(locn) of the first
c                       bit of the next value to pack.  (input/output)
c               ia(k) = data to pack (k=1,nxy).  (input)
c                 nxy = dimension of ia( ).  also, the number of values 
c                       to be packed unless all members of a group
c                       have the same value.  (input)
c              nov(k) = the number of values per group (k=1,lx).(input)
c             lbit(k) = the number of bits to pack for each group
c                       (k=1,lx).  (input)
c                  lx = the number of values in nov( ) and lbit( ).
c                       (input)
c              l3264b = integer word length of machine being used.
c                       (input)
c              ncount = the number of values packed. this value is
c                       initially zeroed out in this routine.
c                       (output) 
c                 ier = error return.
c                         2 = ipos not in the range 1-l3264b.
c                         3 = lbit( ) not in the range 0-32.
c                       705 = nd5 is not large enough to accommodate the
c                             bits necessary to pack the values starting
c                             at the values locn and ipos.
c                       (output)
c                   * = alternate return when ier ne 0.
c
c        local variables
c                   k = the index of the current data item we are 
c                       packing.
c                 l,m = loop indexing variable.
c               lbitl = the number of bits required to pack the
c                       largest value in each group.
c       lbitl1,lbitl2 = used in packing the data using mvbits. they
c                       keep track of temporary bit positions.
c             newipos = used to keep track of the bit
c                       position to start packing at.
c
c        non system subroutines called 
c           none
c
      dimension ipack(nd5)
      dimension ia(nxy)
      dimension nov(lx),lbit(lx)
c
c        set error return.
c
      ier=0
      ncount=0
c
c        check legitimate values of locn and ipos.
c
      if(ipos.le.0.or.ipos.gt.l3264b)then
         ier=2
c        write(kfildo,101)ipos,ier
c101     format(/' ipos = 'i6,' not in the range 1 to l3264b in pk_c7.',
c    1           ' return from pk_c7 with ier = 'i4)
         go to 900 
      endif
c
      k=0
c
      do 300 l=1,lx
c
      lbitl=lbit(l)
c
c        check legitimate values of lbit(l) and whether
c        nd5 is sufficient for all nov(l) values.
c
      if(lbitl.lt.0.or.lbitl.gt.32)then
         ier=3
c        write(kfildo,102)lbitl,l,ier
c102     format(/' lbit(l) = 'i6,' for l ='i5,'not in the range',
c    1           ' 0 to 32 in pk_c7.',
c    1           ' return from pk_c7 with ier = 'i4)
         go to 900
      endif
c
      if(lbitl*nov(l).gt.(l3264b+1-ipos)+(nd5-locn)*l3264b)then
         ier=705
c        write(kfildo,103)nov(l),lbitl,l,locn,ipos,nd5,ier
c103     format(/' nov(l) = 'i9,' and lbit(l) = 'i6,' for l ='i5,
c    1           ' require more bits than are available in ipack( ),',
c    2           ' with locn ='i8,', ipos ='i4,', and nd5 ='i8,'.'/
c    3           ' return from pk_c7 with ier ='i4)
         go to 900
      endif
c
      do 290 m=1,nov(l)
      k=k+1
      if(lbitl.eq.0)go to 290
c        a group with all values the same is omitted.
c
      ncount=ncount+1
      newipos=ipos+lbitl
c
      if(newipos.le.l3264b+1)then
         call mvbits(ia(k),0,lbitl,ipack(locn),l3264b+1-newipos)
c
         if(newipos.le.l3264b)then
            ipos=newipos
         else
            ipos=1
            locn=locn+1
         endif
c
      else
         lbitl1=l3264b+1-ipos
         lbitl2=lbitl-lbitl1
         call mvbits(ia(k),lbitl2,lbitl1,ipack(locn),0)
         locn=locn+1
         call mvbits(ia(k),0,lbitl2,ipack(locn),l3264b-lbitl2)
         ipos=lbitl2+1
      endif
c
 290  continue
c
 300  continue
c
 900  if(ier.ne.0)return1
c
      return
      end
