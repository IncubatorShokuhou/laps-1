      subroutine prep_sect2_real(as2,js2,nval,ib,id,rmina)
c
c        may     2000   lawrence   gsc/tdl   original coding
c        january 2001   glahn      removed maxa; comments
c        january 2001   glahn/lawrence changed 10**ib to 2**ib;
c                                  changed js2(k)*je to js2(k)/je
c        janaury 2001   lawrence   changed je=2**ib to je=2**(-ib)
c                                  and now multiply the values in
c                                  js2( ) by je instead of dividing
c                                  by it. 
c        november 2001   glahn     comment about ib 
c
c        purpose
c            prepares the local data to be
c            packed using the simple packing method.  the
c            preparation of the data encompasses multiplying
c            all of the data elements in the field by
c            the decimal scaling factor, removing the minimum
c            value, and dividing all of the data elements
c            in the field by the binary scaling factor.
c
c        data set use
c           none
c
c        variables
c              as2(k) = contains the local data to be packed using
c                       the simple packing method in the local use
c                       section of the grib2 message (k=1,nval).
c                       (input)
c              js2(k) = once the data in as2( ) has been scaled and the
c                       reference value has been removed, it is copied
c                       into this integer array and returned to the
c                       caller of this routine.  (output) 
c                nval = the dimension of js2( ) and as2( ).  (input)
c                  ib = the binary scaling factor.  this is
c                       0 in the calling program pk_sect2, and
c                       has no effect on the result.  (input)
c                  id = the decimal scaling factor.  (input)
c               rmina = the field minimum value.  (output)
c            
c             local variables
c                   d = the decimal multiplication factor, 10.**id.
c                   e = the binary multiplication factor, 2.**ie.
c                   k = loop indexing variable.
c
c        non system subroutines called
c           none
c
      dimension as2(nval),js2(nval)
c
c     scale the data field by the decimal scaling factor.
c
      if(id.ne.0)then
         d=10.**(int(id))
c
         do 10 k=1,nval
            as2(k)=as2(k)*d
 10      continue
c
      endif
c
c        find the minimum of the data field.
      rmina=as2(1)
c
      do 20 k=2,nval
c
         if(as2(k).lt.rmina)then
           rmina=as2(k)
         endif
c
 20   continue
c
c        subtract out the minimum from the data
c        field.
c
      do 30 k=1,nval
         as2(k)=as2(k)-rmina
 30   continue
c
c        multiply the data field by the binary
c        scaling factor.
c
      if(ib.ne.0)then
         e=2.**(int(-ib))
c
         do 40 k=1,nval
            as2(k)=as2(k)*e
 40      continue
c
      endif
c
c        copy the floating point values into the
c        integer output array. 
      do 50 k=1,nval
         js2(k)=nint(as2(k))  
 50   continue
c
      return
      end
