      subroutine prep_sect2_int(js2,nval,ib,id,mina)
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
c              js2(k) = contains the local data to be packed using
c                       the simple packing method in the local
c                       use section of the grib2 message (k=1,nval).
c                       (input)
c                nval = the number of values in js2( ).  (input)
c                  ib = the binary scaling factor.  this is
c                       0 in the calling program pk_sect2, and
c                       has no effect on the result.  (input)
c                  id = the decimal scaling factor.  (input)
c                mina = the field minimum value.  (output)
c            
c             local variables
c                  jd = the decimal multiplication factor, 10**id.
c                  je = the binary multiplication factor, 2**ie.
c                   k = loop indexing variable.
c
c        non system subroutines called
c           none
c
      dimension js2(nval)
c
c     scale the data field by the decimal scaling factor.
c
      if(id.ne.0)then
         jd=10**id
c
         do 10 k=1,nval
            js2(k)=js2(k)*jd
 10      continue
c
      endif
c
c        find the minimum of the data field.
      mina=js2(1)
c
      do 20 k=2,nval
c
         if(js2(k).lt.mina)then
            mina=js2(k)
         endif
c
 20   continue
c
c        subtract out the minimum from the data
c        field.
c
      do 30 k=1,nval
         js2(k)=js2(k)-mina
 30   continue
c
c        multiply the data field by the binary
c        scaling factor.
c
      if(ib.ne.0)then
         je=2**(-ib)
c
         do 40 k=1,nval
            js2(k)=js2(k)*je
 40      continue
c
      endif
c
      return
      end
