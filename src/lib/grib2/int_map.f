      subroutine int_map(ia,ib,nxy,is5,ns5,iclean,ibitmap,
     1                   missp,jer,ndjer,kjer)
c
c        may      2000   lawrence original coding
c        january  2001   glahn    comments
c        november 2001   glahn    added jer, ndfer, and kjer to
c                                 call and called trace
c
c        purpose
c            processes the bit-map associated with a grid of integer
c            data values.  exactly how a bit-map is processed depends
c            primarily upon which packing method is being used to
c            compress the gridded data.  the bit-map indicates where
c            missing values are located in the data field being packed.
c            a 0 indicates a bad value while 1 indicates a good value.
c            there is a one-to-one correspondence between the values
c            in the bit-map and the data in the data field before 
c            missing values are removed.
c
c            a bit-map is required with simple packing to
c            indicate where missing values are located in the
c            data field.  if a field to be packed using the
c            simple method isn't accompanied by a bit-map and
c            there are missing values embedded in the data,
c            then this routine will generate a bit-map and
c            remove the missing values from the data field.
c
c            the complex packing method does not need a bit-map
c            to indicate where missing values are located in the
c            data field.  a bit-map can be supplied by the user
c            if they are providing a data field without the missing
c            values in it.  this routine will place the missing
c            values into their proper places in the grid in that
c            case and eliminate the bit map.
c
c        data set use
c           none
c
c        variables
c               ia(k) = ia( ) contains the data to be compressed
c                       (k=1,nxy).  (input/output)
c               ib(k) = the bit map when one is used (k=1,nxy).
c                       it can be input or it can be calculated if
c                       the simple method is used.  complex and spatial
c                       differencing methods do not use a bit map,
c                       but will accept one and insert the missing
c                       values into the field before packing.
c                       (input/output)
c                 nxy = dimension of ia( ), ib( ),and iwork( ).
c                       (input)
c              is5(k) = the values associated with section 5, keyed
c                       to the octet number(k=1,ns5). the element
c                       used in this routine is:
c                       is5(10), template number:
c                         0 = simple
c                         1 = not supported
c                         2 = complex
c                         3 = spatial differencing
c                       (input)
c                 ns5 = the dimension of is5( ).  (input)
c              iclean = 1 when there are no missing values in ia( ).
c                       0 otherwise.  (input/output)
c             ibitmap = 1 when there is a bitmap in ib( , ).
c                       0 otherwise.  (input/output)
c               missp = when missing points can be present in the data,
c                       they will have the value missp.  (input)
c            jer(j,k) = return status codes and severity levels
c                       (j=1,ndjer)(k=1,2).  values can come from
c                       subroutines; otherwise: 0 = good return.
c                       (input/output)
c               ndjer = dimension of jer( ).  (input)
c                kjer = number of values in jer( ).  (input/output)
c
c        local variables
c            iwork(k) = used to contain the data field when we
c                       are expanding it to contain missing
c                       values (only done with the complex
c                       packing method).  (automatic array)
c                   k = a looping/array indexing variable.
c                   m = used to keep track of the position of
c                       real values versus missing values in 
c                       the data field.
c
c        non system subroutines called
c          none 
c
      dimension ia(nxy),ib(nxy),iwork(nxy)
      dimension is5(ns5)
      dimension jer(ndjer,2)
c
      if(is5(10).eq.0)then
c
c           simple packing is being used. if the data field
c           contains missing values, then see if the user has
c           supplied a bitmap. if one hasn't been supplied,then
c           generate one.
c
         if(iclean.eq.0)then
c
            if(ibitmap.eq.0)then
               ibitmap=1
c                 ibitmap being overridden and a bit map being
c                 generated.
               call pk_trace(kfildo,jer,ndjer,kjer,913,0)
c
               m=1
c
               do 10 k=1,nxy
c
                 if(ia(k).eq.missp)then
                    ib(k)=0
                 else
                    ib(k)=1
                    ia(m)=ia(k)
                    m=m+1
                 endif
c
 10            continue
c
               iclean=1
            endif
c
         endif
c
c           if the user doesn't supply a bit-map
c           and the field is indicated as clean,
c           then it is assumed that there are no
c           missing values. in that case the bit-map has
c           all one's and does not need to be packed in
c           section 6 nor supplied by the user.
      elseif(is5(10).ge.2)then
c
c           a form of complex packing is being used.
c           if a bit-map has been supplied and the field
c           is clean, then expand the data field to include
c           the missing values and discard the bit-map.
c
         if(ibitmap.eq.1)then
c
            if(iclean.eq.1)then
               m=1
c
               do 20 k=1,nxy
                  if(ib(k).eq.0)then
                     iwork(k)=missp
                  else
                     iwork(k)=ia(m)
                     m=m+1
                  endif
 20            continue
c
               do 30 k=1,nxy
                  ia(k)=iwork(k)
 30            continue
c
               ibitmap=0
               iclean=0
c                 ibitmap and iclean being overridden and a
c                 the bit map being incorporated into the data.
               call pk_trace(kfildo,jer,ndjer,kjer,914,0)
c
            endif
c
         endif
c
      endif
c
      return
      end
