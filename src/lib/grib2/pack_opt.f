      subroutine pack_opt(kfildo,ia,ib,nxy,nval,iclean,ibitmap,
     1                    is5,ns5,is7,ns7,jmisss,
     2                    missp,minpk,numoctet,ipkopt,jer,
     3                    ndjer,kjer,*)
c
c        may      2000   lawrence  original coding - based on logic
c                                  developed by harry glahn.
c        january  2001   glahn     comments; xmissp changed to
c                                  xmisss in do 130 loop;
c        january  2001   glahn/lawrence removed unused miss, is6( ),
c                                  and ns6 from call; removed
c                                  setting is5(6)
c        november 2001   glahn     several diagnostic calls to pk_trace 
c        january  2002   glahn     inserted ier = 906 and 907
c
c        purpose
c            determines which packing method is to be used
c            to pack the data.  the user has the first say
c            in what type of packing method to use by setting
c            the appropriate value in is5(10).  this routine
c            will then determine if it is beneficial and
c            possible to use the packing method that the user
c            chose.
c
c        data set use
c           kfildo - unit number for output (print) file. (output)
c
c        variable
c             kfildo = unit number for output (print) file.  (input)
c               ia(k) = contains the integer data field.
c                       (k=1,nval).  (input/output)
c               ib(k) = the bit map when one is used.  (input/output)
c                       it can be input or in can be calculated if
c                       the simple method is used (k=1,nxy).
c                       complex and spatial differencing do not
c                       use a bit map, but will accept one and insert
c                       the missing values.
c                 nxy = dimension of ia( ) and ib( ).  (input)
c                nval = the number of values in ia( ).  (input)
c              iclean = 1 when there are no missing values ia( ).
c                       0 otherwise.  (input)
c             ibitmap = 1 when there is a bitmap in ib( , ).
c                       0 otherwise.  (input)
c              is5(k) = the values associated with section 5, keyed
c                       to the octet number (k=1,ns5). (input/output)
c                 ns5 = the dimension of is5( ).  (input)
c              is7(k) = the values associated with section 7, keyed
c                       to the octet number (k=1,ns7).  (output)
c                 ns7 = the dimension of is7( ). (input)
c              jmisss = .true. if there is a secondary missing
c                        value in the data field.  .false. otherwise
c                        (logical).  (input)
c               missp = the primary missing value.  (input)
c               minpk = increment in which ranges will be computed.
c                       (input)
c            numoctet = the minimum number of octets required to
c                       pack the extra descriptors when second
c                       order spatial differencing is used with
c                       complex packing.  (output) 
c              ipkopt = packing indicator:
c                       0 = error, don't pack
c                       1 = pack ia( ), simple
c                       2 = pack ia( ) and ib( ), simple
c                       3 = pack complex or spatial differencing
c                       4 = pack complex.
c                       (output)
c              jer(j) = array of errors (j=1,ndjer), max of nderr.
c                       906 = simple packing, no missing values in
c                             array, bit map provided, and 
c                             nval = nxy.  unusual; notify user.
c                       907 = simple packing, no missing values in
c                             array, no bit map provided, but 
c                             nval ne nxy, unrecoverable error.
c                       908 = is5(23) set = 0 consistent with
c                             iclean = 1 and not simple packing.
c                       910 = is5(23) set = 1 because there are
c                             no secondary missing values in field.
c                       911 = is5(10) set = 2 to indicate complex
c                             because it is more efficient than
c                             second order differencing.
c                       912 = is5(23) set = 2 to indicate secondary
c                             missing values.
c                       915 = is5(10) set = 2 because secondary missing
c                             values are present and second order
c                             differencing not supported.
c                       (input/output)
c               ndjer = dimension of jer( ).  (input)
c                kjer = number of values in jer( ).  (input/output)
c                   * = alternate return when jer(kjer) ge 900.
c
c        local variables
c                 ier = contains any error codes generated from
c                       utilities called by this routine.
c                   k = an array/looping index.
c               kount = an array index that is used when either
c                       expanding ia( ) to place missing values in
c                       it or generating a bit-map from ia( ).
c              second = .true. if it is determined that second
c                       order differences would be efficient to use.
c                       .false. if it is determined that second
c                       order differences would not be efficient
c                       to use.
c
c        non system subroutines called
c           pk_trace, pk_missp, pk_nomiss
c
      logical jmisss, second
c
      dimension ia(nxy),ib(nxy)
      dimension is5(ns5)
      dimension is7(ns7)
      dimension jer(ndjer,2)
c
      data second/.false./
c
      ier=0
c     write(kfildo,100)iclean,ibitmap,is5(10),nval
c100  format(/' in pack_opt--iclean,ibitmap,is5(10),nval',4i8)
c
      if(iclean.eq.1)then
c
c           there are no missing values in ia( ).
         if(ibitmap.eq.1)then
c
c              there is a bit map.
            if(is5(10).eq.0)then
c                 simple packing specified.
c
               if(nval.eq.nxy)then
c                    no missing values in array, bit-map
c                    provided, simple, and nval=nxy.
c                    unusual; notify user.
c                    pack ia( ).
                  ipkopt=1
                  ier=906
                  call pk_trace(kfildo,jer,ndjer,kjer,906,1)
               else
c                    no missing values in array, bit map
c                    provided, simple, and nval ne nxy.
c                    pack ia( ) and ib( ).
                  ipkopt=2
               endif
c
            endif
         else
            if(is5(10).eq.0)then
c                 there is no bit map
c
               if(nval.eq.nxy)then
c                     no missing values in array, no bit 
c                     map provided, and nval=nxy.  simple.
c                     pack ia( )
                  ipkopt=1
               else
c                    no missing values in array, no bit
c                    map provided, and nval ne nxy.simple.
c                    unrecoverable error, do not pack.
                  ipkopt=0
                  ier=907
                  call pk_trace(kfildo,jer,ndjer,kjer,907,2)
                  go to 900
               endif
c
            else if(is5(10).eq.2)then
c                 complex packing specified.
c
c                 no missing values in the array, no bit-map
c                 provided, pack complex without second order
c                 differences.
               ipkopt=4
c
               if(is5(23).ne.0)then
                  is5(23)=0
c                    no missing values in data.
                  call pk_trace(kfildo,jer,ndjer,kjer,908,1)
               endif
c
            else
c                 complex packing with 2nd order spatial differencing.
c
c                 no missing values in the array, no bit-map
c                 provided, not simple. pack complex with
c                 or without spatial differences.
               call pk_nomiss(kfildo,ia,ib,nxy,minpk,
     1                        ifirst,isecond,imina,
     2                        numoctet,second,ier)
c
               if(ier.ne.0)then
                  call pk_trace(kfildo,jer,ndjer,kjer,ier,2)
                  go to 900
               endif
c
               if(second)then
                  ipkopt=3
                  is7(6)=ifirst
                  is7(7)=isecond
                  is7(8)=imina
               else
                  ipkopt=4
c
                  is5(10)=2
c                    pack complex without 2nd order differencing.
                  call pk_trace(kfildo,jer,ndjer,kjer,911,0)
               endif
c
               if(is5(23).ne.0)then
                  is5(23)=0
c                    no missing values in data.
                  call pk_trace(kfildo,jer,ndjer,kjer,908,1)
               endif
c
            endif
c
         endif
c
      else
c
c           at this point, complex or spatial difference packing
c           is to be done, no bit map.  jmissp and jmisss have been 
c           set to indicate missing primary and secondary missing 
c           values, respectively, and there are at least primary
c           missing values.  if secondary missing values are present,
c           packing has to consider primary missing values also.
c           jmissp = 1.  spatial differencing is not done when
c           secondary missing values are present.
c
         if(jmisss)then
c
            if(is5(23).ne.2)then
               is5(23)=2
c                 there are secondary missing values.
               call pk_trace(kfildo,jer,ndjer,kjer,912,1)
            endif
c
         else
c
            if(is5(23).ne.1)then
               is5(23)=1
c                 there are no secondary missing values.
               call pk_trace(kfildo,jer,ndjer,kjer,910,1)
            endif
c
         endif
c
         if((is5(10).eq.2).or.(jmisss))then
c              primary and secondary missing values in array.
c              bit map not provided; complex packing.  pack ia( ).
            ipkopt=4
c
            if(is5(10).eq.3)then
               is5(10)=2
               call pk_trace(kfildo,jer,ndjer,kjer,915,1)
            endif
c
         else
c              primary missing values in array, but no
c              secondary missing values.  bit map not provided,
c              spatial difference or complex.  pack ia( ).
            call pk_missp(kfildo,ia,ib,nxy,missp,minpk,
     1                    ifirst,isecond,imina,numoctet,
     2                    second,ier)
c
            if(ier.ne.0)then
               call pk_trace(kfildo,jer,ndjer,kjer,ier,2)
               go to 900
            endif
c
            if(second)then
               ipkopt=3
               is7(6)=ifirst
               is7(7)=isecond
               is7(8)=imina
            else
               ipkopt=4
c
               if(is5(10).ne.2)then
                  is5(10)=2
c                    pack complex without 2nd order differencing.
                  call pk_trace(kfildo,jer,ndjer,kjer,911,0)
               endif
c
            endif
c
         endif
c
      endif
c
      return
c
 900  return 1
      end
