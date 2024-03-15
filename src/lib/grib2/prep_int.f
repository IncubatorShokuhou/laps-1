      subroutine prep_int(ia,nxy,nval,iclean,id,ie,
     1                    mina,jmissp,jmisss,missp,misss,ier,*)
c
c        may      2000   lawrence  original coding
c        january  2001   glahn     comments; missp changed to
c                                  misss in do 130 loop;  eliminated
c                                  unused is5( ) and ns5 from call
c        january  2001   lawrence  changed je=2**ie to je=2**(-ie)
c                                  and now multiply the values in
c                                  ia( ) by je instead of dividing
c                                  by it.
c        november 2001   glahn     modified purpose to indicate 
c                                  integer data
c        january  2002   glahn     added ier return for data with
c                                  too wide a range
c
c        purpose
c            finds the reference value and subtracts it from the
c            data values when the original data field is integer. 
c            the data are scaled by the decimal and binary
c            scale factors as well.  both primary and secondary
c            missing values are taken into account.
c
c            if a data field composed of all missing values is
c            encountered, then the data field will be packed
c            using the simple packing method if that is what
c            the user specified. if the simple method was not
c            specified, then the complex packing method is used.
c
c            operations within loops are kept to a minimum at the
c            expense of more code for efficiency.
c
c        data set use
c           kfildo - unit number for output (print) file. (output)
c
c        variables
c              kfildo = unit number for output (print) file.  (input)
c               ia(k) = contains the integer data field.
c                       (k=1,nval).  (input/output)
c                 nxy = dimension of ia( ).  (input)
c                nval = the number of values in ia( ).  (input)
c              iclean = 1 when there are no missing values ia( ).
c                       0 otherwise.  (input)
c                  id = the decimal scaling exponent.  (input)
c                  ie = the binary scaling exponent.  (input)
c                mina = the field minimum value when the original data
c                       are integer.  (output)
c              jmisss = .true. if there is a secondary missing
c                       value in the data field.  .false. otherwise.
c                       (logical)
c              jmissp = .true. if there is a primary missing value
c                       in the data field.  .false. otherwise.
c                       (logical)
c               missp = when missing points can be present in the data,
c                       they will have the value missp or misss when
c                       the data are integer.  xmissp is the primary
c                       missing value.  (input)
c               misss = secondary missing value indicator when the data
c                       are integer.  (input)
c                 ier = error return.
c                         0 = good return
c                       920 = a value larger than what can be packed
c                             into 30 bits has been encountered.
c                       (output) 
c
c        local variables
c                ipos = the position in the array of the first
c                       non-missing value when there are missing
c                       values.  (internal)
c                  jd = the decimal scaling factor.
c                  je = the binary scaling factor.
c                   k = a loop index variable.
c
c        non system subroutines called
c           none.
c
      logical jmisss, jmissp
c
      dimension ia(nxy)
c
      ier=0
c
      if(iclean.eq.1)then
c
c           there are no missing values in the array.
c           this is the simplest case ...
c           process the data straight away.
c
         if(id.ne.0)then
            jd=10**id
c
            do 10 k=1,nval
               ia(k)=ia(k)*jd
 10         continue
c
         endif
c
         mina=ia(1)
c
         do 20 k=1,nval
            if(ia(k).lt.mina)mina=ia(k)
 20      continue
c
         if(mina.ne.0)then
c
            do 30 k=1,nval
               ia(k)=ia(k)-mina
               if(ia(k).lt.0)go to 901
c                 ia(k) should not be negative.  this is probably an
c                 overflow because of too large a range of values.
 30         continue
c
         endif
c
         if(ie.ne.0)then
            je=2**(-ie)
c
            do 40 k=1,nval
               ia(k)=ia(k)*je
 40         continue
c
         endif
c
      else if(jmissp.and..not.jmisss)then
c
c           there are primary missing values in this field. 
c
         if(id.ne.0)then
            jd=10**id
c
            do 50 k=1,nval
               if(ia(k).eq.missp)cycle
               ia(k)=ia(k)*jd
 50         continue
c
         endif
c
c           find the first non-missing value
c           and set mina to that value.  ipos
c           will end up containing the 
c           position of the first non-missing
c           value.
         ipos=1
c
         do 55 k=1,nval
            if(ia(k).ne.missp)exit
            ipos=ipos+1
 55      continue
c
         mina=ia(ipos)
c
         do 60 k=ipos+1,nval
            if(ia(k).eq.missp)cycle
            if(ia(k).lt.mina)mina=ia(k)
 60      continue
c
         if(mina.ne.0)then
c
            do 70 k=1,nval
               if(ia(k).eq.missp)cycle
               ia(k)=ia(k)-mina
               if(ia(k).lt.0)go to 901
c                 ia(k) should not be negative.  this is probably an
c                 overflow because of too large a range of values.
 70         continue
c
         endif
c
         if(ie.ne.0)then
            je=2**(-ie)
c
            do 80 k=1,nval
               if(ia(k).eq.missp)cycle
               ia(k)=ia(k)*je
 80         continue
c
         endif
c
      else if(jmisss.and.jmissp)then
c
c           there are both primary and secondary 
c           missing values.
         if(id.ne.0)then
            jd=10**id
c
            do 90 k=1,nval
               if((ia(k).eq.missp).or.(ia(k).eq.misss))cycle
               ia(k)=ia(k)*jd
 90         continue
c
         endif
c
c           find the first non-missing value
c           and set mina equal to that value.  ipos
c           will end up containing the position
c           of the first non-missing value.
         ipos=1
c
         do 95 k=1,nval
            if((ia(k).ne.missp).and.(ia(k).ne.misss))exit
            ipos=ipos+1
 95      continue
c
         mina=ia(ipos)
c
         do 100 k=ipos+1,nval
            if((ia(k).eq.missp).or.(ia(k).eq.misss))cycle
            if(ia(k).lt.mina)mina=ia(k)
 100     continue
c
         if(mina.ne.0)then
c
            do 110 k=1,nval
               if((ia(k).eq.missp).or.(ia(k).eq.misss))cycle
               ia(k)=ia(k)-mina
               if(ia(k).lt.0)go to 901
c                 ia(k) should not be negative.  this is probably an
c                 overflow because of too large a range of values.
 110        continue
c
         endif
c
         if(ie.ne.0)then
            je=2**(-ie)
c
            do 120 k=1,nval
               if((ia(k).eq.missp).or.(ia(k).eq.misss))cycle
               ia(k)=ia(k)*je
 120        continue
c
         endif
c
      else if(jmisss.and..not.jmissp)then
c
c           there are only secondary missing values.
         if(id.ne.0)then
            jd=10**id
c
            do 130 k=1,nval
               if(ia(k).eq.misss)cycle
               ia(k)=ia(k)*jd
 130        continue
c
         endif
c
c           find the first non-missing value
c           and set mina to that value.  ipos
c           will end up containing the 
c           position of the first non-missing
c           value.
         ipos=1
c
         do 135 k=1,nval
            if(ia(k).ne.misss)exit
            ipos=ipos+1
 135     continue
c
         mina=ia(ipos)
c
         do 140 k=ipos+1,nval
            if(ia(k).eq.misss)cycle
            if(ia(k).lt.mina)mina=ia(k)
 140     continue
c
         if(mina.ne.0)then
c
            do 150 k=1,nval
               if(ia(k).eq.misss)cycle
               ia(k)=ia(k)-mina
               if(ia(k).lt.0)go to 901
c                 ia(k) should not be negative.  this is probably an
c                 overflow because of too large a range of values.
 150        continue
c
         endif
c
         if(ie.ne.0)then
            je=2**(-ie)
c
            do 160 k=1,nval
               if(ia(k).eq.misss)cycle
               ia(k)=ia(k)*je
 160        continue
c
         endif
c
      endif
c
 900  return
c
c        alternate return section.
c
 901  ier=920
c     write(12,902)k,ia(k),mina
c902  format(' too large a range of values to pack in prep_int.',
c    1       'k, ia(k), mina are',3i12)
      return1
      end
