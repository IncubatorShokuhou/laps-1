      subroutine prep_flt(a,ia,nxy,nval,iclean,id,ie,
     1                    xmina,jmissp,jmisss,xmissp,xmisss,ier,*)
c
c        may     2000   lawrence   original coding
c        january 2001   glahn      comments; xmissp changed to
c                                  xmisss in do 130 loop;  eliminated
c                                  unused is5( ) and ns5 from call
c        janaury 2001   lawrence   changed e=2.**ie to e=2.**(-ie)
c                                  and now multiply the values in
c                                  a( ) by e instead of dividing
c                                  by it.
c        january 2002   glahn      chanaged int( ) to nint( ); added
c                                  ier and alternate return when a( )
c                                  gt mallow = 2**30  
c
c        purpose
c            finds the reference value and subtracts it from the
c            data values when the original data field is floating
c            point. the data are scaled by the decimal and binary
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
c                a(k) = contains the floating point data field.
c                       (k=1,nval).  it is modified by scaling
c                       and subtracting the minimum value.  (input)
c               ia(k) = once the floating point data values are
c                       scaled and the reference value has been 
c                       removed they are truncated to integers
c                       and returned to the calling routine via
c                       this array.  (output)
c                 nxy = dimension of a( ).  (input)
c                nval = the number of values in a( ).  (input)
c              iclean = 1 when there are no missing values in a( ).
c                       0 otherwise.
c                  id = the decimal scaling exponent.  (input)
c                  ie = the binary scaling exponent.  (input)
c               xmina = the field minimum value when the original data
c                       are floating point.  (output)
c              jmisss = .true. if there is a secondary missing
c                       value in the data field.  .false. otherwise.
c                       (logical)  (input)
c              jmissp = .true. if there is a primary missing value
c                       in the data field.  .false. otherwise.
c                       (logical)  (input)
c              xmissp = when missing points can be present in the data,
c                       they will have the value xmissp or xmisss when
c                       the data are floating point.  xmissp
c                       is the primary missing value.  (input)
c              xmisss = secondary missing value indicator when the data
c                       are floating point.  (input)
c                 ier = error return
c                         0 = good return
c                       920 = value too large to be packed in 30 bits.
c                       (output)
c
c        local variables
c                   d = the decimal scaling factor.
c                   e = the binary scaling factor.
c                ipos = the position in the array of the first 
c                       non-missing value when there are missing 
c                       values.  (internal)
c                   k = a loop index variable.
c
c        non system subroutines called
c           none.
c
      parameter (mallow=2**30)
c
      logical jmisss,jmissp
c
      dimension a(nxy),ia(nxy)
c
      ier=0
c
      if(iclean.eq.1)then
c
c           there are no missing values in the array.
c           this is the simplest case.
c           process the data straightaway.
c
         if(id.ne.0)then
            d=10.**id
c
            do 10 k=1,nval
               a(k)=a(k)*d
 10         continue
c
         endif
c
         xmina=a(1)
c
         do 20 k=1,nval
            if(a(k).lt.xmina)xmina=a(k)
 20      continue
c
         if(xmina.ne.0)then
c
            do 30 k=1,nval
               a(k)=a(k)-xmina
               if(a(k).gt.mallow)go to 901
c                 a(k) should not exceed 2**30.  this is caused by too large
c                 a range of values.
 30         continue
c
         endif
c
         if(ie.ne.0)then
            e=2.**(-ie)
c
            do 40 k=1,nval
               a(k)=a(k)*e
               ia(k)=nint(a(k)) 
 40         continue
c
         else 
c
            do 45 k=1,nval
               ia(k)=nint(a(k))
 45         continue
c
         endif
c
      else if(jmissp.and..not.jmisss)then
c
c           there are primary missing values in this field. 
         if(id.ne.0)then
            d=10.**id
c
            do 50 k=1,nval
               if(a(k).eq.xmissp)cycle
               a(k)=a(k)*d
 50         continue
c
         endif
c
c           find the first non-missing value and set xmina to that 
c           value.  ipos will end up containing the position
c           of the first non-missing value.
         ipos=1
c
         do 55 k=1,nval
            if(a(k).ne.xmissp)exit
            ipos=ipos+1
 55      continue
            
         xmina=a(ipos)
c
         do 60 k=ipos+1,nval
            if(a(k).eq.xmissp)cycle
            if(a(k).lt.xmina)xmina=a(k)
 60      continue
c
         if(xmina.ne.0.)then
c
            do 70 k=1,nval
               if(a(k).eq.xmissp)cycle
               a(k)=a(k)-xmina
               if(a(k).gt.mallow)go to 901
c                 a(k) should not exceed 2**30.  this is caused by too large
c                 a range of values.
 70         continue
c
         endif
c
         if(ie.ne.0)then
            e=2.**(-ie)
c
            do 80 k=1,nval
               if(a(k).ne.xmissp)a(k)=a(k)*e
               ia(k)=nint(a(k)) 
 80         continue
c
         else
c
            do 85 k=1,nval
               ia(k)=nint(a(k))
 85         continue

         endif
c
      elseif(jmisss.and.jmissp)then
c
c           there are both primary and secondary
c           missing values.
c
         if(id.ne.0)then
            d=10.**id
c
            do 90 k=1,nval
               if((a(k).eq.xmissp).or.(a(k).eq.xmisss))cycle
               a(k)=a(k)*d
 90         continue
c
         endif
c
c           find the first non-missing value
c           and set xmina equal to that value.  ipos
c           will end up containing the position
c           of the first non-missing value.
         ipos=1
c
         do 95 k=1,nval
            if((a(k).ne.xmissp).and.(a(k).ne.xmisss))exit
            ipos=ipos+1
 95      continue
            
         xmina=a(ipos)
c
         do 100 k=ipos+1,nval
            if((a(k).eq.xmissp).or.(a(k).eq.xmisss))cycle
            if(a(k).lt.xmina)xmina=a(k)
 100     continue
c
         if(xmina.ne.0.)then
c
            do 110 k=1,nval
               if((a(k).eq.xmissp).or.(a(k).eq.xmisss))cycle
               a(k)=a(k)-xmina
               if(a(k).gt.mallow)go to 901
c                 a(k) should not exceed 2**30.  this is caused by too large
c                 a range of values.
 110        continue
c
         endif
c
         if(ie.ne.0)then
            e=2.**(-ie)
c
            do 120 k=1,nval
               if((a(k).ne.xmissp).and.(a(k).ne.xmisss))a(k)=a(k)*e
               ia(k)=nint(a(k))
 120        continue
c
         else
c
            do 125 k=1,nval
               ia(k)=nint(a(k))
 125        continue
c
         endif
c
      elseif(jmisss.and..not.jmissp)then
c
c           there are only secondary missing values.
         if(id.ne.0)then
            d=10.**id
c
            do 130 k=1,nval
               if(a(k).eq.xmisss)cycle
               a(k)=a(k)*d
 130        continue
c
         endif
c
c
c           find the first non-missing value
c           and set xmina to that value.  ipos
c           will end up containing the position
c           of the first non-missing value.
         ipos=1
c
         do 135 k=1,nval
            if(a(k).ne.xmisss)exit
            ipos=ipos+1
 135     continue
            
         xmina=a(ipos)
c
         do 140 k=ipos+1,nval
            if(a(k).eq.xmisss)cycle
            if(a(k).lt.xmina)xmina=a(k)
 140     continue
c
         if(xmina.ne.0.)then
c
            do 150 k=1,nval
               if(a(k).eq.xmisss)cycle
               a(k)=a(k)-xmina
               if(a(k).gt.mallow)go to 901
c                 a(k) should not exceed 2**30.  this is caused by too large
c                 a range of values.
 150        continue
c
         endif
c
         if(ie.ne.0)then
            e=2.**(-ie)
c
            do 160 k=1,nval
               if(a(k).ne.xmisss)a(k)=a(k)*e
               ia(k)=nint(a(k))
 160        continue
c
         else
c
            do 165 k=1,nval
               ia(k)=nint(a(k))
 165        continue
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
c     write(12,902)k,a(k),xmina
c902  format(' too large a range of values to pack in prep_flt.',
c    1       'k, a(k), xmina are',i12,2f12.1)
      return1
      end
