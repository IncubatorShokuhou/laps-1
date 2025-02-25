      subroutine pk_missp(kfildo,ia,ib,nxy,missp,minpk,ifirst,
     1                    isecond,imina,numoctet,second,ier)
c
c        april    1997   glahn   modified packxx to handle primary missing
c                                values
c        march    2000   glahn   revised for grib2;
c                                changed name from packyy
c        january  2001   glahn   comments; eliminated array ibxx2( );
c                                changed ier = 1100 to 920
c        november 2001   glahn   added ier=0 at beginning
c        november 2001   glahn   changed numoctet=(numbit/8)+1 to
c                                numoctet=(numbit+7)/8 near the end
c        january  2002   glahn   added test for 1st order diff ge
c                                2**15.; modified comment for ier = 920
c        january  2002   glahn   added mallow = 2**30+1 and used instead
c                                of 999999
c        december 2002   taylor/glahn   changed statement 225 from
c                                'if(ivalue.lt.0)numbit=numbit+1' to
c                                numbit=numbit+1 to allow for sign
c
c        purpose
c            determines whether or not to use second order
c            differences or original values to pack.  original values
c            are indicated when the average range of consecutive groups
c            of size minpk of the second order differences is
c            larger than the average range of consecutive groups of
c            size minpk of the original values.  note that this
c            procedure does not in general use the same groups as are
c            used in the actual packing because it does not employ the
c            lookback procedure.  this algorithm is relatively cheap
c            and gives a result good enough for the purpose.  this
c            routine is used when there can be missing primary
c            values (i.e., a missing value indicator, missp, can be
c            in the data).
c
c        data set use
c           kfildo - unit number for output (print) file. (output)
c
c        variables
c              kfildo = unit number for output (print) file.  (input) 
c               ia(k) = holds the nxy original values on input
c                       (k=1,nxy).
c                       holds the nxy second order differences on
c                       output when second order differences are to be
c                       used.  in that case, second is .true.  since
c                       there are only nxy-2 values, the first 2
c                       non-missing
c                       values are dummy.  (input-output)
c               ib(k) = work array (k=1,nxy).  (internal)
c                 nxy = number of values in ia( ) on input.  on return,
c                       nxy will also be the number of values in ia( ).
c                       used as dimension of ia( ), ib( ), and iwork( ).
c                       (input)
c               missp = primary missing value indicator.  (input)
c               minpk = increment in which ranges will be computed.
c                       (input)
c              ifirst = ifirst is the first original value.  (output)
c             isecond = isecond is the second original value. (output)
c               imina = the overall minimum value of the second order
c                       spatial differences. (output)
c            numoctet = the minimum number of octets required to
c                       pack each of the values represented by ifirst,
c                       isecond, and imina. (output)
c              second = true (false) when second order differences
c                       are (are not) to be used.  (output)
c                 ier = contains any error codes generated by this
c                       routine
c                         0 = good return
c                       920 = a value larger than what can be packed
c                             into 30 bits has been encountered. 
c
c             local variables
c                avgr = contains the average range of nxy values
c                       in increments of minpk for the original
c                       data field in ia( ). 
c               avgr2 = contains the average range of nxy values
c                       in increments of minpk for the data field
c                       with second order differences in ib( ).
c               cfeed = contains the character representation
c                       of a printer form feed.
c              iakeep = used in determining second order
c                       differences and maintaining the locations
c                       of missing values.
c              ickeep = used in determining second order
c                       differences and maintaining the locations
c                       of missing values.
c               ifeed = contains the integer value of a printer
c                       form feed.
c              irange = the difference of jmax and jmin.
c             isecpos = the position of the second non-missing value
c                       in the data field.
c              ivalue = the value of ifirst, isecond, and imina
c                       that will need the most octets to store.
c            iwork(k) = work array. this array contains all of the
c                       computations until we are certain that
c                       we want to use second order differences
c                       (k=1,nxy).  (automatic array)
c           jmax,jmin = used to keep track of the max and min values 
c                       of a particular group of minpk values.
c              jfirst = contains the absolute value of the first
c                       non-missing value in the original data 
c                       field.
c            jlargest = contains the largest of jfirst, jsecond, and
c                       jmina for the purpose of determining
c                       how many octets to store all three values
c                       in.
c               jmina = contains the absolute value of the minimum
c                       of the field of second order differences.
c             jsecond = contains the absolute value of the second
c                       non-missing value in the original data
c                       field.
c                   k = loop index variable.
c              kfirst = flag indicating whether the first non-missing
c                       value of the field has been found.
c               kount = a counting/looping variable.
c             ksecond = flag indicating whether the second non-missing
c                       value of the field has been found.
c                sumr = used in computing the average range of the
c                       original nxy values in increments of minpk.
c               large = 2**15, the largest 1st order difference
c                       tolerated.  if the 1st order difference 
c                       exceeds this value, it is likely the 2nd
c                       order difference will overflow the integer
c                       computation and cause an undetected error.
c
c        non system subroutines called
c           none
c
      parameter (mallow=2**30+1)
      parameter (large=2**15)
c
      character*1 cfeed
      logical second,kfirst,ksecond
c
      dimension ia(nxy),ib(nxy),iwork(nxy)
c
      data ifeed/12/
c
      ier=0
c
c        initialize the form feed character.
      cfeed=char(ifeed)
c
c        initialize ifirst, isecond, ib(1), and ickeep.
      ifirst=ia(1)
      isecond=missp
      ib(1)=0
      if(ia(1).eq.missp)ib(1)=missp
      ickeep=ia(1)
c
c        compute first order differences.
c
      do 120 k=2,nxy
c
         if(ickeep.eq.missp)then
            ib(k)=missp
            ickeep=ia(k)
         else
c
            if(ia(k).eq.missp)then
               ib(k)=missp
            else
               ib(k)=ia(k)-ickeep
c         
               if(ib(k).gt.large)then
                  second=.false.
c                 write(kfildo,119)            
c119              format(/' 1st order difference exceeds 2**15.',
c    1                    '  second order differences are not pakced.')
                  go to 900
               endif
c
               ickeep=ia(k)
            endif
c
         endif
c
         if(ifirst.eq.missp)then
c
            if(ia(k).ne.missp)then
               ifirst=ia(k)
               ib(k)=0
            endif
c
         else if(isecond.eq.missp)then
c
            if(ia(k).ne.missp)then
               isecond=ia(k)
            endif
c
         endif
c
 120  continue
c
c***d     write(kfildo,121)(ia(k),k=1,10000)
c***d     write(kfildo,121)(ia(k),k=268351,268500)
c***d121  format(/' original scaled values'/(' '25i4))
c
c***d     write(kfildo,1210)ifirst,ifod
c***d1210 format(/' ifirst ='i11,'     ifod ='i11)
c***d     write(kfildo,122)(ib(k),k=268351,268500)
c***d122  format(/ 'first order differences'/(' '25i4))
c
c        compute second order differences.
c
      iakeep=missp
      kfirst=.true.
      ksecond=.true.
c
      do 130 k=1,nxy
c
         if(kfirst)then
            iwork(k)=ib(k)
            if(ib(k).ne.missp)kfirst=.false.
         else
c
            if(iakeep.eq.missp)then
c
               if(ib(k).ne.missp)then
                  iakeep=ib(k)
                  if(ksecond)then
                     ib(k)=0
                     ksecond=.false.
                     isecpos=k
                  endif
                  iwork(k)=ib(k)
               else
                  iwork(k)=missp
               endif
c
            else
c
               if(ib(k).eq.missp)then
                  iwork(k)=missp
               else
                  iwork(k)=ib(k)-iakeep
                  imina=iwork(k)
                  iakeep=ib(k)
               endif
c
            endif
c
         endif
c
 130  continue
c
c***d     write(kfildo,136)cfeed
c***d136  format(1a/' *********************'
c***d    1          /' 2nd order differences'
c***d    2          /' *********************')
c***d     write(kfildo,137) (iwork(k),k=268351,268500)
c***d137  format(/(' '25i4))
c
c        compute the minimum of the field of second order
c        differences and subtract it from of each of the
c        values in the field.
c
      do 132 k=isecpos+1,nxy
         if(iwork(k).eq.missp)cycle
         if(iwork(k).lt.imina)imina=iwork(k)
 132  continue
c
      do 134 k=isecpos+1,nxy
         if(iwork(k).eq.missp)cycle
         iwork(k)=iwork(k)-imina
 134  continue
c
c***d     write(kfildo,131)(iwork(k),k=268351,268500)
c***d131  format(/' second order differences'/(' '25i4))
c
c        compute average range of nxy original values in increments of 
c        minpk.
c
      sumr=0
      kount=0
c 
      do 140 k=1,nxy,minpk
         jmin=mallow
         jmax=-mallow
         if(k+minpk-1.gt.nxy)go to 140
c           the last group may be very small and not be representative
c           of the range.
c
         do 135 j=k,k+minpk-1
            if(ia(j).eq.missp)cycle
            if(ia(j).gt.jmax)jmax=ia(j)
            if(ia(j).lt.jmin)jmin=ia(j)
 135     continue
c
         kount=kount+1
c
         if(jmin.eq.mallow)then
            irange=0
         else
            irange=jmax-jmin
            sumr=sumr+irange 
         endif   
c
 140  continue
c
      avgr=99999.
      if(kount.ne.0)avgr=sumr/kount
c
c        compute average range of nxy-4 2nd order values in 
c        increments of minpk.   don't use the first 2 values,
c        because they could be based on missing values.
c
      sumr=0
      kount=0
c 
      do 150 k=3,nxy-2,minpk
         jmin=mallow
         jmax=-mallow
         if(k+minpk-1.gt.nxy-2)go to 150
c        the last group may be very small and not be representative of
c        the range.
c
         do 145 j=k,k+minpk-1
            if(ib(j).eq.missp)go to 145
            if(ib(j).gt.jmax)jmax=ib(j)
            if(ib(j).lt.jmin)jmin=ib(j)
 145     continue
c
         kount=kount+1
c
         if(jmin.eq.mallow)then
            irange=0
         else
            irange=jmax-jmin
            sumr=sumr+irange 
         endif   
c 
 150  continue
c
      avgr2=99999.
      if(kount.ne.0)avgr2=sumr/kount
c
      second=.false.
c     write(kfildo,155)avgr,avgr2
c155  format(/' average range of original scaled values =  'f10.2/
c    1        ' average range of second order differences ='f10.2)
      if(avgr2.ge.avgr)go to 900
c        second order differences not to be used.
c
c        transfer iwork( ) to ia( ).
      kount=0
c
      do 180 k=1,nxy
         ia(k)=iwork(k)
 180  continue
c
c***d     write(kfildo,185)(ia(k),k=1,100)
c***d185  format(/' second order differences to be packed'/(' '10i11))
c
c           since we are doing second order differences
c           take the time here to determine the minimum
c           number of bytes (not bits) need to contain 
c           the first original value, the second original
c           value, and the minimum of the field of second
c           order differences.
c
      numbit=0
      jmina=abs(imina)
      jfirst=abs(ifirst)
      jsecond=abs(isecond)
c
c           find the largest absolute value of the 
c           the three variables imina, ifirst, and
c           isecond.
c
      if(jmina.gt.jfirst)then
         jlargest=jmina
         ivalue=imina
      else
         jlargest=jfirst
         ivalue=ifirst
      endif
c
      if(jlargest.lt.jsecond)then
         jlargest=jsecond
         ivalue=isecond
      endif
c
      ibxx2=2
c
      do 220 k=1,30
         numbit=k
         if(jlargest.lt.ibxx2)go to 225
         ibxx2=ibxx2*2
 220  continue
c
c        drop through here means the number of bits
c        exceeds 30 and the value to pack exceeds 2**30.
c        2**31 is too large a value to pack.  2**30 
c        is used as the limit throughout the packer,
c        although a value of 2**31-1 could probably
c        be accommodated.
      ier=920
c        this is a fatal error.
      go to 900
c
 225  numbit=numbit+1
c        allow for sign bit.
c
c        find the smallest number of octets 
c        that will contain the values.  
      numoctet=(numbit+7)/8
      second=.true.
c     write(kfildo,899)jmina,jfirst,jsecond,jlargest,
c    1                 ivalue,numbit,numoctet
c899  format(/' at 899 in pk_missp--jmina,jfirst,jsecond,jlargest,',
c    1        'ivalue,numbit,numoctet'/20x,7i10)
 900  return 
      end

