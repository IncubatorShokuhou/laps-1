      subroutine pack_gp(kfildo,ic,nxy,is523,minpk,inc,missp,misss,
     1                   jmin,jmax,lbit,nov,ndg,lx,ibit,jbit,kbit,
     2                   novref,lbitref,ier)            
c
c        february 1994   glahn   tdl   mos-2000
c        june     1995   glahn   modified for lmiss error.
c        july     1996   glahn   added misss
c        february 1997   glahn   removed 4 redundant tests for
c                                missp.eq.0; inserted a test to better
c                                handle a string of 9999's
c        february 1997   glahn   added loops to eliminate test for 
c                                misss when misss = 0
c        march    1997   glahn   corrected for secondary missing value
c        march    1997   glahn   corrected for use of local value
c                                of minpk
c        march    1997   glahn   corrected for secondary missing value
c        march    1997   glahn   changed calculating number of bits 
c                                through exponents to an array (improved
c                                overall packing performance by about
c                                35 percent!).  allowed 0 bits for
c                                packing jmin( ), lbit( ), and nov( ).
c        may      1997   glahn   a number of changes for efficiency.
c                                mod functions eliminated and one
c                                ifthen added.  jount removed.
c                                recomputation of bits not made unless
c                                necessary after moving points from
c                                one group to another.  nendb adjusted
c                                to eliminate possibility of very
c                                small group at the end. 
c                                about 8 percent improvement in
c                                overall packing.  iskipa removed;
c                                there is always a group b that can
c                                become group a.  control on size 
c                                of group b (statement below 150)
c                                added.  added adda, and use
c                                of ge and le instead of gt and lt
c                                in loops between 150 and 160.
c                                ibitbs added to shorten trips 
c                                through loop.
c        march    2000   glahn   modified for grib2; changed name from 
c                                packgp
c        january  2001   glahn   comments; ier = 706 substituted for
c                                stops; added return1; removed statement
c                                number 110; added ier and * return
c        november 2001   glahn   changed some diagnostic formats to 
c                                allow printing larger numbers
c        november 2001   glahn   added misslx( ) to put maximum value
c                                into jmin( ) when all values missing
c                                to agree with grib standard.
c        november 2001   glahn   changed two tests on missp and misss
c                                eq 0 to tests on is523.  however,
c                                missp and misss cannot in general be
c                                = 0.
c        november 2001   glahn   added call to reduce; defined itest
c                                before loops to reduce computation;
c                                started large group when all same
c                                value
c        december 2001   glahn   modified and added a few comments
c        january  2002   glahn   removed loop before 150 to determine
c                                a group of all same value
c        january  2002   glahn   changed mallow from 9999999 to 2**30+1,
c                                and made it a parameter
c        march    2002   glahn   added non fatal ier = 716, 717;
c                                removed nendb=nxy above 150;
c                                added iersav=0; comments
c
c        purpose
c            determines groups of variable size, but at least of
c            size minpk, the associated max (jmax( )) and min (jmin( )),
c            the number of bits necessary to hold the values in each
c            group (lbit( )), the number of values in each group
c            (nov( )), the number of bits necessary to pack the jmin( )
c            values (ibit), the number of bits necessary to pack the
c            lbit( ) values (jbit), and the number of bits necessary
c            to pack the nov( ) values (kbit).  the routine is designed
c            to determine the groups such that a small number of bits
c            is necessary to pack the data without excessive
c            computations.  if all values in the group are zero, the
c            number of bits to use in packing is defined as zero when
c            there can be no missing values; when there can be missing
c            values, the number of bits must be at least 1 to have
c            the capability to recognize the missing value.  however,
c            if all values in a group are missing, the number of bits
c            needed is 0, and the unpacker recognizes this.
c            all variables are integer.  even though the groups are 
c            initially of size minpk or larger, an adjustment between
c            two groups (the lookback procedure) may make a group 
c            smaller than minpk.  the control on group size is that
c            the sum of the sizes of the two consecutive groups, each of
c            size minpk or larger, is not decreased.  when determining
c            the number of bits necessary for packing, the largest
c            value that can be accommodated in, say, mbits, is
c            2**mbits-1; this largest value (and the next smallest
c            value) is reserved for the missing value indicator (only)
c            when is523 ne 0.  if the dimension ndg
c            is not large enough to hold all the groups, the local value
c            of minpk is increased by 50 percent.  this is repeated
c            until ndg will suffice.  a diagnostic is printed whenever
c            this happens, which should be very rarely.  if it happens
c            often, ndg in subroutine pack should be increased and
c            a corresponding increase in subroutine unpack made. 
c            considerable code is provided so that no more checking
c            for missing values within loops is done than necessary;
c            the added efficiency of this is relatively minor,
c            but does no harm.  for grib2, the reference value for
c            the length of groups in nov( ) and for the number of
c            bits necessary to pack group values are determined,
c            and subtracted before jbit and kbit are determined.
c
c            when 1 or more groups are large compared to the others,
c            the width of all groups must be as large as the largest.
c            a subroutine reduce breaks up large groups into 2 or
c            more to reduce total bits required.  if reduce should
c            abort, pack_gp will be executed again without the call
c            to reduce.
c
c        data set use 
c           kfildo - unit number for output (print) file. (output) 
c
c        variables in call sequence 
c              kfildo = unit number for output (print) file.  (input)
c               ic( ) = array to hold data for packing.  the values
c                       do not have to be positive at this point, but
c                       must be in the range -2**30 to +2**30 (the 
c                       the value of mallow).  these integer values
c                       will be retained exactly through packing and
c                       unpacking.  (input)
c                 nxy = number of values in ic( ).  also treated
c                       as its dimension.  (input)
c              is523  = missing value management
c                       0=data contains no missing values
c                       1=data contains primary missing values
c                       2=data contains primary and secondary missing values
c                       (input)
c               minpk = the minimum size of each group, except possibly
c                       the last one.  (input)
c                 inc = the number of values to add to an already
c                       existing group in determining whether or not
c                       to start a new group.  ideally, this would be
c                       1, but each time inc values are attempted, the
c                       max and min of the next minpk values must be
c                       found.  this is "a loop within a loop," and
c                       a slightly larger value may give about as good
c                       results with slightly less computational time.
c                       if inc is le 0, 1 is used, and a diagnostic is
c                       output.  note:  it is expected that inc will
c                       equal 1.  the code uses inc primarily in the
c                       loops starting at statement 180.  if inc
c                       were 1, there would not need to be loops
c                       as such.  however, kinc (the local value of
c                       inc) is set ge 1 when near the end of the data
c                       to forestall a very small group at the end. 
c                       (input)
c               missp = when missing points can be present in the data,
c                       they will have the value missp or misss.
c                       missp is the primary missing value and  misss
c                       is the secondary missing value .  these must
c                       not be values that would occur with subtracting
c                       the minimum (reference) value or scaling.
c                       for example, missp = 0 would not be advisable.
c                       (input)
c               misss = secondary missing value indicator (see missp).
c                       (input)
c             jmin(j) = the minimum of each group (j=1,lx).  (output)
c             jmax(j) = the maximum of each group (j=1,lx).  this is
c                       not really needed, but since the max of each
c                       group must be found, saving it here is cheap
c                       in case the user wants it.  (output)
c             lbit(j) = the number of bits necessary to pack each group
c                       (j=1,lx).  it is assumed the minimum of each
c                       group will be removed before packing, and the
c                       values to pack will, therefore, all be positive.
c                       however, ic( ) does not necessarily contain
c                       all positive values.  if the overall minimum
c                       has been removed (the usual case), then ic( )
c                       will contain only positive values.  (output)
c              nov(j) = the number of values in each group (j=1,lx).
c                       (output)
c                 ndg = the dimension of jmin( ), jmax( ), lbit( ), and
c                       nov( ).  (input)
c                  lx = the number of groups determined.  (output)
c                ibit = the number of bits necessary to pack the jmin(j)
c                       values, j=1,lx.  (output)
c                jbit = the number of bits necessary to pack the lbit(j)
c                       values, j=1,lx.  (output)
c                kbit = the number of bits necessary to pack the nov(j)
c                       values, j=1,lx.  (output)
c              novref = reference value for nov( ).  (output)
c             lbitref = reference value for lbit( ).  (output)
c                 ier = error return.
c                       706 = value will not pack in 30 bits--fatal
c                       714 = error in reduce--non-fatal
c                       715 = ngp not large enough in reduce--non-fatal
c                       716 = minpk inceased--non-fatal
c                       717 = inc set = 1--non-fatal
c                       (output)
c                   * = alternate return when ier ne 0 and fatal error.
c
c        internal variables 
c               cfeed = contains the character representation
c                       of a printer form feed.
c               ifeed = contains the integer value of a printer
c                       form feed.
c                kinc = working copy of inc.  may be modified.
c                mina = minimum value in group a.
c                maxa = maximum value in group a.
c               nenda = the place in ic( ) where group a ends.
c              kstart = the place in ic( ) where group a starts.
c               ibita = number of bits needed to hold values in group a.
c                minb = minimum value in group b.
c                maxb = maximum value in group b.
c               nendb = the place in ic( ) where group b ends.
c               ibitb = number of bits needed to hold values in group b.
c                minc = minimum value in group c.
c                maxc = maximum value in group c.
c              ktotal = count of number of values in ic( ) processed.
c               nount = number of values added to group a.
c               lmiss = 0 when is523 = 0.  when packing into a 
c                       specific number of bits, say mbits,
c                       the maximum value that can be handled is
c                       2**mbits-1.  when is523 = 1, indicating
c                       primary missing values, this maximum value 
c                       is reserved to hold the primary missing value 
c                       indicator and lmiss = 1.  when is523 = 2,
c                       the value just below the maximum (i.e.,
c                       2**mbits-2) is reserved to hold the secondary 
c                       missing value indicator and lmiss = 2.
c              lminpk = local value of minpk.  this will be adjusted
c                       upward whenever ndg is not large enough to hold
c                       all the groups.
c              mallow = the largest allowable value for packing.
c              mislla = set to 1 when all values in group a are missing.
c                       this is used to distinguish between a real
c                       minimum when all values are not missing
c                       and a minimum that has been set to zero when
c                       all values are missing.  0 otherwise.
c                       note that this does not distinguish between
c                       primary and secondary missings when secondary
c                       missings are present.  this means that 
c                       lbit( ) will not be zero with the resulting
c                       compression efficiency when secondary missings
c                       are present.  also note that a check has been
c                       made earlier to determine that secondary
c                       missings are really there.
c              misllb = set to 1 when all values in group b are missing.
c                       this is used to distinguish between a real
c                       minimum when all values are not missing
c                       and a minimum that has been set to zero when
c                       all values are missing.  0 otherwise.
c              misllc = performs the same function for group c that
c                       mislla and misllb do for groups b and c,
c                       respectively.
c            ibxx2(j) = an array that when this routine is first entered
c                       is set to 2**j, j=0,30. ibxx2(30) = 2**30, which
c                       is the largest value packable, because 2**31
c                       is larger than the integer word size.
c              ifirst = set by data statement to 0.  changed to 1 on
c                       first
c                       entry when ibxx2( ) is filled.
c               minak = keeps track of the location in ic( ) where the 
c                       minimum value in group a is located.
c               maxak = does the same as minak, except for the maximum.
c               minbk = the same as minak for group b.
c               maxbk = the same as maxak for group b.
c               minck = the same as minak for group c.
c               maxck = the same as maxak for group c.
c                adda = keeps track whether or not an attempt to add
c                       points to group a was made.  if so, then adda
c                       keeps from trying to put one back into b.
c                       (logical)
c              ibitbs = keeps current value if ibitb so that loop
c                       ending at 166 doesn't have to start at
c                       ibitb = 0 every time.
c           misslx(j) = mallow except when a group is all one value (and
c                       lbit(j) = 0) and that value is missing.  in
c                       that case, misslx(j) is missp or misss.  this
c                       gets inserted into jmin(j) later as the 
c                       missing indicator; it can't be put in until
c                       the end, because jmin( ) is used to calculate
c                       the maximum number of bits (ibits) needed to
c                       pack jmin( ).
c        1         2         3         4         5         6         7 x
c
c        non system subroutines called 
c           none
c
      parameter (mallow=2**30+1)
c
      character*1 cfeed
      logical adda
c
      dimension ic(nxy)
      dimension jmin(ndg),jmax(ndg),lbit(ndg),nov(ndg)
      dimension misslx(ndg)
c        misslx( ) is an automatic array.
      dimension ibxx2(0:30)
c
      save ibxx2
c
      data ifeed/12/
      data ifirst/0/
c
      ier=0
      iersav=0
c     call timpr(kfildo,kfildo,'start pack_gp        ')
      cfeed=char(ifeed)
c
      ired=0
c        ired is a flag.  when zero, reduce will be called.
c        if reduce aborts, ired = 1 and is not called.  in
c        this case pack_gp executes again except for reduce.
c
      if(inc.le.0)then
         iersav=717
c        write(kfildo,101)inc
c101     format(/' ****inc ='i8,' not correct in pack_gp.  1 is used.')
      endif
c
c        there will be a restart of pack_gp if subroutine reduce
c        aborts.  this should not happen, but if it does, pack_gp
c        will complete without subroutine reduce.  a non fatal
c        diagnostic return is provided.
c
 102  kinc=max(inc,1)
      lminpk=minpk
c
c         calculate the powers of 2 the first time entered.
c
      if(ifirst.eq.0)then
         ifirst=1
         ibxx2(0)=1
c
         do 104 j=1,30
         ibxx2(j)=ibxx2(j-1)*2
 104     continue
c
      endif
c
c        there will be a restart at 105 is ndg is not large enough.
c        a non fatal diagnostic return is provided.
c
 105  kstart=1
      ktotal=0
      lx=0
      adda=.false.
      lmiss=0
      if(is523.eq.1)lmiss=1
      if(is523.eq.2)lmiss=2
c
c        *************************************
c
c        this section computes statistics for group a.  group a is
c        a group of size lminpk.
c
c        *************************************
c
      ibita=0
      mina=mallow
      maxa=-mallow
      minak=mallow
      maxak=-mallow
c
c        find the min and max of group a.  this will initially be of
c        size lminpk (if there are still lminpk values in ic( )), but
c        will increase in size in increments of inc until a new
c        group is started.  the definition of group a is done here
c        only once (upon initial entry), because a group b can always
c        become a new group a after a is packed, except if lminpk 
c        has to be increased because ndg is too small.  therefore,
c        the separate loops for missing and non-missing here buys
c        almost nothing.
c
      nenda=min(kstart+lminpk-1,nxy)
      if(nxy-nenda.le.lminpk/2)nenda=nxy
c        above statement guarantees the last group is gt lminpk/2 by 
c        making the actual group larger.  if a provision like this is 
c        not included, there will many times be a very small group
c        at the end.  use separate loops for missing and no missing
c        values for efficiency.
c
c        determine whether there is a long string of the same value
c        unless nenda = nxy.  this may allow a large group a to
c        start with, as with missing values.   separate loops for
c        missing options.  this section is only executed once,
c        in determining the first group.  it helps for an array
c        of mostly missing values or of one value, such as
c        radar or precip data.
c
      if(nenda.ne.nxy.and.ic(kstart).eq.ic(kstart+1))then
c           no need to execute if first two values are not equal.
c
         if(is523.eq.0)then
c              this loop is for no missing values.
c
            do 111 k=kstart+1,nxy
c
               if(ic(k).ne.ic(kstart))then
                  nenda=max(nenda,k-1)
                  go to 114
               endif
c
 111        continue
c
            nenda=nxy
c              fall through the loop means all values are the same.
c
         elseif(is523.eq.1)then
c              this loop is for primary missing values only.
c
            do 112 k=kstart+1,nxy
c        
               if(ic(k).ne.missp)then
c
                  if(ic(k).ne.ic(kstart))then
                     nenda=max(nenda,k-1)
                     go to 114
                  endif
c
               endif
c
 112        continue
c
            nenda=nxy
c              fall through the loop means all values are the same.
c
         else
c              this loop is for primary and secondary missing values.
c
            do 113 k=kstart+1,nxy
c        
               if(ic(k).ne.missp.and.ic(k).ne.misss)then
c
                  if(ic(k).ne.ic(kstart))then
                     nenda=max(nenda,k-1)
                     go to 114
                  endif
c
               endif
c
 113        continue
c
            nenda=nxy
c              fall through the loop means all values are the same.
         endif
c
      endif
c
 114  if(is523.eq.0)then
c
         do 115 k=kstart,nenda
         if(ic(k).lt.mina)then
            mina=ic(k)
            minak=k
         endif
         if(ic(k).gt.maxa)then
            maxa=ic(k)
            maxak=k
         endif
 115     continue
c
      elseif(is523.eq.1)then
c
         do 117 k=kstart,nenda
         if(ic(k).eq.missp)go to 117
         if(ic(k).lt.mina)then
            mina=ic(k)
            minak=k
         endif
         if(ic(k).gt.maxa)then
            maxa=ic(k)
            maxak=k
         endif
 117     continue
c
      else
c
         do 120 k=kstart,nenda
         if(ic(k).eq.missp.or.ic(k).eq.misss)go to 120
         if(ic(k).lt.mina)then
            mina=ic(k)
            minak=k
         endif
         if(ic(k).gt.maxa)then
            maxa=ic(k)
            maxak=k
         endif
 120     continue
c
      endif
c
      kounta=nenda-kstart+1
c
c        increment ktotal and find the bits needed to pack the a group.
c
      ktotal=ktotal+kounta
      mislla=0
      if(mina.ne.mallow)go to 125
c        all missing values must be accommodated.
      mina=0
      maxa=0
      mislla=1
      ibitb=0
      if(is523.ne.2)go to 130
c        when all values are missing and there are no
c        secondary missing values, ibita = 0.
c        otherwise, ibita must be calculated.
c
 125  itest=maxa-mina+lmiss
c  
      do 126 ibita=0,30
      if(itest.lt.ibxx2(ibita))go to 130
c***        this test is the same as:
c***     if(maxa-mina.lt.ibxx2(ibita)-lmiss)go to 130
 126  continue
c
c     write(kfildo,127)maxa,mina
c127  format(' ****error in pack_gp.  value will not pack in 30 bits.',
c    1       '  maxa ='i13,'  mina ='i13,'.  error at 127.')
      ier=706
      go to 900
c
 130  continue
c
c***d     write(kfildo,131)kounta,ktotal,mina,maxa,ibita,mislla
c***d131  format(' at 130, kounta ='i8,'  ktotal ='i8,'  mina ='i8,
c***d    1       '  maxa ='i8,'  ibita ='i3,'  mislla ='i3) 
c
 133  if(ktotal.ge.nxy)go to 200
c
c        *************************************
c
c        this section computes statistics for group b.  group b is a
c        group of size lminpk immediately following group a.
c
c        *************************************
c
 140  minb=mallow
      maxb=-mallow
      minbk=mallow
      maxbk=-mallow
      ibitbs=0
      mstart=ktotal+1
c
c        determine whether there is a long string of the same value.
c        this works when there are no missing values.
c
      nendb=1
c
      if(mstart.lt.nxy)then
c
         if(is523.eq.0)then
c              this loop is for no missing values.
c
            do 145 k=mstart+1,nxy
c
               if(ic(k).ne.ic(mstart))then
                  nendb=k-1
                  go to 150
               endif
c
 145        continue
c
            nendb=nxy
c              fall through the loop means all remaining values
c              are the same.
         endif
c
      endif
c         
 150  nendb=max(nendb,min(ktotal+lminpk,nxy))
c**** 150  nendb=min(ktotal+lminpk,nxy)
c
      if(nxy-nendb.le.lminpk/2)nendb=nxy
c        above statement guarantees the last group is gt lminpk/2 by 
c        making the actual group larger.  if a provision like this is 
c        not included, there will many times be a very small group
c        at the end.  use separate loops for missing and no missing
c
c        use separate loops for missing and no missing values
c        for efficiency.
c
      if(is523.eq.0)then
c              
         do 155 k=mstart,nendb
         if(ic(k).le.minb)then
            minb=ic(k)
c              note le, not lt.  lt could be used but then a 
c              recompute over the whole group would be needed
c              more often.  same reasoning for ge and other
c              loops below.
            minbk=k
         endif
         if(ic(k).ge.maxb)then
            maxb=ic(k)
            maxbk=k
         endif
 155     continue
c
      elseif(is523.eq.1)then
c
         do 157 k=mstart,nendb
         if(ic(k).eq.missp)go to 157
         if(ic(k).le.minb)then
            minb=ic(k)
            minbk=k
         endif
         if(ic(k).ge.maxb)then
            maxb=ic(k)
            maxbk=k
         endif
 157     continue
c
      else
c
         do 160 k=mstart,nendb
         if(ic(k).eq.missp.or.ic(k).eq.misss)go to 160
         if(ic(k).le.minb)then
            minb=ic(k)
            minbk=k
         endif
         if(ic(k).ge.maxb)then
            maxb=ic(k)
            maxbk=k
         endif
 160     continue
c
      endif
c
      kountb=nendb-ktotal
      misllb=0
      if(minb.ne.mallow)go to 165
c        all missing values must be accommodated.
      minb=0
      maxb=0
      misllb=1
      ibitb=0
c
      if(is523.ne.2)go to 170
c        when all values are missing and there are no secondary
c        missing values, ibitb = 0.  otherwise, ibitb must be
c        calculated.
c
 165  do 166 ibitb=ibitbs,30
         if(maxb-minb.lt.ibxx2(ibitb)-lmiss)go to 170
 166  continue
c
c     write(kfildo,167)maxb,minb
c167  format(' ****error in pack_gp.  value will not pack in 30 bits.',
c    1       '  maxb ='i13,'  minb ='i13,'.  error at 167.')
      ier=706
      go to 900
c
c        compare the bits needed to pack group b with those needed
c        to pack group a.  if ibitb ge ibita, try to add to group a.
c        if not, try to add a's points to b, unless addition to a
c        has been done.  this latter is controlled with adda.
c
 170  continue
c
c***d     write(kfildo,171)kounta,ktotal,mina,maxa,ibita,mislla,
c***d    1                               minb,maxb,ibitb,misllb
c***d171  format(' at 171, kounta ='i8,'  ktotal ='i8,'  mina ='i8,
c***d    1       '  maxa ='i8,'  ibita ='i3,'  mislla ='i3,
c***d    2       '  minb ='i8,'  maxb ='i8,'  ibitb ='i3,'  misllb ='i3)  
c
      if(ibitb.ge.ibita)go to 180
      if(adda)go to 200
c
c        *************************************
c
c        group b requires less bits than group a.  put as many of a's
c        points into b as possible without exceeding the number of
c        bits necessary to pack group b.
c
c        *************************************
c
      kounts=kounta
c        kounta refers to the present group a.
      mintst=minb
      maxtst=maxb
      mintstk=minbk
      maxtstk=maxbk
c
c        use separate loops for missing and no missing values
c        for efficiency.
c
      if(is523.eq.0)then
c 
         do 1715 k=ktotal,kstart,-1
c           start with the end of the group and work backwards.
         if(ic(k).lt.minb)then
            mintst=ic(k)
            mintstk=k
         elseif(ic(k).gt.maxb)then
            maxtst=ic(k)
            maxtstk=k
         endif
         if(maxtst-mintst.ge.ibxx2(ibitb))go to 174
c           note that for this loop, lmiss = 0.
         minb=mintst
         maxb=maxtst
         minbk=mintstk
         maxbk=maxtstk
         kounta=kounta-1
c           there is one less point now in a.
 1715    continue  
c
      elseif(is523.eq.1)then            
c
         do 1719 k=ktotal,kstart,-1
c           start with the end of the group and work backwards.
         if(ic(k).eq.missp)go to 1718
         if(ic(k).lt.minb)then
            mintst=ic(k)
            mintstk=k
         elseif(ic(k).gt.maxb)then
            maxtst=ic(k)
            maxtstk=k
         endif
         if(maxtst-mintst.ge.ibxx2(ibitb)-lmiss)go to 174
c           for this loop, lmiss = 1.
         minb=mintst
         maxb=maxtst
         minbk=mintstk
         maxbk=maxtstk
         misllb=0
c           when the point is non missing, misllb set = 0.
 1718    kounta=kounta-1
c           there is one less point now in a.
 1719    continue  
c
      else             
c
         do 173 k=ktotal,kstart,-1
c           start with the end of the group and work backwards.
         if(ic(k).eq.missp.or.ic(k).eq.misss)go to 1729
         if(ic(k).lt.minb)then
            mintst=ic(k)
            mintstk=k
         elseif(ic(k).gt.maxb)then
            maxtst=ic(k)
            maxtstk=k
         endif
         if(maxtst-mintst.ge.ibxx2(ibitb)-lmiss)go to 174
c           for this loop, lmiss = 2.
         minb=mintst
         maxb=maxtst
         minbk=mintstk
         maxbk=maxtstk
         misllb=0
c           when the point is non missing, misllb set = 0.
 1729    kounta=kounta-1
c           there is one less point now in a.
 173     continue  
c
      endif
c
c        at this point, kounta contains the number of points to close
c        out group a with.  group b now starts with kstart+kounta and
c        ends with nendb.  minb and maxb have been adjusted as
c        necessary to reflect group b (even though the number of bits
c        needed to pack group b have not increased, the end points
c        of the range may have).
c
 174  if(kounta.eq.kounts)go to 200
c        on transfer, group a was not changed.  close it out.
c
c        one or more points were taken out of a.  range and ibita
c        may have to be recomputed; ibita could be less than
c        originally computed.  in fact, group a can now contain
c        only one point and be packed with zero bits
c        (unless misss ne 0).
c
      nouta=kounts-kounta
      ktotal=ktotal-nouta
      kountb=kountb+nouta
      if(nenda-nouta.gt.minak.and.nenda-nouta.gt.maxak)go to 200
c        when the above test is met, the min and max of the 
c        current group a were within the old group a, so the
c        range and ibita do not need to be recomputed.
c        note that minak and maxak are no longer needed.
      ibita=0
      mina=mallow
      maxa=-mallow
c
c        use separate loops for missing and no missing values
c        for efficiency.
c
      if(is523.eq.0)then
c 
         do 1742 k=kstart,nenda-nouta
         if(ic(k).lt.mina)then
            mina=ic(k)
         endif
         if(ic(k).gt.maxa)then
            maxa=ic(k)
         endif
 1742    continue
c
      elseif(is523.eq.1)then 
c
         do 1744 k=kstart,nenda-nouta
         if(ic(k).eq.missp)go to 1744
         if(ic(k).lt.mina)then
            mina=ic(k)
         endif
         if(ic(k).gt.maxa)then
            maxa=ic(k)
         endif
 1744    continue
c
      else 
c
         do 175 k=kstart,nenda-nouta
         if(ic(k).eq.missp.or.ic(k).eq.misss)go to 175
         if(ic(k).lt.mina)then
            mina=ic(k)
         endif
         if(ic(k).gt.maxa)then
            maxa=ic(k)
         endif
 175     continue
c
      endif
c
      mislla=0
      if(mina.ne.mallow)go to 1750
c        all missing values must be accommodated.
      mina=0
      maxa=0
      mislla=1
      if(is523.ne.2)go to 177
c        when all values are missing and there are no secondary
c        missing values ibita = 0 as originally set.  otherwise,
c        ibita must be calculated.
c
 1750 itest=maxa-mina+lmiss
c
      do 176 ibita=0,30
      if(itest.lt.ibxx2(ibita))go to 177
c***        this test is the same as:
c***         if(maxa-mina.lt.ibxx2(ibita)-lmiss)go to 177
 176  continue
c
c     write(kfildo,1760)maxa,mina
c1760 format(' ****error in pack_gp.  value will not pack in 30 bits.',
c    1       '  maxa ='i13,'  mina ='i13,'.  error at 1760.')
      ier=706
      go to 900
c
 177  continue
      go to 200
c
c        *************************************
c
c        at this point, group b requires as many bits to pack as groupa.
c        therefore, try to add inc points to group a without increasing
c        ibita.  this augmented group is called group c.
c
c        *************************************
c
 180  if(mislla.eq.1)then
         minc=mallow
         minck=mallow
         maxc=-mallow
         maxck=-mallow
      else
         minc=mina
         maxc=maxa
         minck=minak
         maxck=minak
      endif
c
      nount=0
      if(nxy-(ktotal+kinc).le.lminpk/2)kinc=nxy-ktotal
c        above statement constrains the last group to be not less than
c        lminpk/2 in size.  if a provision like this is not included,
c        there will many times be a very small group at the end.
c
c        use separate loops for missing and no missing values
c        for efficiency.  since kinc is usually 1, using separate
c        loops here doesn't buy much.  a missing value will always
c        transfer back to group a.
c
      if(is523.eq.0)then
c
         do 185 k=ktotal+1,min(ktotal+kinc,nxy)
         if(ic(k).lt.minc)then
            minc=ic(k)
            minck=k
         endif
         if(ic(k).gt.maxc)then
            maxc=ic(k)
            maxck=k
         endif
         nount=nount+1
 185     continue
c
      elseif(is523.eq.1)then
c
         do 187 k=ktotal+1,min(ktotal+kinc,nxy)
         if(ic(k).eq.missp)go to 186
         if(ic(k).lt.minc)then
            minc=ic(k)
            minck=k
         endif
         if(ic(k).gt.maxc)then
            maxc=ic(k)
            maxck=k
         endif
 186     nount=nount+1
 187     continue
c
      else
c
         do 190 k=ktotal+1,min(ktotal+kinc,nxy)
         if(ic(k).eq.missp.or.ic(k).eq.misss)go to 189
         if(ic(k).lt.minc)then
            minc=ic(k)
            minck=k
         endif
         if(ic(k).gt.maxc)then
            maxc=ic(k)
            maxck=k
         endif
 189     nount=nount+1
 190     continue
c
      endif
c
c***d     write(kfildo,191)kounta,ktotal,mina,maxa,ibita,mislla,
c***d    1   minc,maxc,nount,ic(ktotal),ic(ktotal+1)
c***d191  format(' at 191, kounta ='i8,'  ktotal ='i8,'  mina ='i8,
c***d    1       '  maxa ='i8,'  ibita ='i3,'  mislla ='i3,
c***d    2       '  minc ='i8,'  maxc ='i8,
c***d    3       '  nount ='i5,'  ic(ktotal) ='i9,'  ic(ktotal+1) =',i9) 
c
c        if the number of bits needed for group c is gt ibita,
c        then this group a is a group to pack.
c
      if(minc.eq.mallow)then
         minc=mina
         maxc=maxa
         minck=minak
         maxck=maxak
         misllc=1
         go to 195
c           when the new value(s) are missing, they can always
c           be added.
c
      else
         misllc=0
      endif
c
      if(maxc-minc.ge.ibxx2(ibita)-lmiss) go to 200
c
c        the bits necessary for group c has not increased from the
c        bits necessary for group a.  add this point(s) to group a.
c        compute the next group b, etc., unless all points have been
c        used.
c 
 195  ktotal=ktotal+nount
      kounta=kounta+nount
      mina=minc
      maxa=maxc
      minak=minck
      maxak=maxck
      mislla=misllc
      adda=.true.
      if(ktotal.ge.nxy)go to 200
c
      if(minbk.gt.ktotal.and.maxbk.gt.ktotal)then
         mstart=nendb+1
c           the max and min of group b were not from the points
c           removed, so the whole group does not have to be looked
c           at to determine the new max and min.  rather start
c           just beyond the old nendb.
         ibitbs=ibitb
         nendb=1
         go to 150
      else
         go to 140
      endif
c
c        *************************************
c
c        group a is to be packed.  store values in jmin( ), jmax( ),
c        lbit( ), and nov( ).
c
c        *************************************
c
 200  lx=lx+1
      if(lx.le.ndg)go to 205
      lminpk=lminpk+lminpk/2
c     write(kfildo,201)ndg,lminpk,lx
c201  format(' ****ndg ='i5,' not large enough.',
c    1       '  lminpk is increased to 'i3,' for this field.'/
c    2       '  lx = 'i10)
      iersav=716
      go to 105
c
 205  jmin(lx)=mina
      jmax(lx)=maxa
      lbit(lx)=ibita
      nov(lx)=kounta
      kstart=ktotal+1
c
      if(mislla.eq.0)then
         misslx(lx)=mallow
      else
         misslx(lx)=ic(ktotal)
c           ic(ktotal) was the last value processed.  if mislla ne 0,
c           this must be the missing value for this group.
      endif
c
c***d     write(kfildo,206)mislla,ic(ktotal),ktotal,lx,jmin(lx),jmax(lx),
c***d    1                 lbit(lx),nov(lx),misslx(lx)
c***d206  format(' at 206,  mislla ='i2,'  ic(ktotal) ='i5,'  ktotal ='i8,
c***d    1       '  lx ='i6,'  jmin(lx) ='i8,'  jmax(lx) ='i8,
c***d    2       '  lbit(lx) ='i5,'  nov(lx) ='i8,'  misslx(lx) =',i7) 
c
      if(ktotal.ge.nxy)go to 209
c
c        the new group a will be the previous group b.  set limits, etc.
c
      ibita=ibitb
      mina=minb
      maxa=maxb
      minak=minbk
      maxak=maxbk
      mislla=misllb
      nenda=nendb
      kounta=kountb
      ktotal=ktotal+kounta
      adda=.false.
      go to 133
c
c        *************************************
c
c        calculate ibit, the number of bits needed to hold the group
c        minimum values.
c
c        *************************************
c
 209  ibit=0
c
      do 220 l=1,lx
 210  if(jmin(l).lt.ibxx2(ibit))go to 220
      ibit=ibit+1
      go to 210
 220  continue
c
c        insert the value in jmin( ) to be used for all missing
c        values when lbit( ) = 0.  when secondary missing 
c        values can be present, lbit(l) will not = 0.
c
      if(is523.eq.1)then
c
         do 226 l=1,lx
c   
         if(lbit(l).eq.0)then
c
            if(misslx(l).eq.missp)then
               jmin(l)=ibxx2(ibit)-1
            endif
c
         endif
c
 226     continue
c
      endif
c
c        *************************************
c
c        calculate jbit, the number of bits needed to hold the bits
c        needed to pack the values in the groups.  but find and
c        remove the reference value first.
c
c        *************************************
c
c     write(kfildo,228)cfeed,lx
c228  format(a1,/' *****************************************'
c    1          /' the group widths lbit( ) for ',i8,' groups'
c    2          /' *****************************************')
c     write(kfildo,229) (lbit(j),j=1,min(lx,100))
c229  format(/' '20i6)
c
      lbitref=lbit(1)
c
      do 230 k=1,lx
      if(lbit(k).lt.lbitref)lbitref=lbit(k)
 230  continue
c
      if(lbitref.ne.0)then
c
         do 240 k=1,lx
         lbit(k)=lbit(k)-lbitref
 240     continue
c
      endif
c
c     write(kfildo,241)cfeed,lbitref
c241  format(a1,/' *****************************************'
c    1          /' the group widths lbit( ) after removing reference ',
c    2             i8,
c    3          /' *****************************************')
c     write(kfildo,242) (lbit(j),j=1,min(lx,100))
c242  format(/' '20i6)
c
      jbit=0
c
      do 320 k=1,lx
 310  if(lbit(k).lt.ibxx2(jbit))go to 320
      jbit=jbit+1
      go to 310
 320  continue
c
c        *************************************
c
c        calculate kbit, the number of bits needed to hold the number
c        of values in the groups.  but find and remove the
c        reference first.
c
c        *************************************
c
c     write(kfildo,321)cfeed,lx
c321  format(a1,/' *****************************************'
c    1          /' the group sizes nov( ) for ',i8,' groups'
c    2          /' *****************************************')
c     write(kfildo,322) (nov(j),j=1,min(lx,100))
c322  format(/' '20i6)
c
      novref=nov(1)
c
      do 400 k=1,lx
      if(nov(k).lt.novref)novref=nov(k)
 400  continue
c
      if(novref.gt.0)then
c
         do 405 k=1,lx
         nov(k)=nov(k)-novref
 405     continue
c
      endif
c
c     write(kfildo,406)cfeed,novref
c406  format(a1,/' *****************************************'
c    1          /' the group sizes nov( ) after removing reference ',i8,
c    2          /' *****************************************')
c     write(kfildo,407) (nov(j),j=1,min(lx,100))
c407  format(/' '20i6)
c     write(kfildo,408)cfeed
c408  format(a1,/' *****************************************'
c    1          /' the group references jmin( )'
c    2          /' *****************************************')
c     write(kfildo,409) (jmin(j),j=1,min(lx,100))
c409  format(/' '20i6)
c
      kbit=0
c
      do 420 k=1,lx
 410  if(nov(k).lt.ibxx2(kbit))go to 420
      kbit=kbit+1
      go to 410
 420  continue
c
c        determine whether the group sizes should be reduced
c        for space efficiency.
c
      if(ired.eq.0)then
         call reduce(kfildo,jmin,jmax,lbit,nov,lx,ndg,ibit,jbit,kbit,
     1               novref,ibxx2,ier)
c
         if(ier.eq.714.or.ier.eq.715)then
c              reduce has aborted.  reexecute pack_gp without reduce.
c              provide for a non fatal return from reduce.  
            iersav=ier
            ired=1
            ier=0
            go to 102 
         endif
c
      endif         
c
c     call timpr(kfildo,kfildo,'end   pack_gp        ')
      if(iersav.ne.0)then
         ier=iersav
         return
      endif
c
c 900  if(ier.ne.0)return1
c
 900  return
      end
