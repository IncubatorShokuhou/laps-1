      subroutine reduce(kfildo,jmin,jmax,lbit,nov,lx,ndg,ibit,jbit,kbit,
     1                  novref,ibxx2,ier)            
c
c        november 2001   glahn   tdl   grib2
c        march    2002   glahn   comment ier = 715
c        march    2002   glahn   modified to accommodate lx=1 on entry
c
c        purpose
c            determines whether the number of groups should be
c            increased in order to reduce the size of the large
c            groups, and to make that adjustment.  by reducing the
c            size of the large groups, less bits may be necessary
c            for packing the group sizes and all the information
c            about the groups.
c
c            the reference for nov( ) was removed in the calling
c            routine so that kbit could be determined.  this
c            furnishes a starting point for the iterations in reduce.
c            however, the reference must be considered.
c
c        data set use 
c           kfildo - unit number for output (print) file. (output) 
c
c        variables in call sequence 
c              kfildo = unit number for output (print) file.  (input)
c             jmin(j) = the minimum of each group (j=1,lx).  it is
c                       possible after splitting the groups, jmin( )
c                       will not be the minimum of the new group.
c                       this doesn't matter; jmin( ) is really the
c                       group reference and doesn't have to be the
c                       smallest value.  (input/output)
c             jmax(j) = the maximum of each group (j=1,lx). 
c                       (input/output)
c             lbit(j) = the number of bits necessary to pack each group
c                       (j=1,lx).  (input/output)
c              nov(j) = the number of values in each group (j=1,lx).
c                       (input/output)
c                  lx = the number of groups.  this will be increased
c                       if groups are split.  (input/output)
c                 ndg = the dimension of jmin( ), jmax( ), lbit( ), and
c                       nov( ).  (input)
c                ibit = the number of bits necessary to pack the jmin(j)
c                       values, j=1,lx.  (input)
c                jbit = the number of bits necessary to pack the lbit(j)
c                       values, j=1,lx.  (input)
c                kbit = the number of bits necessary to pack the nov(j)
c                       values, j=1,lx.  if the groups are split, kbit
c                       is reduced.  (input/output)
c              novref = reference value for nov( ).  (input)
c            ibxx2(j) = 2**j (j=0,30).  (input)
c                 ier = error return.  (output)
c                         0 = good return.
c                       714 = problem in algorithm.  reduce aborted.
c                       715 = ngp not large enough.  reduce aborted.
c           ntotbt(j) = the total bits used for the packing bits j
c                       (j=1,30).  (internal)
c            nboxj(j) = new boxes needed for the packing bits j
c                       (j=1,30).  (internal)
c           newbox(l) = number of new boxes (groups) for each original
c                       group (l=1,lx) for the current j.  (automatic)
c                       (internal)
c          newboxp(l) = same as newbox( ) but for the previous j.
c                       this eliminates recomputation.  (automatic)
c                       (internal)
c               cfeed = contains the character representation
c                       of a printer form feed.  (character) (internal)
c               ifeed = contains the integer value of a printer
c                       form feed.  (internal)
c              iorigb = the original number of bits necessary
c                       for the group values.  (internal)
c        1         2         3         4         5         6         7 x
c
c        non system subroutines called 
c           none
c
      character*1 cfeed
c
      dimension jmin(ndg),jmax(ndg),lbit(ndg),nov(ndg)
      dimension newbox(ndg),newboxp(ndg)
c        newbox( ) and newboxp( ) are automatic arrays.
      dimension ntotbt(31),nboxj(31)
      dimension ibxx2(0:30)
c
      data ifeed/12/
c
      ier=0
      if(lx.eq.1)go to 410
c        if there is only one group, return.
c
      cfeed=char(ifeed)
c
c        initialize number of new boxes per group to zero.
c
      do 110 l=1,lx
         newbox(l)=0
 110  continue
c
c        initialize number of total new boxes per j to zero.
c
      do 112 j=1,31
         ntotbt(j)=999999999
         nboxj(j)=0
 112  continue
c
      iorigb=(ibit+jbit+kbit)*lx
c        ibit = bits to pack the jmin( ).
c        jbit = bits to pack the lbit( ).
c        kbit = bits to pack the nov( ).
c        lx = number of groups.
         ntotbt(kbit)=iorigb
c           this is the value of total bits for the original lx
c           groups, which requires kbits to pack the group
c           lenghts.  setting this here makes one less loops
c           necessary below.
c
c        compute bits now used for the parameters defined.
c
c        determine other possibilites by increasing lx and decreasing
c        nov( ) with values greater than thresholds.  assume a group is
c        split into 2 or more groups so that kbit is reduced without
c        changing ibit or jbit.
c
      jj=0
c
      do 200 j=min(30,kbit-1),2,-1
c           values ge kbit will not require splits.  once the total
c           bits start increasing with decreasing j, stop.  also, the
c           number of bits required is known for kbits = ntotbt(kbit).
c
         newboxt=0
c
         do 190 l=1,lx
c
            if(nov(l).lt.ibxx2(j))then
               newbox(l)=0
c                 no splits or new boxes.
               go to 190
            else
               novl=nov(l)
c
               m=(nov(l)-1)/(ibxx2(j)-1)+1
c                 m is found by solving the equation below for m:
c                 (nov(l)+m-1)/m lt ibxx2(j)
c                 m gt (nov(l)-1)/(ibxx2(j)-1)
c                 set m = (nov(l)-1)/(ibxx2(j)-1)+1
 130           novl=(nov(l)+m-1)/m
c                 the +m-1 is necessary.  for instance, 15 will fit
c                 into a box 4 bits wide, but won't divide into
c                 two boxes 3 bits wide each.
c      
               if(novl.lt.ibxx2(j))then
                  go to 185
               else
                  m=m+1
c***                  write(kfildo,135)l,nov(l),novl,m,j,ibxx2(j)
c*** 135              format(/' at 135--l,nov(l),novl,m,j,ibxx2(j)',6i10)               
                  go to 130
               endif
c
c                 the above do loop will never complete.
            endif
c
 185        newbox(l)=m-1
            newboxt=newboxt+m-1
 190     continue
c
         nboxj(j)=newboxt
         ntotpr=ntotbt(j+1)
         ntotbt(j)=(ibit+jbit)*(lx+newboxt)+j*(lx+newboxt)
c
         if(ntotbt(j).ge.ntotpr)then
            jj=j+1
c              the plus is used because j decreases per iteration.
            go to 250
         else
c
c              save the total new boxes and newbox( ) in case this
c              is the j to use.
c
            newboxtp=newboxt
c
            do 195 l=1,lx
               newboxp(l)=newbox(l)
 195        continue
c
c           write(kfildo,197)newboxt,ibxx2(j)
c197        format(/' *****************************************'
c    1             /' the number of newboxes per group of the total',
c    2              i10,' for group maxsize plus 1 ='i10
c    3             /' *****************************************')
c           write(kfildo,198) (newbox(l),l=1,lx)
c198        format(/' '20i6/(' '20i6))
    
         endif
c        
c205     write(kfildo,209)kbit,iorigb
c209     format(/' original bits with kbit of',i5,' =',i10)
c        write(kfildo,210)(n,n=2,10),(ibxx2(n),n=2,10),
c    1                    (ntotbt(n),n=2,10),(nboxj(n),n=2,10),
c    2                    (n,n=11,20),(ibxx2(n),n=11,20),
c    3                    (ntotbt(n),n=11,20),(nboxj(n),n=11,20),
c    4                    (n,n=21,30),(ibxx2(n),n=11,20),
c    5                    (ntotbt(n),n=21,30),(nboxj(n),n=21,30)
c210     format(/' the total bytes for maximum group lengths by row'//
c    1      '   j         = the number of bits per group length'/
c    2      '   ibxx2(j)  = the maximum group length plus 1 for this j'/
c    3      '   ntotbt(j) = the total bits for this j'/
c    4      '   nboxj(j)  = the new groups for this j'/
c    5      4(/10x,9i10)/4(/10i10)/4(/10i10))
c
 200  continue
c
 250  pimp=((iorigb-ntotbt(jj))/float(iorigb))*100.
c     write(kfildo,252)pimp,kbit,jj
c252  format(/' percent improvement =',f6.1,
c    1        ' by decreasing group lengths from',i4,' to',i4,' bits')
      if(pimp.ge.2.)then
c
c        write(kfildo,255)cfeed,newboxtp,ibxx2(jj)
c255     format(a1,/' *****************************************'
c    1             /' the number of newboxes per group of the total',
c    2             i10,' for group maxsize plus 1 ='i10
c    2             /' *****************************************')
c        write(kfildo,256) (newboxp(l),l=1,lx)
c256     format(/' '20i6)
c
c           adjust group lengths for maximum length of jj bits.
c           the min per group and the number of bits required
c           per group are not changed.  this may mean that a
c           group has a min (or reference) that is not zero.
c           this should not matter to the unpacker.
c
         lxnkp=lx+newboxtp
c           lxnkp = the new number of boxes
c  
         if(lxnkp.gt.ndg)then
c              dimensions not large enough.  probably an error
c              of some sort.  abort.
c           write(kfildo,257)ndg,lxnpk
c        1         2         3         4         5         6         7 x
c257        format(/' dimensions of jmin, etc. in reduce =',i8,
c    1              ' not large enough for the expanded number of',
c    2              ' groups =',i8,'.  abort reduce.')
            ier=715
            go to 410
c              an abort causes the calling program to reexecute 
c              without calling reduce.
         endif
c
         lxn=lxnkp
c           lxn is the number of the box in the new series being
c           filled.  it decreases per iteration.
         ibxx2m1=ibxx2(jj)-1
c           ibxx2m1 is the maximum number of values per group.
c
         do 300 l=lx,1,-1
c
c              the values is nov( ) represent those values + novref.
c              when values are moved to another box, each value
c              moved to a new box represents that value + novref.
c              this has to be considered in moving values.
c
            if(newboxp(l)*(ibxx2m1+novref)+novref.gt.nov(l)+novref)then
c                 if the above test is met, then moving ibxx2m1 values
c                 for all new boxes will leave a negative number for
c                 the last box.  not a tolerable situation.
               movmin=(nov(l)-(newboxp(l))*novref)/newboxp(l)
               left=nov(l)
c                 left = the number of values to move from the original
c                 box to each new box except the last.  left is the
c                 number left to move.
            else
               movmin=ibxx2m1
c                 movmin values can be moved for each new box.
               left=nov(l)
c                 left is the number of values left to move.
            endif
c
            if(newboxp(l).gt.0)then
               if((movmin+novref)*newboxp(l)+novref.le.nov(l)+novref.
     1          and.(movmin+novref)*(newboxp(l)+1).ge.nov(l)+novref)then
                  go to 288
               else
c***d                 write(kfildo,287)l,movmin,novref,newboxp(l),nov(l)
c***d287              format(/' at 287 in reduce--l,movmin,novref,',
c***d    1                    'newboxp(l),nov(l)',5i12
c***d    2                    ' reduce aborted.')
c              write(kfildo,2870)
c2870          format(/' an error in reduce algorithm.  abort reduce.')
               ier=714
               go to 410
c                 an abort causes the calling program to reexecute 
c                 without calling reduce.
               endif
c
            endif
c
 288        do 290 j=1,newboxp(l)+1
               move=min(movmin,left)
               jmin(lxn)=jmin(l)
               jmax(lxn)=jmax(l)
               lbit(lxn)=lbit(l)
               nov(lxn)=move
               lxn=lxn-1
               left=left-(move+novref)
c                 the move of move values really represents a move of
c                 move + novref values.
 290        continue
c
            if(left.ne.-novref)then
c***               write(kfildo,292)l,lxn,move,lxnkp,ibxx2(jj),left,nov(l),
c***     1                          movmin
c*** 292           format(' at 292 in reduce--l,lxn,move,lxnkp,',
c***     1                'ibxx2(jj),left,nov(l),movmin'/8i12)
            endif
c     
 300     continue
c
         lx=lxnkp
c           lx is now the new number of groups.
         kbit=jj
c           kbit is now the new number of bits required for packing
c           group lenghts.
      endif
c
c     write(kfildo,406)cfeed,lx
c406  format(a1,/' *****************************************'
c    1          /' the group sizes nov( ) after reduction in size',
c    2           ' for'i10,' groups',
c    3          /' *****************************************')
c     write(kfildo,407) (nov(j),j=1,lx)
c407  format(/' '20i6)
c     write(kfildo,408)cfeed,lx
c408  format(a1,/' *****************************************'
c    1          /' the group minima jmin( ) after reduction in size',
c    2           ' for'i10,' groups',
c    3          /' *****************************************')
c     write(kfildo,409) (jmin(j),j=1,lx)
c409  format(/' '20i6)
c
 410  return
      end
      
