      subroutine prepr(kfildo,a,ia,ib,nx,ny,nval,iclean,ibitmap,
     1                 is5,ns5,is6,ns6,is7,ns7,id,ie,mina,
     2                 xmina,missp,misss,xmissp,xmisss,minpk,
     3                 ipkopt,ier,jer,ndjer,kjer,*)
c
c        march    2000   glahn   called by pk_grib2
c        may      2000   lawrence modified routine to reflect
c                                the latest wmo grib2 changes.
c                                the most significant change
c                                is that the values are first
c                                scaled by the decimal scale
c                                factor, then the minimum is
c                                taken from the data field, and
c                                then the values are multiplied
c                                by the binary scale factor.
c        january  2001   glahn   write(kfildo) made /d; comments
c        january  2001   glahn   comment following call to pack_opt;
c                                eliminated is5( ) and ns5 in call to
c                                prep_int and prep_flt
c        january  2001   glahn/lawrence removed unused miss and xmiss;
c                                eliminated is6( ) and ns6 from call
c                                to pack_opt
c        november 2001   glahn   added jer, ndjer, and kjer to call
c                                to check_int, check_flt, int_map,
c                                and flt_map
c        december 2001   glahn   added kfildo to call to check_int
c                                and check_flt
c        december 2001   glahn   moved test on is5(10) = 0, 2, or 3 
c                                from pk_sect5.
c        january  2002   glahn   added ier and *900 to calls to 
c                                prep_int and prep_flt; changed
c                                nval comment from input to output
c        february 2002   glahn   comments
c
c        purpose
c            finds the reference value and subtracts it from the
c            data values.  either floating or integer data are
c            handled, and when there are or aren't missing values.
c            a bit map is generated if necessary, and values are
c            inserted into the grid for complex and spatial
c            differencing.  operations within loops are kept
c            to a relative minimum at the expense of more code
c            for efficiency.  second order differencing is not
c            done when there are secondary missing values present.
c
c            all computations are done on integer data until
c            they are scaled when incoming data are integer in ia( ),
c            then put back in ia( ).
c            all computations are done on floating point data
c            when incoming data are floating point in a( ), then.
c            put in ia( ).
c
c        data set use
c           kfildo - unit number for output (print) file. (output)
c
c        variables
c              kfildo = unit number for output (print) file.  (input)
c                a(k) = when is5(21) = 0, a( ) contains the data
c                       (k=1,nval).  (input)
c               ia(k) = when is5(21) = 1, ia( ) contains the data
c                       (k=1,nval).  the values to pack are in
c                       ia( ) on output.  (input/output)
c               ib(k) = the bit map when one is used.  (input/output)
c                       it can be input or in can be calculated if
c                       the simple method is used (k=1,nxy).
c                       complex and spatial differencing do not
c                       use a bit map, but will accept one and insert
c                       the missing values.
c               nx,ny = dimensions of the grid.  nx*ny is the dimension
c                       of a( ), ia( ), and ib( ).  (input)
c                nval = the number of values in a( ) or ia( ).  (output)
c              iclean = 1 when there are no missing values in a( ) or
c                         ia( ).
c                       0 otherwise.
c                       (input/output)
c             ibitmap = 1 when there is a bitmap in ib( ).
c                       0 otherwise.
c                       (input)
c              is5(j) = the values associated with section 5, keyed
c                       to the octet number (j=1,ns5).  the elements
c                       used in this routine are:
c                       is5(10), template number:
c                         0 = simple
c                         1 = not supported
c                         2 = complex
c                         3 = spatial differencing
c                       is5(21), type of original field:
c                         0 = floating point
c                         1 = integer
c                 ns5 = the dimension of is5( ).  (input)
c              is6(j) = the values associated with section 6, keyed
c                       to the octet number.  the elements used
c                       in this routine are:
c                       is6(6) = bit map indicator:
c                       0 = bit map included
c                       255 = no bit map
c                       is6(6) is modified only as needed.  a value
c                       between 1 and 254 is not disturbed
c                       (j=1,ns6). (input/output)
c                 ns6 = the dimension of is6( ).  (input)
c                  id = the decimal scaling factor.  (input)
c                  ie = the binary scaling factor.  (input)
c                mina = the field minimum value when the original data
c                       are integer.  (output)
c               xmina = the field minimum value when the original data
c                       are floating point.  (output)
c               missp = when missing points can be present in the data,
c                       they will have the value missp or misss when
c                       the data are integer.  missp is the primary
c                       missing value when the original data are
c                       integer.  (input)
c               misss = secondary missing value indicator when the data
c                       are integer.  (input)
c              xmissp = when missing points can be present in the data,
c                       they will have the value xmissp or xmisss when
c                       the data are floating point.  xmissp
c                       is the primary missing value.  (input)
c              xmisss = secondary missing value indicator when the data
c                       are floating point.  (input)
c               minpk = increment in which ranges will be computed.
c                       (input)
c              ipkopt = packing indicator:
c                       0 = error, don't pack
c                       1 = pack ia( ), simple
c                       2 = pack ia( ) and ib( ), simple
c                       3 = pack complex or spatial differencing
c                       4 = pack complex.
c                       (output)
c                 ier = return status.  (output)
c                       0 = good error return.
c                     902 = there are no "good" values in the grid
c                           and the bit-map indicates that.
c                     903 = there are no values in the grid and the
c                           bit-map indicates that there should be.
c                     904 = there are no values in the grid, and there
c                           is no bit-map.
c                     905 = invalid data type indicated in
c                           is5(21).
c                     906 = no missing values in the array, bit-map
c                           provided, simple, and nval=nxy. (warning)
c                     907 = no missing values in the array, no
c                           bit-map provided, and nval ne nxy while
c                           packing simple.
c                     508 = unsupported packing type in is5(10).
c               ndjer = dimension of jer( ).  (input)
c                kjer = number of values in jer( ).  (input/output)
c                   * = alternate return when jer ge 900.
c
c        local variables
c               cfeed = contains the character representation
c                       of a printer form feed.
c               ifeed = contains the integer value of a printer
c                       form feed.
c                 nxy = nx*ny.
c        1         2         3         4         5         6         7 x
c
c        non system subroutines called
c           boustro,check_flt,check_int,pk_trace,flt_map,
c           int_map,pack_opt,prep_flt,prep_int,prep_noval
c
      character*1 cfeed
      logical jmisss,jmissp
c
      dimension a(nx*ny),ia(nx*ny),ib(nx*ny)
      dimension is5(ns5),is6(ns6),is7(ns7)
      dimension jer(ndjer,2)
c
      data ifeed/12/
c
      ier=0
      nxy=nx*ny
      ipkopt=0
      cfeed=char(ifeed)
c
      if(is5(10).ne.0.and.is5(10).ne.2.and.is5(10).ne.3)then
c           checks packing type; must be 1, 2, or 3.
         ier=508
         go to 900
      endif
c
c        process the bit-map, if one was supplied.
      if(is5(21).eq.1)then
c           this is integer data.
         call int_map(ia,ib,nxy,is5,ns5,iclean,
     1                ibitmap,missp,jer,ndjer,kjer)
      elseif(is5(21).eq.0)then
c           this is floating point data.
         call flt_map(a,ib,nxy,is5,ns5,iclean,
     1                ibitmap,xmissp,jer,ndjer,kjer)
      else
c
c           the data are neither integer nor floating point.
c           the type of data is not legitimate.
c           errors may result if packing continues.
         call pk_trace(kfildo,jer,ndjer,kjer,905,2)
         ipkopt=0
         go to 900
      endif
c
c        process the missing values, determining if there
c        are no missing values, primary missing values
c        or secondary missing values in the field.
c
      if(is5(21).eq.1)then
         call check_int(kfildo,ia,ib,nval,nxy,is5,ns5,iclean,ibitmap,
     1                  missp,misss,jmissp,jmisss,ier,
     2                  jer,ndjer,kjer,*900)
c           this is integer data.
      elseif(is5(21).eq.0)then
         call check_flt(kfildo,a,ib,nval,nxy,is5,ns5,iclean,ibitmap,
     1                  xmissp,xmisss,jmissp,jmisss,ier,
     2                  jer,ndjer,kjer,*900)
c           this is floating point data.
      else
         ier=909
         go to 900
      endif
c
c        are there any values in this data field?
      if(nval.eq.0)then
         call prep_noval(ib,nxy,ibitmap,ipkopt,ier)
c
c           was a fatal error encountered?
         if(ier.eq.902)then
            call pk_trace(kfildo,jer,ndjer,kjer,ier,1)
c              ier = 902 is not fatal.
         elseif(ier.ne.0)then
            go to 900
         endif
      else
c
c           determine the type of the data.
         if(is5(21).eq.1)then
c
c              the data are integer.
c              if the complex or complex with spatial differences
c              packing method is being used, then scan the
c              data boustrophedonically.
c
cwdt        if((is5(10).eq.2).or.(is5(10).eq.3))then
cwdt           call boustro_int(ia,nx,ny)
cwdt        endif
c
            call prep_int(ia,nxy,nval,iclean,id,ie,
     1                    mina,jmissp,jmisss,missp,misss,ier,*900)
         else
c
c              the data are floating point.
c              if the complex or complex with spatial differences
c              packing method is being used, then scan the
c              data boustrophedonically.
c
cwdt        if((is5(10).eq.2).or.(is5(10).eq.3))then
cwdt           call boustro_flt(a,nx,ny)
cwdt        endif
c
            call prep_flt(a,ia,nxy,nval,iclean,id,ie,
     1                    xmina,jmissp,jmisss,xmissp,xmisss,ier,*900)
c           write(kfildo,10)cfeed
c10         format(a1,/' **********************'
c    1                /' original scaled values'
c    2                /' **********************')
c           write(kfildo,20) (ia(j),j=1,200)
c20         format(/' '20i6)
c
         endif
c
      endif
c
c        set ipkopt to a proper value if it has not already been
c        done so in prep_noval.
      if(ipkopt.eq.0)then
         call pack_opt(kfildo,ia,ib,nxy,nval,iclean,ibitmap,
     1                 is5,ns5,is7,ns7,jmisss,missp,
     2                 minpk,numoctet,ipkopt,jer,ndjer,
     3                 kjer,*900)
c           pack_opt may have a normal return with a non fatal
c           ier.  this will have been inserted into jer( , ).
c           prepr will return to calling program with that ier.
c           this is a little dangerous, but the calling program
c           pk_grib2 just calls another routine that sets ier = 0.
c        write(kfildo,30)cfeed,is7(8)
c30      format(a1,/' ***************************'
c    1             /' 2nd order differences after'
c    2             /' the removal of the field'
c    3             /' minimum ',i6,
c    4             /' ***************************')
c        write(kfildo,40) (ia(j),j=1,200)
c40      format(/' '20i6)
      endif
c
c        initialize the pertinent values in is5( ) depending
c        on the packing option.
 700  if(ipkopt.eq.1)then
         is5(1)=21
         is5(6)=nval
         is5(10)=0
         is6(6)=255
      elseif(ipkopt.eq.2)then
         is5(1)=21
         is5(6)=nval
         is5(10)=0
         is6(6)=0
      elseif(ipkopt.eq.3)then
         is5(1)=49
         is5(6)=nxy
         is5(10)=3
         is5(48)=2
         is5(49)=numoctet
         is6(6)=255
      elseif(ipkopt.eq.4)then
         is5(1)=47
         is5(6)=nxy
         is5(10)=2
         is6(6)=255
      endif
      return
c
 900  return 1 
      end
