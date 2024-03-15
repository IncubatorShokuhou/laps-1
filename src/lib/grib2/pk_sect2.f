      subroutine pk_sect2(kfildo,ipack,nd5,rdat,nrdat,idat,nidat,
     1                    l3264b,locn,ipos,exists,ier,isevere,*) 
c
c        february 2001   lawrence    gsc/mdl      original coding
c        november 2001   glahn   rearranged dimensions
c        january  2002   glahn   changed int( ) to nint( ) 
c
c        purpose
c            packs section 2, the local use section, of a grib2
c            message.  section 2 is optional. 
c
c            this routine allows the user to pack integer
c            and/or floating point groups of local use data
c            into section 2 of the grib2 message.  each group of
c            data is packed into section 2 using the simple packing
c            method.  the user must specify the decimal scale factor
c            to use in packing each group of data.  for simplicity,
c            the binary scale factor is not used when packing data
c            into the local use section.
c
c            the floating point local use section data to be packed
c            into the grib2 message is passed into this routine
c            through the rdat( ) array calling argument.  likewise,
c            the integer local use section data is passed into this
c            routine through the idat( ) array calling argument. 
c            each group of local use data stored into the rdat( )
c            and idat( ) arrays must be preceded by the number of 
c            values it contains and the decimal scale factor to use
c            in packing the group's data.  the end of the local use
c            data in the rdat( ) and idat( ) arrays is signaled by
c            placing a value of "0" in the array element immediately
c            following the last data value of the last group of local
c            use data.  the data in the rdat( ) and idat( ) arrays
c            must be arranged by the caller of this routine as
c            follows:
c
c            for 1 to k groups of data: 
c
c            rdat(1)        = number of values in the first group
c                             of local use data (n1)
c            rdat(2)        = the decimal scale factor to use in
c                             packing the first group of local use
c                             data (must be a whole number)
c            rdat(3)
c            -rdat(n1+2)    = first group of local use data values
c            rdat(n1+3))    = number of values in the second group of
c                             local use data (n2)
c            rdat(n1+4)     = the decimal scale factor to use in packing
c                             the second group of local use data (must 
c                             be a whole number)
c            rdat(n1+5)
c            -rdat(n1+n2+4) = second group of local use data
c                             values 
c
c                               ........
c
c            rdat((k-1)*2+1+n1+n2+...+n(k-1))  = number of values in
c                                                the kth group of data
c                                                (nk)
c            rdat((k-1)*2+2+n1+n2+...+n(k-1))  = the decimal scale
c                                                factor to use in
c                                                packing the kth group
c                                                of data 
c            rdat((k-1)*2+3+n1+n2+...+n(k-1)) - 
c            rdat((k-1)*2+n1+n2+...+n(k-1)+nk) = the number of
c                                                values in the kth
c                                                group of data   
c            rdat((k-1)*2+1+n1+n2+...+nk)      = 0  no more data
c
c            if the caller is not supplying any local use 
c            data to pack into the grib2 message, then he must 
c            make sure that rdat(1) and idat(1) are both equal 
c            to "0".  in that case, a section 2 will not be packed
c            into the grib2 message.
c
c            the local use data is packed into the message using the
c            following format:
c
c            section 2 octet(s)   description
c            1-4                  total length of section 2
c            5                    section number (2)
c            6                    section 2 format version number
c                                 the current version is 1.
c            7-8                  total number of data groups in
c                                 the local use section (n)
c            9-12                 number of values in first
c                                 group of local use data (n1) 
c            13-16                reference value of first group of
c                                 data 
c            17-18                decimal scale factor
c            19                   number of bits to pack each value
c                                 of the first group of data with
c            20                   type of data in first group
c                                 ("0" = floating point, "1" =
c                                 integer)
c            21-nn                the first group of data packed
c                                 using the simple packing method
c            (nn+1)-(nn+4)        the number of values in the second
c                                 group of data
c            (nn+5)-(nn+8)        the reference value of the second
c                                 group of data
c            (nn+9)-(nn+10)       the decimal scale factor of the
c                                 second group of data
c            (nn+11)              the number of bits to pack each value
c                                 of the second group of data with
c            (nn+12)              type of data in the second group
c                                 ("0" = floating point, "1" = 
c                                 integer)
c            (nn+13) - mm         the second group of data packed
c                                 using the simple packing method
c       
c            this pattern repeats itself for each of the n groups 
c            of local use data specified in octets 7-8 of this 
c            section.
c
c        data set use
c           kfildo - unit number for output (print) file.  (output)
c
c        variables
c              kfildo = unit number for output (print) file.  (input)
c            ipack(j) = the array that holds the actual packed 
c                       message (j=1,nd5).  (input/output)
c                 nd5 = the size of the array ipack( ).  (input)
c             rdat(j) = the array containing the local use groups
c                       consisting of floating point data (j=1,nrdat).
c                       (input)
c               nrdat = the dimension of the rdat( ) array.  (input)
c             idat(j) = the array containing the local use groups
c                       consisting of integer data (j=1,nidat).
c                       (input)
c               nidat = the dimension of the idat( ) array.  (input)
c              l3264b = the integer word length in bits of the
c                       machine being used.  values of 32 and 64 are
c                       accomodated.  (input)
c                locn = the word position to place the next value.
c                       (input/output)
c                ipos = the bit position in locn to start placing
c                       the next value.  (input/output)
c              exists = indicates to the calling routine whether or not
c                       section 2 exists (logical).  (output)
c                 ier = return status code. (output)
c                        0 = good return.
c                      1-4 = error codes generated by pkbg. see the
c                            documentation in the pkbg routine.
c                      5,6 = error codes generated by length function.
c                            see the documentation for the length
c                            function.
c                      202 = the idat( ) or rdat( ) array was not
c                            dimensioned large enough to contain
c                            the local use data. 
c             isevere = the severity level of the error.  the only
c                       value retuned is:
c                       2 = a fatal error  (output)
c                   * = alternate error return.
c
c             local variables
c                ibit = the number of bits required to pack each 
c                       value in one local data group using the
c                       simple packing method.
c                  id = the decimal scale factor for one local
c                       data group.  this is specified by the user
c                       for each group in the idat( ) and rdat( )
c                       arrays. 
c              igroup = keeps a count of the number of groups of 
c                       local use data packed into section 2 of the
c                       grib2 message.
c               index = used to keep track of which data value is
c                       currently being processed in the idat( )
c                       or rdat( ) array. 
c              intdat = flag indicating if there are any local integer 
c                       data groups to be packed.  (logical) 
c             ipos2_1 = saves the bit position in locn2_1 to store the  
c                       length of section 2.
c             ipos2_7 = saves the bit position in locn2_7 to store the
c                       total number of local use data groups packed
c                       into section 2.
c            ipos2_19 = saves the bit position in locn2_19 to store the
c                       number of bits required to pack each value of
c                       a local use data group.
c            itemp(j) = array to contain a local use data group
c                       consisting of integer values.  this array is
c                       used to pass the data to the simple packing
c                       routine (j=1,nidat).   
c              ivalue = this is equivalenced to rvalue.  it is used to
c                       pack the bit pattern of a floating point 
c                       value into the grib2 message.
c            iversion = the version number of the section 2 format.
c                       this routine is designed to pack version 1.
c               izero = contains a value of "0" to be packed into the
c                       grib2 message.
c                  ix = this is a loop indexing variable.
c             locn2_1 = saves the word position in ipack to store
c                       the length of section 2. 
c             locn2_7 = saves the word position in ipack to store
c                       the total number of local use data groups
c                       packed into section 2. 
c            locn2_19 = saves the word position in ipack to store the
c                       number of bits required to pack each value
c                       of a local use data group.
c                mina = the minimum value of a group of integer local
c                       use data. 
c                   n = a short-hand representation of l3264b.
c             realdat = a flag indicating if there are any local
c                       use data groups of floating point values
c                       to pack.  (logical) 
c               rmina = the minimum value of a group of floating 
c                       point local use data.
c            rtemp(j) = an array to temporarily contain one local use
c                       group of floating point data (j=1,nrdat). 
c              rvalue = this is equivalenced to ivalue.  it is used
c                       to pack the bit pattern of a floating point
c                       value into the grib2 message. 
c
c        non system subroutines called
c           length,pkbg,pk_smple,prep_sect2_int,prep_sect2_real
c           
c
      logical exists,intdat,realdat
c
      dimension idat(nidat),rdat(nrdat)
      dimension ipack(nd5)
      dimension itemp(nidat),rtemp(nrdat)
c        itemp( ) and rtemp( ) are automatic arrays.
c
      data iversion/1/
      data izero/0/
c
      equivalence(rvalue,ivalue)
c
      exists=.false.
      intdat=.false.
      realdat=.false.
c
      n=l3264b
      ier=0
      igroup=0
c
c        all errors generated by this routine are fatal.
      isevere=2
c
c        check to determine if data exists for this section.
      if(idat(1).ne.0)intdat=.true.
      if(rdat(1).ne.0)realdat=.true.
c
      if(intdat.or.realdat)then
c
c           there is local use data to be packed into the local
c           use section of the grib2 message.
         exists=.true.
c
c           bytes 1-4 of section 2 must be filled in later with
c           the record length in bytes.  loc2_1 and ipos2_1 hold the
c           location. 
         locn2_1=locn
         ipos2_1=ipos
         call pkbg(kfildo,ipack,nd5,locn,ipos,izero,32,n,ier,*900)
c
c           pack the number of the section.
         call pkbg(kfildo,ipack,nd5,locn,ipos,2,8,n,ier,*900)
c
c           pack the version number of the format of the local
c           use data.
         call pkbg(kfildo,ipack,nd5,locn,ipos,iversion,8,n,ier,*900)
c
c           save the position of octets 7-8 of section 2.
c           these octets will be filled in later with
c           the total number of local use data groups packed
c           into this section. 
         locn2_7=locn
         ipos2_7=ipos
         call pkbg(kfildo,ipack,nd5,locn,ipos,izero,16,n,ier,*900)
c
         if(intdat)then
c
c             pack the integer local use data into 
c             section 2 first.
           isize=idat(1)
           index=1
c
           do while(isize.gt.0)
c
c                has idat( ) been dimensioned large enough?
              if(nidat.lt.(isize+index+2))then
                 ier=202
                 goto 900
              endif
c
              index=index+1
c
c                retrieve the decimal scale factor.
              id=idat(index)
              index=index+1
c
c                copy the local use data into the itemp( ) array.
              do 10 ix=1,isize
                 itemp(ix)=idat(index)
                 index=index+1
 10           continue
c
c                prepare the data in the itemp( ) array
c                to be packed.
              call prep_sect2_int(itemp,isize,0,id,mina)
c
c                pack the number of values in the group.
              call pkbg(kfildo,ipack,nd5,locn,ipos,isize,32,n,ier,*900)
c
c                pack the reference value of the group.
              rmina=float(mina)
              rvalue=fmkieee(rmina)
              call pkbg(kfildo,ipack,nd5,locn,ipos,ivalue,32,n,
     1                  ier,*900)
c
c                pack the decimal scale factor.
              call pkbg(kfildo,ipack,nd5,locn,ipos,id,16,n,ier,*900)
c
c                save the location of the octet to contain
c                the number of bits required to pack each 
c                value of the group of local use data.
              locn2_19 = locn
              ipos2_19 = ipos
              call pkbg(kfildo,ipack,nd5,locn,ipos,izero,8,n,ier,*900)
c
c                pack the type of the local use data. 
              call pkbg(kfildo,ipack,nd5,locn,ipos,1,8,n,ier,*900)
c
c                pack the data using the simple packing method.
              call pk_smple(kfildo,itemp,isize,ipack,nd5,locn,ipos,
     1                      ibit,l3264b,ier,*900)
c
c                pack the number of bits necessary to contain the
c                largest value in the data field.
              call pkbg(kfildo,ipack,nd5,locn2_19,ipos2_19,ibit,8,n,
     1                  ier,*900)
c
c                increment the number of groups and retrieve the size
c                of the next group of local use data to be packed.
              igroup=igroup+1
              isize=idat(index)
           enddo
         endif
c
         if(realdat)then
c
c             pack the floating point local use data.
           isize=nint(rdat(1))
           index=1
c
           do while(isize.gt.0)
c
c                has rdat( ) been dimensioned large enough?
              if(nrdat.lt.(isize+index+2))then
                 ier=202
                 goto 900
              endif
c
              index=index+1
c
c                retrieve the decimal scaling factor.
              id=nint(rdat(index))
              index=index+1
c
c                copy the data into the rtemp( ) array.
              do 20 ix=1,isize
                 rtemp(ix)=rdat(index)
                 index=index+1
 20           continue
c
c                prepare the data in the rtemp( ) array.  the
c                scaled data is returned in the integer array
c                itemp( ).  
              call prep_sect2_real(rtemp,itemp,isize,0,id,rmina)
c
c                pack the number of values in the group.
              call pkbg(kfildo,ipack,nd5,locn,ipos,isize,32,n,ier,*900)
c
c                pack the reference value of the group.
              rvalue=fmkieee(rmina)
              call pkbg(kfildo,ipack,nd5,locn,ipos,ivalue,32,n,
     1                  ier,*900)
c
c                pack the decimal scale factor.
              call pkbg(kfildo,ipack,nd5,locn,ipos,id,16,n,ier,*900)
c 
c                save the location of the octet to contain
c                the number of bits required to pack each 
c                value of the group of local use data.
              locn2_19 = locn
              ipos2_19 = ipos
              call pkbg(kfildo,ipack,nd5,locn,ipos,izero,8,n,ier,*900)
c
c                pack the type of the local use data. 
              call pkbg(kfildo,ipack,nd5,locn,ipos,0,8,n,ier,*900)
c
c                pack the data using the simple packing method.
              call pk_smple(kfildo,itemp,isize,ipack,nd5,locn,ipos,
     1                      ibit,l3264b,ier,*900)
c
c                pack the number of bits necessary to contain the
c                largest value in the data field.
              call pkbg(kfildo,ipack,nd5,locn2_19,ipos2_19,ibit,8,n,
     1                  ier,*900)
c
c                increment the group count and retrieve the size of
c                the next group of local use data to be packed.
              igroup=igroup+1
              isize=nint(rdat(index))
c
           enddo
         endif
c
c           pack the total number of groups into octets 6-7 of 
c           section 2. 
         call pkbg(kfildo,ipack,nd5,locn2_7,ipos2_7,igroup,16,n,
     1             ier,*900)
c
c           compute the length of the section and pack it.  locn2_1 and
c           ipos2_1 represent the length of the record before
c           section 2.  8 is the number of bits in a byte, and each
c           section ends at the end of a byte.
         isize=length(kfildo,ipack,nd5,l3264b,locn2_1,ipos2_1,locn,
     1                 ipos,ier)
      endif          
c
c        error return section
 900  if(ier.ne.0)return 1
c
      return
      end
