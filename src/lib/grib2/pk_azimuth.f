      subroutine pk_azimuth(kfildo,ipack,nd5,is3,ns3,ipkopt,l3264b,
     1                      locn,ipos,ier,*)
c
c        march    2000   lawrence  gsc/tdl   original coding
c        january  2001   glahn     comments; ier = 303 changed to 304;
c                                  added check for size of is3( );
c        february 2001   glahn/lawrence changed j=39+4*(k-1)
c                                  to j=40+4*(k-1)
c
c        purpose
c            packs template 3.120, an azimuth-range projection
c            template, in section 3 of a grib2 message. it is
c            the responsibility of the calling routine to pack 
c            the first 13 octets in section 3.
c
c        data set use
c           kfildo - unit number for output (print) file. (output)
c
c        variables
c              kfildo = unit number for output (print) file. (input)
c            ipack(j) = the array that holds the actual packed message
c                       (j=1,nd5). (input/output)
c                 nd5 = the size of the array ipack( ). (input)
c              is3(j) = contains the azimuth-range projection
c                       information (in elements 15 through 39+4(nr-1)
c                       where nr is the number of radials in the
c                       projection) that will be packed into ipack( )
c                       (j=1,ns3). (input)
c                 ns3 = size of is3( ). (input) 
c              ipkopt = packing indicator:
c                       0 = error, don't pack
c                       1 = pack ia( ), simple
c                       2 = pack ia( ) and ib( ), simple
c                       3 = pack complex or spatial differencing
c                       4 = pack complex.
c                       (input)
c              l3264b = the integer word length in bits of the machine
c                       being used. values of 32 and 64 are
c                       accommodated. (input)
c                locn = the word position to place the next value.
c                       (input/output)
c                ipos = the bit position in locn to start placing
c                       the next value. (input/output)
c                 ier = return status code. (output)
c                         0 = good return.
c                       1-4 = error codes generated by pkbg. see the 
c                             documentation in the pkbg routine.
c                       302 = is3( ) has not been dimensioned large
c                             enough to contain the entire template. 
c                       304 = is3(13) does not indicate polar
c                             stereographic map projection; 
c                             incorrectly called subroutine. 
c                   * = alternate return when ier .ne. 0.
c
c             local variables
c               isign = a flag indicating whether a value to be
c                       packed is positive or negative. 
c                   j = array indexing variable.
c                   k = looping variable.
c                   n = l3264b = the integer word length in bits of 
c                       the machine being used. values of 32 and
c                       64 are accommodated. 
c             minsize = the smallest allowable dimension for is3( ).
c
c        non system subroutines called
c           pkbg
c
      parameter(minsize=39)
c
      dimension ipack(nd5),is3(ns3)
c
      n=l3264b
      ier=0
c
c        check to make sure that the user wants to process
c        an azimuth-range projection.  check for ns3 ge 13
c        assures a value for is3(13). 
c
      if(ns3.ge.13)then
c
         if(is3(13).ne.120)then
c           write(kfildo,10)is3(13)
c10         format(/' map projection code ',i4,' indicated by is3(13)'/
c    1              ' is not azimuth-range. please refer to'/
c    2              ' the grib2 documentation to determine the'/
c    3              ' correct map projection packer to call.'/) 
            ier=304
            go to 900
         endif
c
      endif
c
c        check the dimensions of is3( ).  check for ns3 ge 19
c        assures a value for is3(19).
c
      if(ns3.ge.19)then
c
         if(ns3.lt.(minsize+4*(is3(19)-1)))then
            ier=302
            go to 900
         endif
c
      else
         ier=302
         go to 900
      endif
c
c        pack the number of data bins along the radials (nb).
      call pkbg(kfildo,ipack,nd5,locn,ipos,is3(15),32,n,ier,*900)
c
c        pack the number of radials in the product (nr).
      call pkbg(kfildo,ipack,nd5,locn,ipos,is3(19),32,n,ier,*900)
c
c        pack the latitude & longitude of the center point
      isign=0
      if(is3(23).lt.0)isign=1
      call pkbg(kfildo,ipack,nd5,locn,ipos,isign,1,n,ier,*900)
      call pkbg(kfildo,ipack,nd5,locn,ipos,abs(is3(23)),
     1          31,n,ier,*900)
      isign=0
      if(is3(27).lt.0) isign=1
      call pkbg(kfildo,ipack,nd5,locn,ipos,isign,1,n,ier,*900)
      call pkbg(kfildo,ipack,nd5,locn,ipos,abs(is3(27)),
     1          31,n,ier,*900)
c
c        pack the spacing of bins along the radials
      call pkbg(kfildo,ipack,nd5,locn,ipos,is3(31),32,n,ier,*900)
c
c        pack the offset from the origin to the inner bound.
      call pkbg(kfildo,ipack,nd5,locn,ipos,is3(35),32,n,ier,*900)
c
c        pack the scanning mode. check to see if the data points
c        have been ordered boustrophedonically. if so, set bit 4
c        to 1 in the scanning mode octet.
c
cwdt  if((ipkopt.eq.3).or.(ipkopt.eq.4))then
cwdt     is3(39)=ior(is3(39),6)
cwdt  else
         is3(39)=iand(is3(39),239)
cwdt  endif
c
      call pkbg(kfildo,ipack,nd5,locn,ipos,is3(39),8,n,ier,*900)
c
c        pack the starting azimuth and azimuthal width
c        for each of the radials (the is3(19) value)
c
      do k=1,is3(19)
         j = 40 + 4*(k-1)
         call pkbg(kfildo,ipack,nd5,locn,ipos,is3(j),16,n,ier,*900)
         isign=0
         if(is3(j+2).lt.0)isign=1
         call pkbg(kfildo,ipack,nd5,locn,ipos,isign,1,n,ier,*900)
         call pkbg(kfildo,ipack,nd5,locn,ipos,abs(is3(j+2)),
     1             15,n,ier,*900)
      enddo 
c
c        error return section
c
 900  if(ier.ne.0)return 1
c
      return
      end
