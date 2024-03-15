      subroutine pk_sect0(kfildo,ipack,nd5,is0,ns0,l3264b,
     1                    new,locn,ipos,iedition,locn0_9,
     2                    ipos0_9,ier,isevere,*)
c
c        march   2000   lawrence   gsc/tdl   original coding
c        january 2001   glahn      comments; added test for size
c                                  of is0( )
c        january 2001   glahn/lawrence initialized is0(1) and is0(8)
c        march 2001     lawrence   changed how this routine packs 
c                                  the numeric equivalent of the
c                                  string "grib".
c
c        purpose
c            packs section 0, the indicator section, of a grib2
c            message. section 0 contains information pertaining
c            to the size of the message and the edition of
c            the grib standards to be followed in the packing
c            process.
c
c            if this is a new grib2 message, then the contents of
c            is0( ) are packed into the message. if this is a
c            product that is getting appended onto a preexisting
c            message, then section 0 is not packed again. instead,
c            the end of message indicator of the previously packed
c            product is overwritten, and the size of the grib2
c            message (as stated by octets 9-16 of section 0 of
c            the previous packed product) is zeroed out. the total
c            size of the grib2 message is recalculated in another
c            routine every time new product is packed onto the message.
c
c        data set use
c           kfildo - unit number for output (print) file. (output)
c
c        variables
c              kfildo = unit number for output (print) file. (input)
c            ipack(j) = the array that holds the actual packed message
c                       (j=1,nd5). (input/output)
c                 nd5 = the size of the array ipack( ). (input)
c              is0(j) = contains the indicator data that
c                       will be packed into ipack( ) (j=1,ns0).
c                       (input)
c                 ns0 = size of is0( ). (input)
c              l3264b = the integer word length in bits of the machine
c                       being used. values of 32 and 64 are
c                       accommodated. (input)
c                 new = when new = 1, this is a new product. when
c                       new = 0, this is another grid to put into
c                       the same product.  (input)
c                locn = the word position to place the next value.
c                       (input/output)
c                ipos = the bit position in locn to start placing
c                       the next value. (input/output)
c            iedition = the edition number of the grib2 encoder used.
c                       (input) 
c             locn0_9 = contains the word position in ipack to
c                       pack the total length of the packed grib2
c                       message once this value can be determined,
c                       i.e. when the message has been completely
c                       packed.  (output)
c             ipos0_9 = contains the bit position in word locn0_9 of
c                       ipack to pack the total length of the packed
c                       grib2 message once this value can be
c                       determined, i.e. when the message has been
c                       completely packed. (output)
c                 ier = return status code. (output)
c                          0 = good return.
c                        1-4 = fatal error codes returned from pkbg.
c                       1002 = is0( ) has not been dimensioned large
c                              enough.
c             isevere = the severity level of the error.  the only
c                       value retuned is:
c                       2 = a fatal error  (output)
c                   * = alternate error return. (output)
c
c             local variables
c               igrib = contains the hexadecimal representation
c                       of the string "grib" as it would appear
c                       when the string is equivalenced to
c                       an integer*4 variable on a "big endian"
c                       platform.
c              igrib1 = contains the value of igrib if a 64-bit machine
c                       is being used.
c               izero = contains the value '0'.  this is used in the
c                       code simply to emphasize that a certain 
c                       group of octets in the message are being 
c                       zeroed out.
c                   n = l3264b = the integer word length in bits of
c                       the machine being used. values of 32 and
c                       64 are accommodated.
c             minsize = the minimum size for is0( ).  is0(7) is
c                       filled in pk_sect0, and is0(9) is filled in
c                       pk_grib.  this only applies when this is 
c                       a "new" message.      
c
c        non system subroutines called
c           pkbg
c
      parameter(minsize=9)
c
      dimension ipack(nd5),is0(ns0)
c
      data igrib/'47524942'x/
      data izero/0/
c
      n=l3264b
      ier=0
c
c        all errors generated by this routine are fatal.
      isevere=2
c
      if(new.eq.1)then
c
c           check minimum size of is0( ).
c
         if(ns0.lt.minsize)then
c              is0(7) filled in pk_sect0, and is0(9) filled in
c              pk_grib.         
            ier=1002
            go to 900
         endif
c
c           this is the first grid to be packed.
c           locn = word position in ipack( ) to start packing.
c           pkbg updates it.
         locn=1
c
c           ipos = bit position in ipack(locn) to start putting value.
c           pkbg updates it.
         ipos=1
c
c           ipack is zeroed if this is the first grid to be packed
c
         do k=1,nd5
            ipack(k)=0
         enddo
c
c           pack the word "grib".
c           accommodate a 64-bit word, if need be, by
c           moving the 4 characters to the right half of the word for
c           packing.
         igrib1=igrib
         if(l3264b.eq.64)igrib1=ishft(igrib,-32)
         call pkbg(kfildo,ipack,nd5,locn,ipos,igrib1,32,n,ier,*900)
         is0(1)=igrib
c
c           skip over bytes 5-6, which are reserved
         call pkbg(kfildo,ipack,nd5,locn,ipos,izero,16,n,ier,*900)
         is0(5)=0
         is0(6)=0
c
c           pack the discipline - master table number
         call pkbg(kfildo,ipack,nd5,locn,ipos,is0(7),8,n,ier,*900)
c
c           pack grib edition number.
         call pkbg(kfildo,ipack,nd5,locn,ipos,iedition,8,n,ier,*900)
         is0(8)=iedition
c
c           bytes 9-16 must be filled in later with the grib record
c           length in bytes; locn0_9 and ipos0_9 hold the location, and
c           the below statements hold the place.
         locn0_9=locn
         ipos0_9=ipos
         call pkbg(kfildo,ipack,nd5,locn,ipos,izero,32,n,ier,*900)
         call pkbg(kfildo,ipack,nd5,locn,ipos,izero,32,n,ier,*900)
c
c           note that only a 32-bit size is supported, which
c           accommodates over 4 gigabytes.
         do k=9,ns0
            is0(k)=0
         enddo
c
      else
c
c           when this is an addition to the product, locn and 
c           ipos have values on entry that represent the 7777
c           end of message.  this 7777 should be overwritten.
c           do it by writing a zero, then setting locn = locn-1.
c           also need to zero out the product size, so that it can be
c           filled in at the end.
c
         locn=locn-1
         call pkbg(kfildo,ipack,nd5,locn,ipos,izero,32,n,ier,*900)
         locn=locn-1
         call pkbg(kfildo,ipack,nd5,locn0_9,ipos0_9,izero,32,n,ier,*900)
         locn0_9=locn0_9-1
      endif
c
 900  if(ier.ne.0)return 1
c
      return
      end
