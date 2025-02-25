      subroutine pk_sect3(kfildo,ipack,nd5,is3,ns3,ipkopt,
     1                    l3264b,locn,ipos,ier,isevere,*)
c
c        march    2000   lawrence  gsc/tdl   original coding
c        january  2001   glahn     comments; added check on size of
c                                  is3( ); ier numbers altered
c        january  2001   glahn/lawrence removed new from call
c        february 2001   glahn/lawrence added template names
c        february 2002   glahn     removed test on ier=301 at end
c
c        purpose
c            packs section 3, the grid definition section, of a grib2
c            message.
c
c            this routine is capable of packing the following
c            grid definition templates:
c               template 3.0   equidistant cylindrical latitude/longitude
c               template 3.10  mercator
c               template 3.20  polar stereographic
c               template 3.30  lambert
c               template 3.90  orthographic space view
c               template 3.110 equatorial azimuthal equidistant
c               template 3.120 azimuth-range (radar)
c
c        data set use
c           kfildo - unit number for output (print) file. (output)
c
c        variables
c              kfildo = unit number for output (print) file. (input)
c            ipack(j) = the array that holds the actual packed message
c                       (j=1,nd5). (input/output)
c                 nd5 = the size of ipack( ). (input)
c              is3(j) = contains the grid definition data that
c                       will be packed into ipack( ) (j=1,ns3).
c                       (input/output)
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
c                        0 = good return.
c                       1-4 = error codes generated by pkbg. see the 
c                             documentation in the pkbg routine.
c                       5,6 = error codes generated by the length
c                             function. see the documentation for the
c                             length function.
c                       301 = is3(5) does not indicate section 5.
c                       302 = is3( ) has not been dimensioned large
c                             enough.
c                       303 = map projection in is3(13) is not
c                             supported.
c                       304 = returned by a routine like pk_polster
c                             which indicates it was incorrectly called.
c                             this should not happen.
c                       310 = unrecognized or unsupported shape of
c                             earth code in is3(15) returned by
c                             subroutine earth.
c             isevere = the severity level of the error.  the only
c                       value retuned is:
c                       2 = a fatal error  (output)
c                   * = alternate return when ier ne 0.
c
c             local variables
c               ipos3 = saves the bit position in loc3
c                       to store the length of section 3
c                       after the routine is done packing
c                       data into the section.
c                loc3 = saves the word position in ipack
c                       to store the length of section 3
c                       after the routine is done packing 
c                       data into the section.
c               izero = contains the value '0'.  this is used in the
c                       code simply to emphasize that a certain 
c                       group of octets in the message are being 
c                       zeroed out.
c                   k = a looping index variable.
c                   n = l3264b = the integer word length in bits of
c                       the machine being used. values of 32 and
c                       64 are accommodated.
c
c        non system subroutines called
c           length, pk_azimuth, pk_cylinder, pk_equator, pk_lambert,
c           pk_mercator, pk_polster, pk_orthographic, pkbg,
c
      parameter(minsize=13)
c
      dimension ipack(nd5),is3(ns3)
c
      data izero/0/
c
      n=l3264b
      ier=0
c
      loc3=locn
      ipos3=ipos
c
c        all errors generated by this routine are fatal.
      isevere=2
c
c        check minimum size of is3( ).  template sizes checked
c        in template subroutines.
c
      if(ns3.lt.minsize)then
         ier=302
         go to 900
      endif
c
c        check for correct section number.
c  
      if(is3(5).ne.3)then
         ier=301
         go to 900
      endif
c
c        bytes 1-4 must be filled in later with the record length
c        in bytes; below statement holds the place.  loc3 and ipos3
c        hold the location.
      call pkbg(kfildo,ipack,nd5,locn,ipos,izero,32,n,ier,*900)
c
c        pack the section number
      call pkbg(kfildo,ipack,nd5,locn,ipos,3,8,n,ier,*900)
c
c        pack source of grid definition
      call pkbg(kfildo,ipack,nd5,locn,ipos,is3(6),8,n,ier,*900)
c
c        pack number of data points
      call pkbg(kfildo,ipack,nd5,locn,ipos,is3(7),32,n,ier,*900)
c
c        pack the number of octets for optional list of numbers
c        defining number of points.
      call pkbg(kfildo,ipack,nd5,locn,ipos,is3(11),8,n,ier,*900)
c
c        pack the interpretation of list of numbers defining
c        number of points. 
      call pkbg(kfildo,ipack,nd5,locn,ipos,is3(12),8,n,ier,*900)
c
c        pack grid definition template number
      call pkbg(kfildo,ipack,nd5,locn,ipos,is3(13),16,n,ier,*900)
c
c        pack the values for the type of grid definition template
c
      select case (is3(13))
c
         case (0)
c
c              this is a latitude/longitude (or equidistant
c              cylindrical) projection.
            call pk_cylinder(kfildo,ipack,nd5,is3,ns3,ipkopt,n,
     1                       locn,ipos,ier,*900)
c
         case (10)
c
c              this is a mercator projection
            call pk_mercator(kfildo,ipack,nd5,is3,ns3,ipkopt,n,
     1                       locn,ipos,ier,*900)
c
         case (20)
c
c              polar stereographic map projection
            call pk_polster(kfildo,ipack,nd5,is3,ns3,ipkopt,n,
     1                      locn,ipos,ier,*900)
c
         case (30)
c
c              this is a lambert conformal projection
            call pk_lambert(kfildo,ipack,nd5,is3,ns3,ipkopt,n,
     1                      locn,ipos,ier,*900)
c
         case (90)
c
c              this is a space view perspective or
c              orthographic projection
            call pk_orthographic(kfildo,ipack,nd5,is3,ns3,
     1                           ipkopt,n,locn,ipos,ier,
     2                           *900)
c
         case (110)
c
c              equatorial azimuthal equidistant projection
            call pk_equator(kfildo,ipack,nd5,is3,ns3,ipkopt,n,
     1                      locn,ipos,ier,*900)
c
         case (120)
c
c              azimuthal range projection
            call pk_azimuth(kfildo,ipack,nd5,is3,ns3,ipkopt,n,
     1                      locn,ipos,ier,*900)
c
         case default
c
c              the map projection template is not supported.
            ier=303 
            go to 900
      end select
c
c        compute the length of the section and pack it. loc3 and
c        ipos3 represent the length of the record before section 3.
c        8 is the number of bits in a byte, and each section ends
c        at the end of a byte.
      is3(1)=length(kfildo,ipack,nd5,l3264b,loc3,ipos3,
     1              locn,ipos,ier)
c
c       error return section
 900  if(ier.ne.0)return 1
c
      return
      end
