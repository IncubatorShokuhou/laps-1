      subroutine pk_grib2(kfildo,ain,iain,nx,ny,idat,nidat,rdat,nrdat,
     1                    is0,ns0,is1,ns1,is3,ns3,is4,ns4,is5,
     2                    ns5,is6,ns6,is7,ns7,ib,ibitmap,ipack,nd5,
     3                    missp,xmissp,misss,xmisss,new,minpk,iclean,
     4                    l3264b,jer,ndjer,kjer)
c
c        august   1999   glahn/lawrence  gsc/tdl     original coding
c        february 2001   lawrence    updated this routine's 
c                                    prologue.
c        november 2001   glahn       added optional writing of message
c                                    length to kfildo
c        february 2002   glahn       added is6( ), ns6( ) to call to 
c                                    pk_sect7; changed call method and
c                                    sequence for pk_sect3
c
c        purpose
c           this utility packs a gridded data field according to the
c           rules and structure put forth in the grib version 2 
c           (grib2) documentation.  grib2 is a data representation
c           form for general regularly distributed information
c           in binary.  it provides a means of compressing and
c           encrypting data so as to make its transmission more
c           efficient.  it is recommended that the user of this
c           routine review the external documentation that 
c           explains what grib2 is and how to use it.
c
c           there are eight mandatory and one optional sections
c           contained within a grib2 message.  these sections are
c
c              section 0 - the indicator section
c              section 1 - the identification section
c              section 2 - the local use section (optional)
c              section 3 - the grid definition section
c              section 4 - the product definition section
c              section 5 - the data representation section
c              section 6 - the bit-map section
c              section 7 - the data section
c              section 8 - the end section
c
c           section 0, the indicator section, contains grib2
c           discipline, edition, and message length information.
c
c           section 1, the identification section, contains data
c           that describes characteristics that apply to all of the
c           processed data in the grib2 message.
c
c           section 2, the local use section, contains supplemental
c           information for local use by the originating centers.
c           this section is optional; it does not need to be included
c           in the grib2 message.
c
c           section 3, the grid definition section, describes the
c           geometry of the values within a plane described by
c           two fixed coordinates.
c
c           section 4, the product definition section, provides a
c           description of the data packed within section 7.
c
c           section 5, the data representation section, describes
c           what method is used to compress the data in section 7.
c
c           section 6, the bit-map section, contains a bit-map which
c           indicates the presence or absence of data at each grid
c           point.  a bit-map is only applicable to the simple
c           and complex packing methods.  a bit-map does not
c           apply to complex packing with spatial differences.
c
c           section 7, the data section, contains the packed data
c           values.
c
c           section 8, the end section, contains the string "7777"
c           indicating the end of the grib2 message.
c
c           sections 3, 4, 5, and 7 provide a variety of templates
c           that describe and define the gridded product
c           contained within the grib2 message.  this grib2 packer
c           does not support all of the templates provided in grib2.
c           the more commonly used templates are supported, and
c           this software can be easily extended in the future to
c           accommodate additional templates as the need arises.
c
c           the section 3 grid definition templates currently
c           recognized by this packer are:
c              template 3.0   equidistant cylindrical latitude/
c                             longitude
c              template 3.10  mercator
c              template 3.20  polar stereographic
c              template 3.30  lambert
c              template 3.90  orthographic space view
c              template 3.110 equatorial azimuthal equidistant
c              template 3.120 azimuth-range (radar)
c
c           the section 4 product definition templates currently
c           recognized by this packer are:
c              template 4.0  analysis or forecast at a level and point
c              template 4.1  individual emsemble
c              template 4.2  derived forecast based on ensembles
c              template 4.8  average, accumulation, extremes
c              template 4.20 radar
c              template 4.30 satellite
c
c           the section 5 data representation templates currently
c           recognized by this packer are:
c              template 5.0  simple packing
c              template 5.2  complex packing
c              template 5.3  complex packing with spatial differencing
c
c           the section 7 data templates currently recognized by this
c           packer are:
c              template 7.0  simple packing
c              template 7.2  complex packing
c              template 7.3  complex packing and spatial differencing
c
c           the user supplies the data for the sections of the grib2
c           message through arrays is0( ), is1( ), is3( ),
c           is4( ), is5( ), is6( ), and is7( ) which correspond 
c           to sections 0, 1, 2, 3, 4, 5, 6, and 7, respectively.
c           the user should refer to the external documentation
c           that accompanies this routine to determine which values
c           he needs to supply in each of these arrays.
c
c           the grid of data is passed into the packer using either the
c           ain( ) or the iain( ) arrays.  the data is passed into this 
c           routine using ain( ) if the data in the grid are floating
c           point values.  if the data in the grid contains integer  
c           values, then the data must be passed into this routine
c           using the iain( ) array.
c
c           the representation of missing values in the data field
c           plays an integral role in how the data is packed into
c           the grib2 message.  the handling of missing values in 
c           grib2 depends upon the packing method used to compress
c           the data into section 7 of the grib2 message.
c         
c           when using the simple packing method, a bit-map must be
c           used to indicate the locations of missing values in 
c           the data field.  this bit-map is ultimately packed into
c           section 6 of the grib2 message.  the corresponding data
c           field is packed into section 7 of the message with the
c           missing values removed.  when using the simple packing
c           method, the caller of this routine is the given the 
c           option of supplying a bit-map along with the data field
c           with the missing values removed from it.  the caller 
c           may also supply the data field with the missing values
c           in it and request that the grib2 packer generate a 
c           bit-map from it.  if the user has a data field that has no
c           missing data in it, then a bit-map is not needed and
c           should not be packed into section 6 of the grib2 message.  
c
c           the complex and complex with spatial differencing
c           packing methods allow for the use of primary
c           and secondary missing values.  these packing methods
c           utilize an efficient technique of packing the missing
c           values along with the data field into section 7 of the
c           grib2 message which eliminates the need for a bit-map
c           in section 6.  the missing values are left in the data 
c           field as it is packed.  in order for this to work, the
c           user must specify the missing value management value
c           in is5(23).  a value of "0" indicates that there are no
c           missing values in the data field.  a value of "1" indicates
c           that there may be primary missing values.  a value of "2"
c           indicates that there may be both primary and secondary
c           missing values.  the missp and misss (see below)
c           calling arguments represent the primary and secondary
c           missing values, respectively, when packing an integer 
c           data field.  the xmissp and xmisss (see below) calling
c           arguments represent the primary and secondary missing
c           values when packing a floating point data field.
c           like the simple packing method, the user has the option
c           to provide a bit-map and a data field with all of the
c           primary missing values removed.  however, instead of
c           packing this bit-map into section 6 of the grib2 message,
c           this packer will use the bit-map to place the primary
c           missing values back into the data field before packing
c           it.
c
c           a grib2 message may contain one or more gridded data
c           fields.  when packing multiple grids into a grib2 message
c           the first grid must provide information for sections 0
c           through 7 (with, of course, section 2 being optional).
c           subsequent grids should not provide data for sections
c           0 and 1.  they need only provide data for sections
c           2, 3, 4, 5, 6, 7, or 3, 4, 5, 6, 7, or 4, 5, 6, 7.
c           the final data grid packed into the grib2 message
c           must be followed by a section 8 indicating the end
c           of the grib2 message.  this routine must be called
c           once for each grid packed in the grib2 message.  the
c           "new" calling argument (see below) is used when
c           unpacking grib2 messages containing multiple
c           data grids.
c
c           there are a number of different error conditions that
c           that can be generated while packing a grib2 message.
c           error tracking is handled by the jer( , ) array (see below).
c           each row of this array can contain an error code followed
c           by its severity level.  the severity levels are
c           "0" for "ok", "1" for "warning", and "2" for "fatal".
c           the packer will immediately abort execution when
c           it encounters a fatal error.  however, execution is not
c           impeded by a status code of "ok" or "warning".  see the
c           argument list below for a complete list of the error
c           codes that can be returned by this routine.  unless 
c           otherwise stated, all of the error codes are fatal.
c
c        data set use
c           kfildo - unit number for output (print) file.  (output)
c
c        variables
c              kfildo = unit number of the output (print) file.
c                       (input)
c          ain(ix,jy) = array used to pass the grid of data into the 
c                       packer when packing a field consisting of
c                       floating point values (ix=1,nx) (jy=1,ny).
c                       (input)
c         iain(ix,jy) = array used to pass the grid of data into the
c                       packer when packing a field consisting of
c                       integer values (ix=1,nx) (jy=1,ny).  (input)
c               nx,ny = size of ain( , ), iain( , ), ib( , ).  (input) 
c             idat(j) = the array containing the local use groups
c                       consisting of integer data unpacked from
c                       section 2 (j=1,nidat).  (input)
c               nidat = the dimension of the idat( ) array.  (input)
c             rdat(j) = the array containing the local use groups
c                       consisting of floating point data unpacked
c                       from section 2 (j=1,nrdat).  (input)
c               nrdat = the dimension of the rdat( ) array.  (input)
c              is0(l) = holds the values for the grib indicator
c                       section, section 0 (l=1,ns0). (input/output)
c                 ns0 = size of is0( ).  (input)
c              is1(l) = holds the values for the grib identification
c                       section, section 1 (l=1,ns1).  (input/output)
c                 ns1 = size of is1( ).  (input)
c              is3(l) = holds the values for the grib grid definition
c                       section, section 3 (l=1,ns3).  (input/output)
c                 ns3 = size of is3( ).  (input)
c              is4(l) = holds the values for the grib product 
c                       definition section, section 4 (l=1,ns4).
c                       (input/output)
c                 ns4 = size of is4( ).  (input)
c              is5(l) = holds the values for the grib data 
c                       representation section, section 5 (l=1,ns5).
c                       (input/output)
c                 ns5 = size of is5( ).  (input)
c              is6(l) = holds the values for the bit-map section,
c                       section 6. (this is optional).
c                       (l=1,ns6). (input/output)
c                 ns6 = size of is6( ). must be nx*ny + 6 if a
c                       bit-map is to be used. (input)
c              is7(l) = holds the values for the data section,
c                       section 7 (l=1,ns7).  (input/output)
c                 ns7 = size of is7( ).  (input)
c           ib(ix,jy) = contains the primary bit-map.  the bit-map
c                       is used to represent the locations of missing
c                       values in the data field being packed into the
c                       grib2 message.  it has a one-to-one
c                       correspondence with the data points in the
c                       data grid.   a value of "1" in the bit-map
c                       indicates that the corresponding value in the
c                       data grid is valid.  a value of "0" in the
c                       bit-map indicates that the corresponding value
c                       in the data grid is missing.
c
c                       a bit-map accompanied by a data field with 
c                       the primary missing values removed from it
c                       may be supplied to this packer when using
c                       either the simple or complex packing methods.
c                       however, the bit-map is handled differently
c                       depending on which packing method is being 
c                       used.  for the simple packing method, 
c                       the bit-map is packed into section 6 of the
c                       grib2 message.  for the complex packing 
c                       method, the bit-map is used to insert the
c                       missing values into the data field which is
c                       then packed.  when using this packing method
c                       the bit-map is not packed into section 6 of
c                       the grib2 message.
c
c                       with the simple packing method, the user 
c                       also has the option of passing in a data field
c                       with missing values embedded in it and having
c                       the packer create the bit-map.  care must
c                       be exercised in setting the ibitmap and 
c                       iclean calling arguments (see below) as 
c                       these calling arguments inform the packer as
c                       to whether or not the user is supplying a 
c                       bit-map and whether or not the data in
c                       the gridded field contains missing values. 
c
c                       note that if a data field has no missing values
c                       in it, then a bit-map does not need to be
c                       packed into section 6 of the grib2 message
c                       when using the simple packing method. 
c                       (input/output)
c             ibitmap = flag indicating whether or not a bit-map
c                       is being passed into this packer.  if
c                       a bit-map is being passed into this 
c                       routine, then this argument should be set to
c                       "1".  if a bit-map is not being passed
c                       into this routine, then this argument should
c                       be set to "0".  the user must take care
c                       to set is6(6), the bit-map indicator,
c                       to the proper value.  it is recommended that
c                       if the user is passing in a bit-map then 
c                       he should remove the missing values from the
c                       data field and set the iclean flag to "1".
c                       if the user is using the simple packing 
c                       method and wants the packer to generate 
c                       a bit-map, then he should pass in the data
c                       field with the missing values embedded 
c                       in it along with the ibitmap and iclean flags
c                       set to "0",
c            ipack(j) = the array holding the actual packed message
c                       (j=1,max of nd5).  (output)
c                 nd5 = dimension of ipack( ).  (input)
c               missp = primary missing value, used when integer
c                       data is being packed.  (input)
c              xmissp = primary missing value, used when floating 
c                       point data is being packed.  (input)
c               misss = secondary missing value, used when integer
c                       data is being packed.  (input)
c              xmisss = secondary missing value, used when floating 
c                       point data is being packed.  (input)
c                 new = flag indicating whether or not this is the 
c                       first data grid to be packed into a 
c                       grib2 message.  a value of "1" indicates
c                       that this is the first grid to be packed into
c                       the message.  a value of "0" indicates that
c                       this is not the first data grid to be 
c                       packed into the message.  (input)
c               minpk = this calling argument applies only to the
c                       complex packing method.  it determines
c                       the minimum size of the groups of values that
c                       the data field is broken down into.  (input)  
c              iclean = flag indicating whether or not the
c                       data field (in iain( , ) or ain( , )) contains
c                       missing values.  a value of "0" means that
c                       the data field does contain missing values.
c                       a value of "1" means that the data field does
c                       not contain any missing values.  (input)
c              l3264b = integer word length of machine being used.
c                       values of 32 and 64 are accommodated.  (input)
c            jer(j,k) = return status codes and severity levels
c                       (j=1,ndjer)(k=1,2). values can come from
c                       subroutines; otherwise: 0 = good return.
c                       (output)
c               ndjer = the maximum number of error codes jer( , )
c                       may contain. (input)
c                kjer = the actual number of error messages contained
c                       in jer( , ). (output)
c
c        local variables
c            a(ix,jy) = grid of real values to be packed (ix=1,nx)
c                       (jy=1,ny).  (automatic array)
c              exists = boolean flag indicating whether or not a
c                       grib2 section exists (logical).
c           ia(ix,jy) = grid of integer values to be packed (ix=1,nx)
c                       (jy=1,ny).  (automatic array)
c                  id = the decimal scaling factor.
c                  ie = the binary scaling factor.
c            iedition = the version number of the grib2 encoder.
c                 ier = return code from most subroutines.  an error
c                       will be transferred to jer(kjer, ) upon return.
c                        0 - good return
c                      1-4 - error codes returned from pkbg.
c                      5,6 - error codes generated by the length
c                            function.
c                        8 - is5(21) contains an invalid data type
c                            value. 
c                  101,102 - error codes generated by pk_sect1.
c              301-304,310 - error codes generated by pk_sect3. 
c                  401-403 - error codes generated by pk_sect4. 
c              501,502,508 - error codes generated by pk_sect5.
c                  601,602 - error codes generated by pk_sect6.
c          701-703,711-713 - error codes generated by pk_sect7.
c                      999 - error codes generated by pk_trace.
c                     1002 - error code returned from pk_sect0.
c                 inc = number of values to add to the group to be
c                       packed at a time.  set to 1 in data statement.
c              ioctet = the length of the message in bytes (octets).
c             ipos0_9 = saves the beginning bit position in ipack( )
c                       that the total message length will be stored
c                       at once the grib2 message has been completely
c                       packed. 
c             isevere = the severity level of the error. acceptable
c                       severity levels are:
c                          0 = not a problem ... just a diagnostic
c                          1 = a warning ... something the user should
c                              know about
c                          2 = a fatal error
c            locn0_9  = saves the beginning word position in ipack( ) 
c                       that the total message length will be stored
c                       at once the grib2 message has been completely
c                       packed.
c                mina = the reference value when the original data
c                       are integer.
c               xmina = the reference value when the original data
c                       are floating point.
c        1         2         3         4         5         6         7 x
c
c        non-system subroutines called
c           pk_sect0, pk_sect1, pk_sect2, pk_sect3, pk_sect4,
c           pk_sect5, pk_sect6, pk_sect7, pk_sect8, pk_trace, pkbg,
c           prepr
c
      logical exists
c
      dimension a(nx,ny)
      dimension ia(nx,ny)
      dimension ain(nx,ny)
      dimension iain(nx,ny)
      dimension ib(nx,ny)
      dimension ipack(nd5)
      dimension is0(ns0),is1(ns1),is3(ns3),
     1          is4(ns4),is5(ns5),is6(ns6),is7(ns7)
      dimension jer(ndjer,2)
      dimension idat(nidat),rdat(nrdat)
c
      data exists/.false./
      data iedition/2/
      data inc/1/
c
      save locn, ipos,
     1     locn0_9, ipos0_9
c
      ier=0
      isevere=2
      kjer=0
c
c        initialize jer( , ).
c
      do k=1,ndjer
         jer(k,1)=0
         jer(k,2)=0
      enddo
c
c        determine what type of data the user is packing.
c        it is necessary to make local copies of iain( , ) for
c        integer data or ain( , ) for floating point data to
c        prevent the user-supplied data in these arrays from
c        getting altered by the packing routines.
c
      if(is5(21).eq.0)then
c
         do 20 iy=1,ny
            do 10 ix=1,nx
               a(ix,iy)=ain(ix,iy)
  10        continue
  20     continue
c
      else if(is5(21).eq.1)then
c
         do 40 iy=1,ny
            do 30 ix=1,nx
               ia(ix,iy)=iain(ix,iy)
  30        continue
  40     continue
c
      else
c
c           invalid option in is5(21).
         ier=8
         goto 900
      endif
c
c        *************************************
c
c        pack section 0 of the message into ipack( ).
c
c        *************************************
c
      call pk_trace(kfildo,jer,ndjer,kjer,0,0)
      call pk_sect0(kfildo,ipack,nd5,is0,ns0,l3264b,new,locn,ipos,
     1              iedition,locn0_9,ipos0_9,ier,isevere,*900)
c
c        *************************************
c
c        pack section 1, identification.
c
c        *************************************
c
      call pk_trace(kfildo,jer,ndjer,kjer,100,0)
      call pk_sect1(kfildo,ipack,nd5,is1,ns1,new,l3264b,locn,ipos,
     1              ier,isevere,*900)
c
c        ********************************************************
c
c        pack section 2, the local use section, if it is present.
c
c        ********************************************************
c
      call pk_trace(kfildo,jer,ndjer,kjer,200,0)
      call pk_sect2(kfildo,ipack,nd5,rdat,nrdat,idat,nidat,
     1              l3264b,locn,ipos,exists,ier,isevere,*900)
c
c        *************************************
c
c        prepare ia( ) for packing and ib( ) if necessary.
c        make decisions about packing, etc. since these
c        decisions affect section 3 and section 5, the
c        prepr routine is called before either of those
c        sections are processed.
c
c        *************************************
c
      id=is5(18)
      ie=is5(16)
      call pk_trace(kfildo,jer,ndjer,kjer,900,0)
      call prepr(kfildo,a,ia,ib,nx,ny,nval,iclean,ibitmap,
     1           is5,ns5,is6,ns6,is7,ns7,id,ie,mina,xmina,
     2           missp,misss,xmissp,xmisss,minpk,ipkopt,ier,jer,
     3           ndjer,kjer,*900)
c
c        *************************************
c
c        pack section 3, grid definition.
c
c        *************************************
c
c        if this grid is being packed into an existing message,
c        check to see if there is a section 3 ... note that there
c        must be a section 3 if there was a section 2.
c
      call pk_trace(kfildo,jer,ndjer,kjer,300,0)
c
      if(is3(5).ne.3)then
         if(exists)then
            ier=305
            isevere=2
            go to 900
         endif
c
      endif
         call pk_sect3(kfildo,ipack,nd5,is3,ns3,ipkopt,l3264b,
     1                 locn,ipos,ier,isevere,*900)
c 
c        *************************************
c
c        pack section 4, product definition.
c
c        *************************************
c
      call pk_trace(kfildo,jer,ndjer,kjer,400,0)
      call pk_sect4(kfildo,ipack,nd5,is4,ns4,l3264b,locn,ipos,
     1              ier,isevere,*900)
c
c        *********************************************
c
c        pack section 5, data representation.
c
c        *********************************************
c
      call pk_trace(kfildo,jer,ndjer,kjer,500,0)
      call pk_sect5(kfildo,ipack,nd5,is5,ns5,mina,xmina,
     1              missp,xmissp,misss,xmisss,l3264b,
     2              locn,ipos,locn5_20,ipos5_20,locn5_32,ipos5_32,
     3              ier,isevere,*900)
c     call timpr(kfildo,kfildo,'end section 5        ')
c
c        *******************************************
c
c        pack section 6, bit map.
c
c        *******************************************
c
      call pk_trace(kfildo,jer,ndjer,kjer,600,0)
      call pk_sect6(kfildo,ib,nx*ny,ipack,nd5,locn,ipos,
     1              is6,ns6,ipkopt,l3264b,ier,isevere,*900)
c     call timpr(kfildo,kfildo,'end section 6        ')
c
c        ***********************************************
c
c        pack section 7, binary data.
c
c        ***********************************************
c
      call pk_trace(kfildo,jer,ndjer,kjer,700,0)
      call pk_sect7(kfildo,ia,nx*ny,ipack,nd5,locn,ipos,
     1              is5,ns5,is6,ns6,is7,ns7,inc,minpk,missp,misss,
     2              ipkopt,locn5_20,ipos5_20,
     3              locn5_32,ipos5_32,l3264b,ier,isevere,*900)
c     
      if(ier.ne.0)then
         call pk_trace(kfildo,jer,ndjer,kjer,ier,1)
c           this call to trace is for non fatal errors that
c           can be carried back.  fatal errors are returned
c           to *900.
      endif
c     call timpr(kfildo,kfildo,'end section 7        ')
c
c        *************************************
c
c        pack section 8, end of message.
c
c        *************************************
c
      call pk_trace(kfildo,jer,ndjer,kjer,800,0)
      call pk_sect8(kfildo,ipack,nd5,locn,ipos,l3264b,ier,
     1              isevere,*900)
c
c        fill bytes 9-12 with the total message length in bytes.
      ioctet=locn*l3264b/8-(l3264b+1-ipos)/8
      is0(9)=ioctet
c
      call pkbg(kfildo,ipack,nd5,locn0_9+1,ipos0_9,is0(9),32,l3264b,
     1          ier,*900)
c        this total length is limited to a 4-byte word.  the length is
c        put into is0(9), but is packed into octets 13-16, which is the
c        right half of a 64-bit word.
c
c        write the total message length.
c
c     write(kfildo,800)ioctet
c800  format(/' total message length in bytes =',i12)
c
      return
c
c        error return section.
c
 900  call pk_trace(kfildo,jer,ndjer,kjer,ier,isevere)
c
 901  return
      end
