      subroutine getfield(cgrib,lcgrib,ifldnum,igds,igdstmpl,igdslen,
     &                    ideflist,idefnum,ipdsnum,ipdstmpl,ipdslen,
     &                    coordlist,numcoord,ndpts,idrsnum,idrstmpl,
     &                    idrslen,ibmap,bmap,fld,ierr)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    getfield 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-26
!
! abstract: this subroutine returns the grid definition, product definition,
!   bit-map ( if applicable ), and the unpacked data for a given data
!   field.  since there can be multiple data fields packed into a grib2
!   message, the calling routine indicates which field is being requested
!   with the ifldnum argument.
!
! program history log:
! 2000-05-26  gilbert
!
! usage:    call getfield(cgrib,lcgrib,ifldnum,igds,igdstmpl,igdslen,
!    &                    ideflist,idefnum,ipdsnum,ipdstmpl,ipdslen,
!    &                    coordlist,numcoord,ndpts,idrsnum,idrstmpl,
!    &                    idrslen,ibmap,bmap,fld,ierr)
!   input argument list:
!     cgrib    - character array that contains the grib2 message
!     lcgrib   - length (in bytes) of grib message array cgrib.
!     ifldnum  - specifies which field in the grib2 message to return.
!
!   output argument list:      
!     igds     - contains information read from the appropriate grib grid 
!                definition section 3 for the field being returned.
!                must be dimensioned >= 5.
!                igds(1)=source of grid definition (see code table 3.0)
!                igds(2)=number of grid points in the defined grid.
!                igds(3)=number of octets needed for each 
!                            additional grid points definition.  
!                            used to define number of
!                            points in each row ( or column ) for
!                            non-regular grids.  
!                            = 0, if using regular grid.
!                igds(4)=interpretation of list for optional points 
!                            definition.  (code table 3.11)
!                igds(5)=grid definition template number (code table 3.1)
!     igdstmpl - contains the data values for the specified grid definition
!                template ( nn=igds(5) ).  each element of this integer 
!                array contains an entry (in the order specified) of grid
!                defintion template 3.nn
!                a safe dimension for this array can be obtained in advance
!                from maxvals(2), which is returned from subroutine gribinfo.
!     igdslen  - number of elements in igdstmpl().  i.e. number of entries
!                in grid defintion template 3.nn  ( nn=igds(5) ).
!     ideflist - (used if igds(3) .ne. 0)  this array contains the
!                number of grid points contained in each row ( or column ).
!                (part of section 3)
!                a safe dimension for this array can be obtained in advance
!                from maxvals(3), which is returned from subroutine gribinfo.
!     idefnum  - (used if igds(3) .ne. 0)  the number of entries
!                in array ideflist.  i.e. number of rows ( or columns )
!                for which optional grid points are defined.
!     ipdsnum  - product definition template number ( see code table 4.0)
!     ipdstmpl - contains the data values for the specified product definition
!                template ( n=ipdsnum ).  each element of this integer
!                array contains an entry (in the order specified) of product
!                defintion template 4.n
!                a safe dimension for this array can be obtained in advance
!                from maxvals(4), which is returned from subroutine gribinfo.
!     ipdslen  - number of elements in ipdstmpl().  i.e. number of entries
!                in product defintion template 4.n  ( n=ipdsnum ).
!     coordlist- array containg floating point values intended to document
!                the vertical discretisation associated to model data
!                on hybrid coordinate vertical levels.  (part of section 4)
!                the dimension of this array can be obtained in advance
!                from maxvals(5), which is returned from subroutine gribinfo.
!     numcoord - number of values in array coordlist.
!     ndpts    - number of data points unpacked and returned.
!     idrsnum  - data representation template number ( see code table 5.0)
!     idrstmpl - contains the data values for the specified data representation
!                template ( n=idrsnum ).  each element of this integer
!                array contains an entry (in the order specified) of product
!                defintion template 5.n
!                a safe dimension for this array can be obtained in advance
!                from maxvals(6), which is returned from subroutine gribinfo.
!     idrslen  - number of elements in idrstmpl().  i.e. number of entries
!                in data representation template 5.n  ( n=idrsnum ).
!     ibmap    - bitmap indicator ( see code table 6.0 )
!                0 = bitmap applies and is included in section 6.
!                1-253 = predefined bitmap applies
!                254 = previously defined bitmap applies to this field
!                255 = bit map does not apply to this product.
!     bmap()   - logical*1 array containing decoded bitmap. ( if ibmap=0 )
!                the dimension of this array can be obtained in advance
!                from maxvals(7), which is returned from subroutine gribinfo.
!     fld()    - array of ndpts unpacked data points.
!                a safe dimension for this array can be obtained in advance
!                from maxvals(7), which is returned from subroutine gribinfo.
!     ierr     - error return code.
!                0 = no error
!                1 = beginning characters "grib" not found.
!                2 = grib message is not edition 2.
!                3 = the data field request number was not positive.
!                4 = end string "7777" found, but not where expected.
!                6 = grib message did not contain the requested number of
!                    data fields.
!                7 = end string "7777" not found at end of message.
!                9 = data representation template 5.nn not yet implemented.
!               10 = error unpacking section 3.
!               11 = error unpacking section 4.
!               12 = error unpacking section 5.
!               13 = error unpacking section 6.
!               14 = error unpacking section 7.
!
! remarks: note that subroutine gribinfo can be used to first determine
!          how many data fields exist in a given grib message.
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$

      character(len=1),intent(in) :: cgrib(lcgrib)
      integer,intent(in) :: lcgrib,ifldnum
      integer,intent(out) :: igds(*),igdstmpl(*),ideflist(*)
      integer,intent(out) :: ipdsnum,ipdstmpl(*)
      integer,intent(out) :: idrsnum,idrstmpl(*)
      integer,intent(out) :: ndpts,ibmap,idefnum,numcoord
      integer,intent(out) :: ierr
      logical*1,intent(out) :: bmap(*)
      real,intent(out) :: fld(*),coordlist(*)
      
      character(len=4),parameter :: grib='grib',c7777='7777'
      character(len=4) :: ctemp
      integer:: listsec0(2)
      integer iofst,ibeg,istart
      integer(4) :: ieee
      logical have3,have4,have5,have6,have7

      have3=.false.
      have4=.false.
      have5=.false.
      have6=.false.
      have7=.false.
      ierr=0
      numfld=0
!
!  check for valid request number
!  
      if (ifldnum.le.0) then
        print *,'getfield: request for field number must be positive.'
        ierr=3
        return
      endif
!
!  check for beginning of grib message in the first 100 bytes
!
      istart=0
      do j=1,100
        ctemp=cgrib(j)//cgrib(j+1)//cgrib(j+2)//cgrib(j+3)
        if (ctemp.eq.grib ) then
          istart=j
          exit
        endif
      enddo
      if (istart.eq.0) then
        print *,'getfield:  beginning characters grib not found.'
        ierr=1
        return
      endif
!
!  unpack section 0 - indicator section 
!
      iofst=8*(istart+5)
      call gbyte(cgrib,listsec0(1),iofst,8)     ! discipline
      iofst=iofst+8
      call gbyte(cgrib,listsec0(2),iofst,8)     ! grib edition number
      iofst=iofst+8
      iofst=iofst+32
      call gbyte(cgrib,lengrib,iofst,32)        ! length of grib message
      iofst=iofst+32
      lensec0=16
      ipos=istart+lensec0
!
!  currently handles only grib edition 2.
!  
      if (listsec0(2).ne.2) then
        print *,'getfield: can only decode grib edition 2.'
        ierr=2
        return
      endif
!
!  loop through the remaining sections keeping track of the 
!  length of each.  also keep the latest grid definition section info.
!  unpack the requested field number.
!
      do
        !    check to see if we are at end of grib message
        ctemp=cgrib(ipos)//cgrib(ipos+1)//cgrib(ipos+2)//cgrib(ipos+3)
        if (ctemp.eq.c7777 ) then
          ipos=ipos+4
          !    if end of grib message not where expected, issue error
          if (ipos.ne.(istart+lengrib)) then
            print *,'getfield: "7777" found, but not where expected.'
            ierr=4
            return
          endif
          exit
        endif
        !     get length of section and section number
        iofst=(ipos-1)*8
        call gbyte(cgrib,lensec,iofst,32)        ! get length of section
        iofst=iofst+32
        call gbyte(cgrib,isecnum,iofst,8)         ! get section number
        iofst=iofst+8
        !print *,' lensec= ',lensec,'    secnum= ',isecnum
        !
        !   if found section 3, unpack the gds info using the 
        !   appropriate template.  save in case this is the latest
        !   grid before the requested field.
        !
        if (isecnum.eq.3) then
          iofst=iofst-40       ! reset offset to beginning of section
          call unpack3(cgrib,lcgrib,iofst,igds,igdstmpl,igdslen,
     &                 ideflist,idefnum,jerr)
          if (jerr.eq.0) then
            have3=.true.
          else
            ierr=10
            return
          endif
        endif
        !
        !   if found section 4, check to see if this field is the
        !   one requested.
        !
        if (isecnum.eq.4) then
          numfld=numfld+1
          if (numfld.eq.ifldnum) then
            iofst=iofst-40       ! reset offset to beginning of section
            call unpack4(cgrib,lcgrib,iofst,ipdsnum,ipdstmpl,ipdslen,
     &                   coordlist,numcoord,jerr)
            if (jerr.eq.0) then
              have4=.true.
            else
              ierr=11
              return
            endif
          endif
        endif
        !
        !   if found section 5, check to see if this field is the
        !   one requested.
        !
        if ((isecnum.eq.5).and.(numfld.eq.ifldnum)) then
          iofst=iofst-40       ! reset offset to beginning of section
          call unpack5(cgrib,lcgrib,iofst,ndpts,idrsnum,idrstmpl,
     &                 idrslen,jerr)
          if (jerr.eq.0) then
            have5=.true.
          else
            ierr=12
            return
          endif
        endif
        !
        !   if found section 6, unpack bitmap.
        !   save in case this is the latest
        !   bitmap before the requested field.
        !
        if (isecnum.eq.6) then
          iofst=iofst-40       ! reset offset to beginning of section
          call unpack6(cgrib,lcgrib,iofst,igds(2),ibmap,bmap,jerr)
          if (jerr.eq.0) then
            have6=.true.
          else
            ierr=13
            return
          endif
        endif
        !
        !   if found section 7, check to see if this field is the
        !   one requested.
        !
        if ((isecnum.eq.7).and.(numfld.eq.ifldnum)) then
          if (idrsnum.eq.0) then
            call simunpack(cgrib(ipos+5),lensec-6,idrstmpl,ndpts,fld)
            have7=.true.
          elseif (idrsnum.eq.2.or.idrsnum.eq.3) then
            call comunpack(cgrib(ipos+5),lensec-6,lensec,idrsnum,
     &                     idrstmpl,ndpts,fld,ier)
            if ( ier .ne. 0 ) then
                ierr=14
                return
            endif
            have7=.true.
          elseif (idrsnum.eq.50) then
            call simunpack(cgrib(ipos+5),lensec-6,idrstmpl,ndpts-1,
     &                     fld(2))
            ieee=idrstmpl(5)
            call rdieee(ieee,fld(1),1)
            have7=.true.
          else
            print *,'getfield: data representation template ',idrsnum,
     &              ' not yet implemented.'
            ierr=9
            return
          endif
        endif
        !
        !   check to see if we read pass the end of the grib
        !   message and missed the terminator string '7777'.
        !
        ipos=ipos+lensec                 ! update beginning of section pointer
        if (ipos.gt.(istart+lengrib)) then
          print *,'getfield: "7777"  not found at end of grib message.'
          ierr=7
          return
        endif

        if (have3.and.have4.and.have5.and.have6.and.have7) return
        
      enddo

!
!  if exited from above loop, the end of the grib message was reached
!  before the requested field was found.
!
      print *,'getfield: grib message contained ',numlocal,
     &        ' different fields.'
      print *,'getfield: the request was for the ',ifldnum,
     &        ' field.'
      ierr=6

      return
      end


      subroutine unpack3(cgrib,lcgrib,iofst,igds,igdstmpl,
     &                   mapgridlen,ideflist,idefnum,ierr)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    unpack3 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-26
!
! abstract: this subroutine unpacks section 3 (grid definition section)
!   starting at octet 6 of that section.  
!
! program history log:
! 2000-05-26  gilbert
!
! usage:    call unpack3(cgrib,lcgrib,lensec,iofst,igds,igdstmpl,
!    &                   mapgridlen,ideflist,idefnum,ierr)
!   input argument list:
!     cgrib    - character array that contains the grib2 message
!     lcgrib   - length (in bytes) of grib message array cgrib.
!     iofst    - bit offset of the beginning of section 3.
!
!   output argument list:      
!     iofst    - bit offset at the end of section 3, returned.
!     igds     - contains information read from the appropriate grib grid 
!                definition section 3 for the field being returned.
!                must be dimensioned >= 5.
!                igds(1)=source of grid definition (see code table 3.0)
!                igds(2)=number of grid points in the defined grid.
!                igds(3)=number of octets needed for each 
!                            additional grid points definition.  
!                            used to define number of
!                            points in each row ( or column ) for
!                            non-regular grids.  
!                            = 0, if using regular grid.
!                igds(4)=interpretation of list for optional points 
!                            definition.  (code table 3.11)
!                igds(5)=grid definition template number (code table 3.1)
!     igdstmpl - contains the data values for the specified grid definition
!                template ( nn=igds(5) ).  each element of this integer 
!                array contains an entry (in the order specified) of grid
!                defintion template 3.nn
!     mapgridlen- number of elements in igdstmpl().  i.e. number of entries
!                in grid defintion template 3.nn  ( nn=igds(5) ).
!     ideflist - (used if igds(3) .ne. 0)  this array contains the
!                number of grid points contained in each row ( or column ).
!                (part of section 3)
!     idefnum  - (used if igds(3) .ne. 0)  the number of entries
!                in array ideflist.  i.e. number of rows ( or columns )
!                for which optional grid points are defined.
!     ierr     - error return code.
!                0 = no error
!                5 = "grib" message contains an undefined grid definition
!                    template.
!
! remarks: uses fortran 90 module gridtemplates.
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$

      use gridtemplates

      character(len=1),intent(in) :: cgrib(lcgrib)
      integer,intent(in) :: lcgrib
      integer,intent(inout) :: iofst
      integer,intent(out) :: igds(*),igdstmpl(*),ideflist(*)
      integer,intent(out) :: ierr,idefnum

      integer,allocatable :: mapgrid(:)
      integer :: mapgridlen,ibyttem
      logical needext

      ierr=0

      call gbyte(cgrib,lensec,iofst,32)        ! get length of section
      iofst=iofst+32
      iofst=iofst+8     ! skip section number

      call gbyte(cgrib,igds(1),iofst,8)     ! get source of grid def.
      iofst=iofst+8
      call gbyte(cgrib,igds(2),iofst,32)    ! get number of grid pts.
      iofst=iofst+32
      call gbyte(cgrib,igds(3),iofst,8)     ! get num octets for opt. list
      iofst=iofst+8
      call gbyte(cgrib,igds(4),iofst,8)     ! get interpret. for opt. list
      iofst=iofst+8
      call gbyte(cgrib,igds(5),iofst,16)    ! get grid def template num.
      iofst=iofst+16
      if (igds(1).eq.0) then
!      if (igds(1).eq.0.or.igds(1).eq.255) then  ! for ecmwf test only
        allocate(mapgrid(lensec))
        !   get grid definition template
        call getgridtemplate(igds(5),mapgridlen,mapgrid,needext,
     &                       iret)
        if (iret.ne.0) then
          ierr=5
          return
        endif
      else
!        igdstmpl=-1
        mapgridlen=0
        needext=.false.
      endif
      !
      !   unpack each value into array igdstmpl from the
      !   the appropriate number of octets, which are specified in
      !   corresponding entries in array mapgrid.
      !
      ibyttem=0
      do i=1,mapgridlen
        nbits=iabs(mapgrid(i))*8
        if ( mapgrid(i).ge.0 ) then
          call gbyte(cgrib,igdstmpl(i),iofst,nbits)
        else
          call gbyte(cgrib,isign,iofst,1)
          call gbyte(cgrib,igdstmpl(i),iofst+1,nbits-1)
          if (isign.eq.1) igdstmpl(i)=-igdstmpl(i)
        endif
        iofst=iofst+nbits
        ibyttem=ibyttem+iabs(mapgrid(i))
      enddo
      !
      !   check to see if the grid definition template needs to be
      !   extended.
      !   the number of values in a specific template may vary
      !   depending on data specified in the "static" part of the
      !   template.
      !
      if ( needext ) then
        call extgridtemplate(igds(5),igdstmpl,newmapgridlen,mapgrid)
        !   unpack the rest of the grid definition template
        do i=mapgridlen+1,newmapgridlen
          nbits=iabs(mapgrid(i))*8
          if ( mapgrid(i).ge.0 ) then
            call gbyte(cgrib,igdstmpl(i),iofst,nbits)
          else
            call gbyte(cgrib,isign,iofst,1)
            call gbyte(cgrib,igdstmpl(i),iofst+1,nbits-1)
            if (isign.eq.1) igdstmpl(i)=-igdstmpl(i)
          endif
          iofst=iofst+nbits
          ibyttem=ibyttem+iabs(mapgrid(i))
        enddo
        mapgridlen=newmapgridlen
      endif
      !
      !   unpack optional list of numbers defining number of points
      !   in each row or column, if included.  this is used for non regular
      !   grids.
      !
      if ( igds(3).ne.0 ) then
         nbits=igds(3)*8
         idefnum=(lensec-14-ibyttem)/igds(3)
         call gbytes(cgrib,ideflist,iofst,nbits,0,idefnum)
         iofst=iofst+(nbits*idefnum)
      else
         idefnum=0
      endif
      if( allocated(mapgrid) ) deallocate(mapgrid)
      return    ! end of section 3 processing
      end


      subroutine unpack4(cgrib,lcgrib,iofst,ipdsnum,ipdstmpl,mappdslen,
     &                   coordlist,numcoord,ierr)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    unpack4 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-26
!
! abstract: this subroutine unpacks section 4 (product definition section)
!   starting at octet 6 of that section.  
!
! program history log:
! 2000-05-26  gilbert
!
! usage:    call unpack4(cgrib,lcgrib,iofst,ipdsnum,ipdstmpl,mappdslen,
!    &                   coordlist,numcoord,ierr)
!   input argument list:
!     cgrib    - character array that contains the grib2 message
!     lcgrib   - length (in bytes) of grib message array cgrib.
!     iofst    - bit offset of the beginning of section 4.
!
!   output argument list:      
!     iofst    - bit offset of the end of section 4, returned.
!     ipdsnum  - product definition template number ( see code table 4.0)
!     ipdstmpl - contains the data values for the specified product definition
!                template ( n=ipdsnum ).  each element of this integer
!                array contains an entry (in the order specified) of product
!                defintion template 4.n
!     mappdslen- number of elements in ipdstmpl().  i.e. number of entries
!                in product defintion template 4.n  ( n=ipdsnum ).
!     coordlist- array containg floating point values intended to document
!                the vertical discretisation associated to model data
!                on hybrid coordinate vertical levels.  (part of section 4)
!     numcoord - number of values in array coordlist.
!     ierr     - error return code.
!                0 = no error
!                5 = "grib" message contains an undefined product definition
!                    template.
!
! remarks: uses fortran 90 module pdstemplates.
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$

      use pdstemplates

      character(len=1),intent(in) :: cgrib(lcgrib)
      integer,intent(in) :: lcgrib
      integer,intent(inout) :: iofst
      real,intent(out) :: coordlist(*)
      integer,intent(out) :: ipdsnum,ipdstmpl(*)
      integer,intent(out) :: ierr,numcoord

      real(4),allocatable :: coordieee(:)
      integer,allocatable :: mappds(:)
      integer :: mappdslen
      logical needext

      ierr=0

      call gbyte(cgrib,lensec,iofst,32)        ! get length of section
      iofst=iofst+32
      iofst=iofst+8     ! skip section number
      allocate(mappds(lensec))

      call gbyte(cgrib,numcoord,iofst,16)    ! get num of coordinate values
      iofst=iofst+16
      call gbyte(cgrib,ipdsnum,iofst,16)    ! get prod. def template num.
      iofst=iofst+16
      !   get product definition template
      call getpdstemplate(ipdsnum,mappdslen,mappds,needext,iret)
      if (iret.ne.0) then
        ierr=5
        return
      endif
      !
      !   unpack each value into array ipdstmpl from the
      !   the appropriate number of octets, which are specified in
      !   corresponding entries in array mappds.
      !
      do i=1,mappdslen
        nbits=iabs(mappds(i))*8
        if ( mappds(i).ge.0 ) then
          call gbyte(cgrib,ipdstmpl(i),iofst,nbits)
        else
          call gbyte(cgrib,isign,iofst,1)
          call gbyte(cgrib,ipdstmpl(i),iofst+1,nbits-1)
          if (isign.eq.1) ipdstmpl(i)=-ipdstmpl(i)
        endif
        iofst=iofst+nbits
      enddo
      !
      !   check to see if the product definition template needs to be
      !   extended.
      !   the number of values in a specific template may vary
      !   depending on data specified in the "static" part of the
      !   template.
      !
      if ( needext ) then
        call extpdstemplate(ipdsnum,ipdstmpl,newmappdslen,mappds)
        !   unpack the rest of the product definition template
        do i=mappdslen+1,newmappdslen
          nbits=iabs(mappds(i))*8
          if ( mappds(i).ge.0 ) then
            call gbyte(cgrib,ipdstmpl(i),iofst,nbits)
          else
            call gbyte(cgrib,isign,iofst,1)
            call gbyte(cgrib,ipdstmpl(i),iofst+1,nbits-1)
            if (isign.eq.1) ipdstmpl(i)=-ipdstmpl(i)
          endif
          iofst=iofst+nbits
        enddo
        mappdslen=newmappdslen
      endif
      !
      !   get optional list of vertical coordinate values
      !   after the product definition template, if necessary.
      !
      if ( numcoord .ne. 0 ) then
        allocate (coordieee(numcoord))
        call gbytes(cgrib,coordieee,iofst,32,0,numcoord)
        call rdieee(coordieee,coordlist,numcoord)
        deallocate (coordieee)
        iofst=iofst+(32*numcoord)
      endif
      if( allocated(mappds) ) deallocate(mappds)
      return    ! end of section 4 processing
      end


      subroutine unpack5(cgrib,lcgrib,iofst,ndpts,idrsnum,idrstmpl,
     &                   mapdrslen,ierr)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    unpack5 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-26
!
! abstract: this subroutine unpacks section 5 (data representation section)
!   starting at octet 6 of that section.  
!
! program history log:
! 2000-05-26  gilbert
!
! usage:    call unpack5(cgrib,lcgrib,iofst,ndpts,idrsnum,idrstmpl,
!                        mapdrslen,ierr)
!   input argument list:
!     cgrib    - character array that contains the grib2 message
!     lcgrib   - length (in bytes) of grib message array cgrib.
!     iofst    - bit offset of the beginning of section 5.
!
!   output argument list:      
!     iofst    - bit offset at the end of section 5, returned.
!     ndpts    - number of data points unpacked and returned.
!     idrsnum  - data representation template number ( see code table 5.0)
!     idrstmpl - contains the data values for the specified data representation
!                template ( n=idrsnum ).  each element of this integer
!                array contains an entry (in the order specified) of data
!                representation template 5.n
!     mapdrslen- number of elements in idrstmpl().  i.e. number of entries
!                in data representation template 5.n  ( n=idrsnum ).
!     ierr     - error return code.
!                0 = no error
!                7 = "grib" message contains an undefined data
!                    representation template.
!
! remarks: none
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$

      use drstemplates

      character(len=1),intent(in) :: cgrib(lcgrib)
      integer,intent(in) :: lcgrib
      integer,intent(inout) :: iofst
      integer,intent(out) :: ndpts,idrsnum,idrstmpl(*)
      integer,intent(out) :: ierr

c      integer,allocatable :: mapdrs(:)
      integer,allocatable :: mapdrs(:)
      integer :: mapdrslen
      logical needext

      ierr=0

      call gbyte(cgrib,lensec,iofst,32)        ! get length of section
      iofst=iofst+32
      iofst=iofst+8     ! skip section number
      allocate(mapdrs(lensec))

      call gbyte(cgrib,ndpts,iofst,32)    ! get num of data points
      iofst=iofst+32
      call gbyte(cgrib,idrsnum,iofst,16)     ! get data rep template num.
      iofst=iofst+16
      !   gen data representation template
      call getdrstemplate(idrsnum,mapdrslen,mapdrs,needext,iret)
      if (iret.ne.0) then
        ierr=7
        return
      endif
      !
      !   unpack each value into array ipdstmpl from the
      !   the appropriate number of octets, which are specified in
      !   corresponding entries in array mappds.
      !
      do i=1,mapdrslen
        nbits=iabs(mapdrs(i))*8
        if ( mapdrs(i).ge.0 ) then
          call gbyte(cgrib,idrstmpl(i),iofst,nbits)
        else
          call gbyte(cgrib,isign,iofst,1)
          call gbyte(cgrib,idrstmpl(i),iofst+1,nbits-1)
          if (isign.eq.1) idrstmpl(i)=-idrstmpl(i)
        endif
        iofst=iofst+nbits
      enddo
      !
      !   check to see if the data representation template needs to be
      !   extended.
      !   the number of values in a specific template may vary
      !   depending on data specified in the "static" part of the
      !   template.
      !
      if ( needext ) then
        call extdrstemplate(idrsnum,idrstmpl,newmapdrslen,mapdrs)
        !   unpack the rest of the data representation template
        do i=mapdrslen+1,newmapdrslen
          nbits=iabs(mapdrs(i))*8
          if ( mapdrs(i).ge.0 ) then
            call gbyte(cgrib,idrstmpl(i),iofst,nbits)
          else
            call gbyte(cgrib,isign,iofst,1)
            call gbyte(cgrib,idrstmpl(i),iofst+1,nbits-1)
            if (isign.eq.1) idrstmpl(i)=-idrstmpl(i)
          endif
          iofst=iofst+nbits
        enddo
        mapdrslen=newmapdrslen
      endif
      if( allocated(mapdrs) ) deallocate(mapdrs)
      return    ! end of section 5 processing
      end


      subroutine unpack6(cgrib,lcgrib,iofst,ngpts,ibmap,bmap,ierr)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    unpack6 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-26
!
! abstract: this subroutine unpacks section 6 (bit-map section)
!   starting at octet 6 of that section.  
!
! program history log:
! 2000-05-26  gilbert
!
! usage:    call unpack6(cgrib,lcgrib,iofst,ngpts,ibmap,bmap,ierr)
!   input argument list:
!     cgrib    - character array that contains the grib2 message
!     lcgrib   - length (in bytes) of grib message array cgrib.
!     iofst    - bit offset of the beginning of section 6.
!     ngpts    - number of grid points specified in the bit-map
!
!   output argument list:      
!     iofst    - bit offset at the end of section 6, returned.
!     ibmap    - bitmap indicator ( see code table 6.0 )
!                0 = bitmap applies and is included in section 6.
!                1-253 = predefined bitmap applies
!                254 = previously defined bitmap applies to this field
!                255 = bit map does not apply to this product.
!     bmap()   - logical*1 array containing decoded bitmap. ( if ibmap=0 )
!     ierr     - error return code.
!                0 = no error
!                4 = unrecognized pre-defined bit-map.
!
! remarks: none
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$

      character(len=1),intent(in) :: cgrib(lcgrib)
      integer,intent(in) :: lcgrib,ngpts
      integer,intent(inout) :: iofst
      integer,intent(out) :: ibmap
      integer,intent(out) :: ierr
      logical*1,intent(out) :: bmap(ngpts)

      integer :: intbmap(ngpts)

      ierr=0

      iofst=iofst+32    ! skip length of section
      iofst=iofst+8     ! skip section number

      call gbyte(cgrib,ibmap,iofst,8)    ! get bit-map indicator
      iofst=iofst+8

      if (ibmap.eq.0) then               ! unpack bitmap
        call gbytes(cgrib,intbmap,iofst,1,0,ngpts)
        iofst=iofst+ngpts
        do j=1,ngpts
          bmap(j)=.true.
          if (intbmap(j).eq.0) bmap(j)=.false.
        enddo
      elseif (ibmap.eq.254) then               ! use previous bitmap
        return
      elseif (ibmap.eq.255) then               ! no bitmap in message
        bmap(1:ngpts)=.true.
      else
        print *,'unpack6: predefined bitmap ',ibmap,' not recognized.'
        ierr=4
      endif
      
      return    ! end of section 6 processing
      end

