      subroutine gettemplates(cgrib,lcgrib,ifldnum,igds,igdstmpl,
     &                    igdslen,ideflist,idefnum,ipdsnum,ipdstmpl,
     &                    ipdslen,coordlist,numcoord,ierr)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    gettemplates 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-26
!
! abstract: this subroutine returns the grid definition, and
!   product definition for a given data
!   field.  since there can be multiple data fields packed into a grib2
!   message, the calling routine indicates which field is being requested
!   with the ifldnum argument.
!
! program history log:
! 2000-05-26  gilbert
!
! usage:    call gettemplates(cgrib,lcgrib,ifldnum,igds,igdstmpl,igdslen,
!    &                    ideflist,idefnum,ipdsnum,ipdstmpl,ipdslen,
!    &                    coordlist,numcoord,ierr)
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
!     ierr     - error return code.
!                0 = no error
!                1 = beginning characters "grib" not found.
!                2 = grib message is not edition 2.
!                3 = the data field request number was not positive.
!                4 = end string "7777" found, but not where expected.
!                6 = grib message did not contain the requested number of
!                    data fields.
!                7 = end string "7777" not found at end of message.
!               10 = error unpacking section 3.
!               11 = error unpacking section 4.
!
! remarks: note that subroutine gribinfo can be used to first determine
!          how many data fields exist in the given grib message.
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
      integer,intent(out) :: idefnum,numcoord
      integer,intent(out) :: ierr
      real,intent(out) :: coordlist(*)
      
      character(len=4),parameter :: grib='grib',c7777='7777'
      character(len=4) :: ctemp
      integer:: listsec0(2)
      integer iofst,ibeg,istart
      logical have3,have4

      have3=.false.
      have4=.false.
      ierr=0
      numfld=0
!
!  check for valid request number
!  
      if (ifldnum.le.0) then
        print *,'gettemplates: request for field number must be ',
     &          'positive.'
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
        print *,'gettemplates:  beginning characters grib not found.'
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
        print *,'gettemplates: can only decode grib edition 2.'
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
            print *,'gettemplates: "7777" found, but not where ',
     &              'expected.'
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
        !   check to see if we read pass the end of the grib
        !   message and missed the terminator string '7777'.
        !
        ipos=ipos+lensec                 ! update beginning of section pointer
        if (ipos.gt.(istart+lengrib)) then
          print *,'gettemplates: "7777"  not found at end of grib ',
     &            'message.'
          ierr=7
          return
        endif

        if (have3.and.have4) return
        
      enddo

!
!  if exited from above loop, the end of the grib message was reached
!  before the requested field was found.
!
      print *,'gettemplates: grib message contained ',numlocal,
     &        ' different fields.'
      print *,'gettemplates: the request was for the ',ifldnum,
     &        ' field.'
      ierr=6

      return
      end
