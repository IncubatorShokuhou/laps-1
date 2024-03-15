      subroutine gb_info(cgrib,lcgrib,listsec0,listsec1,
     &                    numfields,numlocal,maxlocal,ierr)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    gb_info 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-25
!
! abstract: this subroutine searches through a grib2 message and
!   returns the number of gridded fields found in the message and
!   the number (and maximum size) of local use sections.
!   also various checks  are performed
!   to see if the message is a valid grib2 message.
!
! program history log:
! 2000-05-25  gilbert
!
! usage:    call gb_info(cgrib,lcgrib,listsec0,listsec1,
!     &                    numfields,numlocal,maxlocal,ierr)
!   input argument list:
!     cgrib    - character array that contains the grib2 message
!     lcgrib   - length (in bytes) of grib message in array cgrib.
!
!   output argument list:      
!     listsec0 - contains information decoded from grib indicator section 0.
!                must be dimensioned >= 2.
!                listsec0(1)=discipline-grib master table number
!                            (see code table 0.0)
!                listsec0(2)=grib edition number (currently 2)
!                listsec0(3)=length of grib message
!     listsec1 - contains information read from grib identification section 1.
!                must be dimensioned >= 13.
!                listsec1(1)=id of orginating centre (common code table c-1)
!                listsec1(2)=id of orginating sub-centre (local table)
!                listsec1(3)=grib master tables version number (code table 1.0)
!                listsec1(4)=grib local tables version number 
!                listsec1(5)=significance of reference time (code table 1.1)
!                listsec1(6)=reference time - year (4 digits)
!                listsec1(7)=reference time - month
!                listsec1(8)=reference time - day
!                listsec1(9)=reference time - hour
!                listsec1(10)=reference time - minute
!                listsec1(11)=reference time - second
!                listsec1(12)=production status of data (code table 1.2)
!                listsec1(13)=type of processed data (code table 1.3)
!     numfields- the number of gridded fieldse found in the grib message.
!     numlocal - the number of local use sections ( section 2 ) found in 
!                the grib message.
!     maxlocal-  the size of the largest local use section ( section 2 ).
!                can be used to ensure that the return array passed
!                to subroutine getlocal is dimensioned large enough.
!     ierr     - error return code.
!                0 = no error
!                1 = beginning characters "grib" not found.
!                2 = grib message is not edition 2.
!                3 = could not find section 1, where expected.
!                4 = end string "7777" found, but not where expected.
!                5 = end string "7777" not found at end of message.
!                6 = invalid section number found.
!
! remarks: none
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$

      character(len=1),intent(in) :: cgrib(lcgrib)
      integer,intent(in) :: lcgrib
      integer,intent(out) :: listsec0(3),listsec1(13)
      integer,intent(out) :: numlocal,numfields,maxlocal,ierr
      
      character(len=4),parameter :: grib='grib',c7777='7777'
      character(len=4) :: ctemp
      integer,parameter :: zero=0,one=1
      integer,parameter :: mapsec1len=13
      integer,parameter :: 
     &        mapsec1(mapsec1len)=(/ 2,2,1,1,1,2,1,1,1,1,1,1,1 /)
      integer iofst,ibeg,istart

      ierr=0
      numlocal=0
      numfields=0
      maxlocal=0
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
        print *,'gb_info:  beginning characters grib not found.'
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
      listsec0(3)=lengrib
      lensec0=16
      ipos=istart+lensec0
!
!  currently handles only grib edition 2.
!  
      if (listsec0(2).ne.2) then
        print *,'gb_info: can only decode grib edition 2.'
        ierr=2
        return
      endif
!
!  unpack section 1 - identification section
!
      call gbyte(cgrib,lensec1,iofst,32)        ! length of section 1
      iofst=iofst+32
      call gbyte(cgrib,isecnum,iofst,8)         ! section number ( 1 )
      iofst=iofst+8
      if (isecnum.ne.1) then
        print *,'gb_info: could not find section 1.'
        ierr=3
        return
      endif
      !
      !   unpack each input value in array listsec1 into the
      !   the appropriate number of octets, which are specified in
      !   corresponding entries in array mapsec1.
      !
      do i=1,mapsec1len
        nbits=mapsec1(i)*8
        call gbyte(cgrib,listsec1(i),iofst,nbits)
        iofst=iofst+nbits
      enddo
      ipos=ipos+lensec1
!
!  loop through the remaining sections to see if they are valid.
!  also count the number of times section 2
!  and section 4 appear.
!
      do
        ctemp=cgrib(ipos)//cgrib(ipos+1)//cgrib(ipos+2)//cgrib(ipos+3)
        if (ctemp.eq.c7777 ) then
          ipos=ipos+4
          if (ipos.ne.(istart+lengrib)) then
            print *,'gb_info: "7777" found, but not where expected.'
            ierr=4
            return
          endif
          exit
        endif
        iofst=(ipos-1)*8
        call gbyte(cgrib,lensec,iofst,32)        ! get length of section
        iofst=iofst+32
        call gbyte(cgrib,isecnum,iofst,8)         ! get section number
        iofst=iofst+8
        ipos=ipos+lensec                 ! update beginning of section pointer
        if (ipos.gt.(istart+lengrib)) then
          print *,'gb_info: "7777"  not found at end of grib message.'
          ierr=5
          return
        endif
        if ( isecnum.ge.2.and.isecnum.le.7 ) then
           if (isecnum.eq.2) then     ! local section 2
              !   increment counter for total number of local sections found
              numlocal=numlocal+1
              lenposs=lensec-5
              if ( lenposs.gt.maxlocal ) maxlocal=lenposs
           elseif (isecnum.eq.4) then
              !   increment counter for total number of fields found
              numfields=numfields+1
           endif
        else
           print *,'gb_info: invalid section number found in grib',
     &             ' message: ',isecnum
           ierr=6
           return
        endif
        
      enddo

      return
      end

