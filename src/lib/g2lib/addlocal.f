      subroutine addlocal(cgrib,lcgrib,csec2,lcsec2,ierr)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    addlocal 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-01
!
! abstract: this subroutine adds a local use section (section 2) to 
!   a grib2 message.
!   this routine is used with routines "gribcreate", "addgrid", "addfield",
!   and "gribend" to create a complete grib2 message.  subroutine
!   gribcreate must be called first to initialize a new grib2 message.
!
! program history log:
! 2000-05-01  gilbert
!
! usage:    call addlocal(cgrib,lcgrib,csec2,lcsec2,ierr)
!   input argument list:
!     cgrib    - character array to contain the grib2 message
!     lcgrib   - maximum length (bytes) of array cgrib.
!     csec2    - character array containing information to be added to
!                section 2.
!     lcsec2   - number of bytes of character array csec2 to be added to
!                section 2.
!
!   output argument list:      
!     cgrib    - character array to contain the grib2 message
!     ierr     - error return code.
!                0 = no error
!                1 = grib message was not initialized.  need to call
!                    routine gribcreate first.
!                2 = grib message already complete.  cannot add new section.
!                3 = sum of section byte counts doesn't add to total byte count.
!                4 = previous section was not 1 or 7.
!
! remarks: note that the local use section ( section 2 ) can only follow
!          section 1 or section 7 in a grib2 message.
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$

      character(len=1),intent(inout) :: cgrib(lcgrib)
      character(len=1),intent(in) :: csec2(lcsec2)
      integer,intent(in) :: lcgrib,lcsec2
      integer,intent(out) :: ierr
      
      character(len=4),parameter :: grib='grib',c7777='7777'
      character(len=4):: ctemp
      integer,parameter :: two=2
      integer lensec2,iofst,ibeg,lencurr,len
 
      ierr=0
!
!  check to see if beginning of grib message exists
!
      ctemp=cgrib(1)//cgrib(2)//cgrib(3)//cgrib(4)
      if ( ctemp.ne.grib ) then
        print *,'addlocal: grib not found in given message.'
        print *,'addlocal: call to routine gribcreate required',
     &          ' to initialize grib messge.'
        ierr=1
        return
      endif
!
!  get current length of grib message
!  
      call gbyte(cgrib,lencurr,96,32)
!
!  check to see if grib message is already complete
!  
      ctemp=cgrib(lencurr-3)//cgrib(lencurr-2)//cgrib(lencurr-1)
     &      //cgrib(lencurr)
      if ( ctemp.eq.c7777 ) then
        print *,'addlocal: grib message already complete.  cannot',
     &          ' add new section.'
        ierr=2
        return
      endif
!
!  loop through all current sections of the grib message to
!  find the last section number.
!
      len=16    ! length of section 0
      do 
      !    get section number and length of next section
        iofst=len*8
        call gbyte(cgrib,ilen,iofst,32)
        iofst=iofst+32
        call gbyte(cgrib,isecnum,iofst,8)
        len=len+ilen
      !    exit loop if last section reached
        if ( len.eq.lencurr ) exit
      !    if byte count for each section doesn't match current
      !    total length, then there is a problem.
        if ( len.gt.lencurr ) then
          print *,'addlocal: section byte counts don''t add to total.'
          print *,'addlocal: sum of section byte counts = ',len
          print *,'addlocal: total byte count in section 0 = ',lencurr
          ierr=3
          return
        endif
      enddo
!
!  section 2 can only be added after sections 1 and 7.
!
      if ( (isecnum.ne.1) .and. (isecnum.ne.7) ) then
        print *,'addlocal: section 2 can only be added after section',
     &          ' 1 or section 7.'
        print *,'addlocal: section ',isecnum,' was the last found in',
     &          ' given grib message.'
        ierr=4
        return
      endif
!
!  add section 2  - local use section
!
      ibeg=lencurr*8        !   calculate offset for beginning of section 2
      iofst=ibeg+32         !   leave space for length of section
      call sbyte(cgrib,two,iofst,8)     ! store section number ( 2 )
      istart=lencurr+5
      cgrib(istart+1:istart+lcsec2)=csec2(1:lcsec2)
      !
      !   calculate length of section 2 and store it in octets
      !   1-4 of section 2.
      !
      lensec2=lcsec2+5      ! bytes
      call sbyte(cgrib,lensec2,ibeg,32)

!
!  update current byte total of message in section 0
!
      call sbyte(cgrib,lencurr+lensec2,96,32)

      return
      end

