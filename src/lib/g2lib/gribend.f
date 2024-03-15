      subroutine gribend(cgrib,lcgrib,lengrib,ierr)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    gribend 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-02
!
! abstract: this subroutine finalizes a grib message after all grids
!   and fields have been added.  it adds the end section ( "7777" )
!   to the end of the grib message and calculates the length and stores
!   it in the appropriate place in section 0.
!   this routine is used with routines "gribcreate", "addlocal", "addgrid",
!   and "addfield" to create a complete grib2 message.  subroutine
!   gribcreate must be called first to initialize a new grib2 message.
!
! program history log:
! 2000-05-02  gilbert
!
! usage:    call gribend(cgrib,lcgrib,lengrib,ierr)
!   input argument list:
!     cgrib    - character array to contain the grib2 message
!     lcgrib   - maximum length (bytes) of array cgrib.
!
!   output argument list:      
!     cgrib    - character array to contain the grib2 message
!     lengrib  - length of the final grib2 message in octets (bytes)
!     ierr     - error return code.
!                0 = no error
!                1 = grib message was not initialized.  need to call
!                    routine gribcreate first.
!                2 = grib message already complete.  
!                3 = sum of section byte counts doesn't add to total byte count.
!                4 = previous section was not 7.
!
! remarks: this routine is intended for use with routines "gribcreate", 
!          "addlocal", "addgrid", and "addfield" to create a complete 
!          grib2 message.
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$

      character(len=1),intent(inout) :: cgrib(lcgrib)
      integer,intent(in) :: lcgrib
      integer,intent(out) :: lengrib,ierr
      
      character(len=4),parameter :: grib='grib',c7777='7777'
      character(len=4):: ctemp
      integer iofst,ibeg,lencurr,len
 
      ierr=0
!
!  check to see if beginning of grib message exists
!
      ctemp=cgrib(1)//cgrib(2)//cgrib(3)//cgrib(4)
      if ( ctemp.ne.grib ) then
        print *,'gribend: grib not found in given message.'
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
!      ctemp=cgrib(lencurr-3)//cgrib(lencurr-2)//cgrib(lencurr-1)
!     &      //cgrib(lencurr)
!      if ( ctemp.eq.c7777 ) then
!        print *,'gribend: grib message already complete.'
!        ierr=2
!        return
!      endif
!
!  loop through all current sections of the grib message to
!  find the last section number.
!
      len=16    ! length of section 0
      do 
      !    get number and length of next section
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
          print *,'gribend: section byte counts don''t add to total.'
          print *,'gribend: sum of section byte counts = ',len
          print *,'gribend: total byte count in section 0 = ',lencurr
          ierr=3
          return
        endif
      enddo
!
!  can only add end section (section 8) after section 7.
!
      if ( isecnum.ne.7 ) then
        print *,'gribend: section 8 can only be added after section 7.'
        print *,'gribend: section ',isecnum,' was the last found in',
     &          ' given grib message.'
        ierr=4
        return
      endif
!
!  add section 8  - end section
!
      cgrib(lencurr+1:lencurr+4)=c7777

!
!  update current byte total of message in section 0
!
      lengrib=lencurr+4
      call sbyte(cgrib,lengrib,96,32)

      return
      end




