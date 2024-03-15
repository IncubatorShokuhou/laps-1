      subroutine getlocal(cgrib,lcgrib,localnum,csec2,lcsec2,ierr)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    getlocal 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-25
!
! abstract: this subroutine returns the contents of section 2 ( local 
!   use section ) from a grib2 message.  since there can be multiple
!   occurrences of section 2 within a grib message, the calling routine
!   indicates which occurrence is being requested with the localnum argument.
!
! program history log:
! 2000-05-25  gilbert
!
! usage:    call getlocal(cgrib,lcgrib,localnum,csec2,lcsec2,ierr)
!   input argument list:
!     cgrib    - character array that contains the grib2 message
!     lcgrib   - length (in bytes) of grib message in array cgrib.
!     localnum - the nth occurrence of section 2 requested.
!
!   output argument list:      
!     csec2    - character array containing information read from 
!                section 2.
!                the dimension of this array can be obtained in advance
!                from argument maxlocal, which is returned from subroutine 
!                gb_info.
!     lcsec2   - number of bytes of character array csec2 read from
!                section 2.
!     ierr     - error return code.
!                0 = no error
!                1 = beginning characters "grib" not found.
!                2 = grib message is not edition 2.
!                3 = the section 2 request number was not positive.
!                4 = end string "7777" found, but not where expected.
!                5 = end string "7777" not found at end of message.
!                6 = grib message did not contain the requested number of
!                    local use sections.
!
! remarks: note that subroutine gb_info can be used to first determine
!          how many local use sections exist in a given grib message.
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$

      character(len=1),intent(in) :: cgrib(lcgrib)
      integer,intent(in) :: lcgrib,localnum
      character(len=1),intent(out) :: csec2(*)
      integer,intent(out) :: lcsec2,ierr
      
      character(len=4),parameter :: grib='grib',c7777='7777'
      character(len=4) :: ctemp
      integer :: listsec0(2)
      integer iofst,ibeg,istart,numlocal

      ierr=0
      numlocal=0
!
!  check for valid request number
!  
      if (localnum.le.0) then
        print *,'getlocal: request for local section must be positive.'
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
        print *,'getlocal:  beginning characters grib not found.'
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
        print *,'getlocal: can only decode grib edition 2.'
        ierr=2
        return
      endif
!
!  loop through the remaining sections keeping track of the 
!  length of each.  also check to see that if the current occurrence
!  of section 2 is the same as the one requested.
!
      do
        !    check to see if we are at end of grib message
        ctemp=cgrib(ipos)//cgrib(ipos+1)//cgrib(ipos+2)//cgrib(ipos+3)
        if (ctemp.eq.c7777 ) then
          ipos=ipos+4
          !    if end of grib message not where expected, issue error
          if (ipos.ne.(istart+lengrib)) then
            print *,'getlocal: "7777" found, but not where expected.'
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
        !   if found the requested occurrence of section 2,
        !   return the section contents.
        if (isecnum.eq.2) then
          numlocal=numlocal+1
          if (numlocal.eq.localnum) then
            lcsec2=lensec-5
            csec2(1:lcsec2)=cgrib(ipos+5:ipos+lensec-1)
            return
          endif
        endif
        !   check to see if we read pass the end of the grib
        !   message and missed the terminator string '7777'.
        ipos=ipos+lensec                 ! update beginning of section pointer
        if (ipos.gt.(istart+lengrib)) then
          print *,'getlocal: "7777"  not found at end of grib message.'
          ierr=5
          return
        endif
        
      enddo

!
!  if exited from above loop, the end of the grib message was reached
!  before the requested occurrence of section 2 was found.
!
      print *,'getlocal: grib message contained ',numlocal,
     &        ' local sections.'
      print *,'getlocal: the request was for the ',localnum,
     &        ' occurrence.'
      ierr=6

      return
      end







