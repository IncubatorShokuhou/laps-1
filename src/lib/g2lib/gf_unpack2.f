      subroutine gf_unpack2(cgrib,lcgrib,iofst,lencsec2,csec2,ierr)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    gf_unpack2 
!   prgmmr: gilbert         org: w/np11    date: 2002-04-09
!
! abstract: this subroutine unpacks section 2 (local use section)
!           as defined in grib edition 2.
!
! program history log:
! 2002-04-09  gilbert
!
! usage:    call gf_unpack2(cgrib,lcgrib,iofst,lencsec2,csec2,ierr)
!   input argument list:
!     cgrib    - character array containing section 2 of the grib2 message
!     lcgrib   - length (in bytes) of grib message array cgrib.
!     iofst    - bit offset of the beginning of section 2.
!
!   output argument list:      
!     iofst    - bit offset at the end of section 2, returned.
!     lencsec2 - length (in octets) of local use data
!     csec2()  - pointer to a character*1 array containing local use data
!     ierr     - error return code.
!                0 = no error
!                2 = array passed is not section 2
!                6 = memory allocation error
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
      integer,intent(inout) :: iofst
      integer,intent(out) :: lencsec2
      integer,intent(out) :: ierr
      character(len=1),pointer,dimension(:) :: csec2

      ierr=0
      lencsec2=0
      nullify(csec2)

      call gbyte(cgrib,lensec,iofst,32)        ! get length of section
      iofst=iofst+32    
      lencsec2=lensec-5
      call gbyte(cgrib,isecnum,iofst,8)         ! get section number
      iofst=iofst+8     
      ipos=(iofst/8)+1

      if ( isecnum.ne.2 ) then
         ierr=6
         print *,'gf_unpack2: not section 2 data. '
         return
      endif

      allocate(csec2(lencsec2),stat=istat)
      if (istat.ne.0) then
         ierr=6
         nullify(csec2)
         return
      endif
      
      csec2(1:lencsec2)=cgrib(ipos:ipos+lencsec2-1)
      iofst=iofst+(lencsec2*8)

      return    ! end of section 2 processing
      end

