      subroutine gf_unpack1(cgrib,lcgrib,iofst,ids,idslen,ierr)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    gf_unpack1 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-26
!
! abstract: this subroutine unpacks section 1 (identification section)
!   starting at octet 6 of that section.  
!
! program history log:
! 2000-05-26  gilbert
! 2002-01-24  gilbert  - changed to dynamically allocate arrays
!                        and to pass pointers to those arrays through
!                        the argument list.
!
! usage:    call gf_unpack1(cgrib,lcgrib,iofst,ids,idslen,ierr)
!   input argument list:
!     cgrib    - character array containing section 1 of the grib2 message
!     lcgrib   - length (in bytes) of grib message array cgrib.
!     iofst    - bit offset of the beginning of section 1.
!
!   output argument list:      
!     iofst    - bit offset at the end of section 1, returned.
!     ids      - pointer to integer array containing information read from 
!                section 1, the identification section.
!            ids(1)  = identification of originating centre
!                                 ( see common code table c-1 )
!            ids(2)  = identification of originating sub-centre
!            ids(3)  = grib master tables version number
!                                 ( see code table 1.0 )
!            ids(4)  = grib local tables version number
!                                 ( see code table 1.1 )
!            ids(5)  = significance of reference time (code table 1.2)
!            ids(6)  = year ( 4 digits )
!            ids(7)  = month
!            ids(8)  = day
!            ids(9)  = hour
!            ids(10)  = minute
!            ids(11)  = second
!            ids(12)  = production status of processed data
!                                 ( see code table 1.3 )
!            ids(13)  = type of processed data ( see code table 1.4 )
!     idslen   - number of elements in ids().
!     ierr     - error return code.
!                0 = no error
!                6 = memory allocation error
!
! remarks: 
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$

      character(len=1),intent(in) :: cgrib(lcgrib)
      integer,intent(in) :: lcgrib
      integer,intent(inout) :: iofst
      integer,pointer,dimension(:) :: ids
      integer,intent(out) :: ierr,idslen

      integer,dimension(:) :: mapid(13)

      data mapid /2,2,1,1,1,2,1,1,1,1,1,1,1/

      ierr=0
      idslen=13
      nullify(ids)

      call gbyte(cgrib,lensec,iofst,32)        ! get length of section
      iofst=iofst+32
      iofst=iofst+8     ! skip section number
      !
      !   unpack each value into array ids from the
      !   the appropriate number of octets, which are specified in
      !   corresponding entries in array mapid.
      !
      istat=0
      allocate(ids(idslen),stat=istat)
      if (istat.ne.0) then
         ierr=6
         nullify(ids)
         return
      endif
      
      do i=1,idslen
        nbits=mapid(i)*8
        call gbyte(cgrib,ids(i),iofst,nbits)
        iofst=iofst+nbits
      enddo
      
      return    ! end of section 1 processing
      end
