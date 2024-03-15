      subroutine gf_unpack6(cgrib,lcgrib,iofst,ngpts,ibmap,bmap,ierr)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    gf_unpack6 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-26
!
! abstract: this subroutine unpacks section 6 (bit-map section)
!   starting at octet 6 of that section.  
!
! program history log:
! 2000-05-26  gilbert
! 2002-01-24  gilbert  - changed to dynamically allocate arrays
!                        and to pass pointers to those arrays through
!                        the argument list.
!
! usage:    call gf_unpack6(cgrib,lcgrib,iofst,ngpts,ibmap,bmap,ierr)
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
!     bmap()   - pointer to a logical*1 array containing decoded bitmap. 
!                ( if ibmap=0 )
!     ierr     - error return code.
!                0 = no error
!                4 = unrecognized pre-defined bit-map.
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
      integer,intent(in) :: lcgrib,ngpts
      integer,intent(inout) :: iofst
      integer,intent(out) :: ibmap
      integer,intent(out) :: ierr
      logical*1,pointer,dimension(:) :: bmap

      integer :: intbmap(ngpts)

      ierr=0
      nullify(bmap)

      iofst=iofst+32    ! skip length of section
      iofst=iofst+8     ! skip section number

      call gbyte(cgrib,ibmap,iofst,8)    ! get bit-map indicator
      iofst=iofst+8

      if (ibmap.eq.0) then               ! unpack bitmap
         istat=0
         if (ngpts.gt.0) allocate(bmap(ngpts),stat=istat)
         if (istat.ne.0) then
            ierr=6
            nullify(bmap)
            return
         endif
         call gbytes(cgrib,intbmap,iofst,1,0,ngpts)
         iofst=iofst+ngpts
         do j=1,ngpts
           bmap(j)=.true.
           if (intbmap(j).eq.0) bmap(j)=.false.
         enddo
!      elseif (ibmap.eq.254) then               ! use previous bitmap
!        return
!      elseif (ibmap.eq.255) then               ! no bitmap in message
!        bmap(1:ngpts)=.true.
!      else
!        print *,'gf_unpack6: predefined bitmap ',ibmap,' not recognized.'
!        ierr=4
      endif
      
      return    ! end of section 6 processing
      end

