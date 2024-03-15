      subroutine simunpack(cpack,len,idrstmpl,ndpts,fld)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    simunpack
!   prgmmr: gilbert          org: w/np11    date: 2000-06-21
!
! abstract: this subroutine unpacks a data field that was packed using a 
!   simple packing algorithm as defined in the grib2 documention,
!   using info from the grib2 data representation template 5.0.
!
! program history log:
! 2000-06-21  gilbert
!
! usage:    call simunpack(cpack,len,idrstmpl,ndpts,fld)
!   input argument list:
!     cpack    - the packed data field (character*1 array)
!     len      - length of packed field cpack().
!     idrstmpl - contains the array of values for data representation
!                template 5.0
!     ndpts    - the number of data values to unpack
!
!   output argument list:
!     fld()    - contains the unpacked data values
!
! remarks: none
!
! attributes:
!   language: xl fortran 90
!   machine:  ibm sp
!
!$$$

      character(len=1),intent(in) :: cpack(len)
      integer,intent(in) :: ndpts,len
      integer,intent(in) :: idrstmpl(*)
      real,intent(out) :: fld(ndpts)

      integer :: ifld(ndpts)
      integer(4) :: ieee
      real :: ref,bscale,dscale

      ieee = idrstmpl(1)
      call rdieee(ieee,ref,1)
      bscale = 2.0**real(idrstmpl(2))
      dscale = 10.0**real(-idrstmpl(3))
      nbits = idrstmpl(4)
      itype = idrstmpl(5)
!
!  if nbits equals 0, we have a constant field where the reference value
!  is the data value at each gridpoint
!
      if (nbits.ne.0) then
         call gbytes(cpack,ifld,0,nbits,0,ndpts)
         do j=1,ndpts
           fld(j)=((real(ifld(j))*bscale)+ref)*dscale
         enddo
      else
         do j=1,ndpts
           fld(j)=ref
         enddo
      endif


      return
      end
