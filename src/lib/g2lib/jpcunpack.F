      subroutine jpcunpack(cpack,len,idrstmpl,ndpts,fld)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    jpcunpack
!   prgmmr: gilbert          org: w/np11    date: 2002-12-17
!
! abstract: this subroutine unpacks a data field that was packed into a 
!   jpeg2000 code stream
!   using info from the grib2 data representation template 5.40 or 5.40000.
!
! program history log:
! 2002-12-17  gilbert
!
! usage:    call jpcunpack(cpack,len,idrstmpl,ndpts,fld)
!   input argument list:
!     cpack    - the packed data field (character*1 array)
!     len      - length of packed field cpack().
!     idrstmpl - contains the array of values for data representation
!                template 5.40 or 5.40000
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
      integer :: dec_jpeg2000
#ifdef use_jpeg2000

      ieee = idrstmpl(1)
      call rdieee(ieee,ref,1)
      bscale = 2.0**real(idrstmpl(2))
      dscale = 10.0**real(-idrstmpl(3))
      nbits = idrstmpl(4)
!
!  if nbits equals 0, we have a constant field where the reference value
!  is the data value at each gridpoint
!
      if (nbits.ne.0) then
         call gbytes(cpack,ifld,0,nbits,0,ndpts)
         iret=dec_jpeg2000(cpack,len,ifld)
!        write(6,*)
!    1     ' error in jpcunpack.f: attempt to use unavailable routines'
!        stop
         do j=1,ndpts
           fld(j)=((real(ifld(j))*bscale)+ref)*dscale
         enddo
      else
         do j=1,ndpts
           fld(j)=ref
         enddo
      endif


#endif /* use_jpeg2000 */
      return
      end
