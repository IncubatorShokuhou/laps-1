      subroutine specunpack(cpack,len,idrstmpl,ndpts,jj,kk,mm,fld)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    specunpack
!   prgmmr: gilbert          org: w/np11    date: 2002-12-19
!
! abstract: this subroutine unpacks a spectral data field that was packed 
!   using the complex packing algorithm for spherical harmonic data as 
!   defined in the grib2 documention,
!   using info from the grib2 data representation template 5.51.
!
! program history log:
! 2002-12-19  gilbert
!
! usage:    call specunpack(cpack,len,idrstmpl,ndpts,jj,kk,mm,fld)
!   input argument list:
!     cpack    - the packed data field (character*1 array)
!     len      - length of packed field cpack().
!     idrstmpl - contains the array of values for data representation
!                template 5.51
!     ndpts    - the number of data values to unpack
!     jj       - j - pentagonal resolution parameter
!     kk       - k - pentagonal resolution parameter
!     mm       - m - pentagonal resolution parameter
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
      integer,intent(in) :: ndpts,len,jj,kk,mm
      integer,intent(in) :: idrstmpl(*)
      real,intent(out) :: fld(ndpts)

      integer :: ifld(ndpts),ts
      integer(4) :: ieee
      real :: ref,bscale,dscale,unpk(ndpts)
      real,allocatable :: pscale(:)

      ieee = idrstmpl(1)
      call rdieee(ieee,ref,1)
      bscale = 2.0**real(idrstmpl(2))
      dscale = 10.0**real(-idrstmpl(3))
      nbits = idrstmpl(4)
      js=idrstmpl(6)
      ks=idrstmpl(7)
      ms=idrstmpl(8)
      ts=idrstmpl(9)

      if (idrstmpl(10).eq.1) then           ! unpacked floats are 32-bit ieee
         !call gbytes(cpack,ifld,0,32,0,ts)
         call rdieee(cpack,unpk,ts)          ! read ieee unpacked floats
         iofst=32*ts
         call gbytes(cpack,ifld,iofst,nbits,0,ndpts-ts)  ! unpack scaled data
!
!   calculate laplacian scaling factors for each possible wave number.
!
         allocate(pscale(jj+mm))
         tscale=real(idrstmpl(5))*1e-6
         do n=js,jj+mm
            pscale(n)=real(n*(n+1))**(-tscale)
         enddo
!
!   assemble spectral coeffs back to original order.
!
         inc=1
         incu=1
         incp=1
         do m=0,mm
            nm=jj      ! triangular or trapezoidal
            if ( kk .eq. jj+mm ) nm=jj+m          ! rhombodial
            ns=js      ! triangular or trapezoidal
            if ( ks .eq. js+ms ) ns=js+m          ! rhombodial
            do n=m,nm
               if (n.le.ns .and. m.le.ms) then    ! grab unpacked value
                  fld(inc)=unpk(incu)         ! real part
                  fld(inc+1)=unpk(incu+1)     ! imaginary part
                  inc=inc+2
                  incu=incu+2
               else                         ! calc coeff from packed value
                  fld(inc)=((real(ifld(incp))*bscale)+ref)*
     &                      dscale*pscale(n)           ! real part
                  fld(inc+1)=((real(ifld(incp+1))*bscale)+ref)*
     &                      dscale*pscale(n)           ! imaginary part
                  inc=inc+2
                  incp=incp+2
               endif
            enddo
         enddo

         deallocate(pscale)

      else
         print *,'specunpack: cannot handle 64 or 128-bit floats.'
         fld=0.0
         return
      endif

      return
      end
