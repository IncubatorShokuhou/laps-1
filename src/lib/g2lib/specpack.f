      subroutine specpack(fld,ndpts,jj,kk,mm,idrstmpl,cpack,lcpack)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    specpack
!   prgmmr: gilbert          org: w/np11    date: 2002-12-19
!
! abstract: this subroutine packs a spectral data field using the complex
!   packing algorithm for spherical harmonic data as 
!   defined in the grib2 data representation template 5.51.
!
! program history log:
! 2002-12-19  gilbert
!
! usage:    call specpack(fld,ndpts,jj,kk,mm,idrstmpl,cpack,lcpack)
!   input argument list:
!     fld()    - contains the packed data values
!     ndpts    - the number of data values to pack
!     jj       - j - pentagonal resolution parameter
!     kk       - k - pentagonal resolution parameter
!     mm       - m - pentagonal resolution parameter
!     idrstmpl - contains the array of values for data representation
!                template 5.51
!
!   output argument list:
!     cpack    - the packed data field (character*1 array)
!     lcpack   - length of packed field cpack().
!
! remarks: none
!
! attributes:
!   language: xl fortran 90
!   machine:  ibm sp
!
!$$$

      real,intent(in) :: fld(ndpts)
      integer,intent(in) :: ndpts,jj,kk,mm
      integer,intent(inout) :: idrstmpl(*)
      character(len=1),intent(out) :: cpack(*)
      integer,intent(out) :: lcpack

      integer :: ifld(ndpts),ts,tmplsim(5)
      real :: bscale,dscale,unpk(ndpts),tfld(ndpts)
      real,allocatable :: pscale(:)

      bscale = 2.0**real(-idrstmpl(2))
      dscale = 10.0**real(idrstmpl(3))
      nbits = idrstmpl(4)
      js=idrstmpl(6)
      ks=idrstmpl(7)
      ms=idrstmpl(8)
      ts=idrstmpl(9)

!
!   calculate laplacian scaling factors for each possible wave number.
!
      allocate(pscale(jj+mm))
      tscale=real(idrstmpl(5))*1e-6
      do n=js,jj+mm
         pscale(n)=real(n*(n+1))**(tscale)
      enddo
!
!   separate spectral coeffs into two lists; one to contain unpacked
!   values within the sub-spectrum js, ks, ms, and the other with values 
!   outside of the sub-spectrum to be packed.
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
            if (n.le.ns .and. m.le.ms) then    ! save unpacked value
               unpk(incu)=fld(inc)         ! real part
               unpk(incu+1)=fld(inc+1)     ! imaginary part
               inc=inc+2
               incu=incu+2
            else                         ! save value to be packed and scale
                                         ! laplacian scale factor
               tfld(incp)=fld(inc)*pscale(n)         ! real part
               tfld(incp+1)=fld(inc+1)*pscale(n)     ! imaginary part
               inc=inc+2
               incp=incp+2
            endif
         enddo
      enddo

      deallocate(pscale)

      incu=incu-1
      if (incu .ne. ts) then
         print *,'specpack: incorrect number of unpacked values ',
     &           'given:',ts     
         print *,'specpack: resetting idrstmpl(9) to ',incu
         ts=incu
      endif
!
!  add unpacked values to the packed data array in 32-bit ieee format
!
      call mkieee(unpk,cpack,ts)
      ipos=4*ts
!
!  scale and pack the rest of the coefficients
! 
      tmplsim(2)=idrstmpl(2)
      tmplsim(3)=idrstmpl(3)
      tmplsim(4)=idrstmpl(4)
      call simpack(tfld,ndpts-ts,tmplsim,cpack(ipos+1),lcpack)
      lcpack=lcpack+ipos
!
!  fill in template 5.51
!
      idrstmpl(1)=tmplsim(1)
      idrstmpl(2)=tmplsim(2)
      idrstmpl(3)=tmplsim(3)
      idrstmpl(4)=tmplsim(4)
      idrstmpl(9)=ts
      idrstmpl(10)=1         ! unpacked spectral data is 32-bit ieee

      return
      end
