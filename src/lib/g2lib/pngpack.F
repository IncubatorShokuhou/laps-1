      subroutine pngpack(fld,width,height,idrstmpl,cpack,lcpack)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    pngpack
!   prgmmr: gilbert          org: w/np11    date: 2002-12-21
!
! abstract: this subroutine packs up a data field into png image format.
!   after the data field is scaled, and the reference value is subtracted out,
!   it is treated as a grayscale image and passed to a png encoder.
!   it also fills in grib2 data representation template 5.41 or 5.40010 with the
!   appropriate values.
!
! program history log:
! 2002-12-21  gilbert
!
! usage:    call pngpack(fld,width,height,idrstmpl,cpack,lcpack)
!   input argument list:
!     fld()    - contains the data values to pack
!     width    - number of points in the x direction
!     height   - number of points in the y direction
!     idrstmpl - contains the array of values for data representation
!                template 5.41 or 5.40010
!                (1) = reference value - ignored on input
!                (2) = binary scale factor
!                (3) = decimal scale factor
!                (4) = number of bits for each data value - ignored on input
!                (5) = original field type - currently ignored on input
!                      data values assumed to be reals.
!
!   output argument list: 
!     idrstmpl - contains the array of values for data representation
!                template 5.41 or 5.40010
!                (1) = reference value - set by pngpack routine.
!                (2) = binary scale factor - unchanged from input
!                (3) = decimal scale factor - unchanged from input
!                (4) = number of bits containing each grayscale pixel value
!                (5) = original field type - currently set = 0 on output.
!                      data values assumed to be reals.
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

      integer,intent(in) :: width,height
      real,intent(in) :: fld(width*height)
      character(len=1),intent(out) :: cpack(*)
      integer,intent(inout) :: idrstmpl(*)
      integer,intent(out) :: lcpack

      real(4) :: ref
      integer(4) :: iref
      integer :: ifld(width*height)
      integer,parameter :: zero=0
      integer :: enc_png
      character(len=1),allocatable :: ctemp(:)
#ifdef use_png
      
      ndpts=width*height
      bscale=2.0**real(-idrstmpl(2))
      dscale=10.0**real(idrstmpl(3))
!
!  find max and min values in the data
!
      rmax=fld(1)
      rmin=fld(1)
      do j=2,ndpts
        if (fld(j).gt.rmax) rmax=fld(j)
        if (fld(j).lt.rmin) rmin=fld(j)
      enddo
      maxdif=nint((rmax-rmin)*dscale*bscale)
!
!  if max and min values are not equal, pack up field.
!  if they are equal, we have a constant field, and the reference
!  value (rmin) is the value for each point in the field and
!  set nbits to 0.
!
      if (rmin.ne.rmax .and. maxdif.ne.0) then
        !
        !  determine which algorithm to use based on user-supplied 
        !  binary scale factor and number of bits.
        !
        if (idrstmpl(2).eq.0) then
           !
           !  no binary scaling and calculate minimum number of 
           !  bits in which the data will fit.
           !
           imin=nint(rmin*dscale)
           imax=nint(rmax*dscale)
           maxdif=imax-imin
           temp=alog(real(maxdif+1))/alog(2.0)
           nbits=ceiling(temp)
           rmin=real(imin)
           !   scale data
           do j=1,ndpts
             ifld(j)=nint(fld(j)*dscale)-imin
           enddo
        else
           !
           !  use binary scaling factor and calculate minimum number of 
           !  bits in which the data will fit.
           !
           rmin=rmin*dscale
           rmax=rmax*dscale
           maxdif=nint((rmax-rmin)*bscale)
           temp=alog(real(maxdif+1))/alog(2.0)
           nbits=ceiling(temp)
           !   scale data
           do j=1,ndpts
             ifld(j)=nint(((fld(j)*dscale)-rmin)*bscale)
           enddo
        endif
        !
        !  pack data into full octets, then do png encode.
        !  and calculate the length of the packed data in bytes
        !
        if (nbits.le.8) then
            nbits=8
        elseif (nbits.le.16) then
            nbits=16
        elseif (nbits.le.24) then
            nbits=24
        else 
            nbits=32
        endif
        nbytes=(nbits/8)*ndpts
        allocate(ctemp(nbytes))
        call sbytes(ctemp,ifld,0,nbits,0,ndpts)
        !
        !  encode data into png format.
        !
        lcpack=enc_png(ctemp,width,height,nbits,cpack)
        if (lcpack.le.0) then
           print *,'pngpack: error encoding png = ',lcpack
        endif
        deallocate(ctemp)

      else
        nbits=0
        lcpack=0
      endif

!
!  fill in ref value and number of bits in template 5.0
!
      call mkieee(rmin,ref,1)   ! ensure reference value is ieee format
!      call gbyte(ref,idrstmpl(1),0,32)
      iref=transfer(ref,iref)
      idrstmpl(1)=iref
      idrstmpl(4)=nbits
      idrstmpl(5)=0         ! original data were reals


#endif /* use_png */
      return
      end
